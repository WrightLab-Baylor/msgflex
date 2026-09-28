#!/usr/bin/env python3
import pandas as pd
import numpy as np
from typing import Sequence
import copy
from scipy.stats import median_abs_deviation, pearsonr
from xgboost import XGBClassifier

#Helper functions
def _accepted_targets(df, qcol, is_decoy_col, q_thresh):
    m = (
        (df[is_decoy_col] == 0) &
        (df[qcol] <= q_thresh)
    )
    return df.loc[m]

def _compute_sample_weights(df, *, fold_auc = None):
    n = len(df)
    w = np.ones(n, dtype=float)
    pos = (df['is_decoy'] == 0).values

    # a noisy target looks bad by every metric but is still a true positive.
    # Reduce aggressiveness of downweighting proportionally to label noise.
    if fold_auc is not None and fold_auc < 0.75:
        noise_scale = fold_auc / 0.75 
    else:
        noise_scale = 1.0

    # 1) Non-top-ranked targets
    if 'is_top_ranked' in df.columns:
        is_not_top = (df['is_top_ranked'].fillna(0).values == 0)
        penalty = 0.5 + (1.0 - 0.5) * (1.0 - noise_scale)   # 0.5 clean → 0.9 noisy
        w[pos & is_not_top] *= penalty

    # 2) Delta score — same logic, reduced penalty under noise
    if 'delta_score' in df.columns:
        ds = pd.to_numeric(df['delta_score'], errors='coerce').fillna(0.0).values
        scale = max(np.percentile(ds[pos], 90) if pos.sum() > 0 else 1.0, 1e-6)
        ds_scaled = (np.tanh(ds / scale) + 1.0) / 2.0
        floor = 0.7 + 0.3 * (1.0 - noise_scale)   # 0.7 clean → 1.0 noisy (no penalty)
        w[pos] *= floor + (1.0 - floor) * ds_scaled[pos]

    # 3) Mass error — same
    if 'mass_error_zscore' in df.columns:
        me = pd.to_numeric(df['mass_error_zscore'], errors='coerce').fillna(0.0).values
        decay = 0.15 * noise_scale   # 0.15 clean → 0.09 noisy (softer penalty)
        w[pos] *= np.exp(-decay * np.clip(me[pos], 0, 10))

    np.clip(w, 0.3, 2.0, out=w)
    return w

def _infer_is_decoy(df):
    if 'is_decoy' in df.columns:
        return df['is_decoy'].astype(int)
    if 'Label' in df.columns:
        return (df['Label'] == -1).astype(int)
    if 'Protein' in df.columns:
        # Common decoy markers
        return df['Protein'].astype(str).str.contains(r'DECOY|REV_|XXX_', case=False, regex=True, na=False).astype(int)
    raise ValueError("Cannot infer decoys: need 'is_decoy' or 'Label' or 'Protein' column.")

# Leakage-free peptide-level features
def _pep_level_features(train_df, apply_df):
    """
    Leakage-free peptide-level aggregated features.
 
    Statistics are derived from the training fold only and mapped onto
    apply_df.  Targets and decoys are handled separately to prevent
    fallback-value inversion artifacts.
    """
    apply_df = apply_df.copy()
 
    _DEFAULTS = {
        'pep_best_specev': 0.0,
        'pep_n_observations': 0.0,
        'pep_score_consistency': 0.0,
        'pep_charge_diversity': 1.0,
        'pep_score_relative_to_best': 0.0,
    }
 
    # Safe guard: If essential columns are absent, fill with safe defaults and return
    if 'neg_log_specEvalue' not in train_df.columns or 'Peptide' not in train_df.columns:
        for col, val in _DEFAULTS.items():
            apply_df[col] = val
        return apply_df
 
    # -------------------------------------------------------------------------
    # Inner helper: compute per-peptide stats from a subset of the training fold
    # -------------------------------------------------------------------------
    def _compute_pep_stats(subset_df):
        if subset_df.empty:
            return pd.DataFrame(
                columns=[
                    'pep_best_specev', 'pep_n_observations',
                    'pep_score_consistency', 'pep_charge_diversity',
                ]
            )
        stats = (
            subset_df
            .groupby('Peptide', sort=False)
            .agg(
                pep_best_specev =('neg_log_specEvalue', 'max'),
                pep_n_observations =('neg_log_specEvalue', 'count'),
                pep_score_consistency =('neg_log_specEvalue', 'std'),
                pep_charge_diversity  =('Charge', 'nunique'),
            )
        )
        # std is NaN for singletons — treat as perfectly consistent (0)
        stats['pep_score_consistency'] = stats['pep_score_consistency'].fillna(0.0)
        return stats
 
    # -------------------------------------------------------------------------
    # Compute stats separately for targets and decoys (leakage-safe)
    # -------------------------------------------------------------------------
    train_targets = train_df[train_df['is_decoy'] == 0]
    train_decoys = train_df[train_df['is_decoy'] == 1]
 
    target_stats = _compute_pep_stats(train_targets)
    decoy_stats = _compute_pep_stats(train_decoys)
 
    target_fallback_best = float(
        train_targets['neg_log_specEvalue'].quantile(0.25)
    ) if len(train_targets) > 0 else 0.0
 
    decoy_fallback_best = float(
        train_decoys['neg_log_specEvalue'].quantile(0.25)
    ) if len(train_decoys) > 0 else 0.0
 
    # -------------------------------------------------------------------------
    # Map statistics onto apply_df, respecting target / decoy identity
    # -------------------------------------------------------------------------
    is_decoy_mask  = apply_df['is_decoy'] == 1
    is_target_mask = ~is_decoy_mask
 
    stat_cols = [
        'pep_best_specev',
        'pep_n_observations',
        'pep_score_consistency',
        'pep_charge_diversity',
    ]
 
    for col in stat_cols:
        apply_df[col] = 0.0  # initialise to a neutral value
 
        # --- target PSMs ---
        if not target_stats.empty and col in target_stats.columns:
            fallback = (
                target_fallback_best if col == 'pep_best_specev'
                else _DEFAULTS.get(col, 0.0)
            )
            apply_df.loc[is_target_mask, col] = (
                apply_df.loc[is_target_mask, 'Peptide']
                .map(target_stats[col])
                .fillna(fallback)
                .values
            )
 
        # --- decoy PSMs ---
        if not decoy_stats.empty and col in decoy_stats.columns:
            fallback = (
                decoy_fallback_best if col == 'pep_best_specev'
                else _DEFAULTS.get(col, 0.0)
            )
            apply_df.loc[is_decoy_mask, col] = (
                apply_df.loc[is_decoy_mask, 'Peptide']
                .map(decoy_stats[col])
                .fillna(fallback)
                .values
            )
 
    apply_df['pep_n_observations'] = np.log1p(apply_df['pep_n_observations'])
 
    # -------------------------------------------------------------------------
    # Derived feature: how does this PSM compare to the peptide's best training score?
    # -------------------------------------------------------------------------
    if 'neg_log_specEvalue' in apply_df.columns:
        apply_df['pep_score_relative_to_best'] = (
            apply_df['neg_log_specEvalue'] - apply_df['pep_best_specev']
        )
    else:
        apply_df['pep_score_relative_to_best'] = 0.0
 
    return apply_df


def _xgb_gpu_cpu(model_params):
    """
    Initialize XGBClassifier and safely fall back to CPU if GPU is not available.
    """
    safe_params = copy.deepcopy(model_params)
    safe_params["early_stopping_rounds"] = None

    try:
        clf = XGBClassifier(**safe_params)

        # Synthetic 2-class dummy data
        X_dummy = np.zeros((2, 2), dtype=np.float32)
        y_dummy = np.array([0, 1], dtype=np.int32)

        clf.fit(X_dummy, y_dummy, verbose=False)

        # Recreate the classifier with full real parameters
        return XGBClassifier(**model_params)

    except Exception as e:
        if "cuda" in str(e).lower() or "gpu" in str(e).lower():
            print("[WARN] GPU unavailable — falling back to CPU")

            model_params["device"] = "cpu"

            # Recreate real classifier (CPU version)
            return XGBClassifier(**model_params)

        raise

# --- Multiprocessing helpers ---
def _worker_init(spectra_dict, extractor_obj):
    """Store shared read-only state in each worker process once."""
    global _WORKER_SPECTRA, _WORKER_EXTRACTOR
    _WORKER_SPECTRA   = spectra_dict
    _WORKER_EXTRACTOR = extractor_obj

def _worker_fn(args):
    """Called once per PSM in each worker process."""
    scan_num, peptide, charge = args
    if scan_num is None or int(scan_num) not in _WORKER_SPECTRA:
        return {"ScanNum": scan_num, "spectral_feature_error": "scan_not_found"}
    return _WORKER_EXTRACTOR._psm_features(
        spectrum = _WORKER_SPECTRA[int(scan_num)],
        peptide  = peptide,
        charge   = int(charge),
        scan_num = int(scan_num),
    )

def _rt_stats(df: pd.DataFrame, group: str) -> dict:
    """
    Compute RT summary statistics for a PSM group (common / gain / lost).
    Returns a flat dict with keys prefixed by group name for the qc_row.
    """
    prefix = f"rt_{group}"

    # Handle empty or missing DataFrame safely
    if df is None or df.empty:
        return {
            f"rt_median_{prefix}": float("nan"),
            f"rt_mad_{prefix}": float("nan"),
            f"rt_pearson_{prefix}": float("nan"),
        }
    # Ensure targeted columns exist in DataFrame
    if "observed_rt" not in df.columns or "predicted_rt" not in df.columns:
        return {
            f"rt_median_{prefix}": float("nan"),
            f"rt_mad_{prefix}": float("nan"),
            f"rt_pearson_{prefix}": float("nan"),
        }

    obs = df["observed_rt"].dropna()
    pred = df["predicted_rt"].dropna()

    # Align on shared index after dropna
    idx  = obs.index.intersection(pred.index)
    obs  = obs.loc[idx].astype(float).to_numpy()
    pred = pred.loc[idx].astype(float).to_numpy()

    # pearsonr requires at least 2 data points
    if len(obs) < 2:
        return {
            f"rt_median_{prefix}": float("nan"),
            f"rt_mad_{prefix}": float("nan"),
            f"rt_pearson_{prefix}": float("nan"),
        }
    # Calculate metrics
    residuals = obs - pred
    median_val = np.median(residuals)
    mad_val = median_abs_deviation(residuals)
    r_val, _ = pearsonr(obs, pred)

    return {
        f"rt_median_{prefix}": round(float(median_val), 4),
        f"rt_mad_{prefix}": round(float(mad_val), 4),
        f"rt_pearson_{prefix}": round(float(r_val), 4),
    }

def _aggregate_feature_importance(fold_importances):
    """
    Aggregate XGBoost feature importances across CV folds.
    """
    if not fold_importances:
        return pd.DataFrame(columns=['feature', 'importance_mean',
                                     'importance_std', 'importance_cv',
                                     'n_folds_present'])

    imp_df = pd.DataFrame(fold_importances).fillna(0.0)

    summary = pd.DataFrame({
        'feature': imp_df.columns,
        'importance_mean': imp_df.mean(),
        'importance_std': imp_df.std().fillna(0.0),
        'importance_cv': (imp_df.std() / (imp_df.mean() + 1e-9)).fillna(0.0),
        'n_folds_present': (imp_df > 0).sum(),
    }).sort_values('importance_mean', ascending=False).reset_index(drop=True)

    return summary