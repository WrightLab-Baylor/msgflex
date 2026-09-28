#!/usr/bin/env python3
import pandas as pd
import numpy as np
from typing import Dict, List
from .utils import _accepted_targets, _pep_level_features

def _validate_fdr_inputs(
    df, score_col, decoy_col
    ):
    """
    Raise useful errors if FDR inputs are malformed
    """
    missing = {score_col, decoy_col} - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if not pd.api.types.is_numeric_dtype(df[score_col]):
        raise TypeError(
            f"'{score_col}' must be numeric, got {df[score_col].dtype}."
        )

    if df[score_col].isna().any():
        raise ValueError(
            f"'{score_col}' has {df[score_col].isna().sum()} NaN(s). "
            "Drop or impute before calling estimate_fdr()."
        )

    invalid = ~df[decoy_col].isin([0, 1])
    if invalid.any():
        raise ValueError(
            f"{invalid.sum()} row(s) have '{decoy_col}' ∉ {{0, 1}}: "
            f"{df.loc[invalid, decoy_col].unique().tolist()}"
        )

    n_targets = int((df[decoy_col] == 0).sum())
    n_decoys = int((df[decoy_col] == 1).sum())
    if n_targets == 0:
        raise ValueError(
            "No target PSMs found. Check is_decoy encoding (0=target, 1=decoy)."
        )
    if n_decoys == 0:
        raise ValueError("No decoy PSMs found. TDA-FDR requires decoys.")

def _build_sort_keys(
    df,
    score_col,
    massdiff_col= None,
    ptm_col = None,
    ):
    """
    Build multi-key sort specification for FDR ranking
    """
    sort_keys = [score_col]
    ascending  = [False]
    if massdiff_col and massdiff_col in df.columns:
        sort_keys.append(massdiff_col)
        ascending.append(True)
    if ptm_col and ptm_col in df.columns:
        sort_keys.append(ptm_col)
        ascending.append(True)
    return sort_keys, ascending

# PSM-level FDR
def estimate_fdr(
    df,
    *,
    score_col = "rescore_prob",
    decoy_col = "is_decoy",
    massdiff_col = "DelM",
    ptm_col = "n_mods",
):
    """
    Estimate PSM-level FDR and q-values (WinnowNet formula).
    """
    _validate_fdr_inputs(df, score_col, decoy_col)
    sort_keys, ascending = _build_sort_keys(df, score_col, massdiff_col, ptm_col)

    sorted_df = df.sort_values(sort_keys, ascending=ascending, kind="stable")
    sort_idx = sorted_df.index

    is_decoy_s = sorted_df[decoy_col].to_numpy(dtype=np.int32)
    is_target_s = 1 - is_decoy_s

    cum_targets = np.cumsum(is_target_s)
    cum_decoys = np.cumsum(is_decoy_s)

    fdr_sorted = np.where(
        cum_targets > 0,
        cum_decoys / cum_targets,
        1.0,
    )
    fdr_sorted = np.clip(fdr_sorted, 0.0, 1.0)
    qval_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]

    out = df.copy()
    out.loc[sort_idx, "xpec_q"] = qval_sorted
    return out

# Peptide-level FDR
def estimate_pep_fdr(
    df,
    *,
    score_col = "rescore_prob",
    decoy_col = "is_decoy",
    massdiff_col = "DelM",
    ptm_col = "n_mods"
    ):
    """
    Estimate peptide-level FDR and q-values.
    """
    _validate_fdr_inputs(df, score_col, decoy_col)

    if 'Peptide' not in df.columns:
        raise ValueError("'Peptide' column required for peptide-level FDR.")

    sort_keys, ascending = _build_sort_keys(df, score_col, massdiff_col, ptm_col)
   
    # preservation across pandas versions.
    best_idx = df.groupby('Peptide')[score_col].idxmax()
    pep_df = df.loc[best_idx].copy()
    pep_df = pep_df.sort_values(sort_keys, ascending=ascending, kind="stable")

    is_target = (pep_df[decoy_col] == 0).astype(int).to_numpy()
    is_decoy_np = (pep_df[decoy_col] == 1).astype(int).to_numpy()

    cum_targets = np.cumsum(is_target)
    cum_decoys = np.cumsum(is_decoy_np)

    fdr_sorted = np.where(
        cum_targets > 0,
        cum_decoys / cum_targets,
        1.0,
    )
    fdr_sorted = np.clip(fdr_sorted, 0.0, 1.0)
    qval_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]

    pep_df["xpec_pepq"] = qval_sorted

    # Map peptide-level q-value back to every PSM row
    pepq_map = pep_df.set_index('Peptide')["xpec_pepq"]
    out = df.copy()
    out["xpec_pepq"] = out["Peptide"].map(pepq_map)
    return out

# Gain comparison
def compare_gains(
    df,
    *,
    msgf_q_col = None,
    resc_q_col = None,
    thresholds = (0.01, 0.05, 0.10),
    verbose = True,
):
    msgf_q_col = msgf_q_col or "QValue"
    resc_q_col = resc_q_col or "xpec_q"

    if resc_q_col not in df.columns:
        raise KeyError(
            "xpec_q column not found — run estimate_fdr() first."
        )

    DECOY_COL = "is_decoy"
    df = df.copy()
    df[DECOY_COL] = df[DECOY_COL].astype(int)
    df = df.dropna(subset=[msgf_q_col])

    improvements = {}
    for q in thresholds:
        msgf_acc = _accepted_targets(df, msgf_q_col, DECOY_COL, q)
        resc_acc = _accepted_targets(df, resc_q_col, DECOY_COL, q)

        n_msgf = len(msgf_acc)
        n_resc = len(resc_acc)
        gain = n_resc - n_msgf
        pct_gain = (
            gain / n_msgf * 100.0 if n_msgf > 0
            else (float("inf") if gain > 0 else 0.0)
        )

        improvements[q] = {
            "n_msgf_targets": n_msgf,
            "n_rescore_targets": n_resc,
            "gain": gain,
            "pct_gain_vs_msgf": pct_gain,
        }

        if verbose:
            pct_str = (
                f"{pct_gain:+.1f}%" if n_msgf > 0
                else ("+∞%" if gain > 0 else "0.0%")
            )
            print(
                f"{int(q*100):>2d}% FDR:  "
                f"MS-GF+={n_msgf:,}  "
                f"XPECTRA={n_resc:,}  "
                f"Gain={gain:+,} ({pct_str})"
            )

    return improvements

# Balanced scan folds
def make_balanced_scan_folds(df, *, n_folds=3, random_state=42):
    """
    Assign each unique OptimalScanNumber to one of n_folds so folds are
    balanced by scan composition (target-only, decoy-only, mixed).
    """
    rng = np.random.default_rng(random_state)

    g = df.groupby('OptimalScanNumber')['is_decoy']
    has_decoy  = g.max() > 0
    has_target = g.min() == 0

    scan_class = pd.Series('other', index=g.size().index)
    scan_class[(has_target) & (~has_decoy)] = 'target_only'
    scan_class[(~has_target) & (has_decoy)] = 'decoy_only'
    scan_class[(has_target) & (has_decoy)]  = 'mixed'

    scan_ids_by_class = {
        c: scan_class.index[scan_class == c].to_numpy()
        for c in ['target_only', 'decoy_only', 'mixed']
    }
    for c in scan_ids_by_class:
        rng.shuffle(scan_ids_by_class[c])

    fold_bins: Dict[int, List[int]] = {k: [] for k in range(n_folds)}
    for c in ['mixed', 'target_only', 'decoy_only']:
        for i, sid in enumerate(scan_ids_by_class[c]):
            fold_bins[i % n_folds].append(int(sid))

    scan_to_fold: Dict[int, int] = {
        int(sid): f
        for f, sids in fold_bins.items()
        for sid in sids
    }
    return scan_to_fold

# Per-fold PSM features 
def compute_fold_psm_features(
    train_df, apply_df
    ):
    """
    Leakage-free peptide-level features.
    Statistics derived from training fold only; applied to apply_df.
    """
    return _pep_level_features(train_df, apply_df)

# Top-1 per-spectrum metrics
def compute_top1_metrics(df, score_col, spectrum_col="OptimalScanNumber"):
    """
    Extract top-1 PSM per spectrum and return labels and scores.

    """
    top1_idx = df.groupby(spectrum_col)[score_col].idxmax()
    top1 = df.loc[top1_idx].copy()

    if top1["is_decoy"].dtype == bool:
        y_true = (~top1["is_decoy"]).astype(int).values
    else:
        try:
            dec = top1["is_decoy"].astype(int).values
        except Exception:
            dec = (
                top1["is_decoy"].astype(str).str.lower()
                .map({"decoy": 1, "target": 0})
                .values.astype(int)
            )
        y_true = 1 - dec

    y_score = top1[score_col].values
    return y_true, y_score

# RT metrics
def compute_rt_metrics(df, rt_pred):
    rt_true = df["observed_rt"].astype(float)
    error = np.abs(rt_true - rt_pred)
    mae = np.mean(error)
    mad = np.median(error)
    rt_range = rt_true.max() - rt_true.min()
    norm_mae = mae / rt_range if rt_range > 0 else np.nan
    return {"MAE": mae, "MAD": mad, "normalized_MAE": norm_mae}

# Noisy-OR evidence fusion
def _calibrated_evalue_prob(x, scale = None):
    if scale is None:
        median = np.median(x)
        scale  = np.log(2.0) / median if median > 0 else 1.0
    return np.clip(1.0 - np.exp(-scale * x), 0.0, 1.0)

def compute_noisy_or_score(
    df,
    *,
    ml_col = 'rescore_prob',
    evalue_col = 'neg_log_specEvalue',
    alpha = 1.0,
    beta = 1.0,
    rescue_gamma = 0.80,
    rescue_threshold = 0.70,
    evalue_scale = None
    ):
    """
    Probabilistic Evidence Fusion via Noisy-OR with an unconditional rescue floor.
    """
    p = df[ml_col].to_numpy(dtype=float).clip(0.0, 1.0)
    x = df[evalue_col].to_numpy(dtype=float)

    e = _calibrated_evalue_prob(x, evalue_scale)
    alpha = min(float(alpha), 1.0)
    beta  = min(float(beta),  1.0)

    # Base noisy-OR
    base = 1.0 - (1.0 - alpha * p) * (1.0 - beta * e)
    floor = np.where(e >= rescue_threshold, rescue_gamma * e, 0.0)

    return np.maximum(base, floor).clip(0.0, 1.0)

# Alpha tuning for noisy-OR
def tune_noisy_or_weights(
    df_scored,
    *,
    alpha_grid = None,
    beta = 1.0,
    verbose = True
    ):
    """
    Grid search over alpha (ML score weight) in the noisy-OR fusion,
    maximising PSM identifications at 1% FDR.

    beta (spectral E-value weight) is held fixed at 1.0.
    """
    if alpha_grid is None:
        alpha_grid = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]

    if verbose:
        print(f'[Tuning] beta={beta}  grid={alpha_grid}')

    records = []
    for alpha in alpha_grid:
        trial = df_scored.copy()
        trial['_trial_score'] = compute_noisy_or_score(
            trial,
            ml_col ='rescore_prob',
            evalue_col ='neg_log_specEvalue',
            alpha =alpha,
            beta =beta
        )
        trial = estimate_fdr(
            trial,
            score_col ='_trial_score',
            decoy_col ='is_decoy',
            massdiff_col ='DelM',
            ptm_col ='n_mods',
        )
        q_col = 'xpec_q' if 'xpec_q' in trial.columns else None
        n_ids = int((trial[q_col] < 0.01).sum()) if q_col else 0
        records.append({'alpha': alpha, 'ids_at_1pct_fdr': n_ids})

        if verbose:
            print(f' alpha={alpha:.2f} => {n_ids:,} PSMs at 1% FDR')

    results = pd.DataFrame(records)
    best_alpha = float(results.loc[results['ids_at_1pct_fdr'].idxmax(), 'alpha'])

    if verbose:
        best_n = int(results['ids_at_1pct_fdr'].max())
        print(f' Best: alpha = {best_alpha}  ({best_n:,} IDs)')

    return best_alpha, results