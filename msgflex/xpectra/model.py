#!/usr/bin/env python3
"""
XPECTRA - Cross-validated PSM rescoring for SPARX outputs
"""

from __future__ import annotations
import re
from typing import Dict, Iterable, Tuple
from functools import lru_cache
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, train_test_split
from scipy.spatial import cKDTree
import joblib, pathlib
import warnings
warnings.filterwarnings('ignore', category=UserWarning,  module='xgboost')
warnings.filterwarnings('ignore', category=FutureWarning, module='sklearn')
try:
    import xgboost as xgb
except Exception as e:
    raise RuntimeError("xgboost is required. pip install xgboost") from e

from .metrics import (
    estimate_fdr,
    estimate_pep_fdr,
    make_balanced_scan_folds,
    compute_fold_psm_features,             
    compute_noisy_or_score,
    tune_noisy_or_weights
)
from .utils import _compute_sample_weights, _xgb_gpu_cpu, _aggregate_feature_importance
from .xrtpred import RTCalibrator

# =============================================================================
# Feature Engineering
# =============================================================================
class FeatureEngineer:
    """
    Leakage-free feature engineering for PSM rescoring.
    """

    # Kyte-Doolittle — used for sequence composition features only
    HYDROPHOBICITY_SCALE = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }

    # Krokhin 300A RC coefficients — used for RT prediction
    KROKHIN_300A = {
        'A': 0.35, 'R': -1.5, 'N': -0.99, 'D': -1.0, 'C': 1.9,
        'Q': -0.7, 'E': -1.0, 'G': 0.0,   'H': -0.5, 'I': 4.8,
        'L': 4.8,  'K': -1.8, 'M': 3.4,   'F': 5.5,  'P': 0.0,
        'S': -0.5, 'T': -0.4, 'W': 6.5,   'Y': 5.2,  'V': 4.2
    }
    _KROKHIN_ARR = np.zeros(26, dtype=np.float32) 

    # Known modification RT offsets 
    MOD_RT_DELTA = {
        'Oxidation': -3.5,   # M+16
        'Phospho': -5.0,   # S/T/Y+80
        'Deamidation': -0.5,   # N/Q+1
        'Acetyl': +3.0,   # protein N-term
        'AcNoTMT': -2.5,   # K-187.15 - acetyl blocked TMTPro;
    }

    def __init__(self, *, verbose = False, rt_coelution_radius = None):
        self.verbose = verbose
        self.is_fitted = False

        # RT statistics -- learned from training data
        self.peptide_rt_median = {}
        self.peptide_rt_std = {}
        self.charge_rt_median = {}
        self.rt_prediction_models = {}
        self.rt_density_tree = None
        self.global_rt_median = 0.0
        self.rt_min_valid = 0.5

        # Mass error statistics -- learned per charge state
        self.charge_mass_error_stats = {}
        # Peptide-level score statistics -- leakage-free
        self.peptide_observation_count = {}
        self.peptide_best_score = {}
        self.peptide_median_score = {}
        self.peptide_score_std = {}
        self.peptide_reliability = {}
        self._rt_coelution_radius = 0.5   # default
        self._rt_coelution_radius_override = rt_coelution_radius

    # -------------------------------------------------------------------------
    # Static helpers
    # -------------------------------------------------------------------------
    @staticmethod
    def _calc_hydrophobicity(seq) :
        """
        Mean Krokhin RC value -- used for sequence composition features.
        """
        if not isinstance(seq, str) or not seq:
            return 0.0
        return float(
            sum(FeatureEngineer.KROKHIN_300A.get(a, 0.0) for a in seq)
            / max(len(seq), 1)
        )
    
    @staticmethod
    def _calc_hydrophobicity_vec(seqs: pd.Series) -> np.ndarray:
        """Vectorized mean Krokhin RC — replaces per-row apply."""
        arr = FeatureEngineer._KROKHIN_ARR
        result = np.zeros(len(seqs), dtype=np.float32)
        for i, seq in enumerate(seqs):
            if seq:
                idx = np.frombuffer(seq.encode('ascii'), dtype=np.uint8) - 65
                valid = (idx >= 0) & (idx < 26)
                if valid.any():
                    result[i] = arr[idx[valid]].mean()
        return result

    @staticmethod
    def _terminal_hydrophobicity(seq):
        if not seq or len(seq) < 2:
            return 0.0, 0.0
        n_term = FeatureEngineer.HYDROPHOBICITY_SCALE.get(seq[0], 0.0)
        c_term = FeatureEngineer.HYDROPHOBICITY_SCALE.get(seq[-1], 0.0)
        return n_term, c_term

    _PHRP_PAT = re.compile(r'^[A-Z\-]\.(.+)\.[A-Z\-]$')

    @staticmethod
    @lru_cache(maxsize=500_000)
    def _peptide_core(peptide):
        s = str(peptide)
        m = FeatureEngineer._PHRP_PAT.match(s)
        core = m.group(1) if m else s
        return re.sub(r'[^A-Z]', '', core.upper())

    @staticmethod
    def _mod_rt_correction(raw_peptide):
        """
        Cumulative RT correction (minutes) from modification annotations in the
        raw PHRP peptide string.
        """
        correction = 0.0
        s = str(raw_peptide)
        correction += len(re.findall(r'M\+15\.9|M\+16', s)) * FeatureEngineer.MOD_RT_DELTA['Oxidation']
        correction += len(re.findall(r'[STY]\+79\.|[STY]\+80', s)) * FeatureEngineer.MOD_RT_DELTA['Phospho']
        correction += len(re.findall(r'[NQ]\+0\.9|[NQ]\+1\.0', s)) * FeatureEngineer.MOD_RT_DELTA['Deamidation']
        correction += len(re.findall(r'(?:^|\.\s*)\+42\.0', s)) * FeatureEngineer.MOD_RT_DELTA['Acetyl']
        correction += len(re.findall(r'K-187\.1|K\-187\.15', s)) * FeatureEngineer.MOD_RT_DELTA['AcNoTMT']

        return correction
    
    @staticmethod
    @lru_cache(maxsize=500_000)
    def _achrom_predict( core_seq):
        """
        SSRCalc-style RT score using Krokhin 300A RC coefficients with
        length correction and N-terminal bonus (Krokhin 2004).
        Returns 0.0 for non-standard sequences.
        """
        if not core_seq:
            return 0.0
        if not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]+', core_seq):
            return 0.0
        rc_sum = sum(FeatureEngineer.KROKHIN_300A.get(a, 0.0) for a in core_seq)
        length_correction = -0.35 * np.log1p(len(core_seq))
        nterm_bonus = FeatureEngineer.KROKHIN_300A.get(core_seq[0], 0.0) * 0.42
        return rc_sum + length_correction + nterm_bonus

    # -------------------------------------------------------------------------
    # Learning methods (called only on training data)
    # -------------------------------------------------------------------------
    def _learn_rt_statistics(self, df):
        """Learn RT statistics from training targets with valid RT values only."""
        rt_col = 'ElutionTime' if 'ElutionTime' in df.columns else 'RetentionTime'
        if rt_col not in df.columns:
            return

        df = df.copy()
        df['observed_rt'] = pd.to_numeric(df[rt_col], errors='coerce')

        # Determine valid RT threshold from the data distribution
        rt_nonzero = df['observed_rt'].dropna()
        rt_nonzero = rt_nonzero[rt_nonzero > 0]
        self.rt_min_valid = float(np.nanpercentile(rt_nonzero, 2)) if len(rt_nonzero) > 0 else 0.5
        valid_rt_mask = (df['observed_rt'] >= self.rt_min_valid) & df['observed_rt'].notna()

        self.global_rt_median = float(df.loc[valid_rt_mask, 'observed_rt'].median())

        clean_peptide = df['Peptide'].apply(self._peptide_core)
        df['_core_seq'] = clean_peptide

        # Per-peptide and per-charge statistics
        pep_stats = df[valid_rt_mask].groupby('Peptide')['observed_rt'].agg(['median', 'std'])
        self.peptide_rt_median = pep_stats['median'].to_dict()
        self.peptide_rt_std  = pep_stats['std'].fillna(0.0).to_dict()

        chg_stats = df[valid_rt_mask].groupby('Charge')['observed_rt'].median()
        self.charge_rt_median = chg_stats.to_dict()

        #using XGB for bossting rt correlation
        target_mask = valid_rt_mask.copy()
        if 'is_decoy' in df.columns:
            target_mask = target_mask & (df['is_decoy'] == 0)
        
        if 'missed_cleavages' in df.columns:
            mc = pd.to_numeric(df['missed_cleavages'], errors='coerce').fillna(0)
            target_mask = target_mask & (mc <= 2)
        reg_df = df[target_mask].copy()

        if self.verbose:
            print(f"[RT] Calibration set after MC filter (≤2): n={len(reg_df):,}")

        reg_df['_achrom_rt'] = reg_df['_core_seq'].apply(self._achrom_predict)
        self.rt_calibrator = RTCalibrator(frac=0.25, min_points=30, verbose=self.verbose)

        try:
            self.rt_calibrator.fit(reg_df)
            if self.verbose:
                print(f"[INFO] Unified RT calibrator fitted (n={len(reg_df):,})")
        except ValueError as e:
            # Not enough data — calibrator stays None, fallback to charge median
            self.rt_calibrator = None
            if self.verbose:
                print(f'[WARN] RT calibration skipped: {e}')

    def _learn_peptide_statistics(self,df):
        """
        Learn peptide-level score statistics from training data.
        """
        if 'Peptide' not in df.columns or 'neg_log_specEvalue' not in df.columns:
            return
        fit_df = df[df['is_decoy'] == 0] if 'is_decoy' in df.columns else df
        pep_stats = (
            fit_df
            .groupby('Peptide')['neg_log_specEvalue']
            .agg(
                peptide_psm_count='count',
                peptide_best_score='max',
                peptide_median_score='median',
                peptide_score_std='std',
            )
        )

        self.peptide_observation_count = pep_stats['peptide_psm_count'].to_dict()
        self.peptide_best_score = pep_stats['peptide_best_score'].to_dict()
        self.peptide_median_score= pep_stats['peptide_median_score'].to_dict()
        self.peptide_score_std = pep_stats['peptide_score_std'].fillna(0.0).to_dict()
        self._global_speceval_median = float(fit_df['neg_log_specEvalue'].median())
        self._global_speceval_best = float(fit_df['neg_log_specEvalue'].quantile(0.75))

        # Reliability: best score normalised by observation count — rewards
        # peptides seen many times with consistently high scores
        reliability = (
            pep_stats['peptide_best_score']
            / (pep_stats['peptide_psm_count'].clip(lower=1) + 1e-6)
        )
        self.peptide_reliability = reliability.to_dict()

        # Build sorted reference arrays for leakage-free pctrank (Fix 2B)
        if 'peptide_length' in df.columns and 'Charge' in df.columns:
            tmp = df[['Charge', 'peptide_length', 'neg_log_specEvalue']].copy()
            tmp['_len_bin'] = pd.cut(
                tmp['peptide_length'], bins=[0, 7, 11, 16, 100], labels=[0, 1, 2, 3]
            ).astype(str)
            tmp['_key'] = tmp['Charge'].astype(str) + '_' + tmp['_len_bin']
            self._speceval_bin_ref: Dict[str, np.ndarray] = {
                key: grp['neg_log_specEvalue'].dropna().sort_values().to_numpy()
                for key, grp in tmp.groupby('_key')
                if len(grp) >= 5
            }
        else:
            self._speceval_bin_ref = {}

  
    def _apply_rt_features(self,df):
        """Apply learned RT statistics. Invalid-RT rows receive neutral (0) values."""
        rt_col = 'ElutionTime' if 'ElutionTime' in df.columns else 'RetentionTime'
        if rt_col not in df.columns:
            return df

        df = df.copy()
        df['observed_rt'] = pd.to_numeric(df[rt_col], errors='coerce')

        rt_min = getattr(self, 'rt_min_valid', 0.5)
        rt_valid = (df['observed_rt'] >= rt_min) & df['observed_rt'].notna()
        df['rt_is_valid'] = rt_valid.astype(int)
        if '_core_seq' in df.columns:
            clean_peptide = df['_core_seq']
        else:
            clean_peptide = df['Peptide'].apply(self._peptide_core)

        for _aa, _val in FeatureEngineer.KROKHIN_300A.items():
            FeatureEngineer._KROKHIN_ARR[ord(_aa) - 65] = _val
        
        df['peptide_hydrophobicity'] = self._calc_hydrophobicity_vec(clean_peptide)

        # Per-peptide RT lookup
        df['peptide_rt_median'] = df['Peptide'].map(self.peptide_rt_median).fillna(self.global_rt_median)
        df['peptide_rt_std'] = df['Peptide'].map(self.peptide_rt_std).fillna(0.0)
        df['rt_deviation_from_peptide'] = (
            (df['observed_rt'] - df['peptide_rt_median']).abs()
        ).where(rt_valid, 0.0)

        # Per-charge mass error z-score
        if 'abs_mass_error_ppm' in df.columns and self.charge_mass_error_stats:
            medians = df['Charge'].map({k: v['median'] for k, v in self.charge_mass_error_stats.items()}).fillna(0.0)
            stds = df['Charge'].map({k: v['std']    for k, v in self.charge_mass_error_stats.items()}).fillna(1.0)
            df['mass_error_zscore'] = ((df['abs_mass_error_ppm'] - medians) / (stds + 1e-6)).abs()

        # Per-charge RT deviation
        df['rt_expected_for_charge'] = df['Charge'].map(self.charge_rt_median).fillna(self.global_rt_median)
        df['rt_deviation_from_charge'] = (
            (df['observed_rt'] - df['rt_expected_for_charge']).abs()
        ).where(rt_valid, 0.0)

        # Determine missed-cleavage group per row for calibrator routing
        achrom_preds = clean_peptide.apply(self._achrom_predict).to_numpy(dtype=float)

        # Build prediction dataframe regardless of calibrator availability
        pred_df = df[['Peptide', 'Charge']].copy()
        pred_df['_achrom_rt'] = achrom_preds
        pred_df['_core_seq']  = clean_peptide
        pred_df['observed_rt'] = df['observed_rt']

        cal = getattr(self, 'rt_calibrator', None)

        if cal is not None:
            # Predict RT and match probability
            res = cal.predict(pred_df)
            df["predicted_rt"] = res["rt_pred"].to_numpy(dtype=float)
            df["rt_match_prob"] = res["rt_score"].to_numpy(dtype=float)

        else:
            # Fallback: charge‑median RT prediction
            df["predicted_rt"] = (
                df['Charge'].map(self.charge_rt_median)
                            .fillna(self.global_rt_median)
                            .to_numpy(dtype=float)
            )
            df["rt_match_prob"] = 0.5

        # RT residuals — zero out for invalid rows
        raw_delta = (df['observed_rt'] - df['predicted_rt']).clip(-100, 100)
        df['rt_delta'] = raw_delta.where(rt_valid, 0.0)
        df['raw_rt_delta'] = raw_delta.abs().where(rt_valid, np.nan)

        # Winsorise abs_rt_delta at 90th percentile of valid rows
        target_mask = rt_valid & (df.get('is_decoy', pd.Series(0, index=df.index)) == 0)
        if target_mask.sum() > 10:
            winsor_cap = float(raw_delta.abs().where(target_mask).quantile(0.90))
        elif rt_valid.sum() > 10:
            winsor_cap = float(raw_delta.abs().where(rt_valid).quantile(0.90))
        else:
            winsor_cap = 40.0
        
        df['abs_rt_delta'] = raw_delta.abs().clip(upper=winsor_cap).where(rt_valid, 0.0)

        # RT per residue
        df['peptide_length'] = clean_peptide.str.len()
        df['rt_per_residue'] = (
            (df['observed_rt'] / df['peptide_length'].clip(lower=1))
        ).where(rt_valid, 0.0)

        # Interaction features
        df['hydrophobicity_x_charge'] = df['peptide_hydrophobicity'] * df['Charge']
        df['rt_x_hydrophobicity'] = (df['observed_rt'] * df['peptide_hydrophobicity']).where(rt_valid, 0.0)

        # Co-elution density from training-fold tree (leakage-free)
        if self.rt_density_tree is not None:
            rt_vals = df['observed_rt'].fillna(self.global_rt_median).to_numpy().reshape(-1, 1)
            counts = self.rt_density_tree.query_ball_point(
                rt_vals, r=self._rt_coelution_radius, return_length=True
            ).astype(np.float32)
            df['rt_coelution_density'] = counts
            df['log_rt_coelution_density'] = np.log1p(counts)

        # Precursor m/z deviation
        if 'ParentIonMZ' in df.columns and hasattr(self, 'charge_precursor_mz_stats'):
            parent_mz = pd.to_numeric(df['ParentIonMZ'], errors='coerce')
            mz_medians = df['Charge'].map({k: v['median'] for k, v in self.charge_precursor_mz_stats.items()}).fillna(parent_mz.median())
            mz_stds = df['Charge'].map({k: v['std'] for k, v in self.charge_precursor_mz_stats.items()}).fillna(1.0)
            df['precursor_mz_deviation'] = (parent_mz - mz_medians) / (mz_stds + 1e-6)

        # Relative parent ion intensity
        if 'ParentIonIntensity' in df.columns and hasattr(self, 'train_log_intensity_median'):
            intensity = pd.to_numeric(df['ParentIonIntensity'], errors='coerce').clip(lower=0)
            df['relative_parent_intensity'] = np.log1p(intensity) - self.train_log_intensity_median

        # Mask invalid observed RT so zeros
        df.loc[~rt_valid, 'observed_rt'] = np.nan

        return df

    def _apply_peptide_statistics(self, df) :
        """
        Apply learned peptide-level score statistics. 
        Used by both fit_transform and transform.
        """
        if 'neg_log_specEvalue' not in df.columns:
            for col in ['peptide_observation_count', 'peptide_best_score',
                        'peptide_median_score', 'peptide_score_std',
                        'peptide_reliability', 'is_best_for_peptide']:
                df[col] = 0.0
            return df

        df['peptide_observation_count'] = df['Peptide'].map(
            self.peptide_observation_count
        ).fillna(1).clip(upper=5)

        _best_prior   = getattr(self, '_global_speceval_best',   0.0)
        _median_prior = getattr(self, '_global_speceval_median', 0.0)

        df['peptide_best_score'] = df['Peptide'].map(
            self.peptide_best_score
        ).fillna(_best_prior)

        df['peptide_median_score'] = df['Peptide'].map(
            self.peptide_median_score
        ).fillna(_median_prior)

        df['peptide_score_std'] = df['Peptide'].map(
            self.peptide_score_std
        ).fillna(0.0)

        df['peptide_reliability'] = df['Peptide'].map(
            self.peptide_reliability
        ).fillna(0.0)

        df['is_best_for_peptide'] = (
            (df['neg_log_specEvalue'] >= df['peptide_best_score'] - 1e-6) &
            (df['peptide_observation_count'] <= 3)
        ).astype(int)

        return df

    def _apply_speceval_rank(self,df):
        """
        Assign each PSM a percent-rank within its charge×length bin, anchored
        to the training-fold score distribution (not the current fold's).
        Prevents distribution shift on this feature across CV folds.
        """
        if not getattr(self, '_speceval_bin_ref', None) or 'neg_log_specEvalue' not in df.columns:
            return df
        df = df.copy()
        df['_len_bin'] = pd.cut(
            df['peptide_length'], bins=[0, 7, 11, 16, 100], labels=[0, 1, 2, 3]
        ).astype(str)
        df['_key'] = df['Charge'].astype(str) + '_' + df['_len_bin']

        pctrank = np.full(len(df), 0.5, dtype=float)
        for key, row_idx in df.groupby('_key').groups.items():
            ref = self._speceval_bin_ref.get(key)
            if ref is None or len(ref) == 0:
                continue   # unseen bin → leave at 0.5 neutral
            vals = df.loc[row_idx, 'neg_log_specEvalue'].to_numpy()
            pctrank[df.index.get_indexer(row_idx)] = np.searchsorted(ref, vals) / len(ref)

        df['specEval_pctrank_in_class'] = pctrank
        df.drop(columns=['_len_bin', '_key'], inplace=True)
        return df

    def save(self, path):
        """
        fe.save('feature_engineer.joblib')
        fe2 = FeatureEngineer.load('feature_engineer.joblib')
        """
        if not self.is_fitted:
            raise RuntimeError("FeatureEngineer.save() called before fit().")
        pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, path, compress=3)

    @classmethod
    def load(cls, path) -> "FeatureEngineer":
        """Load a FeatureEngineer previously saved with save()."""
        obj = joblib.load(path)
        if not isinstance(obj, cls):
            raise TypeError(
                f"File contains {type(obj).__name__}, expected FeatureEngineer."
            )
        return obj

    # Feature engineering
    def engineer_features(self,df):
        """
        Create all sequence- and spectral-derived features that do not require
        training statistics. Safe to call on any fold without leakage.
        """
        df = df.copy()
        clean_peptide = df['Peptide'].apply(self._peptide_core)
        df['_core_seq'] = clean_peptide 

        # MS-GF+ score transforms
        if 'MSGFScore' in df.columns:
            df['log_msgf_score'] = np.log1p(pd.to_numeric(df['MSGFScore'], errors='coerce').clip(lower=0))

        if 'MSGFDB_SpecEValue' in df.columns:
            spec_eval = pd.to_numeric(df['MSGFDB_SpecEValue'], errors='coerce').replace(0, np.nan)
            min_nonzero = spec_eval[spec_eval > 0].min() if (spec_eval > 0).any() else 1e-100
            spec_eval = spec_eval.fillna(min_nonzero / 10)
            df['neg_log_specEvalue'] = -np.log10(spec_eval)
            df['specevalue_rank'] = spec_eval.rank(method='dense')
            
        if 'DeNovoScore' in df.columns:
            df['log_denovo_score'] = np.log1p(pd.to_numeric(df['DeNovoScore'], errors='coerce').clip(lower=0))

        if 'DelM_PPM' in df.columns:
            df['abs_mass_error_ppm'] = df['DelM_PPM'].abs()
            df['signed_mass_error']  = df['DelM_PPM']

        if 'DelM' in df.columns:
            df['abs_delm_da'] = pd.to_numeric(df['DelM'], errors='coerce').fillna(0).abs()

        if 'FWHMInScans' in df.columns:
            df['log_fwhm'] = np.log1p(pd.to_numeric(df['FWHMInScans'], errors='coerce').clip(lower=0))

        if 'PeakArea' in df.columns:
            df['log_peak_area'] = np.log1p(pd.to_numeric(df['PeakArea'], errors='coerce').clip(lower=0))

        # Sequence-derived features
        clean_peptide = df['_core_seq'] 
        df['peptide_length'] = clean_peptide.str.len()
        df['charge_to_length_ratio'] = df['Charge'] / df['peptide_length'].clip(lower=1)
        df['num_basic'] = clean_peptide.str.count(r'[KRH]')
        df['num_acidic'] = clean_peptide.str.count(r'[DE]')
        df['hydrophobic_ratio'] = clean_peptide.str.count(r'[AILMFWV]') / df['peptide_length'].clip(lower=1)

        _HYDRO_MAP = FeatureEngineer.HYDROPHOBICITY_SCALE   # alias for clarity
        df['n_term_hydro'] = clean_peptide.str[0].map(_HYDRO_MAP).fillna(0.0)
        df['c_term_hydro'] = clean_peptide.str[-1].map(_HYDRO_MAP).fillna(0.0)

        # c_half_hydro: build a class-level frozenset once, use str.count per aa
        _HYDRO_POS = frozenset(aa for aa, v in _HYDRO_MAP.items() if v > 0)
        lengths = clean_peptide.str.len().clip(lower=2)
        
        # Correct per-row slicing needs a vectorized approach:
        half_starts = (lengths // 2).to_numpy()
        seqs = clean_peptide.to_numpy()
        c_half_hydro = np.array([
            sum(1 for c in s[h:] if c in _HYDRO_POS) / max(len(s) - h, 1)
            for s, h in zip(seqs, half_starts)
        ], dtype=float)
        df['c_half_hydro'] = c_half_hydro
        
        df['missed_cleavages'] = (
            clean_peptide.str.count(r'[KR](?!P)')
            - clean_peptide.str.endswith(('K', 'R')).astype(int)
        ).clip(lower=0)

        df['nterm_is_basic'] = clean_peptide.str[0].isin(list('KRH')).astype(int)
        df['cterm_is_KR'] = clean_peptide.str[-1].isin(['K', 'R']).astype(int)
        df['has_KP_or_RP'] = clean_peptide.str.contains(r'[KR]P', regex=True).astype(int)

        if 'peptide_length' in df.columns:
            df['charge_x_length'] = df['Charge'] * df['peptide_length']

        if 'abs_mass_error_ppm' in df.columns:
            df['charge_x_mass_error'] = df['Charge'] * df['abs_mass_error_ppm']

        # Within-scan competitive features
        if 'OptimalScanNumber' in df.columns and 'neg_log_specEvalue' in df.columns:
            scan_groups = df.groupby('OptimalScanNumber')['neg_log_specEvalue']
            scan_sorted = df[['OptimalScanNumber', 'neg_log_specEvalue']].copy()
            scan_sorted['_rank'] = scan_groups.rank(method='first', ascending=False)

            r1_map = (
                scan_sorted[scan_sorted['_rank'] == 1]
                .set_index('OptimalScanNumber')['neg_log_specEvalue']
            )
            r2_map = (
                scan_sorted[scan_sorted['_rank'] == 2]
                .set_index('OptimalScanNumber')['neg_log_specEvalue']
            )

            df['_r1'] = df['OptimalScanNumber'].map(r1_map)
            df['_r2'] = df['OptimalScanNumber'].map(r2_map)

            df['delta_score'] = (df['_r1'] - df['_r2']).fillna(-1.0)
            df['is_singleton_scan'] = df['delta_score'].eq(-1.0).astype(int)
            df['delta_score_norm'] = (df['delta_score'] / (df['_r1'].abs() + 1e-8)).clip(0, 1)
            df['n_psms_in_scan'] = scan_groups.transform('size')
            df['log_n_psms_in_scan'] = np.log1p(df['n_psms_in_scan'])
            df.drop(columns=['_r1', '_r2'], inplace=True)

        # Peptide-level interaction features
        if 'pep_is_singleton' in df.columns and 'neg_log_specEvalue' in df.columns:
            df['singleton_x_specEvalue'] = df['pep_is_singleton'] * df['neg_log_specEvalue']

        if 'pep_charge_diversity' in df.columns and 'pep_n_observations' in df.columns:
            df['pep_confidence'] = df['pep_charge_diversity'] * np.log1p(df['pep_n_observations'])

        # Spectral features
        if 'signal_to_noise' in df.columns and 'total_ion_coverage' in df.columns:
            df['snr_x_coverage'] = df['signal_to_noise'] * df['total_ion_coverage']

        if 'total_ion_current' in df.columns:
            df['log_total_ion_current'] = np.log1p(df['total_ion_current'])

        if 'base_peak_intensity' in df.columns:
            df['log_base_peak_intensity'] = np.log1p(df['base_peak_intensity'])

        if 'mean_matched_ion_intensity' in df.columns:
            df['log_mean_matched_intensity'] = np.log1p(df['mean_matched_ion_intensity'])

        if 'total_ion_coverage' in df.columns and 'matched_intensity_fraction' in df.columns:
            df['coverage_x_intensity'] = df['total_ion_coverage'] * df['matched_intensity_fraction']

        if 'b_ion_coverage' in df.columns and 'y_ion_coverage' in df.columns:
            df['min_by_coverage'] = df[['b_ion_coverage', 'y_ion_coverage']].min(axis=1)
            df['max_by_coverage'] = df[['b_ion_coverage', 'y_ion_coverage']].max(axis=1)
            df['by_coverage_ratio'] = df['min_by_coverage'] / (df['max_by_coverage'] + 1e-6)
            df['by_ion_balance'] = (df['b_ion_coverage'] - df['y_ion_coverage']).abs()

        if 'longest_b_series' in df.columns and 'longest_y_series' in df.columns:
            df['max_series_length']   = df[['longest_b_series', 'longest_y_series']].max(axis=1)
            df['series_length_ratio'] = df['max_series_length'] / df['peptide_length'].clip(lower=1)

        if 'mean_fragment_mass_error_ppm' in df.columns and 'abs_mass_error_ppm' in df.columns:
            df['precursor_fragment_error_ratio'] = (
                df['abs_mass_error_ppm'] / (df['mean_fragment_mass_error_ppm'] + 1)
            )

        # Replace inf/nan with column median -Vectorised
        num_cols = df.select_dtypes(include=[np.number]).columns
        df[num_cols] = df[num_cols].replace([np.inf, -np.inf], np.nan)
        medians = df[num_cols].median() 
        df[num_cols] = df[num_cols].fillna(medians)

        return df

    # -------------------------------------------------------------------------
    # Public API configuration
    # -------------------------------------------------------------------------
    def fit(self, df):
        """Learn all statistics from training data (after engineer_features has run)."""
        self._learn_rt_statistics(df)
        self._learn_peptide_statistics(df)

        if 'abs_mass_error_ppm' in df.columns:
            stats = (
                df.groupby('Charge')['abs_mass_error_ppm']
                .agg(median='median', std=lambda x: x.std(ddof=0))
            )
            self.charge_mass_error_stats = stats.to_dict('index')

        # RT density tree — valid rows only to prevent zero-RT contamination
        rt_col = 'ElutionTime' if 'ElutionTime' in df.columns else 'RetentionTime'
        if rt_col in df.columns:
            rt = pd.to_numeric(df[rt_col], errors='coerce')
            rt_valid = rt >= self.rt_min_valid
            target_mask = rt_valid & (df['is_decoy'] == 0) if 'is_decoy' in df.columns else rt_valid
            rt_target = rt.where(target_mask).dropna().to_numpy().reshape(-1, 1)

            if self._rt_coelution_radius_override is not None:
                self._rt_coelution_radius = self._rt_coelution_radius_override
            else:
                rt_median = float(np.median(rt_target)) if len(rt_target) > 0 else 30.0
                self._rt_coelution_radius = 0.5 if rt_median < 300 else 30.0

            if self.verbose:
                unit = 'min' if self._rt_coelution_radius < 5 else 's'
                print(f'[RT] co-elution radius = {self._rt_coelution_radius} {unit} '
                    f'(RT median={np.median(rt_target) if len(rt_target) > 0 else 0.0:.1f})')

            self.rt_density_tree = cKDTree(rt_target) if len(rt_target) >= 10 else None
            
        if 'ParentIonMZ' in df.columns:
            stats = (
                df.groupby('Charge')['ParentIonMZ']
                .agg(median='median', std=lambda x: x.std(ddof=0))
            )
            self.charge_precursor_mz_stats = stats.to_dict('index')

        if 'ParentIonIntensity' in df.columns:
            intensity = pd.to_numeric(df['ParentIonIntensity'], errors='coerce').clip(lower=0)
            self.train_log_intensity_median = float(np.log1p(intensity).median())

        self.is_fitted = True
        return self

    def fit_transform(self, df):
        """
        Fit on and transform training data.
        Applies engineer_features --> fit --> _apply_rt_features --> _apply_peptide_statistics
        """
        df = self.engineer_features(df)
        self.fit(df)
        df = self._apply_rt_features(df)
        df = self._apply_peptide_statistics(df)
        df = self._apply_speceval_rank(df)
        return df

    def transform(self, df):
        """
        Apply learned statistics to test data.
        """
        if not self.is_fitted:
            warnings.warn("FeatureEngineer.transform() called before fit().", UserWarning)
        df = self.engineer_features(df)
        df = self._apply_rt_features(df)
        df = self._apply_peptide_statistics(df)
        df = self._apply_speceval_rank(df)
        return df

# =============================================================================
# Model training utilities
# =============================================================================
def get_xgb_params(n_train, scale_pos_weight, random_state):
    """
    Adaptive XGBoost parameters based on training set size.
    Tiers:  large >= 100k,  medium >= 50k,  small < 50k
    """
    base = dict(
        objective = 'binary:logistic',
        eval_metric = 'auc',
        tree_method = 'hist',
        device = 'cuda',
        scale_pos_weight = scale_pos_weight,
        random_state = random_state,
    )

    if n_train >= 100_000:
        tier = 'large'
        params = dict(n_estimators=3000, max_depth=7, learning_rate=0.01,
                      subsample=0.80, colsample_bytree=0.75, min_child_weight=8,
                      gamma=0.10, reg_alpha=0.05, reg_lambda=2.0,
                      early_stopping_rounds=150)
    elif n_train >= 50_000:
        tier = 'medium'
        params = dict(n_estimators=2000, max_depth=6, learning_rate=0.02,
                      subsample=0.80, colsample_bytree=0.70, min_child_weight=6,
                      gamma=0.15, reg_alpha=0.10, reg_lambda=2.5,
                      early_stopping_rounds=100)
    else:
        tier = 'small'
        params = dict(n_estimators=1000, max_depth=5,  learning_rate=0.05,
                      subsample=0.75, colsample_bytree=0.65, min_child_weight=15,
                      gamma=0.20, reg_alpha=0.20, reg_lambda=3.0,
                      early_stopping_rounds=50)

    base.update(params)
    return base, tier

def train_xgb_classifier(
    X_train,
    y_train,
    *,
    sample_weights=None,
    random_state=42,
    scale_pos_weight=1.0,
    params=None,
):
    n_train = len(X_train)
    model_params, _ = get_xgb_params(n_train, scale_pos_weight, random_state)
    if params:
        model_params.update(params)

    clf = _xgb_gpu_cpu(model_params)

    try:
        tr_idx, val_idx = train_test_split(
            np.arange(n_train),
            test_size=0.1,
            stratify=y_train,
            random_state=random_state,
        )
    except ValueError:
        tr_idx = np.arange(n_train)
        val_idx = tr_idx[:max(int(0.1 * n_train), 1)]

    X_tr, y_tr = X_train[tr_idx],  y_train[tr_idx]
    X_val, y_val = X_train[val_idx], y_train[val_idx]

    if np.unique(y_tr).size == 2 and np.unique(y_val).size == 2:
        fit_kwargs = {'eval_set': [(X_val, y_val)], 'verbose': False}
        if sample_weights is not None:
            fit_kwargs['sample_weight'] = sample_weights[tr_idx]
        clf.fit(X_tr, y_tr, **fit_kwargs)
    else:
        print('[WARN] Single-class split detected, disabling early stopping')
        kw = {'verbose': False}
        if sample_weights is not None:
            kw['sample_weight'] = sample_weights
        clf.fit(X_train, y_train, **kw)

    return clf

# ----------------------------------------------------
# CV splitting and per-fold metrics
# ----------------------------------------------------
def make_grouped_splits(
    df,
    *,
    n_splits=5,
    random_state=42
    ):
    """
    Yield (train_idx, test_idx) splits grouped by OptimalScanNumber and
    stratified by is_decoy. Falls back to balanced scan folds or modulo split.
    """
    try:
        if 'OptimalScanNumber' in df.columns:
            sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
            for tr, te in sgkf.split(df, df['is_decoy'], groups=df['OptimalScanNumber']):
                yield tr, te
            return
    except Exception:
        pass

    scan_to_fold = (
        make_balanced_scan_folds(df, n_folds=n_splits, random_state=random_state)
        if 'OptimalScanNumber' in df.columns else None
    )
    fold_ids = (
        df['OptimalScanNumber'].map(scan_to_fold).to_numpy()
        if scan_to_fold is not None
        else np.arange(len(df)) % n_splits
    )
    for f in range(n_splits):
        yield np.where(fold_ids != f)[0], np.where(fold_ids == f)[0]

def top1_auc_for_fold(test_feat, preds) -> float | None:
    """
    Top-1-per-scan ROC AUC: for each spectrum keep only the highest-scoring
    PSM, then compute AUC against the target/decoy label.
    """
    if 'OptimalScanNumber' not in test_feat.columns or 'is_decoy' not in test_feat.columns:
        return None
    tmp = test_feat[['OptimalScanNumber', 'is_decoy']].copy()
    tmp['pred'] = preds
    top1 = tmp.loc[tmp.groupby('OptimalScanNumber')['pred'].idxmax()]
    y = (1 - top1['is_decoy']).astype(int)
    if y.nunique() < 2:
        return None
    return float(roc_auc_score(y, top1['pred'].to_numpy()))

# --------------------------------------------------------------
# Feature selection
# --------------------------------------------------------------
def select_model_features(df):
    """
    Select numeric columns suitable for XGBoost. Drops identifiers and labels.
    """
    numeric_cols = df.select_dtypes(include=['number', 'bool']).columns.tolist()
    drop_cols = {
        # Labels / outputs
        'is_decoy', 'Label', 'xpec_q', 'xpec_pepq','_core_seq',
        # Scan / spectrum identifiers
        'OptimalScanNumber', 'ScanNum', 'SpectrumID', 'SpectrumIndex',
        'Index', 'SpecIndex', 'Scan', 'PeakScanStart', 'PeakScanEnd',
        'FileID', 'RawFileID', 'RawFileIndex', 'TitleID', 'ResultID', 'ScanType',
        # Raw scores replaced by log transforms
        'MSGFScore', 'MSGFDB_SpecEValue', 'Rank_MSGFDB_SpecEValue', 'specevalue_rank',
        # Raw RT — replaced by residual/deviation features
        'ElutionTime', 'RetentionTime', 'observed_rt', 'predicted_rt', 'raw_rt_delta',
        'peptide_hydrophobicity', 'rt_is_valid',
        # Redundant (log version kept)
        'n_psms_in_scan', 'StatMomentsArea', 'n_mods',
        'total_ion_current', 'base_peak_intensity', 'mean_matched_ion_intensity',
        # Raw ion series lengths replaced by series_length_ratio
        'longest_b_series', 'longest_y_series',
        # Other raw columns
        'PrecursorMZ', 'DelM_PPM', 'EValue', 'QValue', 'PepQValue', 'DelM',
        'ParentIonIntensity',
        # Peptide-level internals 
        'peptide_observation_count', 'peptide_best_score', 'peptide_median_score',
        'peptide_score_std', 'peptide_reliability', 'is_best_for_peptide',
        'is_singleton_scan', 'pep_is_singleton',
        # Raw MS-GF+ columns
        'DeNovoScore', 'BasePeakMZ', 'MH', 'PeakArea',
        'FWHMInScans', 'BasePeakIntensity', 'TotalIonIntensity',
        'PeakMaxIntensity', 'PeakWidthMinutes', 'PeakSignalToNoiseRatio', 'ParentIonMZ',
    }
    candidate_cols = [c for c in df.columns if c not in drop_cols]
    # A feature that never varies across the training fold cannot contribute to any split
    # Removing it before training avoids wasting colsample slots.
    candidate_cols = [
        c for c in df.columns
        if c not in drop_cols
        and pd.api.types.is_numeric_dtype(df[c])
    ]
    non_constant = [
        c for c in candidate_cols
        if df[c].nunique(dropna=False) > 1
    ]

    if len(non_constant) < len(candidate_cols):
        dropped = set(candidate_cols) - set(non_constant)
        warnings.warn(
            f"[Features] Dropping {len(dropped)} constant features: {sorted(dropped)}",
            UserWarning, stacklevel=2,
        )
    return non_constant

# =============================================================================
# Rescoring workflow
# =============================================================================
def cross_validated_rescoring(
    df,
    *,
    n_folds=5,
    n_ensemble=5,
    verbose: bool = True,
    random_state=42,
):
    FeatureEngineer._peptide_core.cache_clear()
    FeatureEngineer._achrom_predict.cache_clear()
    df_output = df.copy()

    if 'log_msgf_score' not in df_output.columns and 'MSGFScore' in df_output.columns:
        df_output['log_msgf_score'] = np.log1p(
            pd.to_numeric(df_output['MSGFScore'], errors='coerce').clip(lower=0)
        )

    if 'neg_log_specEvalue' not in df_output.columns and 'MSGFDB_SpecEValue' in df_output.columns:
        spec_eval = pd.to_numeric(df_output['MSGFDB_SpecEValue'], errors='coerce').replace(0, np.nan)
        min_nonzero = spec_eval[spec_eval > 0].min() if (spec_eval > 0).any() else 1e-100
        spec_eval = spec_eval.fillna(min_nonzero / 10)
        df_output['neg_log_specEvalue'] = -np.log10(spec_eval)
    
    num_psms = len(df_output)
    cv_predictions = np.zeros(num_psms, dtype=float)
    cv_observed_rt = np.full(num_psms, np.nan, dtype=float)
    cv_predicted_rt = np.full(num_psms, np.nan, dtype=float)
    cv_abs_rt_delta = np.full(num_psms, np.nan, dtype=float)
    cv_rt_match_prob = np.full(num_psms, np.nan, dtype=float)
    cv_missed_cleavages = np.full(num_psms, np.nan, dtype=float)
    cv_delta_score = np.full(num_psms, np.nan, dtype=float)
    cv_raw_predicted_rt  = np.full(num_psms, np.nan, dtype=float)
    #entrapment analysis de-bugging
    # cv_rescore_prob_std = np.full(num_psms, np.nan, dtype=float)
    # cv_rescore_prob_range = np.full(num_psms, np.nan, dtype=float)
    # cv_shap_values = None
    # cv_feature_matrix = None

    if verbose:
        print(f'Running {n_folds}-Fold Cross-Validation')

    fold_summaries = []
    fold_importances:  list[dict] = [] 

    for fold_number, (train_idx, test_idx) in enumerate(
        make_grouped_splits(df_output, n_splits=n_folds, random_state=random_state), start=1
    ):

        train_df = df_output.iloc[train_idx].copy()
        test_df  = df_output.iloc[test_idx].copy()

        # Leakage-free fold-level PSM features (computed from training rows only)
        train_df = compute_fold_psm_features(train_df, train_df)
        test_df = compute_fold_psm_features(train_df, test_df)

        # Feature engineering (leakage-free)
        fe = FeatureEngineer(verbose=verbose)
        train_feat = fe.fit_transform(train_df)
        test_feat  = fe.transform(test_df)

        feature_cols = select_model_features(train_feat)
        if fold_number == 1:
            frozen_feature_cols = feature_cols.copy()
            if verbose:
                print(f'[INFO] Using {len(feature_cols)} features')
        else:
            missing = set(frozen_feature_cols) - set(test_feat.columns)
            if missing and verbose:
                print(f'[WARN] Fold {fold_number}: {len(missing)} frozen features missing from test fold')
            feature_cols = frozen_feature_cols

        X_train = train_feat[feature_cols].to_numpy()
        y_train = (1 - train_feat['is_decoy']).astype(int).to_numpy()
        X_test  = test_feat[feature_cols].to_numpy()

        #feature correlation analysis de-bugging
        # if cv_feature_matrix is None and feature_cols:
        #     cv_feature_matrix = np.full((num_psms, len(feature_cols)), np.nan, dtype=np.float32)
        # if cv_feature_matrix is not None and len(feature_cols) == cv_feature_matrix.shape[1]:
        #     cv_feature_matrix[test_idx] = X_test

        y_test  = (1 - test_feat['is_decoy']).astype(int).to_numpy()

        if np.unique(y_train).size < 2 or np.unique(y_test).size < 2:
            if verbose:
                print('[WARN] Single-class fold; assigning 0.5 to test predictions')
            cv_predictions[test_idx] = 0.5
            fold_summaries.append({'fold': fold_number, 'auc': None, 'top1_auc': None})
            continue

        pos = int((y_train == 1).sum())
        neg = int((y_train == 0).sum())
        spw_data = float(neg / max(pos, 1))

        # Estimate label-noise within this fold using the raw spectral score.
        fold_signal_auc = None
        if 'neg_log_specEvalue' in train_feat.columns:
            try:
                fold_signal_auc = float(
                    roc_auc_score(y_train, train_feat['neg_log_specEvalue'].fillna(0).to_numpy())
                )
            except Exception:
                pass

        if fold_signal_auc is not None and fold_signal_auc < 0.70:
            spw = 1.0
            spw_reason = f'noisy label regime (specEvalue train AUC={fold_signal_auc:.3f} < 0.70)'
        elif abs(spw_data - 1.0) < 0.15:
            spw = 1.0
            spw_reason = f'near-balanced classes (computed spw={spw_data:.3f})'
        else:
            spw = spw_data
            spw_reason = f'computed from class ratio'

        # Use the intra-fold AUC signal to soften sample-weight penalties.
        sample_weights = _compute_sample_weights(train_feat, fold_auc=fold_signal_auc)

        # class_weight applies the spw scalar per-sample so XGBoost's internal
        # scale_pos_weight and sample_weight do not double-count the correction.
        class_weight   = np.where(y_train == 1, spw, 1.0)
        sample_weights = sample_weights * class_weight
        sample_weights /= sample_weights.mean()          # keep mean=1 for stable LR

        # # Initialize SHAP Values
        # if cv_shap_values is None and feature_cols:
        #     cv_shap_values = np.zeros((num_psms, len(feature_cols)), dtype=np.float32)

        # fold_shap_sum = np.zeros((len(test_idx), len(feature_cols)), dtype=np.float32)
        # fold_shap_count = 0  # track successful SHAP computations

        ensemble_preds = []
        for m in range(n_ensemble):
            model = train_xgb_classifier(
                X_train, y_train,
                sample_weights=sample_weights,
                random_state=random_state + m,
                scale_pos_weight=1.0,
            )
            dtest = xgb.DMatrix(X_test)
            ensemble_preds.append(model.get_booster().predict(dtest))
            if m == 0:
                fold_importances.append(
                dict(zip(feature_cols, model.feature_importances_))
            )

            # #  SHAP computation
            # if cv_shap_values is not None and len(feature_cols) == cv_shap_values.shape[1]:
            #     try:
            #         explainer = shap.TreeExplainer(model)
            #         fold_shap = explainer.shap_values(X_test)
            #         if isinstance(fold_shap, list):
            #             # For binary classification, shap returns [class0, class1]
            #             fold_shap = fold_shap[1] if len(fold_shap) > 1 else fold_shap[0]
            #         fold_shap = np.asarray(fold_shap, dtype=np.float32)
            #         fold_shap_sum += fold_shap
            #         fold_shap_count += 1
            #     except Exception as e:
            #         if verbose:
            #             print(f'[WARN] SHAP calculation failed on fold {fold_number}, member {m}: {e}')

            if fold_number == 1 and m == 0 and verbose:
                _, tier = get_xgb_params(len(X_train), spw, random_state)
                print(f'[INFO] Best iteration: {model.best_iteration}')
                print(f'[INFO] Best score: {model.best_score:.3f}')
                print(f'[INFO] n_train: {X_train.shape[0]:,} n_features: {X_train.shape[1]}')
                print(f'[INFO] Dataset tier: {tier}')

        preds = np.vstack(ensemble_preds)
        if preds.shape[0] >= 7:
            preds = np.sort(preds, axis=0)[1:-1]   # trim min/max for robustness
        avg_pred = np.clip(preds.mean(axis=0), 1e-4, 1 - 1e-4)

        # # Average SHAP Values
        # if cv_shap_values is not None and fold_shap_count > 0:
        #     cv_shap_values[test_idx] = fold_shap_sum / fold_shap_count

        cv_predictions[test_idx] = avg_pred
        # cv_rescore_prob_std[test_idx]   = preds.std(axis=0)
        # cv_rescore_prob_range[test_idx] = preds.max(axis=0) - preds.min(axis=0)

        # Collect RT columns for post-hoc plotting
        if 'observed_rt' in test_feat.columns:
            cv_observed_rt[test_idx] = test_feat['observed_rt'].to_numpy(dtype=float)
        if 'predicted_rt' in test_feat.columns:
            cv_predicted_rt[test_idx] = test_feat['predicted_rt'].to_numpy(dtype=float)
        if 'abs_rt_delta' in test_feat.columns:
            cv_abs_rt_delta[test_idx] = test_feat['abs_rt_delta'].to_numpy(dtype=float)
        if 'rt_match_prob' in test_feat.columns:
            cv_rt_match_prob[test_idx] = test_feat['rt_match_prob'].to_numpy(dtype=float)
        if 'missed_cleavages' in test_feat.columns:
            cv_missed_cleavages[test_idx] = test_feat['missed_cleavages'].to_numpy(dtype=int)
        if 'delta_score' in test_feat.columns:
            cv_delta_score[test_idx] = test_feat['delta_score'].to_numpy(dtype=float)
        if (fe.rt_calibrator is not None
                and fe.rt_calibrator._raw_achrom_preds is not None
                and '_core_seq' in test_feat.columns):
            raw_preds = test_feat['_core_seq'].map(fe.rt_calibrator._raw_achrom_preds)
            missing = raw_preds.isna() & test_feat['_core_seq'].notna()
            if missing.any():
                raw_preds.loc[missing] = test_feat.loc[missing, '_core_seq'].apply(FeatureEngineer._achrom_predict)
            cv_raw_predicted_rt[test_idx] = raw_preds.to_numpy(dtype=float)
         
        
        # Per-fold diagnostics
        try:
            fold_auc = float(roc_auc_score(y_test, avg_pred)) 
            fold_top1auc = top1_auc_for_fold(test_feat, avg_pred)
            fold_summaries.append({
                'fold': fold_number,
                'auc': fold_auc,
                'top1_auc': fold_top1auc,
            })
        except Exception:
            fold_summaries.append({'fold': fold_number, 'auc': None, 'top1_auc': None})

    if verbose:
        print('Cross-Validation Complete')

    df_output['rescore_prob'] = cv_predictions
    # df_output['rescore_prob_std'] = cv_rescore_prob_std
    # df_output['rescore_prob_range'] = cv_rescore_prob_range
    df_output['observed_rt'] = cv_observed_rt
    df_output['predicted_rt'] = cv_predicted_rt
    df_output['abs_rt_delta'] = cv_abs_rt_delta
    df_output['rt_match_prob'] = cv_rt_match_prob
    df_output['missed_cleavages'] = cv_missed_cleavages
    df_output['delta_score'] = cv_delta_score
    df_output['raw_predicted_rt'] = cv_raw_predicted_rt
    
    # # Entrapment analysis debugging with SHAP
    # if cv_shap_values is not None and frozen_feature_cols is not None:
    #     # Direct assignment — no pd.concat, no attrs, no duplicates
    #     for i, col_name in enumerate(frozen_feature_cols):
    #         df_output[f'shap_{col_name}'] = cv_shap_values[:, i]

    #     if cv_feature_matrix is not None:
    #         shap_feature_matrix = pd.DataFrame(
    #             cv_feature_matrix, columns=frozen_feature_cols, index=df_output.index
    #         )
    #         shap_feature_cols = list(frozen_feature_cols)
    #     else:
    #         shap_feature_matrix = None
    #         shap_feature_cols = None

    #     if verbose:
    #         n_shap = sum(1 for c in df_output.columns if c.startswith('shap_'))
    #         print(f'[INFO] SHAP: {n_shap} columns attached to {len(df_output)} PSMs.')
    # else:
    #     shap_feature_matrix = None
    #     shap_feature_cols = None

    df_output.attrs = {}

    best_alpha, alpha_table = tune_noisy_or_weights(df_output.copy(), verbose=True)
    df_output['xpec_score'] = compute_noisy_or_score(
        df_output,
        ml_col='rescore_prob',
        evalue_col='neg_log_specEvalue',
        alpha=best_alpha,
        beta=1.0,
    )

    df_output = estimate_fdr(df_output, score_col='xpec_score', decoy_col='is_decoy', massdiff_col='DelM', ptm_col='n_mods')
    df_output = estimate_pep_fdr(df_output, score_col='xpec_score')

    df_output = df_output.reset_index(drop=True)
    feature_importance_summary = _aggregate_feature_importance(fold_importances)

    return df_output, feature_importance_summary