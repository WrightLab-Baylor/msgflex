#!/usr/bin/env python3
"""
xrtpred.py --- Two-stage retention-time calibration for PSM rescoring (C18 RP-HPLC).

Pipeline
--------
Stage 1 -- Isotonic warp + LOWESS
Stage 2 -- XGBoost residual regressor

Final prediction
----------------
    rt_pred  = LOWESS(warp(achrom_rt)) + bias_correction + XGB_residual
    rt_score = exp(−|observed_rt − rt_pred| / MAD)   [Laplace likelihood]
"""
import re
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import HuberRegressor
from statsmodels.nonparametric.smoothers_lowess import lowess
from xgboost import XGBRegressor
from functools import lru_cache

class RTCalibrator:
    """
    Two-stage RT calibration for PSM rescoring (C18 RP-HPLC).
    """
    
    # RT offset per modification (minutes, C18 gradient; sign = elution shift)
    _MOD_RT_DELTA_DEFAULTS = {
        "Oxidation": -3.5,   # Met+16:  hydrophilic shift
        "Phospho": -5.0,   # pS/pT/pY+80: strong polar shift
        "Deamidation": -0.5,   # Asn/Gln+1: minor shift
        "Acetyl": +3.0,   # N-term acetylation: slight hydrophobic shift
        "AcNoTMT": -2.5,
    }

    # Plausible physical ranges for empirically fitted mod offsets (minutes)
    _MOD_RT_PLAUSIBLE = {
        "Oxidation": (-6.0, -1.0),
        "Phospho": (-8.0, -1.0),
        "Deamidation": (-2.0,  0.0),
        "Acetyl":  ( 1.0,  6.0),
        "AcNoTMT": (-4.0, -0.5),
    }

    # Monoisotopic masses for total_mod_mass computation
    _MOD_MASSES = {
        "Oxidation": 15.9949,
        "Acetyl": 42.0106,
        "Deamidation": 0.9840,
        "Phospho": 79.9663,
        "AcNoTMT": -187.1524,   # net mass delta relative to TMTPro-K baseline
    }

    # Kyte–Doolittle hydrophobicity scale
    _HYDROPHOBICITY = {
        "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
        "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
        "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
        "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
    }

    _HYDRO_LOOKUP = np.full(128, np.nan)
    for _aa, _val in {
        "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
        "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
        "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
        "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
    }.items():
        _HYDRO_LOOKUP[ord(_aa)] = _val
    del _aa, _val   # keep class namespace clean

    _ALL_AA: str = "ACDEFGHIKLMNPQRSTVWY"

    # Regex patterns for PHRP modification notation
    _MOD_PATTERNS: dict[str, re.Pattern] = {
        "Oxidation": re.compile(r"M\+15\.9"),
        "Phospho": re.compile(r"[STY]\+79\."),
        "Deamidation": re.compile(r"[NQ]\+0\.9"),
        "Acetyl": re.compile(r"(\+42\.0|n\+42\.)"),
        "AcNoTMT": re.compile(r"(K-187\.1|K\-187\.15)"),
    }

    def __init__(
        self,
        frac = 0.25,
        min_points = 30,
        verbose = False,
    ) -> None:
        self.frac = frac
        self.min_points = min_points
        self.verbose = verbose

        # Per-instance copy: empirical fitting must not mutate the class defaults
        self.mod_rt_delta = dict(self._MOD_RT_DELTA_DEFAULTS)

        # Model state (populated by fit())
        self.lowess_model = None
        self.xgb_model: XGBRegressor = None
        self._bias_correction = 0.0
        self._residual_mad = 2.0
        self._achrom_warp: IsotonicRegression | None = None
        self._warp_x: np.ndarray | None = None   # sorted Achrom scores, training set
        self._warp_y: np.ndarray | None = None   # corresponding smoothed RT values
        self._raw_achrom_preds: np.ndarray | None = None  # pre-calibration RT map

    # ------------------------------------------------------------------
    # Public API configuration 
    # ------------------------------------------------------------------

    def fit(self, df) -> None:
        """
        Fit the two-stage RT calibrator.
        """
        df = df.copy()

        # 1. Select high-confidence target PSMs 
        mask = (
            df["_achrom_rt"].notna()
            & df["observed_rt"].notna()
        )
        if "is_decoy" in df.columns:
            mask &= df["is_decoy"] == 0

        # Exclude void-volume peptides (bottom 2 % of observed RT)
        rt_vals = df.loc[mask, "observed_rt"].astype(float)
        positive_rts = rt_vals[rt_vals > 0]
        
        if len(positive_rts) >= 10:
            rt_min = float(np.nanpercentile(positive_rts, 2))
        else:
            rt_min = 0.5
        mask &= df["observed_rt"].astype(float) >= rt_min

        df = df[mask].copy()
        if len(df) < self.min_points:
            raise ValueError(
                f"RTCalibrator.fit: only {len(df)} valid PSMs after filtering; "
                f"need ≥ {self.min_points}."
            )

        # 2. Empirically calibrate per-modification RT offsets 
        self._fit_mod_rt_deltas(df)

        # 3. Mod-corrected observed RT
        mod_corr = df["Peptide"].apply(self._mod_rt_correction).astype(float)
        x_raw = df["_achrom_rt"].astype(float).to_numpy()
        y_raw = (df["observed_rt"].astype(float) - mod_corr).to_numpy()

        if "_core_seq" in df.columns:
            self._raw_achrom_preds = dict(zip(
                df["_core_seq"].to_numpy(),
                x_raw,
            ))
            
        finite = np.isfinite(x_raw) & np.isfinite(y_raw)
        x_raw, y_raw = x_raw[finite], y_raw[finite]

        # 4. Isotonic warp: Achrom score --- observed RT space 
        # Preserves elution ordering while correcting gradient nonlinearity.
        self._fit_achrom_warp(x_raw, y_raw)
        x_warped = self._apply_achrom_warp(x_raw)

        # 5. First-pass LOWESS + IQR outlier removal 
        xs0, ys0 = self._fit_lowess(x_warped, y_raw)
        resid0 = y_raw - np.interp(x_warped, xs0, ys0)
        q1, q3 = np.percentile(resid0, [25, 75])
        iqr = q3 - q1
        keep = (resid0 >= q1 - 1.5 * iqr) & (resid0 <= q3 + 1.5 * iqr)
        x_f, y_f = x_warped[keep], y_raw[keep]

        if len(x_f) < self.min_points:
            raise ValueError("RTCalibrator.fit: too few points after IQR filtering.")

        # 6. Final LOWESS on clean data
        xs, ys = self._fit_lowess(x_f, y_f)
        self.lowess_model = {"x": xs, "y": ys}

        y_lowess_train = np.interp(x_warped, xs, ys)
        self._bias_correction = float(np.median(y_raw - y_lowess_train))
        if self.verbose:
            print(f"[RT] LOWESS bias correction: {self._bias_correction:+.3f} min")

        # 7. XGBoost residual regressor 
        df_xgb = df[finite].copy().reset_index(drop=True)
        cols_to_carry = ['neg_log_specEvalue']
        for col in cols_to_carry:
            if col in df.columns:
                df_xgb[col] = df[col].iloc[
                    np.where(finite)[0]
                ].to_numpy()
    
        x_xgb_warped = self._apply_achrom_warp(
            df_xgb["_achrom_rt"].astype(float).to_numpy()
        )
        df_xgb["_achrom_rt_warped"] = x_xgb_warped

        y_lowess_xgb = np.interp(x_xgb_warped, xs, ys) + self._bias_correction
        mod_corr2 = df_xgb["Peptide"].apply(self._mod_rt_correction).astype(float)
        residual = (df_xgb["observed_rt"].astype(float) - mod_corr2) - y_lowess_xgb
        
        self._residual_mad = float(
            np.median(np.abs(residual - np.median(residual)))
        )
        
        _LEN_BINS = [0, 7, 11, 16, 100]
        _LEN_LABELS = ['s', 'm', 'l', 'xl']

        df_xgb['_residual']  = residual.values
        df_xgb['_len_bin']   = pd.cut(
            df_xgb['_core_seq'].str.len().fillna(8).clip(lower=1),
            bins=_LEN_BINS, labels=_LEN_LABELS,
        ).astype(str)

        self._class_mad = {}
        for (chg, lb), grp in df_xgb.groupby(['Charge', '_len_bin'], observed=True):
            if len(grp) < 10:
                continue
            r = grp['_residual'].to_numpy()
            class_mad = float(np.median(np.abs(r - np.median(r))))
            self._class_mad[f'{int(chg)}_{lb}'] = max(class_mad, 0.3) # 0.3 min floor

        if self.verbose:
            print(f'[RT] Per-class MAD computed for {len(self._class_mad)} charge×length bins')
        
        X = self._build_features(df_xgb)
        n_train = len(X)
        xgb_params, xgb_tier = self._get_xgb_regressor_params(n_train)

        weights = self._compute_sample_weights(df_xgb, residual.to_numpy())
        early_stopping_rounds = xgb_params.pop("early_stopping_rounds")
        val_size = min(max(int(n_train * 0.10), 500), 5000)
        val_idx = np.random.default_rng(42).choice(n_train, size=val_size, replace=False)
        tr_mask = np.ones(n_train, dtype=bool)
        tr_mask[val_idx] = False

        X_arr = X.to_numpy()
        y_arr = residual.to_numpy()

        self.xgb_model = XGBRegressor(**xgb_params, early_stopping_rounds=early_stopping_rounds)
        self.xgb_model.fit(
            X_arr[tr_mask],  y_arr[tr_mask],
            sample_weight=weights[tr_mask],
            eval_set=[(X_arr[val_idx], y_arr[val_idx])],
            verbose=False
        )

    def predict(self, df):
        """
        Generate calibrated RT predictions for a set of PSMs.
        """
        df = df.copy()

        x_warped = self._apply_achrom_warp(df["_achrom_rt"].astype(float).to_numpy())
        df["_achrom_rt_warped"] = x_warped

        lowess_pred = self._predict_lowess(x_warped)
        xgb_delta = self.xgb_model.predict(self._build_features(df))

        rt_pred = lowess_pred + xgb_delta
       
        if "observed_rt" in df.columns:
            rt_error = (df["observed_rt"].astype(float) - rt_pred).abs()
            rt_error = np.minimum(rt_error, 3 * self._residual_mad)
        else:
            rt_error = pd.Series(np.nan, index=df.index)

        if getattr(self, '_class_mad', None):
            len_bins = pd.cut(
                df['_core_seq'].str.len().fillna(8).clip(lower=1),
                bins=[0, 7, 11, 16, 100], labels=['s', 'm', 'l', 'xl'],
            ).astype(str)
            keys = df['Charge'].astype(str) + '_' + len_bins
            per_row_mad = (
                keys.map(self._class_mad)
                    .fillna(self._residual_mad)
                    .clip(lower=0.3)
                    .to_numpy(dtype=float)
            )
        else:
            per_row_mad = np.full(len(df), max(self._residual_mad, 0.3))

        rt_score = np.exp(-rt_error / per_row_mad)

        return {
            "rt_pred": pd.Series(rt_pred, index=df.index),
            "rt_error": rt_error,
            "rt_score": rt_score,
        }

    # ------------------------------------------------------------------
    # Modification RT correction
    # ------------------------------------------------------------------
    @lru_cache(maxsize=500_000)
    def _parse_mod_counts(self, peptide: str) -> dict:
        s = peptide if isinstance(peptide, str) else ""
        return {
            "Oxidation": len(re.findall(self._MOD_PATTERNS["Oxidation"], s)),
            "Phospho": len(re.findall(self._MOD_PATTERNS["Phospho"], s)),
            "Deamidation": len(re.findall(self._MOD_PATTERNS["Deamidation"], s)),
            "Acetyl": len(re.findall(self._MOD_PATTERNS["Acetyl"], s)),
            "AcNoTMT": len(re.findall(self._MOD_PATTERNS["AcNoTMT"], s))
        }

    def _mod_rt_correction(self, peptide: str) -> float:
        counts = self._parse_mod_counts(peptide)
        return (
            counts["Oxidation"] * self.mod_rt_delta["Oxidation"]
            + counts["Phospho"] * self.mod_rt_delta["Phospho"]
            + counts["Deamidation"] * self.mod_rt_delta["Deamidation"]
            + counts["Acetyl"] * self.mod_rt_delta["Acetyl"]
            + counts["AcNoTMT"] * self.mod_rt_delta["AcNoTMT"]
        )

    def _fit_mod_rt_deltas(self, df: pd.DataFrame) -> None:
        """
        Vectorized claculations
        """
        mod_types = list(self._MOD_PATTERNS.keys())
        pep_str = df["Peptide"].astype(str)

        # Step 1 — vectorized modification counts (one pass per modification)
        counts_df = pd.DataFrame(
            {mod: pep_str.str.count(pat.pattern) for mod, pat in self._MOD_PATTERNS.items()},
            index=df.index,
        )
        total_mods = counts_df.sum(axis=1)
        is_unmod = total_mods == 0

        # Step 2 — per-core-seq unmodified RT median (groupby, no Python loop)
        unmod_rt = (
            df.loc[is_unmod, ["_core_seq", "observed_rt"]]
            .groupby("_core_seq")["observed_rt"]
            .median()
            .rename("_rt_unmod")
        )

        # Step 3 — join unmod reference onto modified rows only
        mod_mask = ~is_unmod
        if mod_mask.sum() == 0 or len(unmod_rt) == 0:
            return

        pairs_df = counts_df.loc[mod_mask].copy()
        pairs_df["_rt_unmod"] = df.loc[mod_mask, "_core_seq"].map(unmod_rt)
        pairs_df = pairs_df.dropna(subset=["_rt_unmod"])   # no unmod reference → skip

        if len(pairs_df) < 10:
            return

        pairs_df["rt_delta"] = (
            df.loc[pairs_df.index, "observed_rt"].astype(float).values
            - pairs_df["_rt_unmod"].values
        )

        # Step 4 — fit HuberRegressor on the pairs matrix
        X = pairs_df[mod_types].to_numpy()
        y = pairs_df["rt_delta"].to_numpy()

        model = HuberRegressor(epsilon=1.5, fit_intercept=False)
        model.fit(X, y)
        fitted = dict(zip(mod_types, model.coef_))

        # Step 5 — validate against plausible physical ranges
        for mod, val in fitted.items():
            if mod not in self._MOD_RT_PLAUSIBLE:
                continue
            lo, hi = self._MOD_RT_PLAUSIBLE[mod]
            if lo <= val <= hi:
                self.mod_rt_delta[mod] = val
            elif self.verbose:
                print(
                    f"[RT] Fitted {mod}={val:.2f} outside plausible range "
                    f"{(lo, hi)}; keeping hardcoded default."
                )

    # ------------------------------------------------------------------
    # Achrom isotonic warp
    # ------------------------------------------------------------------
    def _fit_achrom_warp(self, x: np.ndarray, y: np.ndarray) -> None:
        """
        Fit a monotonic warp from Achrom score → observed RT space via
        isotonic regression on a LOWESS-smoothed target.
        """
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0)
        order = np.argsort(x[valid])
        x_sorted = x[valid][order]
        y_sorted = y[valid][order]

        # Smooth target to reduce staircase artefacts in isotonic regression
        frac_warp = float(np.clip(30 / len(x_sorted), 0.05, 0.15))
        y_smooth = lowess(y_sorted, x_sorted, frac=frac_warp, return_sorted=False)

        self._achrom_warp = IsotonicRegression(
            increasing=True,
            out_of_bounds="nan",   # out-of-range points flagged; handled below
        )
        self._achrom_warp.fit(x_sorted, y_smooth)

        # Stored for linear edge extrapolation
        self._warp_x = x_sorted
        self._warp_y = y_smooth

    def _apply_achrom_warp(self, x: np.ndarray) -> np.ndarray:
        """
        Apply the fitted isotonic warp to a vector of Achrom scores.

        Points outside the training range are extrapolated linearly using the
        slope estimated from the outermost 10 % of training data, preventing
        the flat-clip artefact of the default isotonic behaviour.
        """
        x_arr = np.asarray(x, dtype=float)
        if self._achrom_warp is None:
            return x_arr

        warped = self._achrom_warp.predict(x_arr)
        nan_mask = np.isnan(warped)

        if nan_mask.any() and self._warp_x is not None and self._warp_y is not None:
            n = len(self._warp_x)
            tail = max(int(n * 0.1), 2)

            low_mask  = nan_mask & (x_arr < self._warp_x[0])
            high_mask = nan_mask & (x_arr > self._warp_x[-1])

            if low_mask.any():
                slope = (
                    (self._warp_y[tail] - self._warp_y[0])
                    / (self._warp_x[tail] - self._warp_x[0] + 1e-9)
                )
                warped[low_mask] = (
                    self._warp_y[0] + slope * (x_arr[low_mask] - self._warp_x[0])
                )

            if high_mask.any():
                slope = (
                    (self._warp_y[-1] - self._warp_y[-tail])
                    / (self._warp_x[-1] - self._warp_x[-tail] + 1e-9)
                )
                warped[high_mask] = (
                    self._warp_y[-1] + slope * (x_arr[high_mask] - self._warp_x[-1])
                )

        return warped

    # ------------------------------------------------------------------
    # LOWESS helpers
    # ------------------------------------------------------------------
    def _fit_lowess(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Fit LOWESS and return sorted (xs, ys).
        """
        n = len(x)
        frac = float(np.clip(30 / n, 0.08, self.frac))
        result = lowess(y, x, frac=frac, return_sorted=True)
        return result[:, 0], result[:, 1]

    def _predict_lowess(self, x: np.ndarray) -> np.ndarray:
        """
        Evaluate the fitted LOWESS curve at ``x`` and apply the bias correction.

        The bias correction (``self._bias_correction``) is the median residual
        (observed − LOWESS) measured on the full training set after LOWESS was
        fit on IQR-cleaned data.  Adding it here centres final RT predictions
        and ensures the XGB residual model trains on zero-centred targets.
        """
        xs, ys = self.lowess_model["x"], self.lowess_model["y"]
        return np.interp(np.asarray(x, dtype=float), xs, ys) + self._bias_correction

    # ---------------------------------------------------------------
    # Modification feature extraction
    # ---------------------------------------------------------------
    def _mod_features(self, peptide: str) -> dict:
        """
        Parse modification counts and cumulative modification mass from a
        PHRP peptide string.
        """
        s = str(peptide) if isinstance(peptide, str) else ""
        seq_len = max(len(re.sub(r'[^A-Z]', '', s)), 1)

        counts = {
            mod: len(pattern.findall(s))
            for mod, pattern in self._MOD_PATTERNS.items()
        }

        total_mod_mass = sum(
            counts[mod] * self._MOD_MASSES[mod]
            for mod in counts
        )

        return {
            # Aggregate
            'n_mods': sum(counts.values()),
            'total_mod_mass': total_mod_mass,
            'mod_mass_per_res':  total_mod_mass / seq_len,

            # Per-modification counts
            'n_oxidation': counts['Oxidation'],
            'n_acetyl': counts['Acetyl'],
            'n_deamidation': counts['Deamidation'],
            'n_phospho': counts['Phospho'],
            'n_acnotmt': counts['AcNoTMT']
        }

    # ---------------------------------------------------------------
    # Feature engineering
    # ---------------------------------------------------------------

    def _build_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Build the XGBoost feature matrix.

        The (warped) Achrom RT is included as a gradient-position anchor so
        the residual model can capture position-dependent systematic errors.

        Features
        --------
        - achrom_rt : warped Achrom score (gradient anchor)
        - length, length_sq : peptide length and its square
        - charge, charge_len : charge state and charge × length interaction
        - aa_<X>_frac : amino-acid composition as *fractions* (count / length).
          Using fractions instead of raw counts decouples composition from
          length (already captured by ``length`` / ``length_sq``) and removes
          the near-linear dependence between raw AA counts and ``hydro_mean``.
        - nterm_<X>, cterm_<X> : N/C-terminal identity (40 one-hot features)
        - hydro_mean, hydro_std : global Kyte–Doolittle mean and std
        - hydro_nterm3, hydro_cterm3 : mean hydrophobicity of the N/C-terminal
          tripeptide windows (RP-HPLC elution relevant)
        - hydro_nhalf, hydro_chalf : N-half / C-half mean hydrophobicity.
          These four positional features are NOT recoverable from global
          composition, so they are genuinely orthogonal to aa_frac.
          ``hydro_len`` (= hydro_mean × length) is intentionally omitted:
          it is a deterministic product of two features already present and
          adds no orthogonal information.
        - n_mods, mod_mass, mod_mass_per_res, n_oxidation, n_phospho,
          n_deamidation : modification summary
        """
        seq = df["_core_seq"]
        feats = pd.DataFrame(index=df.index)

        # Gradient-position anchor 
        feats["achrom_rt"] = df["_achrom_rt_warped"].astype(float)

        # Sequence length 
        length = seq.str.len().astype(float).clip(lower=1.0)   # avoid /0
        feats["length"] = length
        feats["length_sq"] = length ** 2

        #  Charge 
        feats["charge"] = df["Charge"]
        feats["charge_len"] = feats["charge"] * length

        # Amino-acid composition as fractions 
        seqs_np = seq.to_numpy()
        lens_np = length.to_numpy().astype(int)
        n = len(seqs_np)
        max_len = int(lens_np.max()) if n > 0 else 1

        # (N, max_len) ordinal matrix  — same pattern as _compute_hydro_stats
        ordinal = np.zeros((n, max_len), dtype=np.uint8)
        for i, s in enumerate(seqs_np):
            arr = np.frombuffer(s.encode('ascii'), dtype=np.uint8)
            ordinal[i, :len(arr)] = arr

        aa_codes = np.frombuffer((''.join(self._ALL_AA)).encode('ascii'), dtype=np.uint8)
        counts_mat = (ordinal[:, :, None] == aa_codes[None, None, :]).sum(axis=1)  # (N, 20)
        fracs = counts_mat / lens_np[:, None]

        for j, aa in enumerate(self._ALL_AA):
            feats[f"aa_{aa}_frac"] = fracs[:, j]
        
        #  N/C-terminal identity (one-hot) 
        n_term = seq.str[0]
        c_term = seq.str[-1]
        for aa in self._ALL_AA:
            feats[f"nterm_{aa}"] = (n_term == aa).astype(int)
            feats[f"cterm_{aa}"] = (c_term == aa).astype(int)

        #  Kyte–Doolittle hydrophobicity statistics 
        hydro_feats = self._compute_hydro_stats(seq)
        for col in hydro_feats.columns:
            feats[col] = hydro_feats[col].to_numpy()

        #  Modification features 
        mod_parsed = pd.DataFrame(
            df["Peptide"].apply(self._mod_features).tolist(),
            index=df.index,
        )
        feats["n_mods"] = mod_parsed["n_mods"]
        feats["mod_mass"] = mod_parsed["total_mod_mass"]
        feats["mod_mass_per_res"] = mod_parsed["mod_mass_per_res"]
        feats["n_oxidation"] = mod_parsed["n_oxidation"]
        # If conditions for a few modifications
        if "n_deamidation" in mod_parsed.columns:
            feats["n_deamidation"] = mod_parsed["n_deamidation"]
        if "n_phospho" in mod_parsed.columns:
            feats["n_phospho"] = mod_parsed["n_phospho"]

        if "n_acnotmt" in mod_parsed.columns:
            feats["n_acnotmt"] = mod_parsed["n_acnotmt"]

        return feats.fillna(0.0)

    def _compute_hydro_stats(self, seq: pd.Series):
        """
        Vectorised Kyte–Doolittle hydrophobicity statistics.
        """
        seqs = seq.to_numpy() 
        n = len(seqs)
        lengths = np.fromiter(
            (len(s) for s in seqs), dtype=np.int32, count=n
        )
        max_len = int(lengths.max()) if n > 0 else 1

        #  Step 1: build (N, max_len) ordinal matrix 
        ord_mat = np.zeros((n, max_len), dtype=np.uint8)
        for i, s in enumerate(seqs):
            if s:
                ord_mat[i, :len(s)] = np.frombuffer(
                    s.encode("ascii"), dtype=np.uint8
                )

        # Step 2: lookup  (N, max_len) hydrophobicity matrix
        hydro_mat = self._HYDRO_LOOKUP[ord_mat]   # NaN for unknown AA + padding

        # Mask padding positions explicitly (belt-and-suspenders for ord 0 / 127)
        col_idx = np.arange(max_len, dtype=np.int32)[None, :]  # (1, max_len)
        valid = col_idx < lengths[:, None]                    # (N, max_len)
        hydro_mat = np.where(valid, hydro_mat, np.nan)

        # Step 3: global statistics 

        with np.errstate(all="ignore"):
            hydro_mean = np.nanmean(hydro_mat, axis=1)
            hydro_std  = np.nanstd( hydro_mat, axis=1)

            # N-terminal tripeptide (fixed window — slice is always valid)
            hydro_nterm3 = np.nanmean(hydro_mat[:, :3], axis=1)

         
            rev_col = np.where(
                valid,
                lengths[:, None] - 1 - col_idx,
                0,                    
            )
            rev_mat = hydro_mat[np.arange(n)[:, None], rev_col]
            rev_mat = np.where(valid, rev_mat, np.nan)   # re-mask padding
            hydro_cterm3 = np.nanmean(rev_mat[:, :3], axis=1)

            # N-half / C-half — split at the integer midpoint of each sequence
            mids = np.maximum(lengths // 2, 1)            # (N,)
            nhalf_ok = valid & (col_idx < mids[:, None])
            chalf_ok = valid & (col_idx >= mids[:, None])
            hydro_nhalf = np.nanmean(np.where(nhalf_ok, hydro_mat, np.nan), axis=1)
            hydro_chalf = np.nanmean(np.where(chalf_ok, hydro_mat, np.nan), axis=1)

        hydro_nhalf = np.where(np.isnan(hydro_nhalf), hydro_mean, hydro_nhalf)
        hydro_chalf = np.where(np.isnan(hydro_chalf), hydro_mean, hydro_chalf)

        return pd.DataFrame(
            {
                "hydro_mean": hydro_mean,
                "hydro_std": hydro_std,
                "hydro_nterm3": hydro_nterm3,
                "hydro_cterm3": hydro_cterm3,
                "hydro_nhalf": hydro_nhalf,
                "hydro_chalf": hydro_chalf,
            },
            index=seq.index,
        )
    
    @staticmethod
    def _get_xgb_regressor_params(n_train: int) -> dict:
        """
        Adaptive XGBoost regressor parameters scaled to training-set size.
        Tiers mirror the classifier tiers in best_gain_model.py so that
        RT residual modelling capacity scales with data availability.
        Tiers
        -----
        large  : n >= 50,000   (typical full-run metaproteomics / phospho)
        medium : n >= 20,000   (standard single-run proteomics)
        small  : n <  20,000   (small datasets, fractions, or after MC filter)
        """
        base = dict(
            objective="reg:absoluteerror",
            eval_metric="mae",
            tree_method="hist",
            n_jobs=-1,
            random_state=42,
        )

        if n_train >= 75_000:
            tier = "large"
            params = dict(
                n_estimators=800,
                max_depth=6,
                learning_rate=0.03,
                subsample=0.80,
                colsample_bytree=0.80,
                min_child_weight=10,
                gamma=0.10,
                reg_alpha=0.05,
                reg_lambda=1.0,
                early_stopping_rounds=100,
            )
        elif n_train >= 25_000:
            tier = "medium"
            params = dict(
                n_estimators=600,
                max_depth=5,
                learning_rate=0.04,
                subsample=0.80,
                colsample_bytree=0.80,
                min_child_weight=7,
                gamma=0.10,
                reg_alpha=0.10,
                reg_lambda=1.0,
                early_stopping_rounds=80,
            )    
        else:
            tier = "small"
            params = dict(
                n_estimators=400,
                max_depth=4,
                learning_rate=0.05,
                subsample=0.75,
                colsample_bytree=0.75,
                min_child_weight=10,
                gamma=0.20,
                reg_alpha=0.20,
                reg_lambda=2.0,
                early_stopping_rounds=50,
            )

        base.update(params)
        return base, tier
    
    @staticmethod
    def _compute_sample_weights(
        df,
        residual: np.ndarray,
    ) -> np.ndarray:
        """
        Compute per-PSM sample weights for the XGBoost residual regressor.
        """
        n = len(residual)
        weights = np.ones(n, dtype=np.float32)

        # Component 1 — confidence from spectral score
        if "neg_log_specEvalue" in df.columns:
            spec = pd.to_numeric(
                df["neg_log_specEvalue"], errors="coerce"
            ).fillna(0.0).to_numpy(dtype=np.float32)
            # Rank-based normalisation: avoids sensitivity to score scale
            rank = np.argsort(np.argsort(spec)).astype(np.float32)
            conf_weight = 0.5 + rank / (2.0 * max(n - 1, 1))   # range [0.5, 1.0]
            weights *= conf_weight

        # Component 2 — down-weight extreme residuals
        abs_res = np.abs(residual).astype(np.float32)
        cap = float(np.percentile(abs_res, 95))
        if cap > 0:
            outlier_weight = np.where(
                abs_res > cap,
                cap / abs_res.clip(min=1e-6),   # soft down-weight, never zero
                1.0,
            ).astype(np.float32)
            weights *= outlier_weight

        # Normalise so mean weight = 1.0 (preserves effective sample size)
        mean_w = weights.mean()
        if mean_w > 0:
            weights /= mean_w

        return weights
    
    def save(self, path):
        """
        Persist the fitted calibrator to disk.
        Captures all model state: lowess curve, isotonic warp, XGBoost
        residual model, bias correction, MAD, and modification offsets.
        """
        import joblib, pathlib
        pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, path, compress=3)

    @classmethod
    def load(cls, path: str) -> "RTCalibrator":
        """
        Load a calibrator previously saved with save().
        """
        import joblib
        obj = joblib.load(path)
        if not isinstance(obj, cls):
            raise TypeError(
                f"File contains {type(obj).__name__}, expected RTCalibrator."
            )
        return obj

    def is_fitted(self) -> bool:
        """
        Return True if the calibrator has been fitted.
        """
        return self.lowess_model is not None and self.xgb_model is not None
    