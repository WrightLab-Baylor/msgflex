#!/usr/bin/env python3
"""
Extracts modification-aware, charge-aware spectral features from mzML files
for PSM rescoring downstream.
"""
from __future__ import annotations

from multiprocessing import pool
import re
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm
import multiprocessing as mp
warnings.filterwarnings("ignore")
try:
    from pyteomics import mzml
except ImportError:
    raise ImportError(
        "pyteomics is required. Install with: pip install pyteomics"
    )

# --- Multiprocessing helper functions ---
def _worker_init(spectra_dict, extractor_obj):
    global _WORKER_SPECTRA, _WORKER_EXTRACTOR
    _WORKER_SPECTRA   = spectra_dict
    _WORKER_EXTRACTOR = extractor_obj

def _worker_fn(args):
    scan_num, peptide, charge = args
    if scan_num is None or int(scan_num) not in _WORKER_SPECTRA:
        return {"ScanNum": scan_num, "spectral_feature_error": "scan_not_found"}
    return _WORKER_EXTRACTOR._psm_features(
        spectrum = _WORKER_SPECTRA[int(scan_num)],
        peptide  = peptide,
        charge   = int(charge),
        scan_num = int(scan_num),
    )

# Constants

# Monoisotopic amino acid residue masses (Da)
AA_MASSES = {
    "G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764,
    "V": 99.068414, "T": 101.047670, "C": 103.009185, "L": 113.084064,
    "I": 113.084064, "N": 114.042927, "D": 115.026943, "Q": 128.058578,
    "K": 128.094963, "E": 129.042593, "M": 131.040485, "H": 137.058912,
    "F": 147.068414, "R": 156.101111, "Y": 163.063329, "W": 186.079313,
}

H_MASS = 1.007276   # proton mass
H2O = 18.010565  # water

# MS-GF+ shorthand modification symbols → monoisotopic mass shift
_MSGF_SHORTHAND = {
    "*": 15.994915,  # M* 
    "#": 57.021464,  # C# 
    "@": 79.966331,  # phospho 
    "^": 42.010565,  # acetyl N-term
}

# Isotope spacing constant (12C --> 13C)
_ISO_SPACING = 1.003355

# PeptideParser Class
class PeptideParser:
    """
    Static helpers for parsing MS-GF+ peptide strings.
    """

    _FLANK_RE = re.compile(r'^[A-Za-z\-]\.(.+)\.[A-Za-z\-]$')

    @staticmethod
    def parse(peptide_str):
        """
        Parse a peptide string into a clean sequence and a
        modifications dict {residue_position: mass_shift}.

        """
        # Strip flanking residues  
        # Use regex so dots inside mod masses (+15.995) are not confused with the flank separator dots.
        m = PeptideParser._FLANK_RE.match(peptide_str)
        if m:
            peptide_str = m.group(1)

        clean_seq = ""
        modifications = {}
        i = 0
        pos = 0

        while i < len(peptide_str):
            ch = peptide_str[i]

            if ch.isalpha():
                clean_seq += ch.upper()
                i += 1

                # Explicit mass shift: +15.995 or -17.027
                if i < len(peptide_str) and peptide_str[i] in ("+", "-"):
                    sign = 1.0 if peptide_str[i] == "+" else -1.0
                    j = i + 1
                    while j < len(peptide_str) and (
                        peptide_str[j].isdigit() or peptide_str[j] == "."
                    ):
                        j += 1
                    if j > i + 1:
                        modifications[pos] = sign * float(peptide_str[i + 1 : j])
                        i = j

                # Shorthand symbol: *, #, @, ^
                elif i < len(peptide_str) and peptide_str[i] in _MSGF_SHORTHAND:
                    modifications[pos] = _MSGF_SHORTHAND[peptide_str[i]]
                    i += 1

                pos += 1
            else:
                i += 1  # skip unexpected non-alpha characters

        return clean_seq, modifications

    @staticmethod
    def mod_features(modifications):
        """
        Derive scalar modification features from a modifications dictionary.
        """
        masses = list(modifications.values())
        total = round(sum(masses), 4)
        n = len(masses)

        return {
            "n_mods": n,
            "total_mod_mass": total,
            "has_oxidation": int(any(14.99 < m < 16.01 for m in masses)),
            "has_phospho": int(any(79.95 < m < 79.98 for m in masses)),
            "has_carbamido": int(any(57.01 < m < 57.04 for m in masses)),
            "has_deamidation": int(any(0.980 < m < 0.990 for m in masses)),
            "has_acetyl": int(any(42.00 < m < 42.02 for m in masses)),
            "mod_mass_per_res": 0.0,  # filled after knowing sequence length
        }

# IonCalculator
class IonCalculator:
    """
    Static helpers for computing theoretical fragment ion m/z values.
    """
    @staticmethod
    def fragment_mz(
        sequence,
        modifications,
        charge,
        ion_type,
        ion_number,
    ):
        """
        Monoisotopic m/z for a 'b' or 'y' fragment ion.
        Required:
        sequence, modifications, charge, ion_type, ion_number
        """
        if ion_type == "b":
            frag_seq  = sequence[:ion_number]
            mod_range = range(ion_number)
        elif ion_type == "y":
            frag_seq = sequence[-ion_number:]
            mod_range = range(len(sequence) - ion_number, len(sequence))
        else:
            return None

        mass = sum(AA_MASSES.get(aa, 0.0) for aa in frag_seq)

        for pos, delta in modifications.items():
            if pos in mod_range:
                mass += delta

        # b ions: no terminal addition; y ions: +H2O
        if ion_type == "y":
            mass += H2O

        return (mass + charge * H_MASS) / charge

    @staticmethod
    def build_ion_series(
        sequence ,
        modifications,
        max_charge = 1,
        min_mz = 100.0,
    ):
        """
        Build singly (and optionally doubly) charged b/y ion lists.
        
        """
        b_ions = []
        y_ions = []
        n = len(sequence)

        for i in range(1, n):
            for z in range(1, max_charge + 1):
                b = IonCalculator.fragment_mz(sequence, modifications, z, "b", i)
                y = IonCalculator.fragment_mz(sequence, modifications, z, "y", i)
                if b and b > min_mz:
                    b_ions.append(b)
                if y and y > min_mz:
                    y_ions.append(y)

        return b_ions, y_ions

# SpectralFeatureExtractor
class SpectralFeatureExtractor:
    """
    Extracts spectral features for a set of PSMs from one mzML file.
    """

    def __init__(
        self,
        tolerance = 0.02,
        tolerance_mode = "da",
    ):
        self.tolerance = tolerance
        self.tolerance_mode = tolerance_mode

    # ------------------------------------------------------------------
    # public API configuration
    # ------------------------------------------------------------------
    def extract(self, mzml_file, psm_df, disable_progress):
        """
        Load mzML and compute spectral features for every PSM in psm_df.
        """
        spectra = self._load_spectra(mzml_file)
        if spectra is None:
            return None

        print(f"\n[INFO] Extracting spectral features for {len(psm_df):,} PSMs...")

        records: List[Dict] = []
        failed  = 0

        scan_col = "ScanNum" if "ScanNum" in psm_df.columns else "Scan"
        tasks = [
            (row[scan_col], str(row["Peptide"]), int(row["Charge"]))
            for _, row in psm_df.iterrows()
        ]

        n_workers = min(mp.cpu_count(), 8) # cap at 8
        with mp.Pool(
            processes = n_workers,
            initializer = _worker_init,
            initargs = (spectra, self),
        ) as pool:
            records = list(tqdm(
                pool.imap(_worker_fn, tasks, chunksize=500),
                total = len(tasks),
                desc = "Extracting spectral features",
                unit = "psm",
                leave = True,
                disable=disable_progress,
            ))
        failed = sum(1 for r in records if r.get("spectral_feature_error") == "scan_not_found")

        features_df = pd.DataFrame(records)

        n_ok = len(psm_df) - failed
        n_feat = len(features_df.columns) - 2  # exclude ScanNum + error col
        print(f"[INFO] Feature extraction complete — {n_ok:,}/{len(psm_df):,} PSMs")
        print(f"[INFO] {n_feat} spectral feature columns extracted")

        if "spectral_feature_error" in features_df.columns:
            errs = features_df["spectral_feature_error"].dropna().value_counts()
            if not errs.empty:
                print("[WARNING] Feature extraction errors:")
                for err, cnt in errs.items():
                    print(f"  {err}: {cnt:,} PSMs")

        return features_df

    @staticmethod
    def merge(psm_df, features_df):
        """
        Left-join spectral features onto the PSM table on ScanNum.
        Reports columns with missing values.
        """
        merged = psm_df.merge(
            features_df.drop(columns=["spectral_feature_error"], errors="ignore"),
            on="ScanNum",
            how="left",
            suffixes=("", "_spectral"),
        )

        feat_cols = [
            c for c in features_df.columns
            if c not in ("ScanNum", "spectral_feature_error")
        ]
        missing = merged[feat_cols].isnull().sum()
        missing = missing[missing > 0]
        if not missing.empty:
            print("[WARNING] Missing spectral features after merge:")
            for col, cnt in missing.items():
                print(f"  {col}: {cnt:,} PSMs ({cnt / len(merged) * 100:.1f}%)")

        return merged

    # ------------------------------------------------------------------
    # Internal: mzML loading
    # ------------------------------------------------------------------
    def _load_spectra(self, mzml_file):
        print(f"\n[INFO] Reading mzML file: {mzml_file}")
        print("[INFO] Loading spectra into memory...")
        spectra: Dict[int, dict] = {}
        last_ms1: Optional[dict] = None

        try:
            with mzml.read(mzml_file) as reader:
                for spectrum in reader:
                    raw_id   = spectrum.get("id", "")
                    id_parts = raw_id.split("scan=")
                    if len(id_parts) > 1 and id_parts[-1].strip().isdigit():
                        scan_num = int(id_parts[-1].strip())
                    else:
                        scan_num = spectrum.get("index")

                    if scan_num is None:
                        continue

                    ms_level = int(spectrum.get("ms level",
                                   spectrum.get("msLevel", 2)))

                    if ms_level == 1:
                        last_ms1 = {
                            "m/z array": spectrum["m/z array"],
                            "intensity array": spectrum["intensity array"],
                        }

                    elif ms_level == 2:
                        entry = {
                            "m/z array": spectrum["m/z array"],
                            "intensity array": spectrum["intensity array"],
                            "precursorList": spectrum.get("precursorList", {}),
                        }
                        if last_ms1 is not None:
                            entry["ms1 spectrum"] = last_ms1   # key _isotope_score expects

                        spectra[scan_num] = entry

            n_with_ms1 = sum(1 for s in spectra.values() if "ms1 spectrum" in s)
            print(f"[INFO] Loaded {len(spectra):,} MS2 spectra "
                  f"({n_with_ms1:,} with MS1 survey scan attached)")
            return spectra

        except Exception as exc:
            print(f"[ERROR] Failed to read mzML: {exc}")
            return None
    
    @staticmethod
    def _kmd_variance(matched_mz):
        """
        Calculate Kendrick Mass Defect (KMD) variance for matched ions.
        """
        if not matched_mz or len(matched_mz) < 3:
            return {"kmd_variance": 0.0}

        mz_arr = np.array(matched_mz)
        
        # CH2 base unit for Kendrick Mass calculation
        CH2_EXACT = 14.01565
        CH2_NOMINAL = 14.00000
        
        # Calculate Kendrick Mass and Kendrick Mass Defect
        km = mz_arr * (CH2_NOMINAL / CH2_EXACT)
        kmd = km - np.floor(km)
        # Variance of the defect: true series cluster tightly near 0 variance
        kmd_var = float(np.var(kmd))
        
        return {"kmd_variance": kmd_var}
    
    @staticmethod
    def _classical_ion_weight(sequence, modifications, ion_type, ion_number, charge):
        """
        Rule-based (non-learned) relative intensity weight for a b/y fragment,
        based on well-established physicochemical fragmentation rules:
        - Proline effect: strongly enhanced cleavage N-terminal to proline
        - Enhanced cleavage C-terminal to Asp/Glu at low charge (charge-remote
            fragmentation under the mobile proton model)
        - Very short fragments (1 residue) are rarely intense
        """
        n = len(sequence)
        weight = 1.0
        cleave_idx = ion_number if ion_type == "b" else n - ion_number

        n_term_res = sequence[cleave_idx - 1] if 0 <= cleave_idx - 1 < n else None
        c_term_res = sequence[cleave_idx] if 0 <= cleave_idx < n else None

        if c_term_res == "P":
            weight *= 3.0  # proline effect
        if charge <= 2 and n_term_res in ("D", "E"):
            weight *= 1.8 # acidic-residue enhanced cleavage
        if ion_number <= 1:
            weight *= 0.3  # terminal fragments rarely intense

        return weight
    
    # ------------------------------------------------------------------
    # Internal: per-PSM feature computation
    # ------------------------------------------------------------------
    def _psm_features(
        self, spectrum, peptide, charge, scan_num
        ):
        """
        Compute all features for one PSM. Returns a flat dict.
        Modification features are always populated; spectral features
        fall back to zero on any exception.
        """
        rec: Dict = {"ScanNum": scan_num, "spectral_feature_error": None}

        try:
            #  1. Parse peptide & emit modification features 
            clean_seq, mods = PeptideParser.parse(peptide)
            mod_feat = PeptideParser.mod_features(mods)
            mod_feat["mod_mass_per_res"] = (
                mod_feat["total_mod_mass"] / max(len(clean_seq), 1)
            )
            rec.update(mod_feat)

            if len(clean_seq) < 2:
                rec["spectral_feature_error"] = "peptide_too_short"
                return rec

            # 2. Validate spectrum arrays 
            if "m/z array" not in spectrum or "intensity array" not in spectrum:
                rec["spectral_feature_error"] = "missing_spectrum_data"
                return rec

            mz_arr = np.asarray(spectrum["m/z array"], dtype=float)
            int_arr = np.asarray(spectrum["intensity array"], dtype=float)

            if mz_arr.size == 0:
                rec["spectral_feature_error"] = "empty_spectrum"
                return rec

            # Peaks sorted by intensity descending (used for top-N logic)
            peaks = np.column_stack([mz_arr, int_arr])
            sort_idx = mz_arr.argsort() # bug 1 - sorting based on m/z - not intensity
            peaks = np.column_stack([mz_arr[sort_idx], int_arr[sort_idx]])

            #  3. Basic spectrum-level features
            rec["total_ion_current"] = float(int_arr.sum())
            rec["base_peak_intensity"] = float(int_arr.max())
            rec["num_peaks"] = len(peaks)

            # 4. Theoretical ion generation 
            max_frag_z = 2 if charge >= 3 else 1
            b_ions, y_ions = IonCalculator.build_ion_series(
                clean_seq, mods, max_charge=max_frag_z
            )

            # Singly charged only for coverage counting
            b_c1, y_c1 = IonCalculator.build_ion_series(clean_seq, mods, max_charge=1)
            if max_frag_z > 1:
                b_ions, y_ions = IonCalculator.build_ion_series(clean_seq, mods, max_charge=max_frag_z)
            else:
                b_ions, y_ions = b_c1, y_c1

            # 5. Ion matching 
            matched_b, int_b = self._match_ions(b_c1, peaks)
            matched_y, int_y = self._match_ions(y_c1, peaks)

            n_b = len(matched_b)
            n_y = len(matched_y)

            rec["num_matched_b_ions"] = n_b
            rec["num_matched_y_ions"] = n_y
            rec["num_matched_ions_total"] = n_b + n_y

            rec["b_ion_coverage"] = n_b / max(len(b_c1), 1)
            rec["y_ion_coverage"] = n_y / max(len(y_c1), 1)
            rec["total_ion_coverage"] = (n_b + n_y) / max(len(b_c1) + len(y_c1), 1)

            # 6. Intensity features
            tic = rec["total_ion_current"]
            total_matched  = sum(int_b) + sum(int_y)
            all_intensities = int_b + int_y

            rec["matched_intensity_fraction"]  = total_matched / max(tic, 1.0)
            rec["mean_matched_ion_intensity"]  = (
                float(np.mean(all_intensities)) if all_intensities else 0.0
            )

            total_b_int = sum(int_b) if int_b else 0.0
            total_y_int = sum(int_y) if int_y else 0.0
            rec["b_to_y_intensity_ratio"] = total_b_int / max(total_y_int, 1.0)
            rec["y_to_total_ratio"] = (
                total_y_int / max(total_b_int + total_y_int, 1.0)
            )
            # Spectral contrast: std / mean — measures peak heterogeneity
            if int_arr.size > 1:
                rec["intensity_contrast"] = (
                    float(np.std(int_arr)) / (float(np.mean(int_arr)) + 1.0)
                )
            else:
                rec["intensity_contrast"] = 0.0
            #vectorize for speed gain
            top50_idx = np.argsort(peaks[:, 1])[::-1][:50]
            top50_int = peaks[top50_idx, 1]
            top50_mz = peaks[top50_idx, 0]
            matched_set_arr = np.array(sorted(set(matched_b + matched_y))) if (matched_b or matched_y) else np.array([])
            if len(matched_set_arr):
                diffs = np.abs(top50_mz[:, None] - matched_set_arr[None, :])
                hit_mask = diffs.min(axis=1) <= self.tolerance
            else:
                hit_mask = np.zeros(len(top50_mz), dtype=bool)
            rec['explained_intensity_top50'] = float(top50_int[hit_mask].sum()) / max(float(top50_int.sum()), 1.0)

            #  7. Signal-to-noise
            rec.update(self._signal_to_noise(peaks, matched_b + matched_y))

            #  8. Consecutive ion series length 
            rec["longest_b_series"] = self._longest_series(matched_b, b_c1)
            rec["longest_y_series"] = self._longest_series(matched_y, y_c1)

            # 9. Fragment mass error 
            rec.update(
                self._fragment_mass_error(peaks, matched_b + matched_y)
            )
            rec.update(self._kmd_variance(matched_b + matched_y))

            # 10. Complementary ion pairs & ratio 
            n_complementary_pairs = self._complementary_pairs(
                matched_b, matched_y, b_c1, y_c1, len(clean_seq))
            
            rec["num_complementary_ion_pairs"] = n_complementary_pairs

            # High ratio indicates true fragmentation; low indicates chimeric noise
            rec["complementary_ion_ratio"] = (
                n_complementary_pairs / max(n_b + n_y, 1)
            )

            # 11. Standard spectral similarity (charge-1 ions)
            rec.update(self._spectral_similarity(peaks, b_c1 + y_c1))

            # 12. Charge-weighted spectral angle 
            rec.update(self._charge_weighted_angle(clean_seq, mods, charge, peaks))

            # 13. Isotope pattern & chimeric score 
            precursor_mz = self._precursor_mz(spectrum)
            if precursor_mz is not None:
                rec.update(
                    self._isotope_score(spectrum, precursor_mz, charge)
                )
            else:
                rec["isotope_pattern_score"] = 0.0
                rec["chimeric_indicator"] = 0.0

        except Exception as exc:
            rec["spectral_feature_error"] = str(exc)[:120]

        return rec

    # Internal: ion matching
    def _match_ions(self, theoretical, peaks):
        """
        Match theoretical m/z list against observed peaks.
        """
        matched_mz, matched_int = [], []
        mz_col = peaks[:, 0]
        for theo in theoretical:
            idx = np.searchsorted(mz_col, theo)
            for i in [idx - 1, idx]:
                if 0 <= i < len(mz_col):
                    diff = abs(mz_col[i] - theo)
                    tol_check = (diff / (theo + 1e-10) * 1e6 <= self.tolerance
                             if self.tolerance_mode == "ppm" else diff <= self.tolerance)
                    if tol_check:
                        matched_mz.append(theo)
                        matched_int.append(float(peaks[i, 1]))
                        break
        return matched_mz, matched_int

    # Internal: feature sub-computations
    def _signal_to_noise(self, peaks, matched_mz):
        tol  = self.tolerance
        mode = self.tolerance_mode
        mz_col = peaks[:, 0]

        is_matched = np.zeros(len(peaks), dtype=bool)
        if matched_mz:
            theo_arr = np.array(sorted(set(matched_mz)))
            idxs = np.searchsorted(mz_col, theo_arr)
            for k, i in enumerate(idxs):
                for j in (i - 1, i): # check both neighbours
                    if 0 <= j < len(mz_col):
                        d = abs(mz_col[j] - theo_arr[k])
                        hit = (d / (theo_arr[k] + 1e-10) * 1e6 <= tol) if mode == "ppm" else (d <= tol)
                        if hit:
                            is_matched[j] = True
                            break

        matched_intensities = peaks[is_matched, 1]
        unmatched_intensities = peaks[~is_matched, 1]

        noise = max(float(np.median(unmatched_intensities)) if len(unmatched_intensities) else 1.0, 1.0)
        mean_sig = float(np.mean(matched_intensities)) if len(matched_intensities) else 0.0

        return {
            "signal_to_noise": mean_sig / noise,
            "noise_floor_normalized_tic": float(peaks[:, 1].sum()) / (noise * max(len(peaks), 1)),
        }
    
    @staticmethod
    def _longest_series(
        matched,
        all_ions,
    ):
        """Length of the longest consecutive matched ion run."""
        if not matched:
            return 0
        #vectorize for speed gain
        ion_index = {v: i for i, v in enumerate(all_ions)}
        indices = sorted(ion_index[m] for m in matched if m in ion_index)
        if not indices:
            return 0
        best = cur = 1
        for i in range(1, len(indices)):
            if indices[i] == indices[i - 1] + 1:
                cur += 1
                best = max(best, cur)
            else:
                cur = 1
        return best

    def _fragment_mass_error(
        self,
        peaks: np.ndarray,
        matched_mz: List[float],
    ) -> Dict[str, float]:
        """Mean and std of fragment mass errors in ppm for matched ions."""
        tol  = self.tolerance
        mode = self.tolerance_mode
        errors: List[float] = []
        for theo in matched_mz:
            diff = np.abs(peaks[:, 0] - theo)
            idx  = np.argmin(diff)
            d = diff[idx]
            within = (d / (theo + 1e-10) * 1e6 <= tol) if mode == "ppm" else (d <= tol)
            if within:
                ppm = abs(peaks[idx, 0] - theo) / (theo + 1e-10) * 1e6
                errors.append(ppm)

        if errors:
            return {
                "mean_fragment_mass_error_ppm": float(np.mean(errors)),
                "std_fragment_mass_error_ppm": float(np.std(errors)),
            }
        return {
            "mean_fragment_mass_error_ppm": 0.0,
            "std_fragment_mass_error_ppm": 0.0,
        }

    @staticmethod
    def _complementary_pairs(
        matched_b, matched_y, b_ions, y_ions, pep_len
    ):
        """Count b_i / y_(n-i) complementary matched pairs."""
        b_index_map = {v: i for i, v in enumerate(b_ions)}
        y_index_map = {v: i for i, v in enumerate(y_ions)}
        b_indices = {b_index_map[m] + 1 for m in matched_b if m in b_index_map}
        y_indices = {y_index_map[m] + 1 for m in matched_y if m in y_index_map}
        pairs = 0
        for i in range(1, pep_len):
            if i in b_indices and (pep_len - i) in y_indices:
                pairs += 1
        return pairs

    def _spectral_similarity(
        self,
        peaks,
        all_theo_mz: List[float],
    ) -> Dict[str, float]:
        """
        Cosine similarity between observed and uniform-intensity
        theoretical spectrum.  Also computes Pearson correlation.
        """
        theoretical = np.zeros(len(peaks))
        mz_col = peaks[:, 0]  # sorted ascending searchsort method for speed
        for theo in all_theo_mz:
            i = np.searchsorted(mz_col, theo)
            for j in (i - 1, i):
                if 0 <= j < len(mz_col) and abs(mz_col[j] - theo) <= self.tolerance:
                    theoretical[j] = 1.0
                    break
        observed = peaks[:, 1]
        obs_norm  = observed    / (np.linalg.norm(observed)    + 1e-10)
        theo_norm = theoretical / (np.linalg.norm(theoretical) + 1e-10)

        dot = float(np.clip(np.dot(obs_norm, theo_norm), -1.0, 1.0))
        angle = float(np.arccos(dot))
        pearson = 0.0
        if theoretical.sum() > 0 and observed.size > 2:
            mask = np.isfinite(observed) & np.isfinite(theoretical)
            if mask.sum() > 2:
                o_v = observed[mask]
                t_v = theoretical[mask]
                if o_v.std() > 1e-10 and t_v.std() > 1e-10:
                    try:
                        r, _ = pearsonr(o_v, t_v)
                        pearson = float(r) if np.isfinite(r) else 0.0
                    except Exception:
                        pass

        return {
            "spectral_angle": angle,
            "spectral_sim2": dot ** 2,
            "spectral_pearson": pearson,
        }

    def _charge_weighted_angle(
        self,
        clean_seq, mods, charge, peaks
    ) -> Dict[str, float]:
        theoretical = np.zeros(len(peaks))
        b_c1, y_c1 = IonCalculator.build_ion_series(clean_seq, mods, max_charge=1)
        mz_col = peaks[:, 0]
        for mz in b_c1 + y_c1:
            i = np.searchsorted(mz_col, mz)
            for j in (i - 1, i):
                if 0 <= j < len(mz_col) and abs(mz_col[j] - mz) <= self.tolerance:
                    theoretical[j] = 1.0
                    break

        n_charged_ions = len(b_c1) + len(y_c1)

        if charge >= 3:
            b_c2, y_c2 = IonCalculator.build_ion_series(
                clean_seq, mods, max_charge=2, min_mz=200.0
            )
            b_c2_only = b_c2[len(b_c1):]
            y_c2_only = y_c2[len(y_c1):]
            for mz in b_c2_only + y_c2_only:
                i = np.searchsorted(mz_col, mz)
                for j in (i - 1, i):
                    if 0 <= j < len(mz_col) and abs(mz_col[j] - mz) <= self.tolerance:
                        theoretical[j] = max(theoretical[j], 0.6)
                        break
            n_charged_ions += len(b_c2_only) + len(y_c2_only)

        observed = peaks[:, 1]
        obs_norm = observed / (np.linalg.norm(observed) + 1e-10)
        theo_norm = theoretical / (np.linalg.norm(theoretical) + 1e-10)
        dot = float(np.clip(np.dot(obs_norm, theo_norm), -1.0, 1.0))
        angle = float(np.arccos(dot))

        return {
            "charge_weighted_spectral_angle": angle,
            "charge_weight_sim2": dot ** 2,
            "n_theoretical_ions_charged": n_charged_ions,
        }

    @staticmethod
    def _precursor_mz(spectrum: dict) -> Optional[float]:
        """Extract precursor m/z from spectrum metadata."""
        try:
            prec_list = spectrum.get("precursorList", {}).get("precursor", [])
            if prec_list:
                sel = prec_list[0].get("selectedIonList", {}).get("selectedIon", [])
                if sel:
                    return float(sel[0].get("selected ion m/z", 0)) or None
        except Exception:
            pass
        return None

    @staticmethod
    def _isotope_score(
        spectrum, precursor_mz, charge
    ) -> Dict[str, float]:
        """
        Score the M / M+1 / M+2 isotope pattern against the averagine
        expected ratio. Low score suggests chimeric precursor isolation.
        """
        # Get the MS1 survey scan
        ms1_key = next(
            (k for k in ("ms1 spectrum", "MS1 spectrum", "survey scan") if k in spectrum),
            None,
        )
        if ms1_key is None:
            return {"isotope_pattern_score": 0.0, "chimeric_indicator": 0.0}

        ms1_mz = np.asarray(spectrum[ms1_key]["m/z array"], dtype=float)
        ms1_int = np.asarray(spectrum[ms1_key]["intensity array"], dtype=float)

        spacing = _ISO_SPACING / max(charge, 1)
        peaks_found: List[float] = []

        for k in range(3):
            target = precursor_mz + k * spacing
            d  = np.abs(ms1_mz - target)
            idx    = np.argmin(d)
            peaks_found.append(float(ms1_int[idx]) if d[idx] < 0.02 else 0.0)

        m0, m1, m2 = peaks_found

        if m0 == 0.0:
            return {"isotope_pattern_score": 0.0, "chimeric_indicator": 0.0}

        # Averagine approximation: M+1/M ≈ 0.0011 × peptide_mass
        peptide_mass = precursor_mz * charge
        expected_ratio = min(0.0011 * peptide_mass, 0.8)
        observed_ratio = m1 / (m0 + 1e-10)

        iso_score = float(np.clip(1.0 - abs(observed_ratio - expected_ratio), 0.0, 1.0))

        # M+2 unexpectedly high relative to M+1
        chimeric = float((m2 / (m1 + 1e-10)) > 1.5)

        return {
            "isotope_pattern_score": iso_score,
            "chimeric_indicator": chimeric,
        }

# Module-level convenience functions
def extract_spectral_features(
    mzml_file,
    psm_df,
    tolerance = 0.02,
    tolerance_mode = "da",
    disable_progress=False,
) -> Optional[pd.DataFrame]:
    """
    Extract spectral features for all PSMs.
    """
    extractor = SpectralFeatureExtractor(
        tolerance=tolerance,
        tolerance_mode=tolerance_mode,
    )
    return extractor.extract(mzml_file, psm_df, disable_progress=disable_progress)

def merge_features_with_psms(
    psm_df, features_df,
    ):
    """
    Merge spectral features back into PSM dataframe.
    Convenience wrapper around SpectralFeatureExtractor.merge().
    """
    return SpectralFeatureExtractor.merge(psm_df, features_df)

def extract_spectral_features_mzml(
    msgf_tsv, mzml_file
    ):
    psm_df = pd.read_csv(msgf_tsv, sep="\t")
    if "ScanNum" not in psm_df.columns and "Scan" in psm_df.columns:
        psm_df["ScanNum"] = psm_df["Scan"]

    features_df = extract_spectral_features(mzml_file, psm_df)
    if features_df is None:
        raise RuntimeError("Spectral feature extraction failed")

    return merge_features_with_psms(psm_df, features_df)