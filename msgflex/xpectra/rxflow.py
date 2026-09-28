#!/usr/bin/env python3
"""
XPECTRA End-to-End Workflow
"""
import sys
from pathlib import Path
import logging
import shutil
import pandas as pd
import numpy as np
import importlib.metadata
import contextlib
from contextlib import contextmanager
from tqdm import tqdm as _tqdm
from sklearn.metrics import roc_auc_score, average_precision_score

#from self modules
from .model import cross_validated_rescoring
from .features import merge_features_with_psms, extract_spectral_features
from .metrics import compare_gains
from .utils import _infer_is_decoy, _rt_stats
from .plots import plot_psm_peptide_gains, plot_rt_correlation, plot_target_decoy_distribution, plot_fdr_overlap

def find_matching_pairs(mzml_dir: Path, msgf_dir: Path):
    """
    Match mzML and MSGF+ files by base name.
    mzML: sample.mzML
    MSGF+: sample_fht_PlusSICStats.txt
    """
    mzml_files = {f.stem: f for f in mzml_dir.glob("*.mzML")}
    pairs = []

    for msgf in msgf_dir.glob("*.txt"):
        base = msgf.name.replace("_fht_PlusSICStats.txt", "")
        if base in mzml_files:
            pairs.append((base, mzml_files[base], msgf))

    return pairs

def compute_auc_stats(df, score_col="rescore_prob", decoy_col="is_decoy"):
    """
    Compute ROC AUC and PR AUC from a dataframe with scores and decoy labels.
    """
    y = 1 - df[decoy_col].astype(int).values
    scores = df[score_col].values

    # Guard against degenerate cases
    if len(np.unique(y)) < 2:
        return np.nan, np.nan

    auc = roc_auc_score(y, scores)
    # precision, recall, _ = precision_recall_curve(y, scores)
    pr_auc = average_precision_score(y, scores)
    return auc, pr_auc

def count_targets_at_fdr(df, qcol, thresh):
    """
    Count target PSMs passing FDR threshold.
    """
    return ((df[qcol] <= thresh) & (df["is_decoy"] == 0)).sum()

class StreamToLogger:
    def __init__(self, logger, level, prefix="[MODEL]"):
        self.logger = logger
        self.level = level
        self.prefix = prefix
        self._buffer = ""

    def write(self, message):
        self._buffer += message
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            if line.strip():
                self.logger.log(self.level, f"{self.prefix} {line}")

    def flush(self):
        if self._buffer.strip():
            self.logger.log(self.level, f"{self.prefix} {self._buffer}")
            self._buffer = ""

@contextmanager
def tqdm_to_console():
    import tqdm
    old_tqdm = tqdm.tqdm
    def console_tqdm(*args, **kwargs):
        kwargs["file"] = sys.__stderr__  # force real terminal
        kwargs["leave"] = True           # keep final bar clean
        return _tqdm(*args, **kwargs)
    tqdm.tqdm = console_tqdm
    try:
        yield
    finally:
        tqdm.tqdm = old_tqdm

def rxflow(
    mzml_dir,
    msgf_dir,
    out_dir,
    ensemble=5,
    folds=5,
    log_level="INFO",
    quiet=False,
    ):
    
    mzml_dir = Path(mzml_dir)
    msgf_dir = Path(msgf_dir)
    out_dir  = Path(out_dir)

    features_dir = out_dir / "features"
    rescored_dir = out_dir / "rescored"
    stats_dir    = out_dir / "stats"
    logs_dir     = out_dir / "logs"

    for d in (features_dir, rescored_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    # ---------- ROOT LOGGER ----------
    root_logger = logging.getLogger("rxflow")
    root_logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))

    # Remove any stale terminal StreamHandlers from previous runs.
    for h in root_logger.handlers[:]:
        if (isinstance(h, logging.StreamHandler)
                and not isinstance(h, logging.FileHandler)
                and getattr(h, "stream", None) in (
                    sys.stdout, sys.stderr, sys.__stdout__, sys.__stderr__)):
            root_logger.removeHandler(h)

    # Only add a terminal handler when running standalone
    has_external_handler = any(
        not (isinstance(h, logging.StreamHandler)
             and getattr(h, "stream", None) in (sys.__stdout__, sys.__stderr__))
        for h in root_logger.handlers
    )
    if not has_external_handler and not quiet:
        ch = logging.StreamHandler(sys.__stdout__)
        ch.setFormatter(logging.Formatter(
            "%(asctime)s  %(levelname)-7s  %(name)s  %(message)s",
            "%Y-%m-%d %H:%M:%S",
        ))
        root_logger.addHandler(ch)

    pairs = find_matching_pairs(mzml_dir, msgf_dir)
    if not pairs:
        root_logger.error("No matching mzML / MSGF+ pairs found")
        sys.exit(1)

    try:
        version = importlib.metadata.version("xpectra")
    except importlib.metadata.PackageNotFoundError:
        try:
            from msgflex import __version__ as version
        except ImportError:
            version = "unknown"
    root_logger.info(f"Running XPECTRA v{version}")
    root_logger.info(f"Found {len(pairs)} sample pairs")

    # ---------- Project-level RT plot configuration ----------

    # Scale max PSMs with sample count; cap at 150k to keep render time sane.
    PROJECT_RT_MAX_PSMS  = min(5000 * len(pairs), 150_000)
    PROJECT_RT_PLOT_MODE = "scatter" if len(pairs) <= 50 else "density"

    # ---------- Project-level accumulators (outside the per-sample loop) ----------
    qc_summary_rows = []

    proj_msgf_decoys = []    # for plot_fdr_overlap
    proj_missed_targets = []    # for plot_fdr_overlap

    proj_msgf_counts = {0.01: 0, 0.05: 0}
    proj_resc_counts = {0.01: 0, 0.05: 0}

    proj_pep_msgf_counts = {0.01: 0, 0.05: 0}
    proj_pep_resc_counts = {0.01: 0, 0.05: 0}

    proj_rescored_dfs = []    # for plot_target_decoy_distribution

    proj_rt_gain_dfs = []
    proj_rt_common_dfs = []
    proj_rt_lost_dfs = []

    n_skipped_existing  = 0

    # ======================================================================
    for sample, mzml, msgf in pairs:

        # ---------- SAMPLE LOGGER ----------
        logger = logging.getLogger(f"rxflow.{sample}")
        logger.setLevel(root_logger.level)
        logger.handlers.clear()
        logger.propagate = True

        fh = logging.FileHandler(logs_dir / f"{sample}.log")
        fh.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(fh)

        logger.info(f"Processing sample: {sample}")

        feat_out = features_dir / f"{sample}_spectral_features.tsv"
        resc_out = rescored_dir / f"{sample}_rescored.tsv"

        # skip existing results — check first, before doing any expensive work
        if resc_out.is_file():
            logger.info(f"Skipping: {resc_out.name} already rescored")
            n_skipped_existing += 1
            continue

        # ---------- PART 1 ----------
        psm_df = pd.read_csv(msgf, sep="\t")

        if "ScanNum" not in psm_df.columns and "Scan" in psm_df.columns:
            psm_df["ScanNum"] = psm_df["Scan"]

        logger.info(f"Loaded {len(psm_df):,} PSMs")

        # ---------- PART 2 ----------
        logger.info("Extracting spectral features...")

        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = StreamToLogger(logger, logging.INFO, prefix="[FEATURES]")
        sys.stderr = sys.__stderr__

        try:
            in_gui = any(
                not (isinstance(h, logging.StreamHandler)
                     and getattr(h, "stream", None) in (sys.__stdout__, sys.__stderr__))
                for h in logging.getLogger("rxflow").handlers
            )
            ctx = tqdm_to_console() if not in_gui else contextlib.nullcontext()
            with ctx:
                features_df = extract_spectral_features(
                    str(mzml), psm_df, disable_progress=in_gui
                )

            merged_df = (
                merge_features_with_psms(psm_df, features_df)
                if features_df is not None and not features_df.empty
                else psm_df.copy()
            )
        except Exception as e:
            logger.warning(f"Feature extraction failed: {e}", exc_info=True)
            features_df = None
            merged_df   = psm_df.copy()
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

        merged_df["is_decoy"] = _infer_is_decoy(merged_df)
        merged_df.to_csv(feat_out, sep="\t", index=False)

        # ---------- PART 3 ----------
        logger.info("Rescoring...")

        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = StreamToLogger(logger, logging.INFO)
        sys.stderr = StreamToLogger(logger, logging.ERROR)

        try:
            rescored_df, feat_importance = cross_validated_rescoring(
                merged_df, n_folds=folds, n_ensemble=ensemble, verbose=True, random_state=42,
            )
            feat_importance.to_csv(stats_dir / f"{sample}_feature_importance.tsv", sep="\t", index=False)

            #feature correlation analysis: save feature matrix for downstream analysis
            # if sample_feature_matrix is not None and sample_feature_cols:
            #     pd.DataFrame(sample_feature_matrix, columns=sample_feature_cols).to_csv(
            #         stats_dir / f"{sample}_xgb_matrix.tsv", sep="\t", index=False
            #     )
        except Exception as e:
            logger.error(f"Rescoring failed: {e}", exc_info=True)
            continue
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
        
        if rescored_df is None or rescored_df.empty:
            logger.warning("Empty rescoring output, skipping")
            continue
     
        # ---------- Output columns ----------
        _ALWAYS_DROP = {
            "log_msgf_score", "neg_log_specEvalue", "specevalue_rank",
            "abs_mass_error_ppm", "signed_mass_error",
            "charge_to_length_ratio", "num_basic", "num_acidic",
            "hydrophobic_ratio", "mass_error_zscore",
            "charge_x_length", "charge_x_mass_error",
            "peptide_length",
            "observed_rt", "predicted_rt", "peptide_rt_median",
            "peptide_rt_std", "rt_deviation_from_peptide",
            "rt_expected_for_charge", "rt_deviation_from_charge",
            "rt_per_residue", "peptide_hydrophobicity",
            "hydrophobicity_x_charge", "rt_x_hydrophobicity",
            "n_term_hydro", "c_term_hydro", "c_half_hydro",
            "rt_match_prob",
            "peptide_best_score", "peptide_median_score",
            "peptide_score_std", "peptide_reliability",
            "peptide_observation_count", "is_best_for_peptide",
            "target_cum", "decoy_cum",
            "delta_score_norm", 
            "is_singleton_scan",
            "n_psms_in_scan", "log_n_psms_in_scan",
            "longest_b_series", "longest_y_series",
            "spectral_pearson",
        }

        _ML_SCORES = [
            "rescore_prob", "xpec_q", "xpec_pepq",
            "xpec_score"
        ]

        _DIAG_KEEP = {"neg_log_specEvalue"} 
        base_cols = [c for c in psm_df.columns if c not in _ALWAYS_DROP]
        # base_cols = [c for c in rescored_df.columns if c not in _ALWAYS_DROP]
        # mod_cols = [c for c in _MOD_FEATURES if c in rescored_df.columns]
        score_cols = [c for c in _ML_SCORES if c in rescored_df.columns]
        # shap_cols = [c for c in rescored_df.columns if c.startswith("shap_")]
        
        # keep_cols = list(dict.fromkeys(base_cols + score_cols + list(_DIAG_KEEP)))
        keep_cols = list(dict.fromkeys(base_cols + score_cols))
        
        keep_cols = [c for c in keep_cols if c in rescored_df.columns]

        rescored_df[keep_cols].sort_values(
            "rescore_prob", ascending=False
        ).to_csv(resc_out, sep="\t", index=False)

        logger.info(f"Rescored data saved to {resc_out.name}")

        # ---------- PART 4 ----------
        logger.info("Computing statistics and generating plots...")
    
        qc_row = {"sample": sample}

        # ---------- MISSED PSMs (MS-GF+ yes, XPECTRA no) ----------
        logger.info("\n--- XPECTRA Missed PSMs Analysis ---")

        msgf_decoy_leakage = rescored_df[
            (rescored_df["QValue"] <= 0.01) &
            (rescored_df["is_decoy"] == 1)
        ]
        logger.info(
            f"MS-GF+ decoys that slipped through at 1% FDR: "
            f"{len(msgf_decoy_leakage)}"
        )

        missed_df = rescored_df[
            (rescored_df["is_decoy"] == 0) &
            (rescored_df["QValue"] <= 0.01) &
            (rescored_df["xpec_q"] > 0.01)
        ]
        qc_row["decoy_leak_1"]  = len(msgf_decoy_leakage)
        qc_row["lost_targets_1"] = len(missed_df)

        logger.info(f"Decoy leakage @1% FDR: {qc_row['decoy_leak_1']}")
        logger.info(f"Lost targets @1% FDR: {qc_row['lost_targets_1']}")

        if not msgf_decoy_leakage.empty:
            proj_msgf_decoys.append(msgf_decoy_leakage)
        if not missed_df.empty:
            proj_missed_targets.append(missed_df)

        # ---------- GAIN / LOSS LOGIC ----------
        rx_acc = (rescored_df["is_decoy"] == 0) & (rescored_df["xpec_q"] <= 0.01)
        msgf_acc = (rescored_df["is_decoy"] == 0) & (rescored_df["QValue"] <= 0.01)

        gain_df = rescored_df[rx_acc & ~msgf_acc]
        lost_df = rescored_df[~rx_acc & msgf_acc]
        common_df = rescored_df[rx_acc & msgf_acc]

        qc_row["psm_gained_1"] = len(gain_df)
        qc_row["psm_lost_1"] = len(lost_df)
        qc_row["net_psm_gain_1"] = len(gain_df) - len(lost_df)

        logger.info(f"Net PSM gain @1% FDR: {qc_row['net_psm_gain_1']}")

        zero_imp = feat_importance[feat_importance['importance_mean'] == 0]['feature'].tolist()
        spectral_zeros = [f for f in zero_imp if f in ('chimeric_indicator', 'isotope_pattern_score')]
        if spectral_zeros:
            logger.warning(
            f"Spectral features {spectral_zeros} have zero importance — "
            f"check extract_spectral_features output for these columns"
        )

        # Per-sample RT plot - non-calibrated
        try:
            plot_rt_correlation(
                common_df=common_df,
                gain_df=gain_df,
                lost_df=lost_df,
                observed_col="observed_rt",
                predicted_col="raw_predicted_rt",
                score_col="xpec_score",
                max_psms=5000,
                out_png=stats_dir / f"{sample}_rt_wo_calibration.png",
                title_prefix="RT (1% FDR): Raw Krokhin RC (uncalibrated)",
            )
        except Exception as e:
            logger.warning(f"RT combined correlation plot failed: {e}")

        # Per-sample RT plot - calibrated
        try:
            plot_rt_correlation(
                common_df=common_df,
                gain_df=gain_df,
                lost_df=lost_df,
                observed_col="observed_rt",
                predicted_col="predicted_rt",
                score_col="xpec_score",
                max_psms=5000,
                out_png=stats_dir / f"{sample}_rt_calibrated.png",
                title_prefix="RT (1% FDR): Isotonic warp + LOWESS + XGBoost residual correction",
            )
        except Exception as e:
            logger.warning(f"RT combined correlation plot failed: {e}")

        # RT summary statistics
        qc_row.update(_rt_stats(common_df, "common"))
        qc_row.update(_rt_stats(gain_df,   "gain"))
        qc_row.update(_rt_stats(lost_df,   "lost"))

        # Accumulate for project-level RT plot
        if not common_df.empty:
            proj_rt_common_dfs.append(common_df)
        if not gain_df.empty:
            proj_rt_gain_dfs.append(gain_df)
        if not lost_df.empty:
            proj_rt_lost_dfs.append(lost_df)

        # ---------- PSM & PEPTIDE GAINS (NO PER-SAMPLE PLOTS) ----------
        psm_improvements = compare_gains(
            rescored_df,
            msgf_q_col="QValue",
            resc_q_col="xpec_q",
            thresholds=(0.01, 0.05),
            verbose=False,
        )

        for q, res in psm_improvements.items():
            pct = int(q * 100)
            qc_row[f"psm_msgf_{pct}"] = res["n_msgf_targets"]
            qc_row[f"psm_rescored_{pct}"]  = res["n_rescore_targets"]
            qc_row[f"psm_gain_{pct}"]  = res["gain"]
            qc_row[f"psm_gain_pct_{pct}"]  = res["pct_gain_vs_msgf"]

            proj_msgf_counts[q] += res["n_msgf_targets"]
            proj_resc_counts[q] += res["n_rescore_targets"]

        pep_improvements = compare_gains(
            rescored_df,
            msgf_q_col="PepQValue",
            resc_q_col="xpec_pepq",
            thresholds=(0.01, 0.05),
            verbose=False,
        )

        for q, res in pep_improvements.items():
            pct = int(q * 100)
            qc_row[f"pep_msgf_{pct}"] = res["n_msgf_targets"]
            qc_row[f"pep_rescored_{pct}"]  = res["n_rescore_targets"]
            qc_row[f"pep_gain_{pct}"]  = res["gain"]
            qc_row[f"pep_gain_pct_{pct}"]  = res["pct_gain_vs_msgf"]

            proj_pep_msgf_counts[q] += res["n_msgf_targets"]
            proj_pep_resc_counts[q] += res["n_rescore_targets"]

        # ---------- APPEND TO PROJECT SUMMARY ----------
        qc_summary_rows.append(qc_row)
        proj_df_slice = rescored_df[[
            "rescore_prob", "xpec_score", "is_decoy", "neg_log_specEvalue", "QValue", "PepQValue",
            "xpec_q", "xpec_pepq"
        ]].copy()
        proj_df_slice.attrs = {}
        proj_rescored_dfs.append(proj_df_slice)

        # ---------- SAMPLE SUMMARY ----------
        logger.info("=" * 70)
        logger.info(f"SAMPLE: {sample} - Performance Evaluation")
        logger.info("=" * 70)

        logger.info(" PSM-level:")
        logger.info(
            f"  1% FDR | MS-GF+ = {qc_row['psm_msgf_1']:>4,}  "
            f"Rescored = {qc_row['psm_rescored_1']:>4,}  "
            f"Gain = {qc_row['psm_gain_1']:>+4,}  "
            f"({qc_row['psm_gain_pct_1']:>3.1f}%)"
        )
        logger.info(
            f"  5% FDR | MS-GF+ = {qc_row['psm_msgf_5']:>4,}  "
            f"Rescored = {qc_row['psm_rescored_5']:>4,}  "
            f"Gain = {qc_row['psm_gain_5']:>+4,}  "
            f"({qc_row['psm_gain_pct_5']:>3.1f}%)"
        )

        logger.info(" Peptide-level:")
        logger.info(
            f"  1% FDR | MS-GF+ = {qc_row['pep_msgf_1']:>4,}  "
            f"Rescored = {qc_row['pep_rescored_1']:>4,}  "
            f"Gain = {qc_row['pep_gain_1']:>+4,}  "
            f"({qc_row['pep_gain_pct_1']:>3.1f}%)"
        )
        logger.info(
            f"  5% FDR | MS-GF+ = {qc_row['pep_msgf_5']:>4,}  "
            f"Rescored = {qc_row['pep_rescored_5']:>4,}  "
            f"Gain = {qc_row['pep_gain_5']:>+4,}  "
            f"({qc_row['pep_gain_pct_5']:>3.1f}%)"
        )
        logger.info("=" * 70)
        logger.info(f"SAMPLE {sample} QC done")

    # ======================================================================
    # POST-LOOP: Write QC summary + project-level plots
    # ======================================================================
    if not qc_summary_rows:
        if n_skipped_existing > 0:
            root_logger.info(
                f"No new samples to process — {n_skipped_existing} sample(s) "
                "already had rescored output and were skipped. "
                "Delete the existing *_rescored.tsv files to force a re-run."
            )
            return
        root_logger.error("No samples processed successfully")
        raise RuntimeError("No samples processed successfully")

    # ---------- Write QC summary TSV ----------
    qc_df  = pd.DataFrame(qc_summary_rows)
    qc_out = stats_dir / "XPECTRA_qc_summary.tsv"
    qc_df.to_csv(qc_out, sep="\t", index=False)
    root_logger.info(f"QC summary written to {qc_out}")

    # ---------- Project-level plots ----------
    root_logger.info("Generating project-level QC plots...")

    if proj_msgf_decoys and proj_missed_targets:
        plot_fdr_overlap(
            pd.concat(proj_msgf_decoys,    ignore_index=True),
            pd.concat(proj_missed_targets, ignore_index=True),
            out_png=stats_dir / "project_fdr_overlap.png",
        )

    if proj_rescored_dfs:
        plot_target_decoy_distribution(
            pd.concat(proj_rescored_dfs, ignore_index=True),
            score_col="xpec_score",
            decoy_col="is_decoy",
            out_png=stats_dir / "project_target_decoy_dist.png",
        )

    plot_psm_peptide_gains(
        psm_msgf_counts=proj_msgf_counts,
        psm_resc_counts=proj_resc_counts,
        pep_msgf_counts=proj_pep_msgf_counts,
        pep_resc_counts=proj_pep_resc_counts,
        out_png=stats_dir / "project_psm_peptide_gains.png",
    )

    # ---------- Project-level RT plot ----------
    if proj_rt_common_dfs and proj_rt_gain_dfs:
        try:
            plot_rt_correlation(
                common_df=pd.concat(proj_rt_common_dfs, ignore_index=True),
                gain_df=pd.concat(proj_rt_gain_dfs, ignore_index=True),
                lost_df=(
                    pd.concat(proj_rt_lost_dfs, ignore_index=True)
                    if proj_rt_lost_dfs else None
                ),
                observed_col="observed_rt",
                predicted_col="predicted_rt",
                score_col="xpec_score",
                max_psms=PROJECT_RT_MAX_PSMS,
                plot_mode=PROJECT_RT_PLOT_MODE,
                out_png=stats_dir / "project_rt_combined.png",
                title_prefix="RT (1% FDR): Isotonic warp + LOWESS + XGBoost residual correction",
            )
            root_logger.info(
                f"Project-level RT plot saved  "
                f"(mode={PROJECT_RT_PLOT_MODE}, max_psms={PROJECT_RT_MAX_PSMS:,})"
            )
        except Exception as e:
            root_logger.warning(f"Project RT plot failed: {e}")

    root_logger.info("Rescoring & QC completed successfully")

    shutil.rmtree(features_dir)