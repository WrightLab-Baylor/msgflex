#!/usr/bin/env python3
from __future__ import annotations

import os
import re
import time
import warnings
from typing import Optional

import numpy as np
import pandas as pd

from msgflex.quantix.density_plots import plot_density_msms_vs_is_decoy, plot_density_ppm_vs_is_decoy

warnings.filterwarnings("ignore")

# Score-field routing 

# "Smaller is better" — keep rows with value <= threshold
_SMALLER_IS_BETTER = frozenset({
    "QValue", "PepQValue", "xpec_q", "xpec_pepq",
})

# "Larger is better" — keep rows with value >= threshold
_LARGER_IS_BETTER = frozenset({
    "MSMSScore", "xpec_score",
})

# Sensible default thresholds when the user doesn't override
DEFAULT_THRESHOLDS: dict[str, float] = {
    "QValue": 0.01,
    "PepQValue": 0.01,
    "xpec_q": 0.01,
    "xpec_pepq": 0.01,
    "MSMSScore": 10.0,
    "xpec_score": 0.6,
}

VALID_FILTER_FIELDS = _SMALLER_IS_BETTER | _LARGER_IS_BETTER

# Columns we try to coerce to numeric on every read (both QUANTIX and XPECTRA)
_NUMERIC_COLUMNS = (
    "QValue", "PepQValue", "MSMSScore", "MSGFScore",
    "ParentIonIntensity", "PeakArea", "StatMomentsArea",
    "xpec_q", "xpec_pepq", "xpec_score",
    "absPPM",
)

# Helpers
def ensure_isdecoy(df):
    """
    Ensure is_decoy exists as int 0/1.
    """
    if "is_decoy" not in df.columns:
        df["is_decoy"] = (
            df["Protein"].astype(str).str.startswith("XXX_").astype(int)
        )
        return df

    col = df["is_decoy"]
    if pd.api.types.is_bool_dtype(col):
        df["is_decoy"] = col.astype(int)
    elif pd.api.types.is_numeric_dtype(col):
        df["is_decoy"] = col.fillna(0).astype(int)
    else:
        df["is_decoy"] = (
            col.astype(str).str.lower()
               .isin(["true", "1", "t", "yes"])
               .astype(int)
        )
    return df

def peptide_core(peptide):
    """
    Strip flanks and modification markers from MS-GF+ peptide strings.
    """
    s = str(peptide)
    # Strip flanks first: '<X>.peptide.<Y>' -> 'peptide'
    s = re.sub(r'^[A-Z\-]\.|\.[A-Z\-]$', '', s)
    # Then strip everything that isn't an uppercase letter (mod markers, +12.5, *, @, etc.)
    return re.sub(r'[^A-Z]', '', s.upper())

def add_peptide_core(df) :
    if "PeptideCore" not in df.columns:
        df["PeptideCore"] = df["Peptide"].map(peptide_core)
    return df

#  Filter 
def apply_metric_filter(df, metric, threshold):
    """
    Apply a single direction-aware filter.
    """
    if metric not in VALID_FILTER_FIELDS:
        raise ValueError(
            f"Unknown filter field {metric!r}. "
            f"Valid options: {sorted(VALID_FILTER_FIELDS)}"
        )

    if threshold is None:
        threshold = DEFAULT_THRESHOLDS.get(metric)
        if threshold is None:
            raise ValueError(
                f"No default threshold for {metric!r}. Specify one explicitly."
            )

    if metric not in df.columns:
        print(f"[WARN] Filter column '{metric}' not found — skipping filter.")
        return df

    col = pd.to_numeric(df[metric], errors="coerce")

    if metric in _SMALLER_IS_BETTER:
        return df[col <= threshold]
    else:  # _LARGER_IS_BETTER
        return df[col >= threshold]

def peptide_level_view(df) :
    """
    Collapse to one row per peptide core for peptide-level FDR stats.
    """
    agg = {"is_decoy": "all"}
    if "QValue" in df.columns:
        agg["QValue"] = lambda x: pd.to_numeric(x, errors="coerce").min()

    out = df.groupby("PeptideCore", as_index=False).agg(agg)
    return out.rename(columns={"PeptideCore": "Peptide"})

#  Main processing 
def process_fdrdir(
    input_directory,
    output_directory,
    metric1,
    threshold1 = None,
    metric2 = None,
    threshold2 = None,
    intensity_floor = 1e3
    ):
    """
    Filter *_withsyn.tsv files and write *_filtered.tsv.
    """
   
    start_time = time.time()
    os.makedirs(output_directory, exist_ok=True)

    # Normalize optional second filter
    apply_second = metric2 is not None and threshold2 is not None

    input_files = [f for f in os.listdir(input_directory) if f.endswith("_withsyn.tsv")]
    if not input_files:
        input_files = [f for f in os.listdir(input_directory) if f.endswith(".tsv")]

    # Cumulative counters
    total_files_processed = 0
    rows_processed = 0
    target_peptides_written = 0
    target_proteins_written = 0

    all_assignments = []
    all_peptides_for_stats = []

    for filename in sorted(input_files):
        input_file = os.path.join(input_directory, filename)

        stem = (
            filename[:-len("_withsyn.tsv")] if filename.endswith("_withsyn.tsv")
            else filename[:-len(".tsv")] if filename.endswith(".tsv")
            else filename
        )
        output_file = os.path.join(output_directory, f"{stem}_filtered.tsv")

        df = pd.read_csv(input_file, sep="\t", dtype=str)

        # Ensure MSMSScore exists for plotting
        if "MSMSScore" not in df.columns:
            if "MSGFDB_SpecEValue" in df.columns:
                df = df.copy()
                df["MSGFDB_SpecEValue"] = pd.to_numeric(
                    df["MSGFDB_SpecEValue"], errors="coerce"
                )
                valid_mask = df["MSGFDB_SpecEValue"] > 0
                df.loc[valid_mask, "MSMSScore"] = (
                    -np.log10(df.loc[valid_mask, "MSGFDB_SpecEValue"])
                )
            else:
                raise ValueError(
                    "[WARN] Neither MSMSScore nor MSGFDB_SpecEValue present."
                )

        # Ensure absPPM exists for density plots
        if "absPPM" not in df.columns:
            if "DelM_PPM" in df.columns:
                df["DelM_PPM"] = pd.to_numeric(
                    df["DelM_PPM"], errors="coerce"
                )
                df["absPPM"] = df["DelM_PPM"].abs()
            else:
                print(
                    "[WARN] absPPM not generated: DelM_PPM column missing."
                )

        for col in ("Peptide", "Protein"):
            if col not in df.columns:
                raise ValueError(f"{filename} missing required column: {col}")

        # Coerce numeric columns
        for c in _NUMERIC_COLUMNS:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        df = ensure_isdecoy(df)

        # Drop contaminants
        df = df[
            ~df["Protein"].astype(str).str.contains("Contaminant", case=False, na=False)
        ]

        # Primary filter (required)
        df_filtered = apply_metric_filter(df, metric1, threshold1)

        # Optional secondary filter
        if apply_second:
            df_filtered = apply_metric_filter(df_filtered, metric2, threshold2)

        # Peptide relabeling
        df_filtered = add_peptide_core(df_filtered)
        df_filtered["PeptideFlanked"] = df_filtered["Peptide"]
        df_filtered["Peptide"] = df_filtered["PeptideCore"]

        df_filtered = df_filtered.drop_duplicates()

        # Accumulate for cumulative stats
        stats_cols = ["Peptide", "PeptideCore", "Protein", "is_decoy"]        
        if "QValue" in df_filtered.columns:
            stats_cols.append("QValue")

        if "MSMSScore" in df_filtered.columns:
            stats_cols.append("MSMSScore")

        if "absPPM" in df_filtered.columns:
            stats_cols.append("absPPM")

        all_assignments.append(df_filtered[stats_cols].copy())
        all_peptides_for_stats.append(peptide_level_view(df_filtered))

        total_files_processed += 1
        rows_processed += len(df_filtered)

        # Per-file output: targets only
        out_df = df_filtered[df_filtered["is_decoy"] == 0].copy()

        if "ParentIonIntensity" in out_df.columns:
            out_df = out_df[out_df["ParentIonIntensity"] >= intensity_floor]

        out_df = out_df.drop(columns=["PeptideCore"], errors="ignore")

        target_peptides_written += out_df["Peptide"].nunique()
        target_proteins_written += out_df["Protein"].nunique()

        if "MSMSScore" in out_df.columns and "Protein" in out_df.columns:
            out_df = out_df.sort_values(
                by=["MSMSScore", "Protein"],
                ascending=[False, True],
                na_position="last",
                kind="mergesort",
            )
        elif "Protein" in out_df.columns:
            out_df = out_df.sort_values("Protein", ascending=True, kind="mergesort")

        out_df.to_csv(output_file, sep="\t", index=False)

    # Cumulative stats
    if all_assignments:
        combined_df = pd.concat(all_assignments, ignore_index=True)
        # combined_pep = (
        #     pd.concat(all_peptides_for_stats, ignore_index=True)
        #     .drop_duplicates(subset=["Peptide"])
        # )

        # pep_targets = int((combined_pep["is_decoy"] == 0).sum())
        # pep_decoys = int((combined_pep["is_decoy"] == 1).sum())
        # peptide_fdr_pct = pep_decoys / max(pep_targets + pep_decoys, 1) * 100.0

        # prot_df = (
        #     combined_df.assign(
        #         is_decoy=combined_df["Protein"]
        #         .astype(str)
        #         .str.startswith("XXX_")
        #         .astype(int)
        #     )
        #     .drop_duplicates(subset=["Protein"])[["Protein", "is_decoy"]]
        # )

        # prot_targets = int((prot_df["is_decoy"] == 0).sum())
        # prot_decoys = int((prot_df["is_decoy"] == 1).sum())
        # protein_fdr_pct = prot_decoys / max(prot_targets + prot_decoys, 1) * 100.0

        # Plots
        plot_dir = os.path.join(output_directory, "combined_plots")
        try:
            if "MSMSScore" not in combined_df.columns or combined_df["MSMSScore"].notna().sum() == 0:
                raise ValueError("Skipping MSMS plot: MSMSScore missing or empty")
            if "absPPM" in combined_df.columns and combined_df["absPPM"].notna().sum() == 0:
                print("[WARN] absPPM plot skipped: column empty")
            plot_density_msms_vs_is_decoy(combined_df, plot_dir)
            if "absPPM" in combined_df.columns:
                plot_density_ppm_vs_is_decoy(combined_df, plot_dir)
        except Exception as e:
            print(f"[WARN] Plot generation failed: {e}")
    else:
        peptide_fdr_pct = protein_fdr_pct = 0.0

    elapsed = time.time() - start_time

    # Summary
    filter_desc = (
        f"{metric1} @ "
        f"{threshold1 if threshold1 is not None else DEFAULT_THRESHOLDS.get(metric1)}"
    )
    if apply_second:
        filter_desc += f" + {metric2} @ {threshold2}"
