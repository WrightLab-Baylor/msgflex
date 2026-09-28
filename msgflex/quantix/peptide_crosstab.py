#!/usr/bin/env python3
"""
peptide crosstab generation
Generates a peptide x sample intensity cross-tab from filtered PSM TSV files.
"""

import pandas as pd

def detect_protein_col(df):
    """Detect protein column (Protein or ProteinID)."""
    for col in ("Protein", "ProteinID"):
        if col in df.columns:
            return col
    raise ValueError(f"Protein column not found. Available columns: {list(df.columns)}")

def load_peptides(tsv_path):
    """Load a single filtered PSM file and clean columns."""
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    protein_col = detect_protein_col(df)

    # Standardize and clean columns
    df["Peptide"] = df["Peptide"].astype(str).str.strip()
    if "PeptideFlanked" in df.columns:
        df["PeptideFlanked"] = df["PeptideFlanked"].astype(str).str.strip()
    else:
        df["PeptideFlanked"] = df["Peptide"]  # fallback
    df[protein_col] = df[protein_col].astype(str).str.strip()
    df["ParentIonIntensity"] = pd.to_numeric(df["ParentIonIntensity"], errors="coerce")

    df = df.dropna(subset=["Peptide", protein_col, "ParentIonIntensity"])
    df = df.rename(columns={protein_col: "Protein", "ParentIonIntensity": "Intensity"})

    
    df = df.groupby(["Peptide", "PeptideFlanked", "Protein"], as_index=False)["Intensity"].max() # changed from sum to max for biological relevance

    return df
    
def merge_peptide_ctab(tsvs):
    """
    Merge multiple sample PSM files into a peptide x sample matrix.
    """
    wide = None
    for tsv in tsvs:
        sample_name = tsv.stem.replace("_filtered", "")
        df = load_peptides(tsv)
        df = df.rename(columns={"Intensity": sample_name})

        if wide is None:
            wide = df
        else:
            wide = wide.merge(df, on=["Peptide", "PeptideFlanked", "Protein"], how="outer")

    # Fill missing values with 0
    sample_cols = [c for c in wide.columns if c not in ("Peptide", "PeptideFlanked", "Protein")]
    wide[sample_cols] = wide[sample_cols].fillna(0)

    # Sort for readability
    wide = wide.sort_values(["Protein", "PeptideFlanked", "Peptide"], kind="mergesort").reset_index(drop=True)

    return wide
