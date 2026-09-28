#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import sys
import pandas as pd

"""
Author: Harrison Hall, (PhD)
"""
def _is_decoy_check(protein):
    return isinstance(protein, str) and protein.startswith("XXX_")

def load_safe_tsv(path: Path, sep="\t") :
    try:
        return pd.read_csv(path, sep=sep)
    except Exception as e:
        raise RuntimeError(f"Failed reading {path}: {e}")

def pick_template_rows(fdr) :
    """Pick one template row per peptide (lowest QValue if available)."""
    tmp = fdr.copy()
    if "QValue" in tmp.columns:
        tmp["__Q__"] = pd.to_numeric(tmp["QValue"], errors="coerce")
        tmp = tmp.sort_values(["Peptide", "__Q__"], kind="mergesort")
    else:
        tmp["__Q__"] = pd.NA
        tmp = tmp.sort_values(["Peptide"], kind="mergesort")
    best = tmp.groupby("Peptide", as_index=False).head(1).drop(columns="__Q__", errors="ignore")
    return best

def append_syn_to_fdr(sic_path: Path, syn_path: Path, out_suffix="_withsyn.tsv") -> tuple[int, int]:
    fdr = load_safe_tsv(sic_path, sep="\t")
    syn = load_safe_tsv(syn_path, sep="\t")

    # Required columns
    for col in ("Peptide", "Protein"):
        if col not in fdr.columns:
            raise ValueError(f"{sic_path.name} missing required column: {col}")
        if col not in syn.columns:
            raise ValueError(f"{syn_path.name} missing required column: {col}")

    # Build sets of peptide–protein pairs
    fdr_pairs = fdr[["Peptide", "Protein"]].drop_duplicates()
    syn_pairs = syn[["Scan","Peptide", "Protein"]].drop_duplicates() #added Scan

    # Only consider peptides already present in FDR (so we can template-copy)
    syn_pairs = syn_pairs[syn_pairs["Peptide"].isin(fdr["Peptide"])]

    # New pairs = in SYN but not in FDR
    new_pairs = syn_pairs.merge(
        fdr_pairs.assign(_in_fdr=True),
        on=["Peptide", "Protein"],
        how="left"
    )
    new_pairs = new_pairs[new_pairs["_in_fdr"].isna()][["Peptide", "Protein"]]

    if new_pairs.empty:
        out_path = sic_path.with_name(sic_path.stem.replace("_PlusSICStats", "") + out_suffix)
        fdr.to_csv(out_path, sep="\t", index=False)
        return (0, fdr_pairs.shape[0])

    # Template row per peptide from FDR
    template = pick_template_rows(fdr)

    # Merge new pairs to their template rows
    new_rows = new_pairs.merge(template, on="Peptide", how="left", suffixes=("", "_tmpl"))

    # Replace template protein with the SYN protein (already in 'Protein')
    if "Protein_tmpl" in new_rows.columns:
        new_rows = new_rows.drop(columns=["Protein_tmpl"])

    # Ensure/compute IsDecoy for the appended rows only (existing FDR rows are untouched)
    if "is_decoy" not in new_rows.columns:
        # If the template lacked IsDecoy, create it
        new_rows["is_decoy"] = new_rows["Protein"].astype(str).map(_is_decoy_check).astype(int)
    else:
        # Overwrite for appended rows per your spec
        new_rows["is_decoy"] = new_rows["Protein"].astype(str).map(_is_decoy_check).astype(int)

    # Align appended-row columns to FDR column order
    new_rows = new_rows[[c for c in fdr.columns if c in new_rows.columns]]

    # Append; DO NOT drop duplicates globally (we must not remove any existing FDR rows)
    combined = pd.concat([fdr, new_rows], ignore_index=True)

    # Write output alongside the original FDR file
    out_name = sic_path.stem.replace("_PlusSICStats", "") + out_suffix
    out_path = sic_path.with_name(out_name)
    combined.to_csv(out_path, sep="\t", index=False)

    # Report: how many appended and now how many unique peptide–protein pairs
    return (len(new_rows), combined[["Peptide", "Protein"]].drop_duplicates().shape[0])

def find_pairs(sic_dir: Path, sic_stats_dir="SICs", syn_dir="results/PHRPOut"):
    """
    Find matching FDR/SYN file pairs.
    """
    sic_dir = sic_dir.resolve()
    sic_root = sic_dir / sic_stats_dir
    syn_root = sic_dir.parent / syn_dir  # SYN next to SICdir

    if not sic_root.is_dir():
        raise FileNotFoundError(f"Missing FDR directory: {sic_root}")
    if not syn_root.is_dir():
        raise FileNotFoundError(f"Missing SYN directory: {syn_root}")

    pairs = []
    for sic_file in sic_root.glob("*_PlusSICStats.tsv"):
        base = sic_file.name[:-len("_PlusSICStats.tsv")]
        syn_file = syn_root / f"{base}_syn.txt"
        if syn_file.exists():
            pairs.append((sic_file, syn_file))
        else:
            print(f"[WARN] No matching SYN for {sic_file.name} (looked for {syn_file.name})", file=sys.stderr)
    return pairs
