#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import re

def normalize_sample_name(name) :
    name = str(name).lower().strip()
    name = re.sub(r'\.fasta$|\.fa$|\.faa$|\.fas$', '', name)
    name = re.sub(r'_rescored$', '', name)  
    return name

# ---------- FASTA utilities ----------
def load_fasta_sequences(fasta_path):
    """
    Return {protein_id -> sequence}, where protein_id is the header's first token (up to first whitespace).
    """
    seqs = {}
    acc = None
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                header = line[1:].strip()
                acc = header.split()[0]
                seqs.setdefault(acc, [])
            else:
                if acc is not None:
                    seqs[acc].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}

def load_fasta_for_samples(db_dir, sample_cols, inputfile_path=None):
    """
    Load FASTA(s) for each sample based on inputfile.tsv (must contain 'msfilename' and 'database').
    Automatically updates database column if fasta file has _95percent in its name.
    Returns dict: {sample_col -> {protein_id: sequence}}
    """
    if inputfile_path is None:
        inputfile = Path("inputfile.tsv")
    else:
        inputfile = Path(inputfile_path)
        
    if not inputfile.exists():
        raise FileNotFoundError(f"{inputfile} not found \n This file must contain at least 'msfilename' and 'database' columns.")

    # Read the mapping file
    mapping_df = pd.read_csv(inputfile, sep="\t")
    if "msfilename" not in mapping_df.columns or "database" not in mapping_df.columns:
        raise ValueError("inputfile.tsv must contain 'msfilename' and 'database' columns.")

    # Now load the sequences
    sample_to_fasta = {}
    for _, row in mapping_df.iterrows():
        sample = str(row["msfilename"])
        fasta_file = db_dir / str(row["database"])
        if not fasta_file.exists():
            raise FileNotFoundError(f"FASTA file '{fasta_file}' for sample '{sample}' not found in {db_dir}.")
        sample_to_fasta[sample] = load_fasta_sequences(fasta_file)

    # Check that all required sample columns have mapping
    # Normalize mapping keys
    normalized_map = {
        normalize_sample_name(k): v for k, v in sample_to_fasta.items()
    }

    resolved_sample_to_fasta = {}
    valid_samples = []

    for s in sample_cols:
        norm_s = normalize_sample_name(s)
        if norm_s in normalized_map:
            resolved_sample_to_fasta[s] = normalized_map[norm_s]
            valid_samples.append(s)

    missing = [s for s in sample_cols if s not in resolved_sample_to_fasta]
    if missing:
        print(f"Warning: Samples missing from inputfile.tsv mapping: {missing}")
        print(f"Available samples in mapping: {list(sample_to_fasta.keys())}")
        # Remove missing samples from sample_cols for processing
        sample_cols[:] = valid_samples
        if not sample_cols:
            raise RuntimeError("No valid samples found with FASTA mapping.")

    return resolved_sample_to_fasta, sample_cols

# ---------- Peptide utilities ----------
def parse_peptide_flanked(pflanked):
    parts = str(pflanked).split(".")
    if len(parts) != 3:
        return None, str(pflanked), None
    n, core, c = parts
    n = None if n in ("", "-") else n.upper()
    c = None if c in ("", "-") else c.upper()
    return n, core, c

def find_flanked_positions(protein_seq, core, n_flank, c_flank):
    seq = (protein_seq or "").upper()
    core_u = (core or "").upper()
    L = len(core_u)
    if L == 0:
        return []
    hits, start = [], 0
    while True:
        idx = seq.find(core_u, start)
        if idx == -1:
            break
        n_ok = True if n_flank is None else (idx > 0 and seq[idx - 1] == n_flank)
        c_ok = True if c_flank is None else (idx + L < len(seq) and seq[idx + L] == c_flank)
        if n_ok and c_ok:
            hits.append(idx)
        start = idx + 1
    return hits

def mark_core_coverage(protein_seq, core, starts, cov):
    L = len(core or "")
    if L == 0:
        return
    for s in starts:
        cov[s:s+L] = True

def build_coverage_with_flanks(protein_seq, rows_for_protein):
    cov = np.zeros(len(protein_seq or ""), dtype=bool)
    if rows_for_protein.empty:
        return cov
    for pflanked, core in rows_for_protein[["PeptideFlanked", "Peptide"]].drop_duplicates().itertuples(index=False):
        n_flank, core_seq, c_flank = parse_peptide_flanked(pflanked)
        if not core_seq:
            core_seq = core
        starts = find_flanked_positions(protein_seq, core_seq, n_flank, c_flank)
        mark_core_coverage(protein_seq, core_seq, starts, cov)
    return cov

# ---------- Column utilities ----------
RESERVED_NON_SAMPLE = {
    "Peptide", "PeptideFlanked", "Protein", "Shared",
    "Gene", "Function", "nested"
}

def ensure_required_columns(df):
    required = ["Peptide", "PeptideFlanked", "Protein"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Expected {required}.")
    for c in required:
        df[c] = df[c].astype(str).str.strip()

def detect_sample_columns(df):
    """
    Detect sample columns by excluding known non-sample columns.
    Also excludes columns that look like statistics or metadata.
    """
    potential_samples = []
    
    # Extended list of patterns that shouldn't be treated as samples
    exclude_patterns = [
        r".*_fdrstats$", # FDR statistics columns
        r".*_stats$", # General statistics columns  
        r".*_pvalue$", # P-value columns
        r".*_qvalue$", # Q-value columns
        r".*_score$", # Score columns
        r".*_rank$", # Rank columns
        r".*_count$", # Count columns
        r".*_total$", # Total columns
        r".*_avg$", # Average columns
        r".*_mean$", # Mean columns
        r".*_median$", # Median columns
        r".*_std$", # Standard deviation columns
        r".*_var$", # Variance columns
    ]
    
    for col in df.columns:
        if col in RESERVED_NON_SAMPLE:
            continue
            
        # Check if column matches any exclude pattern
        is_excluded = any(re.match(pattern, col, re.IGNORECASE) for pattern in exclude_patterns)
        
        if not is_excluded:
            # Additional check: see if column contains numeric data that could be sample data
            try:
                numeric_data = pd.to_numeric(df[col], errors='coerce')
                # If more than 50% of values are numeric and not all zero, likely a sample column
                if numeric_data.notna().sum() > len(df) * 0.5:
                    potential_samples.append(col)
            except:
                # If conversion fails completely, probably not a sample column
                continue
    
    return potential_samples

# ---------- Per-sample coverage ----------
def prepare_subtable_for_sample(df, sample_col, mode):
    cols = ["PeptideFlanked", "Peptide", "Protein", sample_col]
    sub = df[cols].copy()
    sub[sample_col] = pd.to_numeric(sub[sample_col], errors="coerce").fillna(0)
    sub = sub[sub[sample_col] > 0]

    if mode == "unique_only":
        if "Unique" not in df.columns:
            raise ValueError("Requested mode 'unique_only' but 'Unique' column is missing.")
        if df["Unique"].dtype != bool:
            df["Unique"] = df["Unique"].astype(bool)
        sub = sub.loc[df.loc[sub.index, "Unique"].values]

    elif mode == "requires_unique":
        if "Unique" not in df.columns:
            raise ValueError("Requested mode 'requires_unique' but 'Unique' column is missing.")
        keep_by_prot = df.groupby("Protein")["Unique"].any()
        sub = sub[sub["Protein"].map(keep_by_prot).fillna(False)]

    if not sub.empty:
        sub = (sub.sort_values(by=sample_col, ascending=False)
                  .drop_duplicates(subset=["Protein", "PeptideFlanked"]))
    return sub

def compute_coverage_for_sample(sub, protein_seqs, min_num_pep, sample_col):
    rows = []
    for prot_id, rows_for_prot in sub.groupby("Protein"):
        seq = protein_seqs.get(str(prot_id))
        if not seq:
            continue
        cov_vec = build_coverage_with_flanks(seq, rows_for_prot[["PeptideFlanked", "Peptide"]])
        coverage_pct = 100.0 * float(cov_vec.sum()) / max(1, len(seq))
        num_pep = int(rows_for_prot["PeptideFlanked"].nunique())
        rows.append((prot_id, num_pep, coverage_pct))

    out = pd.DataFrame(rows, columns=["Protein", "num_pep", "Coverage"])
    if min_num_pep > 1 and not out.empty:
        out = out[out["num_pep"] >= min_num_pep]

    return out.rename(columns={"Coverage": sample_col})