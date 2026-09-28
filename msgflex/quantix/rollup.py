#!/usr/bin/env python3
"""
Protein rollup from annotated peptide data.

Author: Tulasi Rao Relangi, PhD
"""

from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import warnings
import re
warnings.filterwarnings("ignore")

# ---------- Helpers ----------
def _coerce_bool_series(s):
    """
    Convert series to boolean, handling various input formats.
    """
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.fillna(False).apply(lambda x: str(x).strip().lower() in ("true", "1", "t", "yes", "y"))

def _is_boolish_col(col):
    """Check if column contains boolean-like data."""
    if pd.api.types.is_bool_dtype(col):
        return True
    if pd.api.types.is_numeric_dtype(col):
        unique_vals = set(pd.unique(col.dropna()))
        return unique_vals <= {0, 1}
    uniq_vals = set(str(v).strip().lower() for v in pd.unique(col.dropna()))
    return uniq_vals <= {"true", "false", "t", "f", "yes", "no", "y", "n", "0", "1"}

def _infer_sample_cols(df, exclude=()):
    """Infer which columns contain sample intensity data."""
    exclude = set(exclude)
    numeric = [c for c in df.columns 
               if c not in exclude 
               and pd.api.types.is_numeric_dtype(df[c]) 
               and not _is_boolish_col(df[c])]
    if numeric:
        return numeric
    return [c for c in df.columns if c not in exclude and not _is_boolish_col(df[c])]

def _pick_peptide_key(df, prefer_flanked=True):
    """Select the best peptide identifier column."""
    if prefer_flanked and "PeptideFlanked" in df.columns:
        return "PeptideFlanked"
    if "Peptide" in df.columns:
        return "Peptide"
    peptide_cols = [col for col in df.columns if 'peptide' in col.lower()]
    if peptide_cols:
        return peptide_cols[0]
    raise ValueError("Input must contain a peptide identifier column (Peptide, PeptideFlanked, or similar)")

def _filter_to_coverage_ids(rolled_df, coverage_tsv):
    """Filter results to proteins present in grouped coverage file."""
    if not coverage_tsv.is_file():
        raise FileNotFoundError(f"Grouped coverage file not found: {coverage_tsv}")
    try:
        cov = pd.read_csv(coverage_tsv, sep="\t")
    except Exception as e:
        raise SystemExit(f"[ERROR] Failed reading coverage file {coverage_tsv}: {e}")
    # FIX: Validate required column exists in coverage file
    if "Protein" not in cov.columns:
        raise ValueError(f"Coverage file {coverage_tsv} must contain a 'Protein' column")
    cov_ids = set(cov["Protein"].astype(str))
    filtered_df = rolled_df[rolled_df["Protein"].astype(str).isin(cov_ids)].copy()
    return filtered_df

def _natural_sort_key(s):
    """Natural sorting key for sample column names."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]

def _validate_input_data(df, sample_cols):
    """Validate input data quality."""
    if df.empty:
        raise ValueError("Input dataframe is empty")
    if not sample_cols:
        raise ValueError("No sample columns found in the data")
    if "Protein" not in df.columns:
        raise ValueError("Required column 'Protein' is missing")

def _log2_transform_safe(df, sample_cols, zero_to_nan=True):
    """
    Log2 transform sample intensities.
    """
    df = df.copy()
    for col in sample_cols:
        vals = df[col].astype(float)
        if zero_to_nan:
            vals = np.where(vals > 0, vals, np.nan)
        else:
            vals = np.where(vals > 0, vals, 1.0)
        df[col] = np.log2(vals)
    return df

# ---------- RRollup internals ----------
def _pick_reference(peptab, sample_cols):
    """Select reference peptide with least missing values and highest median intensity."""
    if peptab.empty:
        raise ValueError("Empty peptide table for reference selection")
    miss = peptab[sample_cols].isna().sum(axis=1)
    # FIX: Use median instead of sum for robustness against outliers
    medians = peptab[sample_cols].median(axis=1, skipna=True)
    order = pd.DataFrame({
        'idx': peptab.index, 
        'miss': miss, 
        'median': medians
    })
    best = order.sort_values(['miss', 'median'], ascending=[True, False]).iloc[0]
    return int(best['idx'])

def _median_ratio(ref_vec, pep_vec, log_space=False):
    """Calculate median ratio between reference and peptide vectors."""
    mask = (~np.isnan(ref_vec)) & (~np.isnan(pep_vec))
    if not mask.any():
        return np.nan
    if log_space:
        # FIX: In log space, ratio is subtraction
        ratios = ref_vec[mask] - pep_vec[mask]
    else:
        mask = mask & (pep_vec != 0) & (ref_vec != 0)
        if not mask.any():
            return np.nan
        ratios = ref_vec[mask] / pep_vec[mask]
    ratios = ratios[np.isfinite(ratios)]
    return float(np.median(ratios)) if len(ratios) > 0 else np.nan

def _grubbs_filter(values, alpha=0.05):
    """
    Apply Grubbs test to filter outliers.
    """
    x = values.astype(float)
    mask = ~np.isnan(x)
    n = mask.sum()
    if n < 3:
        return mask
    vals = x[mask]
    mean_val = vals.mean()
    sd = vals.std(ddof=1)
    if sd == 0:
        return mask
    z = np.abs(vals - mean_val) / sd
    i = z.argmax()
    G = z[i]
    t_crit = stats.t.ppf(1 - alpha / (2 * n), n - 2)
    G_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit**2 / (n - 2 + t_crit**2))
    if G > G_crit:
        outlier_idx = np.flatnonzero(mask)[i]
        mask[outlier_idx] = False
    return mask

def _first_non_null(x):
    """Get first non-null value from series."""
    for v in x:
        if pd.notna(v) and str(v).strip() != "":
            return v
    return ""

# ---------- Core Rollup ----------
def rollup_from_annotated(
    pep_annot_tsv,
    rollup="sum",
    mode="all_matches",
    outlier_alpha=None,
    log_transform=None,
    report_log_transform=False,
    qrollup_fraction=1.0 / 3.0
):
    """
    Perform protein rollup from annotated peptide data.
    """
    
    try:
        df = pd.read_csv(pep_annot_tsv, sep="\t")
    except Exception as e:
        raise SystemExit(f"[ERROR] Failed reading {pep_annot_tsv}: {e}")
    
    pep_key = _pick_peptide_key(df, prefer_flanked=True)
    
    if "Unique" not in df.columns:
        df["Unique"] = False
    df["Unique"] = _coerce_bool_series(df["Unique"])
    
    gene_cols = [c for c in ("Gene", "gene") if c in df.columns]
    func_cols = [c for c in ("Function", "function") if c in df.columns]
    
    meta_cols = [pep_key, "Protein", "Unique", "Cluster"] + gene_cols + func_cols
    sample_cols = _infer_sample_cols(df, exclude=meta_cols)
    
    _validate_input_data(df, sample_cols)
    
    if "Peptide" not in df.columns and pep_key != "Peptide":
        df["Peptide"] = df[pep_key]
    
    required_cols = ["Protein", "Unique", pep_key]
    available_cols = [col for col in required_cols if col in df.columns]
    
    working_cols = available_cols.copy()
    if "Peptide" not in working_cols:
        working_cols.append("Peptide")
    working_cols.extend(sample_cols)
    
    if "Cluster" in df.columns:
        working_cols.append("Cluster")
    working_cols.extend(gene_cols + func_cols)
    
    slim = df[working_cols].copy()
    slim = slim.rename(columns={pep_key: "Peptide_pref"})
    if "Peptide" in slim.columns and pep_key != "Peptide":
        slim["Peptide_plain"] = slim["Peptide"].astype(str)
    else:
        slim["Peptide_plain"] = slim["Peptide_pref"].astype(str)


    rollups_log = {"rrollup", "zrollup", "qrollup"}
    true_apply_log = log_transform and rollup.lower() in rollups_log
  
    if true_apply_log:
        slim = _log2_transform_safe(slim, sample_cols)
    
    if mode == "unique_only":
        working = slim[slim["Unique"]]
    elif mode == "requires_unique":
        proteins_with_unique = slim.groupby("Protein")["Unique"].any()
        keep_proteins = proteins_with_unique[proteins_with_unique].index
        working = slim[slim["Protein"].isin(keep_proteins)]
    else:
        working = slim
    
    if working.empty:
        empty_cols = ["Protein"] + sample_cols + [
            "Peptide Number", "Peptides", "Flanked Peptides", 
            "Unique Peptide(s)", "Shared Peptide(s)"
        ]
        if "Cluster" in working_cols:
            empty_cols.insert(1, "Cluster")
        for c in gene_cols:
            if c not in empty_cols:
                empty_cols.append(c)
        for c in func_cols:
            if c not in empty_cols:
                empty_cols.append(c)
        return pd.DataFrame(columns=empty_cols)
    
    cluster_map = working.groupby("Protein")["Cluster"].first() if "Cluster" in working.columns else pd.Series(dtype=object)
    gene_map = {c: working.groupby("Protein")[c].apply(_first_non_null) for c in gene_cols}
    func_map = {c: working.groupby("Protein")[c].apply(_first_non_null) for c in func_cols}
    
    uniq_map = working.groupby("Protein")["Unique"].any()
    shared_map = working.groupby("Protein")["Unique"].apply(lambda s: (~s).any())
    
    flanked_peptide_map = working.groupby("Protein")["Peptide_pref"].apply(
        lambda x: ";".join(sorted(set(x.dropna().astype(str))))
    )
    plain_peptide_map = working.groupby("Protein")["Peptide_plain"].apply(
        lambda x: ";".join(sorted({s for s in x.dropna().astype(str) if s}))
    )
    
    if rollup.lower() == "sum":
        agg_fn = {c: "sum" for c in sample_cols}
        grouped = working.groupby("Protein", as_index=False).agg({
            **agg_fn,
            "Peptide_pref": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
            "Peptide_plain": lambda x: ";".join(sorted({s for s in map(str, x.dropna()) if s})),
        })
    
    else:
        prot_rows = []
        # fix weighted rollup
        if rollup.lower() == "dwrollup":
            global_pep_prot_counts = working.groupby("Peptide_pref")["Protein"].nunique()

        for prot, peptab in working.groupby("Protein"):
            vals = peptab[sample_cols].values.astype(float)
            
            if rollup.lower() == "rrollup":
                if peptab.shape[0] == 1:
                    prot_vals = vals[0]
                else:
                    try:
                        ref_idx = _pick_reference(peptab, sample_cols)
                        ref = peptab.loc[ref_idx, sample_cols].values.astype(float)
                        
                        scaled = []
                        for _, row in peptab.iterrows():
                            v = row[sample_cols].values.astype(float)
                            if row.name == ref_idx:
                                v_scaled = v.copy()
                            else:
                                sf = _median_ratio(ref, v, log_space=log_transform)
                                if np.isfinite(sf):
                                    if log_transform:
                                        v_scaled = v + sf
                                    else:
                                        if sf > 0:
                                            v_scaled = v * sf
                                        else:
                                            v_scaled = np.full_like(v, np.nan)
                                else:
                                    v_scaled = np.full_like(v, np.nan)
                            scaled.append(v_scaled)
                        
                        scaled = np.vstack(scaled)
                        
                        # Apply outlier filtering if requested
                        if outlier_alpha is not None:
                            for j in range(scaled.shape[1]):
                                col_mask = _grubbs_filter(scaled[:, j], alpha=outlier_alpha)
                                scaled[~col_mask, j] = np.nan
                        
                        prot_vals = np.nanmedian(scaled, axis=0)
                        
                    except Exception as e:
                        warnings.warn(f"RRollup failed for protein {prot}: {e}. Using median.")
                        prot_vals = np.nanmedian(vals, axis=0)
            
            elif rollup.lower() == "zrollup":
                if vals.shape[0] == 1:
                    warnings.warn(f"ZRollup not applicable for protein {prot}: only one peptide present")
                    prot_vals = np.full(vals.shape[1], np.nan)
                elif np.all(np.isnan(vals)):
                    warnings.warn(f"ZRollup not applicable for protein {prot}: all peptide values are NaN")
                    prot_vals = np.full(vals.shape[1], np.nan)
                else:
                    means = np.nanmean(vals, axis=1, keepdims=True)   
                    stds = np.nanstd(vals, axis=1, ddof=1, keepdims=True)

                    # Calculate Z-scores, ignoring division by zero warnings
                    with np.errstate(invalid='ignore', divide='ignore'):
                        vals_z = (vals - means) / stds
                    
                    vals_z[~np.isfinite(vals_z)] = np.nan
                    
                    # Check if still have at least one valid peptide left
                    if np.all(np.isnan(vals_z)):
                        warnings.warn(f"ZRollup not applicable for protein {prot}: all peptides had zero/NaN variance")
                        prot_vals = np.full(vals.shape[1], np.nan)
                    else:
                        # Take the median of the remaining valid Z-scores
                        prot_vals = np.nanmedian(vals_z, axis=0)
            
            elif rollup.lower() == "qrollup":
                pep_means = np.nanmean(vals, axis=1)
                if np.all(np.isnan(pep_means)):
                    prot_vals = np.full(vals.shape[1], np.nan)
                else:
                    threshold = np.nanpercentile(pep_means, (1 - qrollup_fraction) * 100)
                    keep = pep_means >= threshold
                    if not keep.any():
                        keep = pep_means == np.nanmax(pep_means)
                    prot_vals = np.nanmean(vals[keep], axis=0)
            
            elif rollup.lower() == "dwrollup":    
                weights = peptab["Peptide_pref"].map(lambda x: 1.0 / global_pep_prot_counts[x])
                # Normalize weights to sum to 1 for true weighted average
                weights = weights / weights.sum()
                weighted_vals = vals * weights.values[:, np.newaxis]
                prot_vals = np.nansum(weighted_vals, axis=0)
            
            else:
                raise ValueError(f"Unsupported rollup method: {rollup}")
            
            prot_rows.append((prot, prot_vals, len(peptab)))
        
        if not prot_rows:
            empty_cols = ["Protein"] + sample_cols + [
                "Peptide Number", "Peptides", "Flanked Peptides",
                "Unique Peptide(s)", "Shared Peptide(s)"
            ]
            if "Cluster" in working_cols:
                empty_cols.insert(1, "Cluster")
            for c in gene_cols:
                if c not in empty_cols:
                    empty_cols.append(c)
            for c in func_cols:
                if c not in empty_cols:
                    empty_cols.append(c)
            return pd.DataFrame(columns=empty_cols)
        
        grouped = pd.DataFrame({"Protein": [p for p, _, _ in prot_rows]})
        for j, col in enumerate(sample_cols):
            grouped[col] = [vals[j] for _, vals, _ in prot_rows]
        grouped["Peptide Number"] = [count for _, _, count in prot_rows]
    
    if not cluster_map.empty:
        grouped["Cluster"] = grouped["Protein"].map(cluster_map)
    
    for c in gene_cols:
        grouped[c] = grouped["Protein"].map(gene_map[c])
    
    for c in func_cols:
        grouped[c] = grouped["Protein"].map(func_map[c])
    
    grouped["Unique Peptide(s)"] = grouped["Protein"].map(uniq_map).fillna(False).astype(bool)
    grouped["Shared Peptide(s)"] = grouped["Protein"].map(shared_map).fillna(False).astype(bool)
    
    grouped["Flanked Peptides"] = grouped["Protein"].map(flanked_peptide_map).fillna("")
    grouped["Peptides"] = grouped["Protein"].map(plain_peptide_map).fillna("")
    
    if "Peptide Number" not in grouped.columns:
        base_list_col = "Flanked Peptides" if grouped["Flanked Peptides"].astype(str).str.len().gt(0).any() else "Peptides"
        grouped["Peptide Number"] = grouped[base_list_col].apply(
            lambda s: 0 if pd.isna(s) or str(s) == "" else len(set(str(s).split(";")))
        )
    
    for col in ["Peptide_pref", "Peptide_plain"]:
        if col in grouped.columns:
            grouped = grouped.drop(columns=[col])

    if true_apply_log and not report_log_transform and rollup.lower() != "zrollup":
            for col in sample_cols:
                grouped[col] = np.power(2, grouped[col])
    
    return grouped

def _finalize_column_order(df):
    """
    Organize columns in a logical order with naturally sorted sample names.
    """
    meta_cols = {"Protein", "Cluster", "Gene", "Function", "Peptide Number", 
                "Unique Peptide(s)", "Shared Peptide(s)", "Peptides", "Flanked Peptides"}
    
    sample_cols = [c for c in df.columns 
                   if c not in meta_cols 
                   and pd.api.types.is_numeric_dtype(df[c]) 
                   and not _is_boolish_col(df[c])]
    sample_cols = sorted(sample_cols, key=_natural_sort_key)
    
    desired_order = (["Protein"] + sample_cols + 
                    ["Cluster", "Gene", "Function", "Peptide Number",
                     "Unique Peptide(s)", "Shared Peptide(s)", "Peptides", "Flanked Peptides"])
    
    final_cols = [c for c in desired_order if c in df.columns]
    remaining_cols = [c for c in df.columns if c not in final_cols]
    
    return df[final_cols + remaining_cols]

def run_rollup(input_tsv, output_protein_tsv=None,
               rollup="sum", mode="all_matches",
               outlier_alpha=None, coverage_tsv=None,
               log_transform=True, report_log_transform=False):
    """
    Programmatic API that runs the same pipeline as main() but returns the rolled DataFrame.
    """
    input_path = Path(input_tsv).resolve()
    input_dir = input_path.parent

    rolled = rollup_from_annotated(
        pep_annot_tsv=input_path,
        rollup=rollup,
        mode=mode,
        outlier_alpha=outlier_alpha,
        log_transform=log_transform,
        report_log_transform=report_log_transform
    )

    cov_path = Path(coverage_tsv).resolve() if coverage_tsv else (input_dir / "map_files" / "grouped_coverage.tsv")
    rolled = _filter_to_coverage_ids(rolled, cov_path)

    rolled = _finalize_column_order(rolled)

    if output_protein_tsv is not None:
        out_path = Path(output_protein_tsv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        rolled.to_csv(out_path, sep="\t", index=False, na_rep="NA")

    return rolled