#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc, r2_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from .metrics import compute_top1_metrics

def plot_roc_pr_curves(
    df,
    score_col_ml = "rescore_prob",
    score_col_baseline = "MSGFScore",
    spectrum_col = "OptimalScanNumber",
    out_png = None,
    show_psm_level = False
):
    """
    Plot ROC and Precision-Recall curves for PSM rescoring.
    """
    sns.set_theme(style="white", rc={
        "axes.facecolor": "white",
        "figure.facecolor": "white"
    })
    required = ["is_decoy", score_col_ml, score_col_baseline, spectrum_col]
    for col in required:
        if col not in df.columns:
            raise KeyError(f"Missing column: {col}")
    
    # Clean data
    df_clean = df[required].dropna()
    
    # TOP-1 METRICS
    y_true_ml, y_score_ml = compute_top1_metrics(df_clean, score_col_ml, spectrum_col)
    y_true_base, y_score_base = compute_top1_metrics(df_clean, score_col_baseline, spectrum_col)
    
    # Validate
    if len(np.unique(y_true_ml)) < 2 or len(np.unique(y_true_base)) < 2:
        raise ValueError("Top-1 evaluation requires both positive and negative examples.")
    
    # ROC curves (Top-1)
    fpr_ml_top1, tpr_ml_top1, _ = roc_curve(y_true_ml, y_score_ml)
    fpr_base_top1, tpr_base_top1, _ = roc_curve(y_true_base, y_score_base)
    auc_ml_top1 = roc_auc_score(y_true_ml, y_score_ml)
    auc_base_top1 = roc_auc_score(y_true_base, y_score_base)
    
    # PR curves (Top-1)
    precision_ml_top1, recall_ml_top1, _ = precision_recall_curve(y_true_ml, y_score_ml)
    precision_base_top1, recall_base_top1, _ = precision_recall_curve(y_true_base, y_score_base)
    # Use AUC of PR curve (more informative than average_precision_score for imbalanced data)
    auprc_ml_top1 = auc(recall_ml_top1, precision_ml_top1)
    auprc_base_top1 = auc(recall_base_top1, precision_base_top1)
    baseline_prec_top1 = y_true_ml.mean()
    
    # PLOTTING
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ROC PLOT
    ax_roc = axes[1]
    
    # Top-1 curves (MAIN)
    ax_roc.plot(fpr_ml_top1, tpr_ml_top1, "b-", linewidth=3,
                label=f"ReScoreX(AUC={auc_ml_top1:.3f})", zorder=10)
    ax_roc.plot(fpr_base_top1, tpr_base_top1, "r-", linewidth=3,
                label=f"MS-GF+ (AUC={auc_base_top1:.3f})", zorder=9)
    
    # Random baseline
    ax_roc.plot([0, 1], [0, 1], color="gray", linestyle=":", linewidth=2,
                label="Random (AUC=0.500)", zorder=0)
    
    ax_roc.set_title("ROC Curve (Best PSMs)", fontsize=14, fontweight="bold")
    ax_roc.set_xlabel("False Positive Rate", fontsize=12)
    ax_roc.set_ylabel("True Positive Rate", fontsize=12)
    ax_roc.legend(loc="lower right", fontsize=10)
    ax_roc.grid(alpha=0.3)
    ax_roc.set_xlim([0, 1])
    ax_roc.set_ylim([0, 1.05])
    
    # PRECISION-RECALL PLOT
    ax_pr = axes[0]
    
    # Top-1 curves (MAIN)
    ax_pr.plot(recall_ml_top1, precision_ml_top1, "b-", linewidth=3,
               label=f"ReScoreX (AUPRC={auprc_ml_top1:.3f})", zorder=10)
    ax_pr.plot(recall_base_top1, precision_base_top1, "r-", linewidth=3,
               label=f"MS-GF+ (AUPRC={auprc_base_top1:.3f})", zorder=9)
        
    # Random baseline
    ax_pr.axhline(baseline_prec_top1, color="gray", linestyle=":", linewidth=2,
                  label=f"Random (Precision={baseline_prec_top1:.3f})", zorder=0)
    
    ax_pr.set_title("Precision-Recall Curve", fontsize=14, fontweight="bold")
    ax_pr.set_xlabel("Recall", fontsize=12)
    ax_pr.set_ylabel("Precision", fontsize=12)
    ax_pr.legend(loc="best", fontsize=10)
    ax_pr.grid(alpha=0.3)
    ax_pr.set_xlim([0, 1])
    ax_pr.set_ylim([0, 1.05])
    
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    
    if out_png:
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
    else:
        plt.show()
    
    plt.close()

def plot_psm_peptide_gains(
    psm_msgf_counts,
    psm_resc_counts,
    pep_msgf_counts,
    pep_resc_counts,
    out_png,
):
    sns.set_theme(style="white", rc={
        "axes.facecolor": "white",
        "figure.facecolor": "white"
    })
    MSGF_COLOR = "#C4C4C4"   # gray (baseline)
    RESCORE_COLOR = "#0072B2" 
    GAIN_POS = "#006B6E"     # same as rescored (consistent signal)
    GAIN_NEG = "#D55E00"    

    thresholds = list(psm_msgf_counts.keys())
    labels = [f"{int(q * 100)}%" for q in thresholds]
    x = np.arange(len(labels)) * 0.6
    width = 0.2
    offset = width * 0.45

    fig, axes = plt.subplots(
        nrows=2, ncols=1, figsize=(7, 7), sharex=True
    )

    def _plot_panel(ax, msgf_vals, resc_vals, ylabel, title):
        ax.bar(x - offset, msgf_vals, width, label="MS-GF+", color=MSGF_COLOR, edgecolor="black", linewidth=0.5)
        ax.bar(x + offset, resc_vals, width, label="XPECTRA", color=RESCORE_COLOR, edgecolor="black", linewidth=0.5)

        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=8, fontweight='bold')

        ymax = max(max(msgf_vals), max(resc_vals))
        pad = max(1, 0.03 * ymax)   # slightly larger buffer
        text_offset = 1.2 * pad 

        for i, (m, r) in enumerate(zip(msgf_vals, resc_vals)):
            gain = r - m
            gain_str = f"+{gain:,}" if gain >= 0 else f"{gain:,}"
            color = GAIN_POS if gain >= 0 else GAIN_NEG

            ax.text(
                x[i],
                max(m, r) + text_offset,
                gain_str,
                ha="center",
                va="bottom",
                fontsize=8,
                color=color,
            )

        ax.set_ylim(0, ymax + 5 * pad)
        ax.legend(fontsize=8, frameon=True)

    # ----- PSM panel -----
    _plot_panel(
        axes[0],
        list(psm_msgf_counts.values()),
        list(psm_resc_counts.values()),
        ylabel="Accepted Target PSMs",
        title="PSM-level gains",
    )

    # ----- Peptide panel -----
    _plot_panel(
        axes[1],
        list(pep_msgf_counts.values()),
        list(pep_resc_counts.values()),
        ylabel="Accepted Target Peptides",
        title="Peptide-level gains",
    )

    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels)
    axes[1].set_xlabel("FDR threshold", fontweight='bold')

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def plot_rt_correlation(
    common_df,
    gain_df,
    lost_df=None,
    *,
    out_png,
    observed_col="observed_rt",
    predicted_col="predicted_rt",
    score_col="neg_log_specEvalue",
    max_psms=5000,
    plot_mode=None,
    title_prefix=None
):
    """
    Combined RT calibration plot showing:
      - Common PSMs (MS-GF+ ∩ rescored)
      - Gained PSMs (rescored only)
      - Optional lost PSMs (MS-GF+ only)
    """

    sns.set_theme(style="white", rc={
        "axes.facecolor": "white",
        "figure.facecolor": "white",
    })

    # ------------------------------ helpers -----------------------------
    def _prep(df):
        if df is None or df.empty:
            return None
        df = df.copy()
        if "is_decoy" in df.columns:
            df = df[df["is_decoy"] == 0]
        df = df[
            df[observed_col].notna() &
            df[predicted_col].notna()
        ]
        return (
            df.sort_values(score_col, ascending=False)
              .head(max_psms)
              .copy()
        )

    def _scatter(ax, df, label, color, alpha, zorder, point_size):
        if df is None or df.empty:
            return
        ax.scatter(
            df[predicted_col],
            df[observed_col],
            s=point_size,
            alpha=alpha,
            color=color,
            edgecolors="none",
            rasterized=True,
            label=label,
            zorder=zorder,
        )

    def _kde(ax, df, label, color, zorder):
        if df is None or df.empty:
            return
        delta = df[observed_col] - df[predicted_col]
        sns.kdeplot(
            delta,
            ax=ax,
            label=label,
            color=color,
            lw=2,
            fill=False,
            alpha=0.9,
            zorder=zorder,
        )

    # ------------------------------ prep data ---------------------------
    common = _prep(common_df)
    gain = _prep(gain_df)
    lost = _prep(lost_df)

    if common is None or len(common) < 5:
        raise ValueError("Not enough common PSMs for RT evaluation.")

    metric_df = pd.concat(
        [d for d in (common, gain) if d is not None],
        axis=0,
    )

    observed  = metric_df[observed_col].astype(float).to_numpy()
    predicted = metric_df[predicted_col].astype(float).to_numpy()

    r, _  = pearsonr(observed, predicted)
    r2 = r2_score(observed, predicted)
    mae = mean_absolute_error(observed, predicted)
    rt_range = observed.max() - observed.min()
    nmae = mae / rt_range if rt_range > 0 else np.nan

    # ------------------------------ resolve plot mode -------------------
    n_total = len(common) + (len(gain) if gain is not None else 0)
    effective_mode = plot_mode or ("density" if n_total > 50_000 else "scatter")

    # Scatter appearance scales with mode
    if effective_mode == "density":
        pt_size = 4
        alpha_common = 0.08
        alpha_gain_lost = 0.20
    else:
        pt_size = 8
        alpha_common = 0.35
        alpha_gain_lost = 0.70

    # Wong (2011) colorblind-safe palette
    COLOR_COMMON = "#C0C0C0"   # grey      — recedes into background
    COLOR_GAIN   = "#0072B2"   # blue      — XPECTRA-gained PSMs
    COLOR_LOST   = "#E69F00"   

    # ------------------------------ plotting ----------------------------
    fig, axes = plt.subplots(2, 1, figsize=(7, 9))

    # -------------- Panel 1: Observed vs Predicted RT ------------------- 
    ax = axes[0]

    # Draw common first (lowest zorder) so gained/lost render on top
    _scatter(ax, common, "Common PSMs", COLOR_COMMON, alpha_common, 1, pt_size)
    _scatter(ax, gain, "Gained PSMs", COLOR_GAIN, alpha_gain_lost, 3, pt_size)
    _scatter(ax, lost, "Lost PSMs", COLOR_LOST, alpha_gain_lost, 4, pt_size)

    lim_lo = min(observed.min(), predicted.min())
    lim_hi = max(observed.max(), predicted.max())
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", lw=1.5, zorder=2)

    ax.set_xlabel("Predicted RT (min)", fontsize=8)
    ax.set_ylabel("Observed RT (min)", fontsize=8)
    ax.set_title(title_prefix, fontsize=8, fontweight='bold')

    ax.text(
        0.02, 0.98,
        f"r\u209A = {r:.2f}\n"
        f"R² = {r2:.2f}\n"
        f"nMAE = {nmae:.2f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox=dict(
            boxstyle="round,pad=0.35",
            fc="white",
            ec="black",
            alpha=1.0,
        ),
    )

    ax.legend(
        loc="lower right",
        fontsize=8,
        frameon=True,
        edgecolor="black",
        framealpha=0.95,
        borderpad=0.3,
        labelspacing=0.3,
    )

    # -------------- Panel 2: (delta)RT KDE ------------------------------------ 
    ax = axes[1]

    _kde(ax, common, "Common PSMs", COLOR_COMMON, 1)
    _kde(ax, gain, "Gained PSMs", COLOR_GAIN, 3)
    _kde(ax, lost, "Lost PSMs", COLOR_LOST, 4)

    ax.axvline(0, color="black", lw=1, ls="--", alpha=0.6)
    ax.set_xlabel("ΔRT = Observed - Predicted (min)", fontsize=8)
    ax.set_ylabel("Density", fontsize=8)
    ax.set_title("RT residual density", fontsize=8, fontweight='bold')

    ax.legend(
        loc="upper left",
        fontsize=8,
        frameon=True,
        borderpad=0.3,
        labelspacing=0.3,
    )

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

def plot_target_decoy_distribution(
    df,
    score_col,
    decoy_col="is_decoy",
    out_png=None
):
    """
    Target vs decoy histogram after global z-score normalization
    of the score column.
    """
    sns.set_theme(style="white", rc={
        "axes.facecolor": "white",
        "figure.facecolor": "white"
    })
    if score_col not in df.columns or decoy_col not in df.columns:
        raise ValueError("Required columns not found in dataframe")
    scaler = StandardScaler()
    df[f"{score_col}_z"] = scaler.fit_transform(df[[score_col]].fillna(0))
    
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.kdeplot(df[df[decoy_col] == 0][f"{score_col}_z"], label="Target", fill=True, alpha=0.5, color="#0072B2")
    sns.kdeplot(df[df[decoy_col] == 1][f"{score_col}_z"], label="Decoy", fill=True, alpha=0.5, color="#E69F00")

    ax.set_xlabel("Xpec_score (Z-scored)", fontsize=8)
    ax.set_ylabel("Density", fontsize=8)
    ax.set_title("Target vs Decoy Score Distribution", fontsize=8)
    ax.legend(fontsize=6, frameon=True)
    ax.tick_params(axis='both', labelsize=6)
    
    plt.tight_layout()
    if out_png:
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    
def plot_fdr_overlap(msgf_decoy_leakage, missed_df, out_png):
    sns.set_theme(style="white", rc={
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "xtick.labelsize": 8, "ytick.labelsize": 8
    })
    fig, ax = plt.subplots(figsize=(6, 4))

    sns.kdeplot(
        msgf_decoy_leakage["neg_log_specEvalue"],
        ax=ax,
        fill=True,
        alpha=0.5,
        color="#0072B2",
        linewidth=0.25,
        label=f"Leakage Decoys ({len(msgf_decoy_leakage)})",
    )

    sns.kdeplot(
        missed_df["neg_log_specEvalue"],
        ax=ax,
        fill=True,
        alpha=0.5,
        color="#E69F00",
        linewidth=0.25,
        label=f"Missed Targets ({len(missed_df)})",
        
    )
    ax.set_xlabel("-log SpecEValue", fontsize=6)
    ax.set_ylabel("Count", fontsize=6)
    ax.set_title("MS-GF+ Score Overlap: Leakage Decoys vs Missed Targets", fontsize=8, fontweight='bold')
    ax.legend(fontsize=6, frameon=True)
    ax.tick_params(axis='both', labelsize=6)
    plt.tight_layout()
    if out_png:
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

# def _logit(p, eps = 1e-6):
#     """
#     Safe logit transform. Clips to avoid inf at p=0 or p=1.
#     """
#     p = np.clip(p, eps, 1 - eps)
#     return np.log(p / (1 - p))


# def _zscore(a: np.ndarray) -> np.ndarray:
#     """Standard z-transform: (x - mean) / std."""
#     mu, sd = np.nanmean(a), np.nanstd(a)
#     if sd == 0 or np.isnan(sd):
#         return a - mu
#     return (a - mu) / sd
