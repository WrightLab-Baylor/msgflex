#!/usr/bin/env python3
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

def plot_density_msms_vs_is_decoy(df, output_prefix, bins = 100):
    df = df.copy()
    df['MSMSScore'] = pd.to_numeric(df['MSMSScore'], errors='coerce')
    df['is_decoy'] = df['is_decoy'].astype('category')
    df['MSMSScore'] = pd.to_numeric(df['MSMSScore'], errors='coerce')
    df['is_decoy'] = df['is_decoy'].astype('category')
    plt.figure(figsize=(8, 6))
    sns.kdeplot(data=df, x='MSMSScore', hue='is_decoy', fill=False, alpha=0.5, palette={False: 'blue', True: 'orange'})
    plt.xlabel('MSMS_Score')
    plt.ylabel('Density')
    plt.title('Density Plot of MSMS Score by Decoy Status')

    legend_elements = [
        Line2D([0], [0], color='blue',   label='Target'),
        Line2D([0], [0], color='orange', label='Decoy'),
    ]
    plt.legend(handles=legend_elements, title='is_decoy', loc='upper right', frameon=True)
    
    #histogram helper
    def _histogram_density(df, bins=10):
        hist, bin_edges = np.histogram(df, bins=bins, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        density_values = hist * np.diff(bin_edges)
        return bin_centers, density_values

    # Compute histogram-based density estimate
    bin_centers, density_values = _histogram_density(df['MSMSScore'], bins=bins)
    # Find peak density
    max_density_index = np.argmax(density_values)
    max_density_x = bin_centers[max_density_index]
    max_density_y = density_values[max_density_index]

    plt.axvline(x=max_density_x, color='cyan', linestyle='--', label=f'Peak Density ({max_density_y:.2f})')

    output_file = f"{output_prefix}_msms_score.pdf"
    plt.savefig(output_file)
    plt.close()  

def plot_density_ppm_vs_is_decoy(df, output_prefix):
    df = df.copy()
    df['absPPM'] = pd.to_numeric(df['absPPM'], errors='coerce')
    df = df.dropna(subset=['absPPM'])
    df['is_decoy'] = df['is_decoy'].astype('category')
    plt.figure(figsize=(8, 6))
    sns.kdeplot(data=df, x='absPPM', hue='is_decoy', fill=False, alpha=0.5, palette={False: 'blue', True: 'orange'})
    plt.xlabel('absParentMassError(ppm)')
    plt.ylabel('Density')
    plt.title('Density Plot of absParentMassError by Decoy Status')

    legend_elements = [
        Line2D([0], [0], color='blue',   label='Target'),
        Line2D([0], [0], color='orange', label='Decoy'),
    ]
    plt.legend(handles=legend_elements, title='is_decoy', loc='upper right', frameon=True)

    output_file = f"{output_prefix}_absppm.pdf"
    plt.savefig(output_file)
    plt.close()  
    
def generate_coverage_heatmap(merged_df, outdir, top_n = 100, cmap="viridis"):
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Extract sample columns
    sample_cols = [c for c in merged_df.columns if c != "Protein"]
    if not sample_cols:
        return  # nothing to plot

    data = merged_df.set_index("Protein")[sample_cols].copy()
    
    # Select top N proteins if there are too many
    if len(data) > top_n:
        top_idx = data.max(axis=1).sort_values(ascending=False).head(top_n).index
        data = data.loc[top_idx]

    # Only plot if we have at least 2 proteins and 2 samples
    if data.shape[0] > 1 and data.shape[1] > 1:
        plt.figure(figsize=(10, max(6, len(data) * 0.15)))
        sns.heatmap(data, cmap=cmap)
        plt.tight_layout()
        pdf_path = outdir / "coverage_heatmap.pdf"
        plt.savefig(pdf_path, bbox_inches="tight")
        plt.close()
        return pdf_path  # return path for logging
    return None