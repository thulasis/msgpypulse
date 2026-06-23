import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from pathlib import Path

def histogram_density(df, bins=10):
    hist, bin_edges = np.histogram(df, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    density_values = hist * np.diff(bin_edges)
    return bin_centers, density_values

def plot_density_msms_vs_is_decoy(df, output_prefix, bins = 100):
    plt.figure(figsize=(8, 6))
    sns.kdeplot(data=df, x='MSMSScore', hue='IsDecoy', fill=False, alpha=0.5, palette={False: 'blue', True: 'orange'})
    plt.xlabel('MSMS_Score')
    plt.ylabel('Density')
    plt.title('Density Plot of MSMS Score by Decoy Status')
    plt.legend(title='IsDecoy', labels=['False', 'True'], loc='upper right', frameon=True)
    # Compute histogram-based density estimate
    bin_centers, density_values = histogram_density(df['MSMSScore'], bins=bins)
    # Find peak density
    max_density_index = np.argmax(density_values)
    max_density_x = bin_centers[max_density_index]
    max_density_y = density_values[max_density_index]

    plt.axvline(x=max_density_x, color='cyan', linestyle='--', label=f'Peak Density ({max_density_y:.2f})')

    output_file = f"{output_prefix}_msms_score.pdf"
    plt.savefig(output_file)

def plot_density_ppm_vs_is_decoy(df, output_prefix):
    plt.figure(figsize=(8, 6))
    sns.kdeplot(data=df, x='absPPM', hue='IsDecoy', fill=False, alpha=0.5, palette={False: 'blue', True: 'orange'})
    plt.xlabel('absParentMassError(ppm)')
    plt.ylabel('Density')
    plt.title('Density Plot of absParentMassError by Decoy Status')
    plt.legend(title='IsDecoy', labels=['False', 'True'], loc='upper right', frameon=True)
    output_file = f"{output_prefix}_absppm.pdf"
    plt.savefig(output_file)

def generate_coverage_heatmap(merged_df: pd.DataFrame, outdir: Path, top_n: int = 100, cmap="viridis"):
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

