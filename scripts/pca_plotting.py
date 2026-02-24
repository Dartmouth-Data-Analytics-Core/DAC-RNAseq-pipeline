#!/usr/bin/env python3

import pandas as pd
import numpy as np
from plotnine import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.optimize import curve_fit
import argparse
import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt

# ----------------------------
# README
#
# This script calculates simple PCA on raw read counts.
# Counts undergo median-of-ratios normalization followed by log2(x+1).
# Data are centered and scaled prior to PCA since no VST or Rlog is done
#
# An automatic variance asymptote detection is used to select
# highly variable genes (HVGs) for PCA.
#
# NOTE:
# This is intended for basic exploratory PCA only.
# vst/rlog remain preferred for downstream analyses.
# ----------------------------

# ----------------------------
# Argument Parsing
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("tsv_file", help="Path to TSV readcounts file")
parser.add_argument("output_path", help="Output folder for plots")
parser.add_argument("-m", "--metadata", help="Sample metadata TSV")
parser.add_argument("-p", "--pca_comp", help="Number of PCA components", type=int, default=10)
args = parser.parse_args()

os.makedirs(args.output_path, exist_ok=True)

# ----------------------------
# Load Data
# ----------------------------
try:
    with open(args.tsv_file) as f:
        first_line = f.readline()
    skiprows = 1 if first_line.startswith("#") else 0
    df = pd.read_csv(
        args.tsv_file,
        sep="\t",
        index_col=0,
        skiprows=skiprows,
        low_memory=False
    )
except Exception as e:
    print(f"Failed to load input TSV file: {e}")
    sys.exit(1)

# Remove featureCounts metadata columns if present
for col_set in [
    ["Gene Name", "Chr", "Start", "End", "Strand", "Length"],
    ["Chr", "Start", "End", "Strand", "Length"]
]:
    df = df.drop(columns=[c for c in col_set if c in df.columns], errors="ignore")

# Clean sample names
sample_names = [s.split("/")[-1].replace(".srt.bam", "") for s in df.columns.tolist()]
df.columns = sample_names

# ----------------------------
# Normalization
# ----------------------------
def median_of_ratios(df):
    df = df.astype(float)
    df_nozero = df[(df != 0).all(axis=1)]
    if df_nozero.shape[0] == 0:
        return df
    gene_geom = np.exp(np.log(df_nozero).mean(axis=1))
    size_factors = np.median(df_nozero.div(gene_geom, axis=0), axis=0)
    return df.div(size_factors, axis=1)

df = median_of_ratios(df)
df = np.log2(df + 1)

# ----------------------------
# Gene Variance & Asymptote Detection
# ----------------------------
gene_var = df.var(axis=1).sort_values(ascending=False)

gene_var_df = pd.DataFrame({
    "Rank": np.arange(1, len(gene_var) + 1),
    "Variance": gene_var.values
})

def decay(x, a, b, c):
    return a * np.exp(-b * x) + c

xdata = np.arange(len(gene_var))
ydata = gene_var.values

try:
    params, _ = curve_fit(
        decay,
        xdata,
        ydata,
        p0=[ydata.max() - ydata.min(), 0.001, ydata.min()]
    )
    a, b, c = params
except Exception:
    c = ydata.min()
    print("Curve fitting failed, using minimum variance as plateau")

tol = 0.01 * (ydata.max() - ydata.min())
plateau_indices = np.where(ydata <= c + tol)[0]
if len(plateau_indices) == 0:
    plateau_idx = len(ydata)
    print(f"No plateau detected, using all {plateau_idx} genes")
else:
    plateau_idx = plateau_indices[0] + 1
    print(f"Auto-detected plateau at rank {plateau_idx}, variance ~ {ydata[plateau_idx-1]:.4f}")

# ----------------------------
# Logging
# ----------------------------
log_file = os.path.join(args.output_path, "pca_hvg_log.txt")
with open(log_file, "w") as log:
    log.write("PCA automatic HVG thresholding log\n")
    log.write("---------------------------------\n")
    log.write(f"Total genes: {df.shape[0]}\n")
    log.write(f"Selected genes (plateau): {plateau_idx}\n")
    log.write(f"Plateau variance: {ydata[plateau_idx-1]:.6f}\n")

# ----------------------------
# Gene Variance Plot
# ----------------------------
p1 = (
    ggplot(gene_var_df, aes(x="Rank", y="Variance")) +
    geom_point(color="red", size=2.5) +
    geom_vline(xintercept=plateau_idx, linetype="dashed", color="blue", size=1) +
    theme_bw(base_size=16) +
    labs(
        title="Gene Expression Variance over Samples",
        x="Gene Rank",
        y="Variance"
    )
)

p1.save(
    f"{args.output_path}/Gene_Variance_Plot.png",
    width=8,
    height=6,
    dpi=300
)

# ----------------------------
# Select Top Genes
# ----------------------------
top_genes = gene_var.index[:plateau_idx]
df_top = df.loc[top_genes]

# ----------------------------
# Scaling
# ----------------------------
X_scaled_all = StandardScaler().fit_transform(df.T)
X_scaled_top = StandardScaler().fit_transform(df_top.T)

# ----------------------------
# PCA Function
# ----------------------------
def run_pca(X_scaled, sample_names, prefix, ncomp=args.pca_comp):
    ncomp = min(ncomp, X_scaled.shape[0], X_scaled.shape[1])
    pca = PCA(n_components=ncomp)
    pcs = pca.fit_transform(X_scaled)

    pcadf = pd.DataFrame(
        pcs,
        columns=[f"PC{i+1}" for i in range(ncomp)]
    )
    pcadf["sample"] = sample_names

    def plot_pc_pair(pc_x, pc_y):
        x = f"PC{pc_x}"
        y = f"PC{pc_y}"

        pcadf["x_offset"] = pcadf[x] + 0.02 * (pcadf[x].max() - pcadf[x].min())
        pcadf["y_offset"] = pcadf[y] + 0.02 * (pcadf[y].max() - pcadf[y].min())

        p = (
            ggplot(pcadf, aes(x=x, y=y)) +
            geom_point(color="black", size=4.5, fill="red") +
            geom_text(
                aes(x="x_offset", y="y_offset", label="sample"),
                ha="left",
                va="bottom",
                size=8
            ) +
            theme_bw(base_size=16) +
            labs(
                x=f"{x}: {pca.explained_variance_ratio_[pc_x-1]*100:.2f}%",
                y=f"{y}: {pca.explained_variance_ratio_[pc_y-1]*100:.2f}%"
            )
        )

        p.save(
            f"{args.output_path}/{prefix}_{x}_vs_{y}.png",
            width=8,
            height=6,
            dpi=300
        )

    plot_pc_pair(1, 2)
    if ncomp >= 3:
        plot_pc_pair(2, 3)
    if ncomp >= 4:
        plot_pc_pair(3, 4)

    pca_var_df = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(ncomp)],
        "Variance": pca.explained_variance_ratio_ * 100
    })

    pca_var_df["PC"] = pd.Categorical(
        pca_var_df["PC"],
        categories=[f"PC{i+1}" for i in range(ncomp)],
        ordered=True
    )

    pca_bar = (
        ggplot(pca_var_df, aes(x="PC", y="Variance")) +
        geom_bar(stat="identity", fill="skyblue", color="black") +
        theme_bw(base_size=16) +
        labs(x="Principal Component", y="Percent Variance Explained")
    )

    pca_bar.save(
        f"{args.output_path}/{prefix}_PCA_variance_bar.png",
        width=8,
        height=6,
        dpi=300
    )

# ----------------------------
# Run PCA
# ----------------------------
run_pca(X_scaled_all, sample_names, prefix="PCA_all")
run_pca(X_scaled_top, sample_names, prefix="PCA_top")

# ----------------------------
# Heatmap of Top Genes
# ----------------------------

# Z-score normalization (row-wise)
row_std = df_top.std(axis=1)
row_mean = df_top.mean(axis=1)

# Only keep rows with non-zero variance to avoid NaN from division by zero
valid_rows = row_std > 0
df_top_filtered = df_top.loc[valid_rows]
row_std = row_std[valid_rows]
row_mean = row_mean[valid_rows]

df_top_scaled = (
    (df_top_filtered - row_mean.values[:, None]) /
    row_std.values[:, None]
)

# Drop any remaining NaN values
df_top_scaled = df_top_scaled.dropna(axis=0, how='any')
df_top_scaled = df_top_scaled.dropna(axis=1, how='any')

# Remove columns with zero variance (e.g., identical samples)
col_variance = df_top_scaled.var(axis=0)
df_top_scaled = df_top_scaled.loc[:, col_variance > 0]

if df_top_scaled.shape[0] >= 2 and df_top_scaled.shape[1] >= 2:
    sns_clustermap = sns.clustermap(
        df_top_scaled,
        cmap="vlag",
        method="average",
        metric="euclidean",
        col_cluster=True,
        row_cluster=True,
        yticklabels=False,
        figsize=(12, 14)
    )
    sns_clustermap.fig.suptitle(
        f"Heatmap of Top {len(df_top_scaled)} Genes",
        y=1.05
    )
    sns_clustermap.savefig(
        f"{args.output_path}/Top_Genes_Heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
else:
    print(
        f"Skipping clustermap: insufficient data after filtering "
        f"({df_top_scaled.shape[0]} genes, {df_top_scaled.shape[1]} samples). "
        f"Need at least 2 of each."
    )
    # Fallback: simple heatmap without hierarchical clustering
    fig, ax = plt.subplots(figsize=(12, 14))
    sns.heatmap(
        df_top_scaled if not df_top_scaled.empty else df_top.fillna(0),
        cmap="vlag",
        yticklabels=False,
        ax=ax
    )
    ax.set_title(f"Heatmap of Top {len(df_top)} Genes (unclustered)")
    fig.savefig(
        f"{args.output_path}/Top_Genes_Heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

