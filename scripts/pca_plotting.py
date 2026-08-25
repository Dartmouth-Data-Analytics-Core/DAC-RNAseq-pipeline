#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.optimize import curve_fit
import argparse
import sys
import os

# Plotting is imported behind a guard. The tabular outputs written below (PC scores,
# variance explained, normalised counts, gene variance) are what the dashboard and any
# downstream analysis consume, so a missing plotting backend must not block them.
try:
    from plotnine import *
    import seaborn as sns
    import matplotlib.pyplot as plt
    HAVE_PLOTTING = True
except ImportError as _e:
    HAVE_PLOTTING = False
    print(f"Plotting libraries unavailable ({_e}); writing tabular outputs only")

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
# Alongside the PNGs the script writes the underlying numbers, so they can be reused
# without recomputing the PCA (the QC dashboard reads these rather than redoing it):
#
#   PCA_top_PCs.csv                 sample x PC scores, high-variance gene set
#   PCA_top_variance_explained.csv  percent variance per PC, high-variance gene set
#   PCA_all_PCs.csv                 sample x PC scores, all genes
#   PCA_all_variance_explained.csv  percent variance per PC, all genes
#   PCA_top_normalized_counts.tsv   the log2 median-of-ratios matrix fed to the top PCA
#   PCA_all_normalized_counts.tsv   the same, for the all-gene PCA
#   PCA_gene_variance.csv           per-gene variance and rank, i.e. the HVG cut
#
# If the sample sheet passed with -m carries an optional `group` column, the PCA points
# are coloured by it and the PCs CSVs gain a `group` column. Without it nothing changes.
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
parser.add_argument("-m", "--metadata",
                    help="Sample sheet or metadata table (CSV or TSV) with a sample_id "
                         "column. An optional 'group' column colours the PCA by group.")
parser.add_argument("--colors",
                    help="Two-column group_name/color table for the group palette "
                         "(default: sample_ref/sample_colors_hex.tsv next to this script)")
parser.add_argument("-p", "--pca_comp", help="Number of PCA components", type=int, default=10)
args = parser.parse_args()

os.makedirs(args.output_path, exist_ok=True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ----------------------------
# Optional sample grouping
#
# `group` is an optional sample-sheet column. When it is present the PCA is coloured
# by it and the grouping travels into PCA_*_PCs.csv, so the dashboard colours its
# interactive PCA the same way. When it is absent nothing changes.
# ----------------------------
def read_table(path):
    """Read a CSV or TSV without caring which it is."""
    with open(path, errors="ignore") as fh:
        head = fh.readline()
    sep = "\t" if head.count("\t") > head.count(",") else ","
    return pd.read_csv(path, sep=sep)


def load_groups(path):
    """sample_id -> group from the sample sheet, or {} when there is no group column."""
    if not path:
        return {}
    if not os.path.exists(path):
        print(f"Sample sheet not found, PCA will not be grouped: {path}")
        return {}
    try:
        meta = read_table(path)
    except Exception as e:
        print(f"Could not read sample sheet ({e}); PCA will not be grouped")
        return {}

    cols = {str(c).strip().lower(): c for c in meta.columns}
    if "sample_id" not in cols or "group" not in cols:
        return {}

    groups = {}
    for sid, grp in zip(meta[cols["sample_id"]], meta[cols["group"]]):
        if pd.isna(grp) or str(grp).strip() == "":
            continue
        groups[str(sid).strip()] = str(grp).strip()
    return groups


def _hex_to_rgb(h):
    h = str(h).strip().lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    try:
        return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))
    except ValueError:
        return (128, 128, 128)


def spread_palette(colors):
    """
    Reorder a palette so consecutive picks are as visually distinct as possible.

    sample_colors_hex.tsv is a lookup keyed by group_name 1-6, not an ordered palette,
    and taken in file order its 3rd and 4th entries are two near-identical purples - a
    four-group run would draw two of them the same. Greedy farthest-point ordering keeps
    the file's own colours but makes the first N of them the most separable N.
    """
    rgb = [_hex_to_rgb(c) for c in colors]
    # Rough perceptual weighting; enough to keep two purples apart without a colour lib.
    w = (0.30, 0.59, 0.11)

    def dist(a, b):
        return sum(w[i] * (a[i] - b[i]) ** 2 for i in range(3))

    order = [0]
    while len(order) < len(colors):
        best = max(
            (i for i in range(len(colors)) if i not in order),
            key=lambda i: min(dist(rgb[i], rgb[j]) for j in order),
        )
        order.append(best)
    return [colors[i] for i in order]


def load_palette(path, levels):
    """
    Colours for the group levels.

    Prefers the repo's sample_ref/sample_colors_hex.tsv: matched by name when the group
    names are the ones it lists, otherwise spread for separability. Falls back to a
    built-in list, and cycles if a run has more groups than the palette has colours.
    """
    fallback = ["#dd0000", "#1037ff", "#26e79a", "#ffa34d", "#834194", "#630599",
                "#2b7a78", "#c23b8c", "#67932a", "#3f8ee0", "#c99700", "#8a6d3b"]
    table = path or os.path.join(SCRIPT_DIR, "..", "sample_ref", "sample_colors_hex.tsv")

    named, ordered = {}, []
    if os.path.exists(table):
        try:
            pal = read_table(table)
            cols = {str(c).strip().lower(): c for c in pal.columns}
            if "color" in cols:
                ordered = [str(c).strip() for c in pal[cols["color"]] if str(c).strip()]
                if "group_name" in cols:
                    named = {str(n).strip(): str(c).strip()
                             for n, c in zip(pal[cols["group_name"]], pal[cols["color"]])}
        except Exception as e:
            print(f"Could not read colour table ({e}); using the built-in palette")

    if named and all(l in named for l in levels):
        return [named[l] for l in levels]
    pool = spread_palette(ordered) if ordered else fallback
    return [pool[i % len(pool)] for i in range(len(levels))]


GROUP_MAP = load_groups(args.metadata)

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

# The normalised matrix behind PCA_all. Written before HVG selection so that anything
# downstream (sample distances, ad-hoc analysis) can use the same numbers the PCA saw.
df.rename_axis("Gene").to_csv(
    os.path.join(args.output_path, "PCA_all_normalized_counts.tsv"),
    sep="\t", float_format="%.6f"
)

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
gene_var_df.assign(Gene=gene_var.index)[["Gene", "Rank", "Variance"]].to_csv(
    os.path.join(args.output_path, "PCA_gene_variance.csv"), index=False
)

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
if HAVE_PLOTTING:
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

# The exact matrix the top PCA is computed from. Emitted so that sample-to-sample
# distances are derived from the same numbers as the PCA rather than a re-derivation.
df_top.rename_axis("Gene").to_csv(
    os.path.join(args.output_path, "PCA_top_normalized_counts.tsv"),
    sep="\t", float_format="%.6f"
)

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

    #----- optional grouping
    levels, palette = [], []
    if GROUP_MAP:
        pcadf["group"] = [GROUP_MAP.get(s, "unassigned") for s in sample_names]
        levels = sorted(set(pcadf["group"]), key=lambda g: (g == "unassigned", str(g)))
        pcadf["group"] = pd.Categorical(pcadf["group"], categories=levels, ordered=True)
        palette = load_palette(args.colors, levels)

    #----- persist the numbers behind the plots
    score_cols = ["sample"] + (["group"] if GROUP_MAP else []) + \
        [f"PC{i+1}" for i in range(ncomp)]
    pcadf[score_cols].to_csv(
        f"{args.output_path}/{prefix}_PCs.csv", index=False, float_format="%.6f"
    )
    pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(ncomp)],
        "Percent_Variance": pca.explained_variance_ratio_[:ncomp] * 100,
    }).to_csv(
        f"{args.output_path}/{prefix}_variance_explained.csv", index=False,
        float_format="%.6f"
    )

    def plot_pc_pair(pc_x, pc_y):
        x = f"PC{pc_x}"
        y = f"PC{pc_y}"

        pcadf["x_offset"] = pcadf[x] + 0.02 * (pcadf[x].max() - pcadf[x].min())
        pcadf["y_offset"] = pcadf[y] + 0.02 * (pcadf[y].max() - pcadf[y].min())

        if GROUP_MAP:
            p = (
                ggplot(pcadf, aes(x=x, y=y, color="group")) +
                geom_point(size=4.5) +
                scale_color_manual(values=palette, name="Group") +
                geom_text(
                    aes(x="x_offset", y="y_offset", label="sample"),
                    ha="left",
                    va="bottom",
                    size=8,
                    show_legend=False
                ) +
                theme_bw(base_size=16) +
                labs(
                    x=f"{x}: {pca.explained_variance_ratio_[pc_x-1]*100:.2f}%",
                    y=f"{y}: {pca.explained_variance_ratio_[pc_y-1]*100:.2f}%"
                )
            )
        else:
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

    if not HAVE_PLOTTING:
        return

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

if not HAVE_PLOTTING:
    print("Skipping heatmap: plotting libraries unavailable")
elif df_top_scaled.shape[0] >= 2 and df_top_scaled.shape[1] >= 2:
    try:
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
    except RecursionError:
        print(
            f"Skipping clustermap: recursion depth exceeded with "
            f"{df_top_scaled.shape[0]} genes. Try reducing the number of top genes."
        )
        plt.close("all")
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
