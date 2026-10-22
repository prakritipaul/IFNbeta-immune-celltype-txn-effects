"""Figs 1–2: composition bars and UMAP, plus the shared Kang 2018 colours and orders.

Contract: audit/contracts/F12_figures.md. Decisions C1 (no integration) and C2 (UMAP recipe) in audit/decisions.md.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# Kang 2018 Fig 3 colours (sampled from the paper's legend), donor order and stacking order (bottom → top), to match the paper
CT_COLORS = {"B cells": "#67b294", "CD14+ Monocytes": "#e48164", "CD4 T cells": "#818eb8", "CD8 T cells": "#d181a8",
             "Dendritic cells": "#95be56", "FCGR3A+ Monocytes": "#ebc62f", "Megakaryocytes": "#d2b283", "NK cells": "#a1a1a1"}
COND_COLORS = {"ctrl": "#1b74ac", "stim": "#dd1b26"}
DONOR_ORDER = ["patient_1039", "patient_107", "patient_101", "patient_1016", "patient_1244", "patient_1256", "patient_1488", "patient_1015"]
STACK_ORDER = ["B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells"]


def composition_counts(obs):
    """Cells per (donor, condition) × cell type.

    In:  obs with columns donor, condition, cell_type.
    Out: DataFrame, rows (donor, condition) in DONOR_ORDER × [ctrl, stim], columns STACK_ORDER.
    """
    t = pd.crosstab([obs["donor"], obs["condition"]], obs["cell_type"])
    return t.loc[pd.MultiIndex.from_product([DONOR_ORDER, list(COND_COLORS)]), STACK_ORDER]  # KeyError if a group is missing


def plot_fig1(counts):
    """Fig 1: per donor, ctrl and stim bars side by side, stacked by cell type; counts and proportions panels.

    In:  counts from composition_counts.
    Out: matplotlib Figure.
    """
    props = counts.div(counts.sum(axis=1), axis=0)
    x = np.arange(len(DONOR_ORDER)) * 3.0
    pos = np.concatenate([[xi, xi + 1] for xi in x])  # ctrl at xi, stim at xi + 1
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for ax, tab, ylabel in [(axes[0], counts, "cells"), (axes[1], props, "proportion of cells")]:
        bottom = np.zeros(len(tab))
        for ct in STACK_ORDER:
            ax.bar(pos, tab[ct].values, bottom=bottom, color=CT_COLORS[ct], label=ct, width=0.9, edgecolor="white", linewidth=0.5)
            bottom += tab[ct].values
        ax.set_ylabel(ylabel)
    axes[0].set_title("Cells per donor, condition and cell type")
    axes[1].set_title("Cell-type proportions per donor and condition")
    axes[1].set_xticks(pos, [c for _ in DONOR_ORDER for c in COND_COLORS], fontsize=7)
    for xi, d in zip(x, DONOR_ORDER):
        axes[1].text(xi + 0.5, -0.12, d, ha="center", va="top", fontsize=8, transform=axes[1].get_xaxis_transform())
    h, l = axes[0].get_legend_handles_labels()
    axes[0].legend(h[::-1], l[::-1], frameon=False, fontsize=8, bbox_to_anchor=(1, 1), loc="upper left")
    fig.tight_layout()
    return fig


def compute_umap(adata):
    """UMAP of all cells, recipe C2, no integration (C1). Every result-affecting parameter is passed explicitly.

    In:  adata with raw counts in X (not modified).
    Out: DataFrame indexed by adata.obs_names with columns umap1, umap2.
    Routine: 1) normalize_total to 1e4, log1p  2) 2,000 HVGs (seurat flavor, all cells pooled)
             3) keep HVGs, scale (clip 10)  4) PCA 50 comps  5) kNN graph k=15 on 30 PCs  6) UMAP min_dist 0.5.
    """
    a = sc.AnnData(X=adata.X.copy(), obs=adata.obs[[]].copy(), var=adata.var[[]].copy())  # fresh copy: no stale obsm/uns
    sc.pp.normalize_total(a, target_sum=1e4, exclude_highly_expressed=False)
    sc.pp.log1p(a, base=None)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat", n_bins=20, batch_key=None, subset=False)
    a = a[:, a.var["highly_variable"].values].copy()  # scale is per gene, so subsetting first gives the same PCA input
    sc.pp.scale(a, zero_center=True, max_value=10)
    sc.pp.pca(a, n_comps=50, zero_center=True, svd_solver="arpack", random_state=0, mask_var=None, dtype="float32")
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, use_rep="X_pca", knn=True, method="umap", metric="euclidean", random_state=0)
    sc.tl.umap(a, min_dist=0.5, spread=1.0, n_components=2, init_pos="spectral", random_state=0)
    return pd.DataFrame(a.obsm["X_umap"], index=adata.obs_names, columns=["umap1", "umap2"])


def plot_fig2(df):
    """Fig 2: one embedding, two panels coloured by condition and by cell type.

    In:  df with columns umap1, umap2, condition, cell_type (one row per cell).
    Out: matplotlib Figure; both panels draw the same rows in the same (seeded random) order.
    """
    o = df.iloc[np.random.default_rng(0).permutation(len(df))]  # random draw order so no group hides another
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, col, colors, title in [(axes[0], "condition", COND_COLORS, "UMAP: ctrl vs stim"),
                                   (axes[1], "cell_type", CT_COLORS, "UMAP: authors' cell types")]:
        ax.scatter(o["umap1"], o["umap2"], s=1, c=o[col].map(colors).values, rasterized=True)
        ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=v, label=k) for k, v in colors.items()],
                  frameon=False, fontsize=8, bbox_to_anchor=(1, 1), loc="upper left")
        ax.set(title=title, xlabel="UMAP1", ylabel="UMAP2", xticks=[], yticks=[])
        ax.grid(False)
    fig.tight_layout()
    return fig
