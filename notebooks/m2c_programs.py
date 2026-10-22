# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # M2c: shared programs (top genes and GSEA)
#
# Contract: `audit/contracts/M2c_programs.md` (decisions G1–G6). Python: `src/m2c_programs.py`; fgsea: `R/m2c_fgsea.R`
# (tested in `tests/test_m2c.py`). This notebook reads their output files (written by `run_m2c(data, data/processed, seed=0)`) and shows them.
#
# Two questions: (1) do cell types share their **strongest** genes? (2) do they run the same **programs** (Hallmark + KEGG gene sets)?

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform

pd.set_option("display.width", 160, "display.max_columns", 20, "display.max_colwidth", 60)
PROC = ROOT / "data/processed"
FILES = ["m2c_top_genes.csv", "m2c_top_sharing.csv", "m2c_jaccard.csv", "m2c_ranks.csv",
         "m2c_gene_set_coverage.csv", "m2c_gsea.csv", "m2c_programs.csv"]
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})
LINKAGE = "average"  # decided (D2 row)

top = pd.read_csv(PROC / "m2c_top_genes.csv")
sharing = pd.read_csv(PROC / "m2c_top_sharing.csv")
jac = pd.read_csv(PROC / "m2c_jaccard.csv", index_col=[0, 1])
gsea = pd.read_csv(PROC / "m2c_gsea.csv")
programs = pd.read_csv(PROC / "m2c_programs.csv")
cov = pd.read_csv(PROC / "m2c_gene_set_coverage.csv")
CELL_TYPES = sorted(top["cell_type"].unique())

# %% [markdown]
# ## 1. Every output file: shape, head and tail

# %%
for f in FILES:
    df = pd.read_csv(PROC / f)
    print(f"{f}: {df.shape}")
    display(df.head()); display(df.tail())

# %% [markdown]
# ## 2. Gene sets: how many of their genes are in our data
# Set genes are matched by exact symbol after the 6 B4 renames. Absent genes are listed, never fuzzily matched.
# GSEA then only tests a set in a cell type if 15–500 of its genes are among that cell type's ranked (baseMean > 1) genes.

# %%
s = cov.groupby("collection").agg(sets=("pathway", "size"), genes=("n_genes", "sum"), found=("n_found", "sum"))
s["pct_found"] = (100 * s["found"] / s["genes"]).round(1)
display(s)
print("sets tested per cell type:"); display(gsea.groupby(["cell_type", "collection"]).size().unstack())

# %% [markdown]
# ## 3. Top 20 up and top 20 down genes per cell type
# Only DE genes (q < 0.05 and |log2FC| ≥ 1), ranked by q, then |log2FC|. Fewer than 20 where fewer are DE.

# %%
for d in ("up", "down"):
    t = top[top["direction"] == d]
    wide = t.assign(g=t["gene"] + " (" + t["log2FC"].round(1).astype(str) + ")").pivot(
        index="rank", columns="cell_type", values="g").reindex(columns=CELL_TYPES)
    print(f"Top 20 {d} (log2FC in brackets)"); display(wide.fillna(""))

# %% [markdown]
# ## 4. Which genes are in the top 20 of several cell types?
# One row per gene in any top list; ↑ / ↓ = in that cell type's top 20 up / down; – = not in its top 20 (it may still be DE).

# %%
for d, col in (("up", "n_up"), ("down", "n_down")):
    print(f"\nTop-20 {d}, by number of cell types:")
    for k, g in sharing[sharing[col] > 0].groupby(col, sort=False)["gene"]:
        print(f"  in {k} cell types ({len(g)}): {', '.join(sorted(g))}")
display(sharing[(sharing["n_up"] + sharing["n_down"]) >= 3])


# %% [markdown]
# ## 5. Which cell types share their full DE lists? (Jaccard)
# Jaccard = shared DE genes / DE genes in either cell type, up and down separately (1 = identical lists, 0 = no overlap).
# Cell types are clustered on 1 − Jaccard; the dendrogram shows which cell types group together.

# %%
def heatmap_with_tree(mat, dist, title, ax_tree, ax_map, fmt="{:.2f}", cmap="viridis", vmin=0, vmax=1, marks=None):
    """Heatmap with rows/columns ordered by a dendrogram of `dist` (condensed); returns the linkage."""
    Z = linkage(dist, method=LINKAGE)
    order = dendrogram(Z, labels=list(mat.index), ax=ax_tree, orientation="left", color_threshold=0,
                       above_threshold_color="#555")["ivl"][::-1]
    ax_tree.set_axis_off()
    m = mat.loc[order]
    if m.shape[0] == m.shape[1]:
        m = m[order]
    im = ax_map.imshow(m.to_numpy(dtype=float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax_map.set_xticks(range(m.shape[1]), m.columns, rotation=60, ha="right"); ax_map.set_yticks(range(m.shape[0]), m.index)
    if fmt:
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                ax_map.text(j, i, fmt.format(m.iat[i, j]), ha="center", va="center", fontsize=7,
                            color="white" if m.iat[i, j] < (vmin + vmax) / 2 else "black")
    if marks is not None:
        mk = marks.loc[m.index, m.columns]
        ys, xs = np.nonzero(mk.to_numpy())
        ax_map.scatter(xs, ys, s=6, c="black", marker="o")
    ax_map.set_title(title)
    plt.colorbar(im, ax=ax_map, fraction=0.04)
    return Z


fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), gridspec_kw={"width_ratios": [1, 3, 1, 3]})
jac_groups = {}
for k, d in enumerate(("up", "down")):
    m = jac.loc[d][CELL_TYPES]
    assert not m.isna().any().any(), "Jaccard undefined (two empty DE sets)"
    Z = heatmap_with_tree(m, squareform(1 - m.to_numpy(), checks=False), f"Jaccard, DE {d}", axes[2 * k], axes[2 * k + 1])
    jac_groups[d] = Z
plt.tight_layout(); plt.show()

for d, Z in jac_groups.items():
    print(f"\nDE {d}: groups when the tree is cut into 2, 3, 4 clusters")
    for k in (2, 3, 4):
        lab = fcluster(Z, k, criterion="maxclust")
        print(f"  k={k}:", " | ".join(", ".join(c for c, l in zip(CELL_TYPES, lab) if l == g) for g in sorted(set(lab))))

# %% [markdown]
# ## 6. Programs: NES heatmap
# NES > 0: the gene set sits towards the up-regulated end of the cell type's ranking (Wald statistic); NES < 0: towards the down end.
# Rows = every program affected (BH-adjusted p < 0.05) in at least one cell type; dots = affected in that cell type.
# Rows and columns are clustered on correlation distance (1 − Pearson r) of NES. Grey hatched **n/t** = set not tested there (fewer than 15 or more than 500 of its genes ranked).

# %%
nes = gsea.pivot(index="pathway", columns="cell_type", values="NES").reindex(columns=CELL_TYPES)
sig = gsea.assign(s=gsea["padj"] < 0.05).pivot(index="pathway", columns="cell_type", values="s").reindex(columns=CELL_TYPES)
rows = programs["pathway"]
nes_aff, sig_aff = nes.loc[rows], sig.loc[rows].fillna(False).astype(bool)
print(f"programs affected anywhere: {len(rows)}; NES missing (set untested) in {int(nes_aff.isna().sum().sum())} cells")
complete = nes_aff.dropna()  # display only: correlation needs complete rows; counts reported above
print(f"programs with NES in all 8 cell types (used for clustering): {len(complete)}")

Zc = linkage(complete.T.to_numpy(), method=LINKAGE, metric="correlation")
Zr = linkage(complete.to_numpy(), method=LINKAGE, metric="correlation")
col_order = [CELL_TYPES[i] for i in dendrogram(Zc, no_plot=True)["leaves"]]
row_order = [complete.index[i] for i in dendrogram(Zr, no_plot=True)["leaves"]]
row_order += [p for p in nes_aff.index if p not in row_order]  # untested-somewhere programs at the bottom
lim = np.nanmax(np.abs(nes_aff.to_numpy()))

H = 0.16 * len(row_order) + 3
fig = plt.figure(figsize=(9, H))
ax_top = fig.add_axes([0.45, 1 - 1.1 / H, 0.4, 0.9 / H])
dendrogram(Zc, ax=ax_top, color_threshold=0, above_threshold_color="#555", no_labels=True)
ax_top.set_xlim(0, 10 * len(col_order)); ax_top.set_axis_off()  # leaves sit at 5, 15, ... = column centres
ax = fig.add_axes([0.45, 1.6 / H, 0.4, 1 - 2.8 / H])
cax = fig.add_axes([0.88, 0.4, 0.02, 0.2])
m = nes_aff.loc[row_order, col_order]
ax.set_facecolor("#bdbdbd")  # n/t (not tested: set size outside 15-500) = grey cells, hatched
im = ax.imshow(m.to_numpy(dtype=float), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
for y, x in zip(*np.nonzero(m.isna().to_numpy())):
    ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor="#bdbdbd", hatch="///", edgecolor="white", lw=0))
    ax.text(x, y, "n/t", ha="center", va="center", fontsize=5)
ys, xs = np.nonzero(sig_aff.loc[row_order, col_order].to_numpy())
ax.scatter(xs, ys, s=5, c="black")
ax.set_xticks(range(len(col_order)), col_order, rotation=60, ha="right", fontsize=8)
ax.set_yticks(range(len(row_order)), [p.replace("HALLMARK_", "H: ").replace("KEGG_", "K: ") for p in row_order], fontsize=6)
fig.colorbar(im, cax=cax, label="NES")
plt.show()

print("Cell-type groups from NES (correlation distance):")
for k in (2, 3, 4):
    lab = fcluster(Zc, k, criterion="maxclust")
    print(f"  k={k}:", " | ".join(", ".join(c for c, l in zip(CELL_TYPES, lab) if l == g) for g in sorted(set(lab))))
print("\nNES correlation between cell types (affected programs with complete NES):")
display(complete.corr().round(2))

# %% [markdown]
# ## 7. Program labels across the 8 cell types
# **shared** = affected in all 8, same direction; **partial** = 2–7, same direction; **specific** = 1;
# **mixed** = affected up in some cell types and down in others. Not affected anywhere → not listed.
# In `m2c_programs.csv`: ↑ / ↓ affected; – tested, not affected; **n/t** not tested (counts as not affected).

# %%
programs["direction"] = np.where(programs["label"] == "mixed", "mixed",
                                 np.where(programs["n_up"] > 0, "up", "down"))
display(pd.crosstab([programs["collection"]], [programs["label"], programs["direction"]], margins=True))
nt = (programs[CELL_TYPES] == "n/t")
print(f"n/t cells: {int(nt.sum().sum())} in {int(nt.any(axis=1).sum())} programs; per cell type:", nt.sum().to_dict())
for (lab, d), g in programs.groupby(["label", "direction"]):
    if lab == "mixed":
        items = [f"{p} (↑{u}/↓{n})" for p, u, n in zip(g["pathway"], g["n_up"], g["n_down"])]
    elif lab == "specific":
        items = [f"{p} [{next(c for c in CELL_TYPES if r[c] in ('↑', '↓'))}]" for p, (_, r) in zip(g["pathway"], g.iterrows())]
    else:
        items = [f"{p} ({u + n})" for p, u, n in zip(g["pathway"], g["n_up"], g["n_down"])]
    print(f"\n{lab} / {d} ({len(g)}): " + "; ".join(sorted(items)))

# %% [markdown]
# ## 8. Sanity check: IFN programs
# Hallmark IFN-α and IFN-γ response should be up in every cell type (M3 scores are up in all 8).

# %%
display(gsea[gsea["pathway"].str.contains("INTERFERON")].pivot(
    index="cell_type", columns="pathway", values=["NES", "padj"]).round(4))
