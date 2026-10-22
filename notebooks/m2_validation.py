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
# # M2b: validation checks
#
# Contract: `audit/contracts/M2_shared_specific.md` (M2b). Definitions: `audit/validation_plan.md`, plus the user's
# "M2b checks" decisions (`audit/decisions.md`, M2 section). Logic: `src/m2b_validation.py`.
# Numbers and plots only; the one pass mark is V9 (> 70% of the paper's DE genes recovered).
# DE calls are re-derived with the current rule (q < 0.05 and |log2FC| ≥ 1, `m1b_de.classify`); the number of calls that
# differ from the files is printed.

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display
from scipy import stats

from src import m2b_validation as v
from src.figures import CT_COLORS, DONOR_ORDER, STACK_ORDER

pd.set_option("display.width", 160, "display.max_columns", 30)
PROC, FIG = ROOT / "data/processed", ROOT / "figures"
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e3", "grid.linewidth": 0.6, "axes.axisbelow": True, "figure.dpi": 110})

de, n_changed = v.rederive(pd.read_csv(PROC / "m1b_de.csv"))
assert (de["run"] == "main").all() and not de.duplicated(["cell_type", "gene"]).any()
print(f"main-run calls changed by the >= rule: {n_changed}")
fits = pd.read_csv(PROC / "m1b_fits.csv")
runs, n_changed_runs = v.read_runs(PROC / "m1b_runs.csv", ["shuffle_", "lodo_", "unpaired", "no_outliers"])
print(f"other-run calls changed by the >= rule: {n_changed_runs}")
ifna = v.ifna_genes(de["gene"].unique())
counts = pd.read_csv(PROC / "m1a_counts.csv", index_col="gene")
samples = pd.read_csv(PROC / "m1a_samples.csv")
print(f"main run: {len(de)} rows; other runs read: {runs['run'].nunique()}; Hallmark IFN-α genes: {len(ifna)}; "
      f"pseudobulk: {counts.shape[0]} genes × {counts.shape[1]} samples")


def save(df, name):
    df.to_csv(PROC / f"m2_validation_{name}.csv", index=False)
    print(f"m2_validation_{name}.csv: {df.shape}")


# %% [markdown]
# ## P1: canonical ISGs DE per cell type

# %%
p1 = v.p1_isgs(de).loc[STACK_ORDER]
display(p1)
save(p1.reset_index(), "p1_isgs")

fig, axes = plt.subplots(1, 2, figsize=(13, 4))
axes[0].barh(p1.index, p1["n_de"], color=[CT_COLORS[c] for c in p1.index])
axes[0].set(xlim=(0, 10), xlabel="canonical ISGs DE (of 10)", title="P1: ISGs DE per cell type"); axes[0].invert_yaxis()
hist = p1["n_de"].value_counts(normalize=True).reindex(range(11), fill_value=0)
axes[1].bar(hist.index, hist.values, color="#7a5aa6")
axes[1].set(xticks=range(11), xlabel="ISGs DE (0–10)", ylabel="proportion of cell types", title="P1: histogram over 8 cell types")
plt.tight_layout(); fig.savefig(FIG / "m2b_p1_isgs.png", dpi=200, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## P2 / P3: Fisher tests on Hallmark IFN-α
# P2: up-DEGs vs all tested genes (baseMean > 1) in that cell type. P3: up-DEGs vs down-DEGs. `inf` odds ratio = no IFN gene among the down-DEGs.
# Headline: one-sided ("greater") p, BH-adjusted across the 8 cell types (`p2_q_bh`). Two-sided p shown for reference.

# %%
fisher = v.p2_p3_fisher(de, ifna).set_index("cell_type").loc[STACK_ORDER].reset_index()
display(fisher)
save(fisher, "p2_p3_fisher")

# %% [markdown]
# ## P4: expression-matched null
# Per cell type: control mean = mean DESeq2-normalised count over the kept control samples (median-of-ratios size
# factors over all kept samples and all genes; equal to R `estimateSizeFactorsForMatrix` on Megakaryocytes to 6 decimals);
# detection = fraction of the cell type's control cells from kept donors with count > 0. Each up-DEG gets the 30 nearest
# tested genes on rank-percentiles of the two (Euclidean), focal excluded; 1,000 draws of one gene per up-DEG, no gene
# reused within a draw (next-nearest unused gene if a pool runs out); statistic = IFN-α genes in the matched set;
# p = (1 + #null ≥ observed) / 1,001; seed 0 (one generator per cell type). Diagnostics: % of null genes that are DEGs,
# % that are IFN-α, % of up-DEGs that are IFN-α.

# %%
adata = ad.read_h5ad(PROC / "m0_cells.h5ad")
assert counts.index.equals(adata.var_names)
p4_rows, p4_draws = [], {}
for ct in STACK_ORDER:
    s_ct = samples[samples["cell_type"] == ct]
    donors = s_ct.loc[s_ct["kept"], "donor"].unique()
    cells = ((adata.obs["cell_type"] == ct) & (adata.obs["condition"] == "ctrl") & adata.obs["donor"].isin(donors)).to_numpy()
    detection = pd.Series(np.asarray((adata.X[cells] > 0).mean(axis=0)).ravel(), index=adata.var_names)
    d_ct = de[de["cell_type"] == ct]
    feats = v.p4_features(counts, s_ct, d_ct.loc[d_ct["tested"], "gene"], detection)
    res, null = v.p4_null(feats, d_ct.loc[d_ct["de"] == "up", "gene"], d_ct.loc[d_ct["de"] != "none", "gene"], ifna,
                          k=30, n_draws=1000, seed=0)
    p4_rows.append({"cell_type": ct, "n_ctrl_cells": int(cells.sum()), "n_tested": len(feats), **res})
    p4_draws[ct] = null
p4 = pd.DataFrame(p4_rows)
display(p4.round(4))
save(p4, "p4_matched_null")

fig, axes = plt.subplots(2, 4, figsize=(15, 6))
for ax, ct in zip(axes.flat, STACK_ORDER):
    r = p4.set_index("cell_type").loc[ct]
    ax.hist(p4_draws[ct], bins=range(int(max(p4_draws[ct].max(), r["observed"])) + 2), color="#b0b0ad")
    ax.axvline(r["observed"], c=CT_COLORS[ct], lw=2)
    ax.set_title(f"{ct}\nobserved {r['observed']}, p = {r['p_empirical']:.4f}", fontsize=9)
    ax.set_xlabel("IFN-α genes in matched set")
fig.suptitle("P4: expression-matched null (grey, 1,000 draws) vs up-DEGs (line)")
plt.tight_layout(); fig.savefig(FIG / "m2b_p4_null.png", dpi=200, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## N1: label shuffles (35 / 10 / 3) next to the real DE count

# %%
n1 = v.n1_counts(runs[runs["run"].str.startswith("shuffle_")], fits, de)
save(n1, "n1_shuffles")
display(n1.pivot_table(index="cell_type", columns="kind", values="n_de", aggfunc=["min", "max", "size"]).loc[STACK_ORDER])

fig, ax = plt.subplots(figsize=(11, 4))
for i, ct in enumerate(STACK_ORDER):
    s = n1[(n1["cell_type"] == ct) & (n1["kind"] == "shuffle")]["n_de"]
    ax.scatter(i + np.random.default_rng(i).uniform(-0.15, 0.15, len(s)) - 0.2, s + 1, s=12, c="#b0b0ad",
               label="shuffle" if i == 0 else None)
    real = n1[(n1["cell_type"] == ct) & (n1["kind"] == "real")]["n_de"].iloc[0]
    ax.scatter(i + 0.2, real + 1, s=60, marker="D", c=CT_COLORS[ct], edgecolor="k", label="real labels (colour = cell type)" if i == 0 else None)
ax.set_yscale("log"); ax.set_xticks(range(8), STACK_ORDER, rotation=25, ha="right")
ax.set(ylabel="DE genes + 1 (log scale)", title="N1: DE genes in label shuffles vs the real labels")
leg = ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1)); leg.legend_handles[1].set_facecolor("white"); plt.tight_layout(); fig.savefig(FIG / "m2b_n1_shuffles.png", dpi=200, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## R1: leave-one-donor-out direction stability (O4: log2FC keeps its sign in every refit; significance not required)

# %%
r1, r1_calls = v.r1_lodo(de, runs[runs["run"].str.startswith("lodo_")])
r1 = r1.set_index("cell_type").loc[STACK_ORDER].reset_index()
display(r1)
save(r1, "r1_lodo")

# %% [markdown]
# ## R2 / R3: overlap of main DE calls with the unpaired (`~ condition`) and no-outlier runs
# Genes, either direction; `sign_flips` = DE in both runs with opposite directions.

# %%
r23 = pd.concat([v.overlap(de, runs[runs["run"] == name], name).rename(columns={f"n_{name}": "n_other"}).assign(other_run=name)
                 for name in ["unpaired", "no_outliers"]], ignore_index=True)
display(r23)
save(r23, "r2_r3_overlap")

fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), sharey=True)
for ax, col, lab in [(axes[0], "jaccard", "Jaccard"), (axes[1], "pct_main_recovered", "% of main DE calls recovered")]:
    w = r23.pivot(index="cell_type", columns="other_run", values=col).loc[STACK_ORDER]
    w.plot.barh(ax=ax, color={"unpaired": "#1b74ac", "no_outliers": "#dd1b26"}, legend=False)
    ax.set(xlabel=lab, ylabel="")
axes[1].legend(*axes[0].get_legend_handles_labels(), title="other run", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
axes[0].invert_yaxis(); plt.tight_layout(); fig.savefig(FIG / "m2b_r2_r3.png", dpi=200, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## V9 / V9b: Kang 2018 Supp Table 1
# V9 overlap is the M1b computation (same symbol mapping), restated. Correlations on two gene sets, leading with
# `paper_tested` = every paper DE gene we tested; `in_both` = DE in both lists. `median_abs_diff` is there because a
# correlation cannot see shrinkage toward 0. Decision: calls stay on unshrunk log2FC; V9b is recorded as favouring
# shrunken values.

# %%
paper = v.paper_supp1()
v9_overlap, v9_cor, v9_genes = v.v9(de, paper)
v9_overlap = v9_overlap.set_index("cell_type").loc[STACK_ORDER].reset_index()
v9_cor = v9_cor.sort_values(["gene_set", "cell_type"], ascending=[False, True], key=lambda c: c if c.name == "gene_set"
                            else c.map(STACK_ORDER.index)).reset_index(drop=True)
display(v9_overlap); display(v9_cor.round(3))
save(v9_overlap, "v9_overlap"); save(v9_cor, "v9_correlation")

fig, axes = plt.subplots(2, 8, figsize=(22, 6), sharex="col")
for j, ct in enumerate(STACK_ORDER):
    g = v9_genes[(v9_genes["cell_type"] == ct) & v9_genes["tested"]]
    for i, (col, lab) in enumerate([("log2FC", "ours, unshrunk"), ("log2FC_shrunk", "ours, shrunken")]):
        ax = axes[i, j]
        ax.scatter(g["paper_log2FC"], g[col], s=4, alpha=0.5, lw=0, c=np.where(g["in_both"], CT_COLORS[ct], "#b0b0ad"))
        lim = np.nanmax(np.abs(g[["paper_log2FC", col]].to_numpy()))
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.6)
        ax.set_title(ct if i == 0 else "", fontsize=9); ax.set_ylabel(lab if j == 0 else "")
        ax.set_xlabel("paper log2FC" if i == 1 else "")
fig.suptitle("V9b: paper log2FC vs ours, every paper DEG we tested (colour = DE in both; grey = paper only); dashed = identity")
plt.tight_layout(); fig.savefig(FIG / "m2b_v9b_scatter.png", dpi=150, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## F3b: Kang Fig 3b recreated
# log2(CPM + 1) per donor × cell type × condition, all 128 samples (incl. the 12 dropped by A1, as the paper's 8 donors
# per group) plus all PBMCs (sum over all cells per donor × condition); genes DE in ≥ 1 cell type; rows z-scored (ddof 0).
# Rows in the authors' file order (their 7 modules); our DE genes absent from it are at the bottom as "no module".
# The authors' file values are not log2 CPM (0–43,088, many ties), so this recreates the layout, not their values.

# %%
authors = v.authors_modules(paper)
display(authors["module"].value_counts().rename("genes").to_frame().T)
print("Label '3' has 401 genes = Supp Table 2 'Cluster III' (401 genes); other label-to-numeral matches are not established.")

de_genes = sorted(de.loc[de["de"] != "none", "gene"].unique())
z = v.f3b_matrix(counts, samples, de_genes, pseudocount=1, use_dropped_samples=True)
amod = authors.set_index("gene")
row_mod = amod["module"].reindex(z.index).where(z.index.isin(amod.index), "no module")
rows = amod["row_order"].reindex(z.index).where(z.index.isin(amod.index), len(amod)).sort_values(kind="stable").index
GROUPS = ["All PBMCs", "NK cells", "CD8 T cells", "CD4 T cells", "B cells", "Dendritic cells", "FCGR3A+ Monocytes",
          "CD14+ Monocytes", "Megakaryocytes"]  # authors' column order
cols = [f"{g}|{d}|{c}" for g in GROUPS for c in ["ctrl", "stim"] for d in DONOR_ORDER]
assert set(cols) == set(z.columns) and len(cols) == 144
zf = z.loc[rows, cols]
print(f"F3b matrix: {zf.shape[0]} genes × {zf.shape[1]} samples; genes per authors' module:")
display(row_mod.value_counts().rename("genes").to_frame().T)
save(zf.rename_axis("gene").reset_index().assign(authors_module=row_mod[rows].to_numpy()), "f3b_zscores")

fig, (axm, ax) = plt.subplots(1, 2, figsize=(16, 10), gridspec_kw={"width_ratios": [0.4, 20], "wspace": 0.01})
mods = sorted(row_mod.unique())
axm.imshow(row_mod[rows].map({m: i for i, m in enumerate(mods)}).to_numpy()[:, None], aspect="auto", cmap="tab10",
           vmin=0, vmax=9, interpolation="nearest")
axm.set(xticks=[], yticks=[], ylabel="authors' module (colour); 'no module' at the bottom")
im = ax.imshow(zf.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3, interpolation="nearest")
ax.set_xticks([i * 16 + 7.5 for i in range(9)], GROUPS, rotation=30, ha="right"); ax.set_yticks([])
for i in range(1, 9):
    ax.axvline(i * 16 - 0.5, c="k", lw=0.5)
ax.grid(False); axm.grid(False)
fig.colorbar(im, ax=ax, fraction=0.02, anchor=(0.0, 1.0), shrink=0.5, label="z-score of log2(CPM + 1)")
handles = [plt.matplotlib.patches.Patch(color=plt.cm.tab10(i), label=f"{m} ({(row_mod == m).sum():,} genes)")
           for i, m in enumerate(mods)]
ax.legend(handles=handles, title="authors' module label\n(3 = Supp Table 2 Cluster III)", frameon=False,
          loc="lower left", bbox_to_anchor=(1.01, 0.0), fontsize=9, title_fontsize=9)
ax.set_title("F3b: DE genes × samples (each group: 8 ctrl then 8 stim donors)")
fig.savefig(FIG / "m2b_f3b_heatmap.png", dpi=150, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## Modules: our k-means (k = 7, random_state 0, n_init 25) vs the authors' 7 clusters ("method inferred")
# k-means on all rows of the F3b matrix; the ARI uses genes that have an authors' module.

# %%
ours, xtab, ari = v.compare_modules(zf, authors, k=7, random_state=0, n_init=25)
xtab_all = pd.crosstab(ours, row_mod[ours.index])
display(xtab_all)
print(f"adjusted Rand index (genes in both): {ari:.3f}  (n = {int(xtab.to_numpy().sum())})")
save(xtab_all.reset_index(), "modules_overlap")
save(pd.DataFrame({"gene": ours.index, "our_module": ours.to_numpy(), "authors_module": row_mod[ours.index].to_numpy()}),
     "modules_genes")

# %% [markdown]
# ## Module biology: KEGG legacy (O5) enrichment of our modules
# One-sided Fisher per KEGG set, background = genes tested in ≥ 1 cell type, BH across sets within each module.
# KEGG symbols get the 6 B4 renames (e.g. RIGI → DDX58).
# KEGG legacy has no Influenza A, Measles, Herpes simplex or Hepatitis C sets (added to KEGG after the legacy freeze),
# so Supp Table 2's Cluster III can only be compared on terms present in both.

# %%
kegg = v.kegg_sets()
universe = set(de.loc[de["tested"], "gene"])
kegg_genes = set().union(*kegg.values())
print(f"KEGG legacy: {len(kegg)} sets, {len(kegg_genes)} genes; {len(kegg_genes & universe)} in the background "
      f"({len(universe)} genes); {len(kegg_genes & set(counts.index))} in the data at all")
kres = pd.concat([v.kegg_enrichment(ours.index[ours == m], universe, kegg).assign(our_module=m) for m in range(7)],
                 ignore_index=True)
save(kres, "kegg_modules")
display(kres[kres["q_bh"] < 0.05].groupby("our_module").head(5)[["our_module", "set", "overlap", "odds_ratio", "p", "q_bh"]])

for m in xtab.index[xtab["3"] >= 0.25 * xtab["3"].sum()]:  # our modules holding >= 25% of the label-3 genes
    print(f"our module {m}: {xtab.loc[m, '3']} of {int(xtab['3'].sum())} label-3 genes in our DE list; "
          f"module size {int((ours == m).sum())} ({100 * xtab.loc[m, '3'] / (ours == m).sum():.0f}% label 3)")
    top = kres[kres["our_module"] == m].reset_index(drop=True)
    top.insert(0, "rank", top.index + 1)
    display(top.head(8)[["rank", "set", "overlap", "odds_ratio", "p", "q_bh"]])
supp = pd.read_excel(ROOT / "data/kang2018_supp/MOESM5.xlsx", sheet_name=None, header=None)
supp3 = next(d for name, d in supp.items() if name.startswith("Cluster III"))
print("Supp Table 2 Cluster III terms (DAVID):", supp3[supp3.eq("KEGG_PATHWAY").any(axis=1)].apply(
    lambda r: r[r.eq("KEGG_PATHWAY").idxmax() + 1], axis=1).tolist())

# %% [markdown]
# ## I5: ribosomal / MT genes among the top 50 shared DE genes
# Shared = `shared` and `shared*` in M2a's `m2_table.csv`; ranked by mean |log2FC| over the 8 cell types (decision).

# %%
m2_table = pd.read_csv(PROC / "m2_table.csv")
shared_genes = m2_table.loc[m2_table["shared_response"].isin(["shared", "shared*"]), "gene"]
lfc = de[de["gene"].isin(shared_genes)].pivot(index="gene", columns="cell_type", values="log2FC")
assert lfc.shape == (len(shared_genes), 8) and lfc.notna().all().all(), "every shared gene needs a log2FC in all 8 cell types"
i5_rank = lfc.abs().mean(axis=1)
i5 = v.i5_ribo(m2_table, i5_rank, n=50)
print(i5)
top50 = i5_rank.sort_values(ascending=False).head(50)
print("top 50 shared genes by mean |log2FC|:", ", ".join(f"{g} ({x:.1f})" for g, x in top50.items()))
save(pd.DataFrame([i5]), "i5_ribo")

# %% [markdown]
# ## I6: DE effect sizes vs n_cells / total_counts
# Per donor × cell type (kept samples): headline = median over that cell type's DE genes of
# sign(call) × log2((stim CPM + 1) / (ctrl CPM + 1)) (direction-aligned, so up and down genes don't cancel), vs log of
# the mean of ctrl and stim `n_cells` and `total_counts`; Spearman per cell type. The unsigned median is shown alongside.

# %%
i6 = v.i6_effects(counts, samples, de)
save(i6, "i6_effects")
display(i6.groupby("cell_type")[["n_de", "median_log2_ratio_aligned", "median_log2_ratio"]].agg(["min", "max"]).loc[STACK_ORDER].round(2))
i6_cor = pd.DataFrame([{"cell_type": ct, "y": y, "x": x, "n_donors": len(g), "spearman": stats.spearmanr(g[x], g[y]).statistic}
                       for ct, g in i6.groupby("cell_type") for y in ["median_log2_ratio_aligned", "median_log2_ratio"]
                       for x in ["log_n_cells", "log_total_counts"]])
display(i6_cor.pivot_table(index="cell_type", columns=["y", "x"], values="spearman").loc[STACK_ORDER].round(2))
save(i6_cor, "i6_spearman")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, x in zip(axes, ["log_n_cells", "log_total_counts"]):
    for ct, g in i6.groupby("cell_type"):
        ax.plot(g[x], g["median_log2_ratio_aligned"], "o", c=CT_COLORS[ct], label=ct, ms=5)
    ax.set(xlabel=x + " (mean of ctrl, stim)")
axes[0].set_ylabel("median direction-aligned\nlog2((stim+1)/(ctrl+1) CPM), DE genes"); axes[1].legend(fontsize=7, frameon=False, bbox_to_anchor=(1, 1))
fig.suptitle("I6: per-donor effect size vs sample size (each dot a donor × cell type)")
plt.tight_layout(); fig.savefig(FIG / "m2b_i6.png", dpi=200, bbox_inches="tight"); plt.show()
