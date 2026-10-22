# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M0: EDA of Kang 2018 (IFN-β vs control PBMCs)
#
# Contract: `audit/contracts/M0_eda.md`. Logic lives in `src/m0_qc.py` (tested); this notebook shows the data.
# Nothing is dropped here except the agreed basic filters (cells ≥ 200 genes, genes ≥ 3 cells). Everything else is **flagged** for the user.

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import warnings
warnings.filterwarnings("ignore", message=".*OldFormatWarning.*")

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np


def neutral_lines(parts):
    for k in ("cmins", "cmaxes", "cbars", "cmedians"):
        if k in parts:
            parts[k].set_color("#5c5c57"); parts[k].set_linewidth(1)
import pandas as pd
from IPython.display import display

from src import m0_qc

pd.set_option("display.width", 160, "display.max_columns", 20)

# Input files, fingerprints and cell-type names are defined once, in src/m0_qc.py
H5AD, GEO = m0_qc.KANG_H5AD, m0_qc.KANG_GEO
H5AD_SHA, GEO_SHA = m0_qc.KANG_H5AD_SHA256, m0_qc.KANG_GEO_SHA256
CELL_TYPES = m0_qc.KANG_CELL_TYPES
# Thresholds (all agreed; passed explicitly everywhere, and drawn on the plots from these values)
MIN_GENES, MIN_CELLS, N_MADS, MIN_CELLS_GROUP, MIN_PAIRED_DONORS = 200, 3, 5.0, 10, 3
KEY_GENES = ["ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3", "RSAD2", "IRF7", "STAT1"]

# Kang 2018 Fig 3 colours, donor order and stacking order: defined once, in src/figures.py
from src.figures import CT_COLORS, COND_COLORS, DONOR_ORDER, STACK_ORDER
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e3", "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "axes.edgecolor": "#8a8a85", "text.color": "#2b2b28", "axes.labelcolor": "#2b2b28",
                     "xtick.color": "#5c5c57", "ytick.color": "#5c5c57", "figure.dpi": 110})

# %% [markdown]
# ## 1. The raw inputs: what's in each table

# %%
raw, geo = m0_qc.load_inputs(H5AD, GEO, expected_sha256=H5AD_SHA, expected_geo_sha256=GEO_SHA)
print(raw)
print("\nh5ad obs: head"); display(raw.obs.head())
print("h5ad obs: tail"); display(raw.obs.tail())
print("h5ad var: head / tail"); display(raw.var.head()); display(raw.var.tail())
print("GEO metadata: head"); display(geo.head())
print("GEO metadata: tail"); display(geo.tail())


# %%
def column_summary(df):
    return pd.DataFrame({
        "dtype": df.dtypes.astype(str),
        "n_NA": df.isna().sum(),
        "n_unique": df.nunique(),
        "example": df.iloc[0].astype(str),
    })

print("h5ad obs columns"); display(column_summary(raw.obs))
print("GEO columns"); display(column_summary(geo))

counts = pd.Series({
    "h5ad cells": raw.n_obs, "h5ad unique barcodes": raw.obs_names.nunique(),
    "h5ad genes": raw.n_vars, "h5ad duplicated gene names": int(raw.var_names.duplicated().sum()),
    "GEO rows (droplets)": len(geo), "GEO unique barcodes": geo.index.nunique(),
    "donors (h5ad replicate)": raw.obs["replicate"].nunique(),
})
display(counts.to_frame("value"))

# %% [markdown]
# ### Every metadata column, with a plot suited to it
# Categorical → bar of counts; numeric → histogram. Columns from the later re-analysis (`nCount_SCT`, `seurat_clusters`, …) are shown but **not used**.

# %%
cols = list(raw.obs.columns)
fig, axes = plt.subplots(int(np.ceil(len(cols) / 4)), 4, figsize=(16, 3 * int(np.ceil(len(cols) / 4))))
for ax, c in zip(axes.ravel(), cols):
    s = raw.obs[c]
    if s.dtype.name in ("category", "object") or s.nunique() <= 15:
        vc = s.astype(str).value_counts().sort_index()
        ax.bar(range(len(vc)), vc.values, color="#2a78d6", width=0.8)
        ax.set_xticks(range(len(vc)), vc.index, rotation=90, fontsize=7)
    else:
        ax.hist(s, bins=60, color="#2a78d6")
    ax.set_title(c, fontsize=9)
for ax in axes.ravel()[len(cols):]:
    ax.axis("off")
fig.tight_layout()

# %% [markdown]
# Same for every column of the authors' GEO table (all 29,065 droplets, including doublets). Most duplicate h5ad columns and are checked cell by cell in the join (§3).

# %%
gcols = list(geo.columns)
fig, axes = plt.subplots(int(np.ceil(len(gcols) / 4)), 4, figsize=(16, 3 * int(np.ceil(len(gcols) / 4))))
for ax, c in zip(axes.ravel(), gcols):
    s = geo[c]
    if s.dtype.name in ("category", "object") or s.nunique() <= 15:
        vc = s.astype(str).value_counts().sort_index()
        ax.bar(range(len(vc)), vc.values, color="#2a78d6", width=0.8)
        ax.set_xticks(range(len(vc)), vc.index, rotation=90, fontsize=7)
    else:
        ax.hist(s, bins=60, color="#2a78d6")
    ax.set_title(f"GEO: {c}", fontsize=9)
for ax in axes.ravel()[len(gcols):]:
    ax.axis("off")
fig.tight_layout()

# %% [markdown]
# ## 2. The counts matrix

# %%
v = raw.X.data
print("X type:", type(raw.X).__name__, raw.X.dtype, "| non-negative integers:", bool((v >= 0).all() and np.array_equal(v, np.round(v))))
print("Max count:", v.max(), "| fraction of zero entries: %.3f" % (1 - raw.X.nnz / np.prod(raw.shape)))
print("MT- genes:", int(raw.var_names.str.startswith("MT-").sum()), "| RPS/RPL-prefix genes:", int(raw.var_names.str.match(r"^RP[SL]").sum()),
      "| counted as ribosomal (minus 14 non-ribosomal, D4):", int(m0_qc.ribosomal_mask(raw.var_names).sum()))
display(pd.DataFrame({"present": [g in raw.var_names for g in KEY_GENES]}, index=KEY_GENES).T)

# %% [markdown]
# ## 3. Join the authors' labels (key = barcode sequence + condition)
# 234 barcode sequences appear in both conditions among these singlets; a raw-barcode join would silently label only the control cells.

# %%
joined = m0_qc.join_author_metadata(raw, geo)
core = raw.obs_names.str.replace(r"-\d+$", "", regex=True)
n_both = int((pd.Series(raw.obs["label"].astype(str).values, index=core).groupby(level=0).nunique() == 2).sum())
print(f"matched {joined.n_obs} of {raw.n_obs} cells | sequences in both conditions: {n_both}")
print("control / treatment:", m0_qc.infer_conditions(joined.obs["condition"]))
agree = (raw.obs["cell_type"].astype(str).values == joined.obs["cell_type"].values).mean()
print(f"re-analysis labels agree with the authors' on {agree:.1%} of cells")
display(pd.crosstab(joined.obs["cell_type"], raw.obs["cell_type"].values, rownames=["authors"], colnames=["re-analysis"], margins=True))
display(joined.obs.head()); display(joined.obs.tail())

# %% [markdown]
# ## 4. QC metrics and the agreed basic filters
# Cutoffs are drawn on the plots: cells need ≥ 200 genes; genes need ≥ 3 cells.

# %%
qc = m0_qc.compute_qc_metrics(joined)
filtered = m0_qc.basic_filter(qc, min_genes=MIN_GENES, min_cells=MIN_CELLS)
print(f"cells {qc.n_obs} → {filtered.n_obs} (removed {qc.n_obs - filtered.n_obs}); genes {qc.n_vars} → {filtered.n_vars} (removed {qc.n_vars - filtered.n_vars})")

fig, axes = plt.subplots(1, 3, figsize=(15, 3.2))
axes[0].hist(qc.obs["n_genes"], bins=100, color="#2a78d6"); axes[0].axvline(MIN_GENES, color="#2b2b28", lw=1, ls="--")
axes[0].set(title="Genes per cell (before filters)", xlabel="genes detected", ylabel="cells")
axes[1].hist(np.log10(qc.var["n_cells"] + 1), bins=100, color="#2a78d6"); axes[1].axvline(np.log10(MIN_CELLS + 1), color="#2b2b28", lw=1, ls="--")
axes[1].set(title="Cells per gene (before filters)", xlabel="log10(cells expressing + 1)", ylabel="genes")
axes[2].hist(qc.obs["n_counts"], bins=100, color="#2a78d6")
axes[2].set(title="Counts per cell (before filters)", xlabel="total counts", ylabel="cells")
fig.tight_layout()
print("lowest genes-per-cell values:", np.sort(qc.obs["n_genes"].values)[:10])

# %% [markdown]
# ## 5. Run the tested pipeline and write the outputs
# `run_m0_kang` repeats the steps above end to end on the real files (fingerprints and cell-type names always checked), validates everything **before** writing, and logs cells in → out.

# %%
flow = m0_qc.run_m0_kang(ROOT / "data/processed", audit_dir=ROOT / "audit", min_genes=MIN_GENES, min_cells=MIN_CELLS,
                         n_mads=N_MADS, min_cells_group=MIN_CELLS_GROUP, min_paired_donors=MIN_PAIRED_DONORS)
display(pd.DataFrame(flow))
m0 = ad.read_h5ad(ROOT / "data/processed/m0_cells.h5ad")
obs = m0.obs
assert m0.n_obs == filtered.n_obs and m0.n_vars == filtered.n_vars, "notebook filtering and run_m0 disagree"
print(m0)
print("uns:", {k: m0.uns[k] for k in ("control", "treatment")})
display(obs.head()); display(obs.tail())
display(m0.var.head()); display(m0.var.tail())

# %% [markdown]
# ## 6. Cells per donor, and genes detected per donor

# %%
donors = [d for d in DONOR_ORDER if d in set(obs["donor"])]  # paper order (Kang Fig 3d)
assert len(donors) == obs["donor"].nunique()
fig, axes = plt.subplots(1, 2, figsize=(15, 4))
by_cond = pd.crosstab(obs["donor"], obs["condition"]).loc[donors]
bottom = np.zeros(len(donors))
for c in ["ctrl", "stim"]:
    axes[0].bar(donors, by_cond[c], bottom=bottom, color=COND_COLORS[c], label=c, width=0.7, edgecolor="white", linewidth=1)
    bottom += by_cond[c].values
axes[0].set(title="Cells per donor, by condition", ylabel="cells"); axes[0].legend(frameon=False)
by_ct = pd.crosstab(obs["donor"], obs["cell_type"]).loc[donors, STACK_ORDER]
bottom = np.zeros(len(donors))
for ct in STACK_ORDER:
    axes[1].bar(donors, by_ct[ct], bottom=bottom, color=CT_COLORS[ct], label=ct, width=0.7, edgecolor="white", linewidth=1)
    bottom += by_ct[ct].values
axes[1].set(title="Cells per donor, by cell type", ylabel="cells")
h, l = axes[1].get_legend_handles_labels(); axes[1].legend(h[::-1], l[::-1], frameon=False, fontsize=7, bbox_to_anchor=(1, 1))
for ax in axes:
    ax.tick_params(axis="x", rotation=45)
fig.tight_layout()
display(by_cond.assign(total=by_cond.sum(1)))

# %%
fig, ax = plt.subplots(figsize=(15, 4))
pos, data, colors = [], [], []
for i, d in enumerate(donors):
    for j, c in enumerate(["ctrl", "stim"]):
        pos.append(i * 3 + j); data.append(obs.loc[(obs.donor == d) & (obs.condition == c), "n_genes"].values); colors.append(COND_COLORS[c])
parts = ax.violinplot(data, positions=pos, showmedians=True, widths=0.9)
neutral_lines(parts)
for body, col in zip(parts["bodies"], colors):
    body.set_facecolor(col); body.set_alpha(0.6)
ax.set_xticks([i * 3 + 0.5 for i in range(len(donors))], donors, rotation=45)
ax.set(title="Genes detected per cell, per donor (blue = ctrl, red = stim)", ylabel="genes detected")
display(obs.groupby(["donor", "condition"], observed=True)["n_genes"].median().unstack().round(0))

# %% [markdown]
# ## 7. Per-cell QC distributions, with outliers flagged (> 5 raw MADs above or below, within donor × condition × cell type; groups < 10 cells not assessed)

# %%
fig, axes = plt.subplots(1, 3, figsize=(16, 4))
for ax, (metric, label) in zip(axes, [("n_counts", "total counts (log10)"), ("n_genes", "genes detected (log10)"), ("pct_ribo", "% ribosomal")]):
    for k, c in enumerate(["ctrl", "stim"]):
        vals = obs.loc[obs.condition == c, metric]
        vals = np.log10(vals) if metric != "pct_ribo" else vals
        part = ax.violinplot(vals, positions=[k], widths=0.8, showmedians=True)
        neutral_lines(part)
        part["bodies"][0].set_facecolor(COND_COLORS[c]); part["bodies"][0].set_alpha(0.6)
        out = obs.loc[(obs.condition == c) & obs.flag_outlier, metric]
        ax.scatter(np.full(len(out), k) + 0.25, np.log10(out) if metric != "pct_ribo" else out, s=10, color="#2b2b28", label="flagged outlier" if k == 0 else None)
    ax.set_xticks([0, 1], ["ctrl", "stim"]); ax.set_title(label)
axes[0].legend(frameon=False)
fig.tight_layout()
display(obs.groupby("condition", observed=True)[["n_counts", "n_genes", "pct_ribo"]].median().round(1))
display(obs.groupby(["donor", "condition"], observed=True)["flag_outlier"].sum().unstack().rename(columns=lambda c: f"outliers_{c}"))
med = obs.groupby(["donor", "condition", "cell_type"], observed=True)["n_counts"].transform("median")
flagged = obs[obs.flag_outlier]
outlier_summary = pd.DataFrame({
    "flagged": flagged.groupby("cell_type", observed=True).size(),
    "cells": obs.groupby("cell_type", observed=True).size(),
}).reindex(CELL_TYPES).fillna(0).astype(int).assign(pct=lambda d: (100 * d.flagged / d.cells).round(1))
print(f"flagged {len(flagged)}: {int((flagged.n_counts > med[flagged.index]).sum())} above their group's median counts, "
      f"{int((flagged.n_counts <= med[flagged.index]).sum())} at or below")
display(outlier_summary)

# %% [markdown]
# ### Outliers within each donor × condition × cell type group (how they are actually computed)
# One panel per cell type; one column per donor × condition. Each column shows its group's **median** (solid) and **± 5 MAD cutoffs** (dashed).
# A cell is flagged if it's outside the cutoffs on log counts (first figure) **or** log genes (second figure).
# **Filled dark** = flagged on the metric shown (outside the dashed lines). **Hollow** = flagged only on the *other* metric, so it can sit inside the lines.
# Groups with < 10 cells are not assessed (grey, no cutoffs).

# %%
# Per-metric flags, recomputed with the pipeline's own function; their union must equal the pipeline's flag_outlier.
per_metric = {m: pd.Series(False, index=obs.index) for m in ("n_counts", "n_genes")}
for _, g in obs.groupby(["donor", "condition", "cell_type"], observed=True):
    if len(g) >= MIN_CELLS_GROUP:
        for m in per_metric:
            per_metric[m].loc[g.index] = m0_qc.mad_outlier(np.log1p(g[m].to_numpy(float)), N_MADS)
assert ((per_metric["n_counts"] | per_metric["n_genes"]) == obs["flag_outlier"]).all(), "per-metric flags don't reproduce flag_outlier"
print("flagged on counts only:", int((per_metric["n_counts"] & ~per_metric["n_genes"]).sum()),
      "| genes only:", int((per_metric["n_genes"] & ~per_metric["n_counts"]).sum()),
      "| both:", int((per_metric["n_counts"] & per_metric["n_genes"]).sum()))


def plot_group_outliers(metric, label):
    other = "n_genes" if metric == "n_counts" else "n_counts"
    fig, axes = plt.subplots(2, 4, figsize=(18, 7), sharey=True)
    rng = np.random.default_rng(0)
    cols = [(d, c) for d in donors for c in ["ctrl", "stim"]]
    for ax, ct in zip(axes.ravel(), CELL_TYPES):
        for x, (d, c) in enumerate(cols):
            g = obs[(obs.cell_type == ct) & (obs.donor == d) & (obs.condition == c)]
            if g.empty:
                continue
            v = np.log1p(g[metric].to_numpy(float))
            jx = x + rng.uniform(-0.3, 0.3, len(v))
            this = per_metric[metric].loc[g.index].to_numpy()
            only_other = per_metric[other].loc[g.index].to_numpy() & ~this
            normal = ~this & ~only_other
            assessed = len(g) >= MIN_CELLS_GROUP
            ax.scatter(jx[normal], v[normal], s=3, color="#86b6ef" if assessed else "#c3c2b7", rasterized=True)
            ax.scatter(jx[only_other], v[only_other], s=12, facecolors="none", edgecolors="#2b2b28", linewidths=0.8)
            ax.scatter(jx[this], v[this], s=12, color="#2b2b28")
            if assessed:
                med = np.median(v); mad = np.median(np.abs(v - med))
                ax.hlines(med, x - 0.4, x + 0.4, color="#2a78d6", lw=1.5)
                ax.hlines([med - N_MADS * mad, med + N_MADS * mad], x - 0.4, x + 0.4, color="#5c5c57", lw=1, ls="--")
        ax.set_title(ct, fontsize=10)
        ax.set_xticks(range(len(cols)), [f"{d.replace('patient_', '')} {c}" for d, c in cols], rotation=90, fontsize=6)
    for ax in axes[:, 0]:
        ax.set_ylabel(label)
    handles = [plt.Line2D([], [], marker="o", ls="", color="#86b6ef", label="not flagged"),
               plt.Line2D([], [], marker="o", ls="", color="#2b2b28", label="flagged on this metric"),
               plt.Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor="#2b2b28", label="flagged on the other metric only"),
               plt.Line2D([], [], marker="o", ls="", color="#c3c2b7", label="group < 10 cells (not assessed)")]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle(f"{label}: each cell vs its own donor × condition × cell type group", y=1.07)
    fig.tight_layout()

plot_group_outliers("n_counts", "log(1 + total counts)")
plot_group_outliers("n_genes", "log(1 + genes detected)")
print(obs.groupby(["cell_type", "donor", "condition"], observed=True)["flag_outlier"].sum()
      .unstack(["donor", "condition"]).reindex(CELL_TYPES).fillna(0).astype(int).to_string())

# %% [markdown]
# ## 8. Design checks: donor × condition × cell type counts
# Groups with < 10 cells are too small for a reliable pseudobulk sample (flagged). This table is the input for decision **D3**.

# %%
table = obs.groupby(["cell_type", "donor", "condition"], observed=True).size().unstack(["donor", "condition"]).reindex(CELL_TYPES).fillna(0).astype(int)
display(table.style.apply(lambda s: ["background-color: #fde2e1" if x < 10 else "" for x in s], axis=1))
pairs = pd.DataFrame({ct: {"donors with both conditions ≥ 10 cells": int(((table.loc[ct].unstack() >= 10).all(axis=1)).sum()),
                           "flag_underpowered_celltype": bool(obs.loc[obs.cell_type == ct, "flag_underpowered_celltype"].iloc[0])}
                      for ct in CELL_TYPES}).T
display(pairs)

# %% [markdown]
# ## 9. Weak donors, composition shifts (Kang Fig 3d), covariate trends: numbers and plots, no automatic flags

# %%
tot = obs.groupby("donor", observed=True).size()
display(pd.DataFrame({"cells": tot, "fraction of median donor": (tot / tot.median()).round(2)}))

prop = (obs.groupby(["donor", "condition", "cell_type"], observed=True).size()
        .groupby(level=["donor", "condition"], observed=True).transform(lambda s: s / s.sum()).unstack("cell_type")[CELL_TYPES])
fig, axes = plt.subplots(1, 2, figsize=(15, 4), sharey=True)
for ax, c in zip(axes, ["ctrl", "stim"]):
    p = prop.xs(c, level="condition").loc[donors]
    bottom = np.zeros(len(donors))
    for ct in STACK_ORDER:
        ax.bar(donors, p[ct], bottom=bottom, color=CT_COLORS[ct], label=ct, width=0.7, edgecolor="white", linewidth=1)
        bottom += p[ct].values
    ax.set_title(f"Cell-type proportions per donor: {c} (recreates Kang Fig 3d)"); ax.tick_params(axis="x", rotation=45)
h, l = axes[1].get_legend_handles_labels(); axes[1].legend(h[::-1], l[::-1], frameon=False, fontsize=7, bbox_to_anchor=(1, 1))
fig.tight_layout()
shift = (prop.xs("stim", level="condition") - prop.xs("ctrl", level="condition")).loc[donors]
print("Proportion stim − ctrl, per donor (percentage points)"); display((100 * shift).round(1))

cov = obs.groupby(["cell_type", "condition"], observed=True)[["n_counts", "n_genes", "pct_ribo"]].median().unstack("condition")
ratio = pd.DataFrame({m: (cov[(m, "stim")] / cov[(m, "ctrl")]).round(2) for m in ["n_counts", "n_genes", "pct_ribo"]})
print("Median stim / ctrl ratio of each technical covariate, per cell type"); display(ratio.reindex(CELL_TYPES))

# Covariates by donor (D5): medians per donor × condition, and each donor's median relative to all donors
by_donor = obs.groupby(["donor", "condition"], observed=True)[["n_counts", "n_genes", "pct_ribo"]].median().round(1)
rel_donor = (obs.groupby("donor", observed=True)[["n_counts", "n_genes", "pct_ribo"]].median()
             / obs[["n_counts", "n_genes", "pct_ribo"]].median()).round(2)
print("Median covariates per donor × condition"); display(by_donor)
print("Each donor's median / overall median"); display(rel_donor)
fig, axes = plt.subplots(1, 3, figsize=(16, 3.5))
for ax, (metric, label) in zip(axes, [("n_counts", "total counts (log10)"), ("n_genes", "genes detected (log10)"), ("pct_ribo", "% ribosomal")]):
    for i, d in enumerate(donors):
        vals = obs.loc[obs.donor == d, metric]
        part = ax.violinplot(np.log10(vals) if metric != "pct_ribo" else vals, positions=[i], widths=0.8, showmedians=True)
        neutral_lines(part); part["bodies"][0].set_facecolor("#2a78d6"); part["bodies"][0].set_alpha(0.6)
    ax.set_xticks(range(len(donors)), donors, rotation=45, fontsize=8); ax.set_title(f"{label} by donor")
fig.tight_layout()

# %% [markdown]
# ## 10. Reading the data correctly: authors' t-SNE (recreates Kang Fig 3a)

# %%
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
order = np.random.default_rng(0).permutation(len(obs))  # random draw order so neither condition hides the other
o = obs.iloc[order]
axes[0].scatter(o.tsne1, o.tsne2, s=1, c=o.condition.map(COND_COLORS).values, rasterized=True)
axes[0].set_title("Authors' t-SNE: ctrl (blue) vs stim (red)")
axes[1].scatter(o.tsne1, o.tsne2, s=1, c=o.cell_type.map(CT_COLORS).values, rasterized=True)
axes[1].legend(handles=[plt.Line2D([], [], marker="o", ls="", color=CT_COLORS[ct], label=ct) for ct in CELL_TYPES],
               frameon=False, fontsize=8, bbox_to_anchor=(1, 1))
axes[0].legend(handles=[plt.Line2D([], [], marker="o", ls="", color=COND_COLORS[c], label=c) for c in COND_COLORS], frameon=False)
axes[1].set_title("Authors' t-SNE: authors' cell types")
for ax in axes:
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
fig.tight_layout()

# %% [markdown]
# ## 11. Write `audit/data_facts.md` and the numbers for `audit/qc_flags.md`

# %%
def md(df):
    return df.to_markdown()

facts = f"""# Data facts (written by notebooks/00_eda.ipynb)

Read this instead of guessing. Source files and fingerprints are in `CLAUDE.md`.

- h5ad: {raw.n_obs} cells × {raw.n_vars} genes; `X` raw integer counts (sparse, float32 storage); no `.raw`, no layers.
- Gene names: symbols in `var_names`, unique. **No `MT-` genes** (so no % mito). {int(raw.var_names.str.match(r"^RP[SL]").sum())} genes with the RPS/RPL prefix; {int(m0_qc.ribosomal_mask(raw.var_names).sum())} counted as ribosomal after excluding 14 non-ribosomal genes (decision D4).
- Key ISGs present: {", ".join(g for g in KEY_GENES if g in raw.var_names)}. Missing: {", ".join(g for g in KEY_GENES if g not in raw.var_names) or "none"}.
- Conditions: `ctrl` (control) / `stim` (treatment); stored as `uns['control']` / `uns['treatment']` in the M0 output.
- Donors: {", ".join(donors)} (column `donor`, from h5ad `replicate`); GEO `ind` holds the bare number.
- Cell types (authors', GEO `cell`): {", ".join(CELL_TYPES)}.
- Join: key (barcode without suffix, condition); {joined.n_obs}/{raw.n_obs} matched; {n_both} sequences occur in both conditions among these singlets (313 across all GEO droplets).
- The h5ad's own `cell_type` (re-analysis) agrees with the authors' on {agree:.1%} of cells; not used.
- Basic filters: cells {qc.n_obs} → {filtered.n_obs}; genes {qc.n_vars} → {filtered.n_vars}.
- M0 output: `data/processed/m0_cells.h5ad`, {m0.n_obs} cells × {m0.n_vars} genes.

## Cells per donor × condition

{md(by_cond.assign(total=by_cond.sum(1)))}

## Cells per cell type × donor × condition

{md(table)}
"""
(ROOT / "audit/data_facts.md").write_text(facts)

small = obs[obs.flag_small_group].groupby(["cell_type", "donor", "condition"], observed=True).size()
flag_numbers = f"""# QC flag numbers (written by notebooks/00_eda.ipynb; proposed actions are added at review)

## Outlier cells (> 5 raw MADs above or below, on log counts or log genes, within donor × condition × cell type; groups < 10 cells not assessed)

{md(obs.groupby(["donor", "condition"], observed=True)["flag_outlier"].sum().unstack())}

By cell type:

{md(outlier_summary)}

## Groups with < 10 cells (donor × condition × cell type)

{md(small.to_frame("cells")) if len(small) else "none"}

## Cell types by donors with both conditions ≥ 10 cells

{md(pairs)}

## Donor sizes

{md(pd.DataFrame({"cells": tot, "fraction of median donor": (tot / tot.median()).round(2)}))}

## Composition shift, stim − ctrl (percentage points)

{md((100 * shift).round(1))}

## Technical covariates: median stim / ctrl ratio per cell type

{md(ratio.reindex(CELL_TYPES))}

## Technical covariates by donor: each donor's median / overall median

{md(rel_donor)}
"""
(ROOT / "audit/qc_flag_numbers.md").write_text(flag_numbers)
print("wrote audit/data_facts.md, audit/qc_flag_numbers.md, audit/cell_flow.md")
