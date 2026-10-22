# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M1a: pseudobulk samples and variance partitioning
#
# Contract: `audit/contracts/M1a_pseudobulk.md`. Pseudobulk logic: `src/m1a_pseudobulk.py` (tested); variance partitioning: `R/m1a_varpart.R`.
# This notebook only reads their output files and shows them. Nothing here drops a sample; outlier samples are shown with distances, **no automatic flag** (decision, Wave 1).

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

pd.set_option("display.width", 160, "display.max_columns", 20)
PROC = ROOT / "data/processed"
FILES = ["m1a_counts.csv", "m1a_samples.csv", "m1a_counts_no_outliers.csv", "m1a_samples_no_outliers.csv",
         "m1a_varpart.csv", "m1a_varpart_vst.csv", "m1a_varpart_fits.csv", "m1a_cancor.csv"]
MIN_BASEMEAN = 1  # the R script was run with this value (DE gene filter, baseMean > 1)
TERMS = ["donor", "condition", "log_n_cells", "log_total_counts", "Residuals"]
# Condition labels come from the M0 output (uns), never hardcoded
_m0 = ad.read_h5ad(ROOT / "data/processed/m0_cells.h5ad", backed="r")
CONTROL, TREATMENT = str(_m0.uns["control"]), str(_m0.uns["treatment"])
_m0.file.close()
CONDS = [CONTROL, TREATMENT]
COND_COLORS = {CONTROL: "#1b74ac", TREATMENT: "#dd1b26"}
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e3", "grid.linewidth": 0.6, "axes.axisbelow": True, "figure.dpi": 110})

samples = pd.read_csv(PROC / "m1a_samples.csv")
varpart = pd.read_csv(PROC / "m1a_varpart.csv")
vst = pd.read_csv(PROC / "m1a_varpart_vst.csv")
fits = pd.read_csv(PROC / "m1a_varpart_fits.csv")
CELL_TYPES = sorted(samples["cell_type"].unique())

# %% [markdown]
# ## 1. Every output file: shape, head and tail

# %%
for f in FILES:
    df = pd.read_csv(PROC / f, index_col=0 if "counts" in f else None)
    print(f"{f}: {df.shape[0]} rows × {df.shape[1]} columns")
    display(df.head()); display(df.tail())

# %% [markdown]
# ## 2. Sample table
# Cells per sample (rows: cell type; columns: donor × condition). A1: a donor is dropped from a cell type (both conditions) if either condition has < 10 cells. Dropped samples are listed below.

# %%
display(samples.pivot_table(index="cell_type", columns=["donor", "condition"], values="n_cells", aggfunc="sum"))
print("Dropped samples (kept = False):")
display(samples[~samples["kept"]])
summary = samples.groupby("cell_type").agg(samples=("sample", "size"), kept=("kept", "sum"),
                                           paired_donors=("donor", lambda d: d[samples.loc[d.index, "kept"]].nunique()),
                                           cells_kept=("n_cells", lambda n: n[samples.loc[n.index, "kept"]].sum()))
display(summary)
print("Total samples:", len(samples), "| kept:", int(samples["kept"].sum()),
      "| n_cells sum:", int(samples["n_cells"].sum()), "| in kept samples:", int(samples.loc[samples["kept"], "n_cells"].sum()))
clean = pd.read_csv(PROC / "m1a_samples_no_outliers.csv")
print("No-outlier variant: samples", len(clean), "| kept", int(clean["kept"].sum()), "| n_cells sum", int(clean["n_cells"].sum()),
      "| samples whose kept status differs from the main run:",
      int((clean.set_index("sample")["kept"] != samples.set_index("sample")["kept"]).sum()))

# %% [markdown]
# ## 3. Technical covariates vs condition (kept samples)
# `n_cells` and `total_counts` per sample, control vs stimulated, one line per donor. These are the two covariates in the variance-partitioning model. (10x run = condition, so it is a caveat, not a model term; % ribo stays out, F6.)

# %%
kept = samples[samples["kept"]]
fig, axes = plt.subplots(2, len(CELL_TYPES), figsize=(2.3 * len(CELL_TYPES), 6), sharey="row")
for row, cov in enumerate(["n_cells", "total_counts"]):
    for ax, ct in zip(axes[row], CELL_TYPES):
        wide = kept[kept["cell_type"] == ct].pivot(index="donor", columns="condition", values=cov)[CONDS]
        for _, r in wide.iterrows():
            ax.plot([0, 1], r.values, color="#8a8a85", lw=0.8)
        for i, c in enumerate(CONDS):
            ax.scatter([i] * len(wide), wide[c], color=COND_COLORS[c], s=18, zorder=3)
        ax.set_yscale("log"); ax.set_xticks([0, 1], CONDS); ax.set_xlim(-0.4, 1.4)
        ax.set_title(ct, fontsize=8)
        log_ratio = np.log2(wide[TREATMENT] / wide[CONTROL])
        ax.text(0.5, 0.02, f"median log2 {TREATMENT}/{CONTROL} {log_ratio.median():+.2f}", transform=ax.transAxes, ha="center", fontsize=6)
    axes[row, 0].set_ylabel(cov + " (log scale)")
fig.tight_layout()

# %% [markdown]
# ## 4. PCA per cell type (kept samples), and distance to the cell type × condition centroid
# Space: the DESeq2 `varianceStabilizingTransformation(blind=TRUE)` values written by the R script (genes with baseMean > 1 in that cell type; the same data variance partitioning was fit on), centred, not scaled (decision L2).
# Distance = Euclidean distance in that full space from each sample to the mean of its cell type × condition; a ratio to the group's median distance is shown for scale. No automatic flag.

# %%
dist_rows = []
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
for ax, ct in zip(axes.ravel(), CELL_TYPES):
    meta = kept[kept["cell_type"] == ct].set_index("sample")
    vs = vst[vst["cell_type"] == ct].pivot(index="gene", columns="sample", values="vst")[meta.index]
    genes = vs.index
    centred = vs.sub(vs.mean(axis=1), axis=0).T.to_numpy()
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    pcs, var_exp = u * s, s**2 / (s**2).sum()
    for c in CONDS:
        m = (meta["condition"] == c).to_numpy()
        ax.scatter(pcs[m, 0], pcs[m, 1], color=COND_COLORS[c], s=30, label=c)
        for x, y, d in zip(pcs[m, 0], pcs[m, 1], meta.loc[m, "donor"]):
            ax.annotate(d.replace("patient_", ""), (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points")
        group = vs.loc[:, m]
        d_ = np.sqrt(((group.sub(group.mean(axis=1), axis=0)) ** 2).sum(axis=0))
        dist_rows += [{"cell_type": ct, "condition": c, "sample": smp, "donor": meta.loc[smp, "donor"],
                       "distance": dist, "ratio_to_group_median": dist / d_.median()} for smp, dist in d_.items()]
    ax.set_title(f"{ct} ({len(genes)} genes)", fontsize=9)
    ax.set_xlabel(f"PC1 ({var_exp[0]:.0%})"); ax.set_ylabel(f"PC2 ({var_exp[1]:.0%})")
axes[0, 0].legend(fontsize=7)
fig.suptitle("PCA of kept pseudobulk samples (DESeq2 VST, blind)", y=1.0)
fig.tight_layout()

# %%
dist = pd.DataFrame(dist_rows)
fig, ax = plt.subplots(figsize=(12, 4))
groups = [(ct, c) for ct in CELL_TYPES for c in CONDS]
for i, (ct, c) in enumerate(groups):
    d = dist[(dist["cell_type"] == ct) & (dist["condition"] == c)]
    ax.scatter(np.full(len(d), i) + np.linspace(-0.15, 0.15, len(d)), d["ratio_to_group_median"], color=COND_COLORS[c], s=18)
ax.set_xticks(range(len(groups)), [f"{ct} · {c}" for ct, c in groups], rotation=90, fontsize=7)
ax.set_ylabel("distance to centroid / group median")
ax.set_title("Distance of each sample to its cell type × condition centroid (VST space)")
fig.tight_layout()
print("Largest 10 ratios:")
display(dist.sort_values("ratio_to_group_median", ascending=False).head(10))

# %% [markdown]
# ## 5. Variance partitioning per cell type
# Fraction of each gene's variance (DESeq2 `varianceStabilizingTransformation(blind=TRUE)`, genes with baseMean > 1) explained by each term of `~ (1|donor) + (1|condition) + log(n_cells) + log(total_counts)`. Each gene's fractions sum to 1.
# A4: no numeric rule; the user judges from these plots whether M1b's `~ donor + condition` needs a covariate.
#
# **How to read these violins.** Each violin is all genes in the cell type, and most genes don't respond to IFN-β,
# so condition's median share is ~0; for the genes M1b calls DE, condition explains the majority (median 57–83%).
# Donor's share of 0 is an artifact: `n_cells` and `total_counts` are almost identical to donor (section 6), so the
# variance donor would explain is assigned to them; donor + covariates together give donor's real share (median
# ~15–34%). The DE model removes donor-level differences, so these covariates don't create stim-vs-ctrl calls.
# DE vs non-DE genes are compared in `notebooks/m1b_de.ipynb`, section 8.

# %%
print("DESeq2 dispersion trend actually used per cell type (parametric, or local when DESeq2 substituted it):")
display(fits)
print("Genes per cell type after baseMean >", MIN_BASEMEAN)
display(varpart.groupby("cell_type").size().to_frame("genes"))
print("Median variance fraction per term:")
display(varpart.groupby("cell_type")[TERMS].median().round(3))
print("Share of genes where a technical covariate explains > 10% of variance:")
display(varpart.groupby("cell_type")[["log_n_cells", "log_total_counts"]].agg(lambda x: (x > 0.1).mean()).round(3))

# %%
fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey=True)
for ax, ct in zip(axes.ravel(), CELL_TYPES):
    v = varpart[varpart["cell_type"] == ct]
    parts = ax.violinplot([v[t] for t in TERMS], showmedians=True, widths=0.8)
    for body in parts["bodies"]:
        body.set_facecolor("#2a78d6"); body.set_alpha(0.5)
    ax.set_xticks(range(1, len(TERMS) + 1), TERMS, rotation=45, ha="right", fontsize=7)
    ax.set_title(f"{ct} ({len(v)} genes)", fontsize=9)
axes[0, 0].set_ylabel("fraction of variance"); axes[1, 0].set_ylabel("fraction of variance")
fig.tight_layout()

# %% [markdown]
# ## 6. Collinearity of the model terms (variancePartition `canCorPairs`, kept samples)
# Canonical correlation between each pair of terms, per cell type (`R/m1a_cancor.R` → `m1a_cancor.csv`). Values near 1 mean two terms carry the same information, so variance partitioning cannot split variance between them.

# %%
cancor = pd.read_csv(PROC / "m1a_cancor.csv")
CC_TERMS = ["donor", "condition", "log_n_cells", "log_total_counts"]
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
for ax, ct in zip(axes.ravel(), CELL_TYPES):
    m = cancor[cancor["cell_type"] == ct].pivot(index="term1", columns="term2", values="cancor").loc[CC_TERMS, CC_TERMS]
    ax.imshow(m.to_numpy(), cmap="Blues", vmin=0, vmax=1)
    for i in range(len(CC_TERMS)):
        for j in range(len(CC_TERMS)):
            ax.text(j, i, f"{m.iat[i, j]:.2f}", ha="center", va="center", fontsize=8, color="white" if m.iat[i, j] > 0.6 else "#2b2b28")
    ax.set_xticks(range(len(CC_TERMS)), CC_TERMS, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(CC_TERMS)), CC_TERMS, fontsize=7)
    ax.set_title(ct, fontsize=9); ax.grid(False)
fig.tight_layout()
display(cancor[cancor["term1"].isin(["donor", "condition"]) & cancor["term2"].isin(["log_n_cells", "log_total_counts"])]
        .pivot_table(index="cell_type", columns=["term1", "term2"], values="cancor").round(3))
