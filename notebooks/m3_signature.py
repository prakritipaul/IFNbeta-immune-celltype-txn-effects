# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M3: IFN signature scores (Hallmark IFN-α, IFN-γ control)
#
# Contract: `audit/contracts/M3_signature.md`. Logic in `src/m3_signature.py` (tested); outputs written by `run_m3`.
# This notebook reads those outputs and shows numbers and plots for P5, P6 and I6. **No pass marks.**

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import Markdown, display

from src import m3_signature as m3

pd.set_option("display.width", 160, "display.max_columns", 20)
SCORES = m3.SCORE_COLUMNS
LABELS = {"ifna_score": "IFN-α", "ifng_score": "IFN-γ (full)", "ifng_noshared_score": "IFN-γ (no shared genes)"}
COND_COLORS = {"ctrl": "#1b74ac", "stim": "#dd1b26"}
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e3", "grid.linewidth": 0.6, "axes.axisbelow": True, "figure.dpi": 110})

adata = ad.read_h5ad(ROOT / "data/processed/m0_cells.h5ad", backed="r")
CONTROL, TREATMENT = adata.uns["control"], adata.uns["treatment"]
cells = pd.read_csv(ROOT / "data/processed/m3_cell_scores.csv")
samples = pd.read_csv(ROOT / "data/processed/m3_sample_scores.csv")
CELL_TYPES = sorted(cells["cell_type"].unique())

# %% [markdown]
# ## 1. Gene-set coverage (after the 6 approved mappings, B4)

# %%
sets = m3.load_gene_sets(ROOT / m3.MSIGDB_GMT, m3.MSIGDB_SHA256)
mapped, coverage = m3.map_gene_sets(sets, adata.var_names)
display(coverage)
display(Markdown((ROOT / "audit/m3_gene_sets.md").read_text()))

# %% [markdown]
# ## 2. Outputs: shape, head and tail

# %%
for name, df in [("m3_cell_scores.csv", cells), ("m3_sample_scores.csv", samples)]:
    print(name, df.shape)
    display(df.head()); display(df.tail())
print("Kept donors per cell type (A1):")
display(samples.groupby("cell_type")["donor"].nunique().to_frame("donors"))

# %% [markdown]
# ## 3. Per-cell scores by cell type × condition

# %%
fig, axes = plt.subplots(len(SCORES), 1, figsize=(12, 3.2 * len(SCORES)), sharex=True)
for ax, score in zip(axes, SCORES):
    for j, cond in enumerate([CONTROL, TREATMENT]):
        data = [cells.loc[(cells.cell_type == ct) & (cells.condition == cond), score] for ct in CELL_TYPES]
        pos = np.arange(len(CELL_TYPES)) + (j - 0.5) * 0.4
        parts = ax.violinplot(data, positions=pos, widths=0.38, showmedians=True, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(COND_COLORS.get(cond, "#888")); b.set_alpha(0.6)
        parts["cmedians"].set_color("#2b2b28")
        ax.plot([], [], color=COND_COLORS.get(cond, "#888"), lw=6, label=cond)
    ax.set_ylabel(LABELS[score])
    ax.legend(loc="upper right", fontsize=8)
axes[-1].set_xticks(range(len(CELL_TYPES)), CELL_TYPES, rotation=30, ha="right")
fig.tight_layout()

# %% [markdown]
# ## 4. P5: donor-paired stim − ctrl of the per-sample mean score
# One point per donor; table = median across donors per cell type.

# %%
diff = m3.paired_differences(samples, control=CONTROL, treatment=TREATMENT)
med = diff.groupby("cell_type")[SCORES].median().rename(columns=LABELS)
med.insert(0, "donors", diff.groupby("cell_type")["donor"].nunique())
display(med.round(3))

fig, ax = plt.subplots(figsize=(10, 3.5))
for j, score in enumerate(SCORES):
    for i, ct in enumerate(CELL_TYPES):
        y = diff.loc[diff.cell_type == ct, score]
        ax.scatter(np.full(len(y), i + (j - 1) * 0.25), y, s=14, color=f"C{j}", label=LABELS[score] if i == 0 else None)
ax.axhline(0, color="#8a8a85", lw=0.8)
ax.set_xticks(range(len(CELL_TYPES)), CELL_TYPES, rotation=30, ha="right")
ax.set_ylabel(f"{TREATMENT} − {CONTROL} (mean score)")
ax.legend(fontsize=8)
fig.tight_layout()

# %% [markdown]
# ## 5. P6: IFN-α increase vs IFN-γ increase (full set and without the shared genes)
# Per donor × cell type; the diagonal is equal increase.

# %%
display((med[LABELS["ifna_score"]] - med[[LABELS["ifng_score"], LABELS["ifng_noshared_score"]]].T).T
        .rename(columns=lambda c: f"IFN-α − {c}").round(3))
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharex=True, sharey=True)
for ax, score in zip(axes, ["ifng_score", "ifng_noshared_score"]):
    for i, ct in enumerate(CELL_TYPES):
        d = diff[diff.cell_type == ct]
        ax.scatter(d["ifna_score"], d[score], s=16, color=f"C{i}", label=ct)
    lim = [min(0, diff[SCORES].min().min()), diff[SCORES].max().max() * 1.05]
    ax.plot(lim, lim, color="#8a8a85", lw=0.8, ls="--")
    ax.set_xlabel("IFN-α increase"); ax.set_ylabel(f"{LABELS[score]} increase")
axes[1].legend(fontsize=7, loc="upper left")
fig.tight_layout()

# %% [markdown]
# ## 6. I6: score vs log library size, per cell
# Pearson and Spearman correlations, pooled and within each cell type × condition (library size differs by condition, so pooled values mix in the IFN effect).

# %%
cells["log_n_counts"] = np.log10(adata.obs.loc[cells["barcode"], "n_counts"].to_numpy(dtype=float))
rows = []
for (ct, cond), g in [(("all", "all"), cells)] + list(cells.groupby(["cell_type", "condition"])):
    for score in SCORES:
        rows.append({"cell_type": ct, "condition": cond, "score": LABELS[score], "n_cells": len(g),
                     "pearson": g[score].corr(g["log_n_counts"], method="pearson"),
                     "spearman": g[score].corr(g["log_n_counts"], method="spearman")})
corr = pd.DataFrame(rows)
display(corr.pivot_table(index=["cell_type", "condition"], columns="score", values="spearman", sort=False).round(3))
display(corr.round(3))

fig, axes = plt.subplots(1, len(SCORES), figsize=(14, 4), sharex=True)
for ax, score in zip(axes, SCORES):
    for cond in [CONTROL, TREATMENT]:
        g = cells[cells.condition == cond]
        ax.scatter(g["log_n_counts"], g[score], s=1, alpha=0.2, color=COND_COLORS.get(cond, "#888"), label=cond, rasterized=True)
    ax.set_xlabel("log10 total counts"); ax.set_ylabel(LABELS[score])
axes[0].legend(markerscale=8, fontsize=8)
fig.tight_layout()
