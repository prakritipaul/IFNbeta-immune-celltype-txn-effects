# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # Figs 1–2: composition and UMAP
#
# Contract: `audit/contracts/F12_figures.md`. Plot and UMAP code lives in `src/figures.py`.
# Input: `data/processed/m0_cells.h5ad` (all 24,562 cells; outliers kept). UMAP recipe C2, no integration (C1).

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
from IPython.display import display

from src import figures

N_CELLS = 24562
(ROOT / "figures").mkdir(exist_ok=True)
adata = ad.read_h5ad(ROOT / "data/processed/m0_cells.h5ad")
assert adata.n_obs == N_CELLS, adata.n_obs
print(adata)

# %% [markdown]
# ## Fig 1: cells per donor, condition and cell type

# %%
counts = figures.composition_counts(adata.obs)

# Expected cells per donor × condition: the table in audit/data_facts.md
lines = (ROOT / "audit/data_facts.md").read_text().split("## Cells per donor × condition")[1].split("##")[0].strip().splitlines()
rows = [[c.strip() for c in l.strip("|").split("|")] for l in lines if l.startswith("| patient_")]
facts = pd.DataFrame(rows, columns=["donor", "ctrl", "stim", "total"]).set_index("donor").astype(int)
expected = facts[["ctrl", "stim"]].stack().rename_axis(["donor", "condition"])
bar_totals = counts.sum(axis=1)
assert (bar_totals.values == expected.loc[bar_totals.index].values).all() and len(expected) == len(bar_totals) == 16
assert bar_totals.sum() == N_CELLS
print("PASS: Fig 1 bar totals equal data_facts.md for all 16 donor × condition bars; sum =", bar_totals.sum())

fig1 = figures.plot_fig1(counts)
for ext in ("png", "pdf"):
    fig1.savefig(ROOT / f"figures/fig1_composition.{ext}", dpi=200, bbox_inches="tight")
display(counts.assign(total=bar_totals))

# %% [markdown]
# ## Fig 2: UMAP (C2 recipe, no integration)

# %%
emb = figures.compute_umap(adata)
df = adata.obs[["barcode", "condition", "cell_type"]].join(emb)
assert len(df) == N_CELLS and df.index.equals(adata.obs_names)
assert np.isfinite(df[["umap1", "umap2"]].to_numpy()).all()

fig2 = figures.plot_fig2(df)
# I3: both panels plot the same cells at identical coordinates
off = [ax.collections[0].get_offsets() for ax in fig2.axes]
assert off[0].shape == off[1].shape == (N_CELLS, 2)
assert np.array_equal(np.asarray(off[0]), np.asarray(off[1]))
assert np.array_equal(np.sort(np.asarray(off[0]), axis=0), np.sort(df[["umap1", "umap2"]].to_numpy(), axis=0))
print(f"PASS I3: both panels draw the same {N_CELLS} cells at identical coordinates")
for ext in ("png", "pdf"):
    fig2.savefig(ROOT / f"figures/fig2_umap.{ext}", dpi=200, bbox_inches="tight")

# %%
out = df[["barcode", "condition", "umap1", "umap2"]]
out.to_csv(ROOT / "data/processed/f2_umap.csv", index=False)
back = pd.read_csv(ROOT / "data/processed/f2_umap.csv")
assert back.shape == (N_CELLS, 4) and back["barcode"].is_unique
print("f2_umap.csv shape:", back.shape)
display(back.head()); display(back.tail())
