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
# # M2a: shared vs cell-type-specific calls
#
# Contract: `audit/contracts/M2_shared_specific.md` (M2a). Logic: `src/m2a_calls.py` (tested in `tests/test_m2a.py`).
# This notebook only reads `m2_calls.csv` and `m2_table.csv` and shows them.
# Labels: `mixed` > `shared` / `shared*` > `partial` > `specific*` / `specific` (decisions O1, O2, specific* rule).

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display

from src.m2a_calls import CELL_TYPES

pd.set_option("display.width", 160, "display.max_columns", 20)
PROC = ROOT / "data/processed"
calls = pd.read_csv(PROC / "m2_calls.csv")
table = pd.read_csv(PROC / "m2_table.csv")
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})

# %% [markdown]
# ## Outputs: shape, head, tail

# %%
for name, df in [("m2_calls.csv", calls), ("m2_table.csv", table)]:
    print(name, df.shape)
    display(df.head(), df.tail())

# %% [markdown]
# ## Label counts and number of responding cell types

# %%
ORDER = ["shared", "shared*", "partial", "specific", "specific*", "mixed"]
display(table["shared_response"].value_counts().reindex(ORDER).rename("genes").to_frame())
display(table["n_cell_types"].value_counts().sort_index().rename("genes").to_frame())
print("megakaryocyte_lower_confidence:", int(table["megakaryocyte_lower_confidence"].sum()), "genes")

# %% [markdown]
# ## Table sorted by label, then number of cell types

# %%
sorted_table = table.assign(_o=table["shared_response"].map(ORDER.index)) \
    .sort_values(["_o", "n_cell_types", "gene"], ascending=[True, False, True]).drop(columns="_o")
display(sorted_table.head(30), sorted_table.tail(10))

# %% [markdown]
# ## Which cell-type combinations occur most (top 25 patterns, upset-style)
# Each column of the dot matrix is one pattern (set of cell types where the gene is DE, either direction);
# bar = number of genes with exactly that pattern.

# %%
pattern = table[CELL_TYPES].ne("–")
counts = pattern.value_counts().head(25)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1.4]})
ax1.bar(range(len(counts)), counts.to_numpy(), color="#4a6fa5")
ax1.set_ylabel("genes")
for i, v in enumerate(counts.to_numpy()):
    ax1.text(i, v, str(v), ha="center", va="bottom", fontsize=7)
for i, combo in enumerate(counts.index):
    on = [j for j, x in enumerate(combo) if x]
    ax2.scatter([i] * len(CELL_TYPES), range(len(CELL_TYPES)), color="#dddddd", s=25)
    ax2.scatter([i] * len(on), on, color="#222222", s=25)
    ax2.plot([i, i], [min(on), max(on)], color="#222222", lw=1)
ax2.set_yticks(range(len(CELL_TYPES)), CELL_TYPES)
ax2.set_xticks([])
ax2.invert_yaxis()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Direction heatmap (genes DE in ≥ 2 cell types, sorted by label and pattern)

# %%
multi = sorted_table[sorted_table["n_cell_types"] >= 2]
val = multi[CELL_TYPES].apply(lambda col: col.map({"↑": 1, "↓": -1, "–": 0}))
fig, ax = plt.subplots(figsize=(6, 8))
ax.imshow(val.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5, interpolation="nearest")
ax.set_xticks(range(len(CELL_TYPES)), CELL_TYPES, rotation=45, ha="right")
ax.set_yticks([])
ax.set_ylabel(f"{len(multi)} genes (red ↑, blue ↓)")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## LODO stability (R1): DE calls keeping their sign in every leave-one-donor-out refit

# %%
de_calls = calls[calls["de"] != "none"]
display(de_calls.groupby("cell_type")["lodo_stable"].agg(n_calls="size", stable="sum", pct=lambda s: 100 * s.mean()))
print(f"overall: {de_calls['lodo_stable'].sum()} / {len(de_calls)} = {100 * de_calls['lodo_stable'].mean():.1f}%")
