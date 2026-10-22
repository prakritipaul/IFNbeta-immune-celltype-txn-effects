# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # Wave 4: Fig 3, Fig 4, Fig 5 and the final table
#
# Contract: `audit/contracts/W4_figures_table.md` (decisions F1–F6, F5). Code: `src/w4_outputs.py`.
# Nothing is recomputed: every number shown is read from M3, M2a or M2c outputs, and each builder asserts that
# what it shows equals its source file.

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import pandas as pd
from IPython.display import display

from src import w4_outputs as w4

pd.set_option("display.width", 200, "display.max_columns", 30)
PROC, FIG, RES = ROOT / "data/processed", ROOT / "figures", ROOT / "results"
RES.mkdir(exist_ok=True)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=200, bbox_inches="tight")


# %% [markdown]
# ## Fig 3: IFN signature score, control vs stimulated, per cell type (F1, F2)
# Violins: per-cell scores of kept samples only (A1: DCs 6 donors, megakaryocytes 4). Dots: donor means from
# `m3_sample_scores.csv`, joined ctrl→stim. `kept_cells` asserts the cell counts and means equal that file.

# %%
cell_scores = pd.read_csv(PROC / "m3_cell_scores.csv")
sample_scores = pd.read_csv(PROC / "m3_sample_scores.csv")
cells = w4.kept_cells(cell_scores, sample_scores)
print(f"cells: {len(cell_scores)} in m3_cell_scores -> {len(cells)} in kept samples "
      f"({len(cell_scores) - len(cells)} from donors dropped by A1); samples: {len(sample_scores)}")
assert len(cells) == sample_scores["n_cells"].sum()
fig3 = w4.plot_fig3(cells, sample_scores)
save(fig3, "fig3_ifn_signature")

# %% [markdown]
# ## Fig 4: top shared and cell-type-specific genes (F3–F5)
# Columns: average-linkage clustering of the up-regulated Jaccard matrix (1 − J). Rows: top 25 shared + shared\*
# by mean |log2FC| over the 8 cell types; then per cell type (in column order) up to 5 specific genes tested in all
# 8 cell types, by |log2FC| (F3 revised). \* in a label = shared\* / specific\*.

# %%
table = pd.read_csv(PROC / "m2_table.csv")
calls = pd.read_csv(PROC / "m2_calls.csv")
jac = pd.read_csv(PROC / "m2c_jaccard.csv", index_col=[0, 1])
col4, Z4 = w4.jaccard_order(jac)
rows4 = w4.fig4_rows(table, calls, col4)
lfc4, de4 = w4.fig4_matrix(rows4, calls, col4)
print("Fig 4 column order:", col4)
print("rows per group:", rows4.groupby("group", sort=False).size().to_dict())
print(f"cells: {lfc4.size}; untested (grey): {int(lfc4.isna().sum().sum())}; DE (dot): {int(de4.sum().sum())}; "
      "every tested cell equals m2_calls log2FC")
display(rows4)
fig4 = w4.plot_fig4(rows4, lfc4, de4, Z4)
save(fig4, "fig4_de_heatmap")

# %% [markdown]
# ## Fig 5: response programs across cell types
# Columns: average-linkage clustering on NES profiles (1 − Pearson r) over affected programs tested in all 8.
# Rows: the 58 affected programs grouped by label; within a label by number of affected cell types, then mean NES.

# %%
gsea = pd.read_csv(PROC / "m2c_gsea.csv")
programs = pd.read_csv(PROC / "m2c_programs.csv")
col5, Z5, n_complete = w4.nes_order(gsea, programs)
rows5, nes5, sig5 = w4.fig5_matrix(gsea, programs, col5)
print("Fig 5 column order:", col5, f"(clustered on {n_complete} programs tested in all 8)")
print("rows per label:", rows5["label"].value_counts().reindex(w4.PROGRAM_LABELS).to_dict())
print(f"cells: {nes5.size}; n/t: {int(nes5.isna().sum().sum())}; significant (dot): {int(sig5.sum().sum())}; "
      "every cell equals m2c_gsea NES; arrows equal m2c_programs")
fig5 = w4.plot_fig5(rows5, nes5, sig5, Z5)
save(fig5, "fig5_programs_heatmap")

# %% [markdown]
# ## Final table (F6)
# `final_table` asserts every arrow equals `m2_table.csv` and the `de` call in `m2_calls.csv`.

# %%
final = w4.final_table(table, calls)
# real-data check (kept here, not in the function, so final_table works on any input)
sp = final["Shared response?"] == "specific"
assert sp.sum() == 241 and (final.loc[sp, "expressed_in_n_cell_types"] == 1).sum() == 229, "expected 229 of 241 specific genes expressed in 1 cell type"
assert len(final) == len(table) == 2801
final.to_csv(RES / "final_table.csv", index=False)
print(final.shape)
print(final["Shared response?"].value_counts().reindex(w4.LABEL_ORDER).to_dict())
display(final.head(10))
display(final.tail(5))
