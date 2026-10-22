# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M0 explained by example
#
# Written by a reviewer who did not write `src/m0_qc.py`. Each section takes one key line, states in one sentence
# what it should do, and runs the **real** function on a tiny hand-checkable input.
# The `expected` value in each code cell is written first, from the wording of `audit/contracts/M0_eda.md`,
# `audit/decisions.md`, `CLAUDE.md` or the M0 reviewer brief at commit `2bd9257` (in git history), never from the code. Then the code runs, and
# each table shows **input | expected | code gives**. An `assert` stops the notebook if they disagree.
#
# Thresholds are passed explicitly. Where the real threshold (200 genes, 10 cells, 3 donors) would need a big input,
# a smaller one is passed so the example stays hand-checkable; the rule (≥ or <) is what is being checked.

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from IPython.display import display

from src import m0_qc


def compare(table, expected_col="expected", actual_col="code gives"):
    """Show the table and stop if any row's expected differs from what the code gives."""
    table["match"] = [e == a for e, a in zip(table[expected_col], table[actual_col])]
    display(table)
    assert table["match"].all(), table[~table["match"]]


# %% [markdown]
# ## 1. Join key and label lookup (`join_author_metadata`)
#
# **Lines:** `src/m0_qc.py:39` (strip `-<digits>`), `src/m0_qc.py:98` (GEO key), `src/m0_qc.py:101` (cell key),
# `src/m0_qc.py:102` (look-up), `src/m0_qc.py:123` (label copied from GEO `cell`).
#
# **Should do:** each h5ad cell gets the authors' cell type from the GEO row with the same barcode sequence **and** the same condition.
#
# Tiny input. `AAA` is a **collision**: the same sequence was captured in both runs. The h5ad calls the stim one `AAA-2`,
# GEO calls it `AAA-11`. `GGG` is stim only, so GEO keeps `-1`.
#
# | h5ad barcode | h5ad label | | GEO barcode | GEO stim | GEO cell |
# |---|---|---|---|---|---|
# | AAA-1 | ctrl | | AAA-1 | ctrl | CD4 T cells |
# | AAA-2 | stim | | AAA-11 | stim | B cells |
# | CCC-1 | ctrl | | CCC-1 | ctrl | NK cells |
# | GGG-2 | stim | | GGG-1 | stim | CD14+ Monocytes |
#
# A raw-barcode join would match only `AAA-1` and `CCC-1` (CLAUDE.md: "a raw-barcode join silently matches only control cells").

# %%
cells = pd.DataFrame(
    {"label": ["ctrl", "stim", "ctrl", "stim"],
     "replicate": ["patient_101", "patient_101", "patient_107", "patient_107"],
     "cluster": [1, 2, 3, 4], "tsne1": [0.1, 0.2, 0.3, 0.4], "tsne2": [1.1, 1.2, 1.3, 1.4]},
    index=["AAA-1", "AAA-2", "CCC-1", "GGG-2"],
)
h5ad_toy = ad.AnnData(X=sp.csr_matrix(np.zeros((4, 2))), obs=cells)
geo_toy = pd.DataFrame(
    {"stim": ["ctrl", "stim", "ctrl", "stim"], "ind": [101, 101, 107, 107],
     "cluster": [1, 2, 3, 4], "cell": ["CD4 T cells", "B cells", "NK cells", "CD14+ Monocytes"],
     "tsne1": [0.1, 0.2, 0.3, 0.4], "tsne2": [1.1, 1.2, 1.3, 1.4], "multiplets": ["singlet"] * 4},
    index=["AAA-1", "AAA-11", "CCC-1", "GGG-1"],
)

# Decision: "Join on (barcode without suffix, condition)". AAA-2 (stim) must get the stim AAA row (B cells),
# not the ctrl AAA row (CD4 T cells); GGG-2 matches GGG-1 because the suffix is ignored.
expected = ["CD4 T cells", "B cells", "NK cells", "CD14+ Monocytes"]

joined = m0_qc.join_author_metadata(h5ad_toy, geo_toy)
compare(pd.DataFrame({
    "input: h5ad barcode": cells.index, "input: label": cells["label"].to_numpy(),
    "expected": expected, "code gives": joined.obs["cell_type"].tolist(),
}))
# Contract step 2: the join keeps the cell count and order.
assert joined.obs["barcode"].tolist() == ["AAA-1", "AAA-2", "CCC-1", "GGG-2"]

# %% [markdown]
# ## 2. QC metrics (`compute_qc_metrics`)
#
# **Lines:** `src/m0_qc.py:152` (n_counts), `src/m0_qc.py:153` (n_genes), `src/m0_qc.py:154–156` (pct_ribo),
# `src/m0_qc.py:140` (`ribosomal_mask`: prefix minus `NON_RIBOSOMAL`, list at `src/m0_qc.py:19`), `src/m0_qc.py:157` (var n_cells).
#
# **Should do:** per cell, total counts, genes with a count > 0, and the % of counts in ribosomal genes (name starts with `RPL`/`RPS`, minus 14 listed non-ribosomal genes); per gene, the number of cells with a count > 0.
#
# Rows are cells, columns are genes:
#
# | | RPL5 | RPS3 | ISG15 | RPS6KA1 |
# |---|---|---|---|---|
# | c1 | 2 | 0 | 6 | 2 |
# | c2 | 0 | 0 | 5 | 0 |
# | c3 | 1 | 3 | 0 | 0 |
#
# `RPS6KA1` is a kinase, not a ribosomal protein. It starts with `RPS` but is on the user's 14-gene exclusion list
# (decision D4), so it must **not** count as ribosomal. Its counts still count towards `n_counts` (it is still a gene).

# %%
genes = ["RPL5", "RPS3", "ISG15", "RPS6KA1"]
X = np.array([[2, 0, 6, 2],
              [0, 0, 5, 0],
              [1, 3, 0, 0]])
qc_toy = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=["c1", "c2", "c3"]), var=pd.DataFrame(index=genes))

# Contract step 3: "total counts, genes detected, % ribosomal" per cell = sum across a ROW.
# decisions.md (D4): '"ribosomal" = RPS/RPL prefix minus 14 non-ribosomal genes (RPS6KA1-KA5, ...)'.
# So ribosomal = RPL5, RPS3 only.
expected_mask = [True, True, False, False]  # RPL5, RPS3, ISG15, RPS6KA1
# c1: 2+0+6+2 = 10, 3 genes > 0, ribo (RPL5+RPS3) = 2+0 = 2 -> 20%   (RPS6KA1's 2 are NOT ribosomal)
# c2: 5, 1 gene, ribo 0 -> 0%
# c3: 1+3 = 4, 2 genes, ribo 4 -> 100%
expected_cells = {"n_counts": [10, 5, 4], "n_genes": [3, 1, 2], "pct_ribo": [20.0, 0.0, 100.0]}
# Contract outputs: var n_cells = cells with the gene detected = count > 0 down a COLUMN.
expected_genes = [2, 1, 2, 1]

compare(pd.DataFrame({"input: gene": genes, "expected": expected_mask,
                      "code gives": m0_qc.ribosomal_mask(genes).tolist()}))
# real list: 98 genes start with RPS/RPL, 84 after the exclusion (D4); check the list has 14 names
assert len(m0_qc.NON_RIBOSOMAL) == 14

qc = m0_qc.compute_qc_metrics(qc_toy)
rows = []
for m in expected_cells:
    for i, c in enumerate(["c1", "c2", "c3"]):
        rows.append({"input": f"{c} row {X[i].tolist()}", "metric": m,
                     "expected": expected_cells[m][i], "code gives": float(qc.obs[m].iloc[i])})
compare(pd.DataFrame(rows))

# %% [markdown]
# ### 2b. Axis check: per cell (row) vs per gene (column)
#
# Same matrix. `n_counts` must be a **row** sum (one number per cell, 3 numbers). If the code summed down **columns**
# instead, it would give one number per gene: RPL5 3, RPS3 3, ISG15 11, RPS6KA1 2. `n_cells` must be a **column**
# count (one number per gene, 4 numbers).

# %%
wrong_axis = X.sum(axis=0).tolist()  # what a column sum would look like: [3, 3, 11, 2]
compare(pd.DataFrame({
    "input": [f"gene {g}: column {X[:, j].tolist()}" for j, g in enumerate(genes)],
    "expected": expected_genes,  # cells with count > 0 in that column
    "code gives": qc.var["n_cells"].tolist(),
    "column SUM (not what n_cells is)": wrong_axis,
}))
compare(pd.DataFrame({
    "input": [f"cell {c}: row {X[i].tolist()}" for i, c in enumerate(["c1", "c2", "c3"])],
    "expected": expected_cells["n_counts"],
    "code gives": qc.obs["n_counts"].tolist(),
}))
assert qc.obs["n_counts"].tolist() != wrong_axis[:3]  # would be [3, 3, 11] if summed down columns

# %% [markdown]
# ## 3. Basic filter (`basic_filter`), cells first, then genes
#
# **Lines:** `src/m0_qc.py:170` (keep cells), `src/m0_qc.py:172` (keep genes, counted on the kept cells).
#
# **Should do:** keep cells with ≥ `min_genes` genes detected, then, among those cells only, keep genes detected in ≥ `min_cells` cells.
#
# Real thresholds are 200 genes / 3 cells; here `min_genes=2`, `min_cells=2`.
#
# | | A | B | C | genes detected |
# |---|---|---|---|---|
# | c1 | 1 | 1 | 0 | 2 |
# | c2 | 1 | 1 | 0 | 2 |
# | c3 | 0 | 0 | 1 | 1 |
# | c4 | 1 | 0 | 1 | 2 |
#
# Gene C is in 2 cells before filtering, but c3 is removed first, so C is left in 1 cell and must go.
# If genes were filtered first, C would wrongly survive.

# %%
Xf = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1]])
f_toy = ad.AnnData(X=sp.csr_matrix(Xf), obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
                   var=pd.DataFrame(index=["A", "B", "C"]))

# CLAUDE.md: "Valid cells: expressing >= 200 genes. Valid genes: expressed in >= 3 cells"; contract step 4:
# "keep cells with >= 200 genes, THEN genes detected in >= 3 cells". c1, c2, c4 have exactly 2 (>= keeps them);
# c3 has 1. On c1, c2, c4: A in 3 cells, B in 2 (kept, boundary), C in 1 (dropped).
expected_cells_kept = ["c1", "c2", "c4"]
expected_genes_kept = ["A", "B"]

filt = m0_qc.basic_filter(f_toy, min_genes=2, min_cells=2)
compare(pd.DataFrame({
    "input": ["cells", "genes"],
    "expected": [expected_cells_kept, expected_genes_kept],
    "code gives": [filt.obs_names.tolist(), filt.var_names.tolist()],
}))

# %% [markdown]
# ## 4. MAD outlier rule (`mad_outlier`), one group
#
# **Lines:** `src/m0_qc.py:182` (median), `src/m0_qc.py:183` (raw MAD), `src/m0_qc.py:185` (flag if > n_mads × MAD).
#
# **Should do:** flag a value whose distance from the group median is **more than** 5 raw (unscaled) MADs.
#
# | values | median | distances from median | MAD | cut (5 × MAD) |
# |---|---|---|---|---|
# | −1, −1, 0, 0, 1, 1, **5** | 0 | 1, 1, 0, 0, 1, 1, 5 | 1 | 5 |
# | −1, −1, 0, 0, 1, 1, **5.01** | 0 | 1, 1, 0, 0, 1, 1, 5.01 | 1 | 5 |
#
# 5 is exactly at the cut (not more than 5): not flagged. 5.01 is just past it: flagged. With a *scaled* MAD
# (× 1.4826) the cut would be 7.41 and 5.01 would not be flagged, so this example also checks "raw".

# %%
at_edge = [-1, -1, 0, 0, 1, 1, 5]
past_edge = [-1, -1, 0, 0, 1, 1, 5.01]
# CLAUDE.md / decisions: "more than 5 raw (unscaled) MADs from the median" -> strictly greater, MAD not scaled.
expected_at = [False] * 7
expected_past = [False] * 6 + [True]

rows = []
for vals, exp in [(at_edge, expected_at), (past_edge, expected_past)]:
    got = m0_qc.mad_outlier(vals, n_mads=5.0).tolist()
    rows.append({"input": vals, "expected": exp, "code gives": got})
compare(pd.DataFrame(rows))

# %% [markdown]
# ## 5. QC flags (`flag_qc`)
#
# ### 5a. Outlier grouping: per donor × condition × cell type, both tails, small groups skipped
#
# **Lines:** `src/m0_qc.py:203` (loop over donor × condition × cell type groups), `src/m0_qc.py:204` (skip groups
# with < `min_cells_group` cells), `src/m0_qc.py:206–207` (log1p, n_counts or n_genes), `src/m0_qc.py:185` (both tails).
#
# **Should do:** within each donor × condition × cell type group of at least `min_cells_group` cells, flag a cell whose log1p(n_counts) or log1p(n_genes) is more than 5 raw MADs above **or below** the group median; cells in smaller groups are never flagged.
#
# Real threshold is 10 cells; here `min_cells_group=4`. `n_counts` is set so that log1p(n_counts) is the round
# number shown. `n_genes` is the same in every cell, so it flags nothing.
#
# | cell | donor | condition | cell type | log1p(n_counts) | group size |
# |---|---|---|---|---|---|
# | c1 | A | ctrl | B cells | 10 | 4 |
# | c2 | A | ctrl | B cells | 11 | 4 |
# | c3 | A | ctrl | B cells | 11.5 | 4 |
# | c4 | A | ctrl | B cells | **5** | 4 |
# | c5 | A | ctrl | NK cells | 5 | 4 |
# | c6 | A | ctrl | NK cells | 5.5 | 4 |
# | c7 | A | ctrl | NK cells | 6 | 4 |
# | c8 | A | ctrl | NK cells | 6.2 | 4 |
# | c9 | A | stim | B cells | 11 | 3 (skipped) |
# | c10 | A | stim | B cells | 11 | 3 (skipped) |
# | c11 | A | stim | B cells | **3** | 3 (skipped) |
#
# By hand:
# - **A/ctrl/B cells:** values 5, 10, 11, 11.5 → median 10.5; distances 0.5, 0.5, 1, 5.5 → MAD 0.75, cut 3.75.
#   c4 is 5.5 **below** the median: flagged (low side).
# - **A/ctrl/NK cells:** values 5, 5.5, 6, 6.2 → median 5.75; distances 0.75, 0.25, 0.25, 0.45 → MAD 0.35, cut 1.75: none flagged.
#   c5 has the same value (5) as c4, but it is normal for an NK cell: each cell is judged against its own cell type.
# - **A/stim/B cells:** 3 cells < 4, not assessed: none flagged. (If it were assessed: median 11, MAD 0, and c11 would be flagged.)
# - **Old rule** (donor × condition, no skipping): A/ctrl pooled = 5, 5, 5.5, 6, 6.2, 10, 11, 11.5 → median 6.1,
#   MAD 1.1, cut 5.5 → nothing flagged, so c4 escapes; A/stim assessed → c11 flagged.

# %%
logs = [10, 11, 11.5, 5, 5, 5.5, 6, 6.2, 11, 11, 3]
obs5 = pd.DataFrame({
    "donor": ["A"] * 11,
    "condition": ["ctrl"] * 8 + ["stim"] * 3,
    "cell_type": ["B cells"] * 4 + ["NK cells"] * 4 + ["B cells"] * 3,
    "n_counts": np.expm1(logs),
    "n_genes": [100] * 11,
}, index=[f"c{i}" for i in range(1, 12)])
fq_toy = ad.AnnData(X=sp.csr_matrix(np.zeros((11, 2))), obs=obs5)

# decisions.md (D2): "> 5 MADs above or below the median on log counts / log genes, computed per
# donor x condition x cell type (groups < 10 cells not assessed)". Hand computation above: only c4 (low side).
expected = [False, False, False, True] + [False] * 4 + [False] * 3
old_rule = [False] * 8 + [False, False, True]  # donor x condition, every group assessed (hand computation above)

fq = m0_qc.flag_qc(fq_toy, n_mads=5.0, min_cells_group=4, min_paired_donors=1)
compare(pd.DataFrame({
    "input": [f"{d}/{c}/{t} log={l}" for d, c, t, l in
              zip(obs5["donor"], obs5["condition"], obs5["cell_type"], logs)],
    "expected": expected, "code gives": fq.obs["flag_outlier"].tolist(),
    "old rule would give": old_rule,
}, index=obs5.index))
assert len(fq) == 11  # "Flags only (no drops)"

# %% [markdown]
# ### 5b. Small groups (`flag_small_group`)
#
# **Lines:** `src/m0_qc.py:210–211`.
#
# **Should do:** flag every cell whose donor × condition × cell type group has fewer than `min_cells_group` cells.
#
# Real threshold is 10; here `min_cells_group=3`.
#
# | cells | donor | condition | cell type | group size |
# |---|---|---|---|---|
# | c1–c3 | A | ctrl | B cells | 3 |
# | c4–c5 | A | stim | B cells | 2 |
# | c6 | B | ctrl | B cells | 1 |
# | c7 | A | ctrl | NK cells | 1 |
#
# c7 shares donor and condition with c1–c3; if cell type were left out of the grouping, A/ctrl would have 4 cells and c7 would escape.

# %%
obs6 = pd.DataFrame({
    "donor": ["A", "A", "A", "A", "A", "B", "A"],
    "condition": ["ctrl", "ctrl", "ctrl", "stim", "stim", "ctrl", "ctrl"],
    "cell_type": ["B cells"] * 6 + ["NK cells"],
    "n_counts": [100.0] * 7, "n_genes": [50.0] * 7,
}, index=[f"c{i}" for i in range(1, 8)])
sg_toy = ad.AnnData(X=sp.csr_matrix(np.zeros((7, 2))), obs=obs6)

# CLAUDE.md: "donor x condition x cell type with < 10 cells". Size 3 is not < 3 (kept); sizes 2 and 1 are flagged.
expected = [False, False, False, True, True, True, True]

sg = m0_qc.flag_qc(sg_toy, n_mads=5.0, min_cells_group=3, min_paired_donors=1)
compare(pd.DataFrame({
    "input": [f"{d}/{c}/{t}" for d, c, t in zip(obs6["donor"], obs6["condition"], obs6["cell_type"])],
    "expected": expected, "code gives": sg.obs["flag_small_group"].tolist(),
}, index=obs6.index))

# %% [markdown]
# ### 5c. Underpowered cell types (`flag_underpowered_celltype`)
#
# **Lines:** `src/m0_qc.py:214–217`.
#
# **Should do**, for each cell type:
# 1. For each donor, check it has **enough cells in both conditions** (≥ `min_cells_group` in ctrl **and** ≥ `min_cells_group` in stim). If so, it's a "paired donor".
# 2. **Count** the paired donors.
# 3. If the count is **at least** `min_paired_donors`, keep the cell type; otherwise flag it as underpowered.
#
# **Real example** (thresholds 10 cells, 3 donors): megakaryocytes have 4 donors with ≥ 10 cells in both conditions → 4 ≥ 3 → **kept** (not flagged). A cell type with only 2 such donors would be flagged.
#
# The toy below uses smaller thresholds so it stays tiny: `min_cells_group=2`, `min_paired_donors=1`.
#
# | cell type | donor | ctrl cells | stim cells | paired donor? |
# |---|---|---|---|---|
# | B cells | A | 2 | 2 | yes |
# | NK cells | A | 2 | **1** | no: stim group is small |
#
# B cells: 1 paired donor, not < 1, so not flagged. NK cells: donor A has both conditions **present**, but the stim
# group (1 cell) fails the small-group rule, so 0 paired donors: flagged. An off-by-one here (counting presence
# instead of surviving groups) would leave NK unflagged.

# %%
obs7 = pd.DataFrame({
    "donor": ["A"] * 7,
    "condition": ["ctrl", "ctrl", "stim", "stim", "ctrl", "ctrl", "stim"],
    "cell_type": ["B cells"] * 4 + ["NK cells"] * 3,
    "n_counts": [100.0] * 7, "n_genes": [50.0] * 7,
}, index=[f"c{i}" for i in range(1, 8)])
up_toy = ad.AnnData(X=sp.csr_matrix(np.zeros((7, 2))), obs=obs7)

# CLAUDE.md: "fewer than 3 donors with both conditions left AFTER the rule above" (the < 10 cells rule);
# brief §2: "a paired donor has both conditions with >= 10 cells each".
expected = [False] * 4 + [True] * 3

up = m0_qc.flag_qc(up_toy, n_mads=5.0, min_cells_group=2, min_paired_donors=1)
compare(pd.DataFrame({
    "input": [f"{t} {d}/{c}" for t, d, c in zip(obs7["cell_type"], obs7["donor"], obs7["condition"])],
    "expected": expected, "code gives": up.obs["flag_underpowered_celltype"].tolist(),
}, index=obs7.index))

# %% [markdown]
# ## 6. Which label is the control (`infer_conditions`)
#
# **Lines:** `src/m0_qc.py:17` (control names), `src/m0_qc.py:227–233`.
#
# **Should do:** with exactly 2 labels, return (control, treatment), where the control is whichever is `ctrl` or `control`; anything else is an error.
#
# | input labels | expected |
# |---|---|
# | ctrl, stim, stim | (ctrl, stim) |
# | stim, control | (control, stim) |
# | a, b | error (no control) |
# | ctrl, control, stim | error (3 labels) |
# | ctrl | error (1 label) |

# %%
def run(labels):
    try:
        return m0_qc.infer_conditions(labels)
    except ValueError:
        return "error"

# decisions.md: "conditions: exactly 2 labels, the control is whichever is `ctrl` or `control` ..., the other is the treatment"
cases = [
    (["ctrl", "stim", "stim"], ("ctrl", "stim")),
    (["stim", "control"], ("control", "stim")),
    (["a", "b"], "error"),
    (["ctrl", "control", "stim"], "error"),
    (["ctrl"], "error"),
]
compare(pd.DataFrame({
    "input": [c for c, _ in cases],
    "expected": [e for _, e in cases],
    "code gives": [run(c) for c, _ in cases],
}))
