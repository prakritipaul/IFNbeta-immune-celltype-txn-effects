# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M3 explained by example
#
# Written by a reviewer who did not write `src/m3_signature.py`. Each section takes one key line, states in one sentence
# what it should do, and runs the **real** function on a tiny hand-checkable input.
# The `expected` value in each code cell is written first, from `audit/contracts/M3_signature.md`, `audit/decisions.md`
# (B1–B4, A1, pinned defaults), `audit/test_plans/M3.md` and the M3 reviewer brief (`audit/reviews/006_…`), never from the code.
# Then the code runs, and each table shows **input | expected | code gives**. A mismatch is marked **POSSIBLE BUG**;
# the notebook keeps going and the last cell fails if any example mismatched.
#
# Scoring parameters are passed explicitly. The real ones are `ctrl_size=50, n_bins=25` (B1); the toys use
# `n_bins=2` so every gene falls in one expression bin, and `ctrl_size=50` is more than the toy has genes, so scanpy
# takes **all** non-set genes as controls. Nothing is random, and every score can be worked out by hand.
#
# Scanpy's score for one cell is **mean(log expression of set genes) − mean(log expression of control genes)**, where
# log expression = `log1p(count / cell total × 10,000)`.

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
from IPython.display import display

from src import m3_signature as m3
from src.m3_signature import IFNA, IFNG

FAILURES = []


def compare(table, name, expected_col="expected", actual_col="code gives", tol=None):
    """Show the table; record (not stop) any row whose expected differs from what the code gives."""
    if tol is None:
        table["match"] = [e == a for e, a in zip(table[expected_col], table[actual_col])]
    else:
        table["match"] = [abs(e - a) <= tol for e, a in zip(table[expected_col], table[actual_col])]
    display(table)
    if not table["match"].all():
        print(f"POSSIBLE BUG in {name}:\n{table[~table['match']]}")
        FAILURES.append(name)


def toy(counts, genes):
    """AnnData of raw counts: rows = cells c1, c2, ...; columns = genes."""
    X = np.array(counts, dtype=float)
    return ad.AnnData(X=X, obs=pd.DataFrame(index=[f"c{i + 1}" for i in range(len(X))]),
                      var=pd.DataFrame(index=genes))


# %% [markdown]
# ## 1. Gene-symbol mapping (`map_gene_sets`)
#
# **Line:** `src/m3_signature.py:56` (brief §4).
#
# **Should do:** rename only the 6 approved symbols (B4) to the names used in our data, then require every set gene to be present; any other missing gene stops the run and names it. No fuzzy matching.
#
# Our data's gene names: `DDX58, PSMA2.1, ISG15, WARS, STAT1, MX1.1`.
#
# | case | set genes (MSigDB names) | expected |
# |---|---|---|
# | a | IFN-α: RIGI, PSMA2, ISG15; IFN-γ: WARS1, STAT1 | IFN-α: DDX58, PSMA2.1, ISG15; IFN-γ: WARS, STAT1; coverage 100% |
# | b | IFN-α: ISG15, FOO | error naming FOO |
# | c | IFN-α: MX1 (data only has MX1.1) | error naming MX1: `MX1 → MX1.1` is not an approved mapping |

# %%
var_names = pd.Index(["DDX58", "PSMA2.1", "ISG15", "WARS", "STAT1", "MX1.1"])

# B4: TENT5A->FAM46A, WARS1->WARS, MARCHF1->MARCH1, RIGI->DDX58, TMT1B->METTL7B, PSMA2->PSMA2.1. Order kept.
expected_a = {IFNA: ["DDX58", "PSMA2.1", "ISG15"], IFNG: ["WARS", "STAT1"]}
expected_cov = {IFNA: (3, 3, 100.0), IFNG: (2, 2, 100.0)}  # n_genes, n_found, pct_found

mapped, cov = m3.map_gene_sets({IFNA: ["RIGI", "PSMA2", "ISG15"], IFNG: ["WARS1", "STAT1"]}, var_names)
cov = cov.set_index("set")
compare(pd.DataFrame({
    "input": ["IFN-α: RIGI, PSMA2, ISG15", "IFN-γ: WARS1, STAT1"] * 2,
    "what": ["mapped genes"] * 2 + ["coverage (n_genes, n_found, %)"] * 2,
    "expected": [expected_a[IFNA], expected_a[IFNG], expected_cov[IFNA], expected_cov[IFNG]],
    "code gives": [mapped[IFNA], mapped[IFNG]]
                  + [tuple(cov.loc[s, ["n_genes", "n_found", "pct_found"]].tolist()) for s in (IFNA, IFNG)],
}), "1a mapping")


def run_map(sets):
    try:
        m3.map_gene_sets(sets, var_names)
        return "no error"
    except ValueError as e:
        return f"error: {e}"


# Contract: "Any other missing gene = stop and report"; test plan: "Raises naming any gene still absent".
cases = [({IFNA: ["ISG15", "FOO"], IFNG: ["STAT1"]}, "FOO"),
         ({IFNA: ["MX1"], IFNG: ["STAT1"]}, "MX1")]
got = [run_map(s) for s, _ in cases]
compare(pd.DataFrame({
    "input": ["IFN-α: ISG15, FOO", "IFN-α: MX1 (data has MX1.1)"],
    "expected": [f"error naming {g}" for _, g in cases],
    "code gives": got,
    "code gives (shortened)": [f"error naming {g}" if g_msg.startswith("error") and g in g_msg else g_msg
                               for (_, g), g_msg in zip(cases, got)],
}), "1b unmapped gene raises", actual_col="code gives (shortened)")

# %% [markdown]
# ## 2. Library-size normalisation (`score_cells`)
#
# **Line:** `src/m3_signature.py:82` (brief §4).
#
# **Should do:** divide each cell's counts by its total and multiply by 10,000, then `log1p`, **before** scoring (B1), so sequencing depth alone cannot change a score.
#
# | | g1 | g2 | g3 | g4 | g5 | g6 | g7 | g8 |
# |---|---|---|---|---|---|---|---|---|
# | c1 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
# | c2 | 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 |
#
# c2 is c1 sequenced 10× deeper. After normalisation both rows are identical, so all three scores must be identical
# for c1 and c2. IFN-α = {g1, g2}, IFN-γ = {g3, g4}. (Scoring raw counts would give c2 a different score: its log
# values are all shifted, and the set genes and control genes shift by different amounts.)

# %%
genes2 = [f"g{i}" for i in range(1, 9)]
ad2 = toy([[1, 2, 3, 4, 5, 6, 7, 8], [10, 20, 30, 40, 50, 60, 70, 80]], genes2)
s2 = m3.score_cells(ad2, {IFNA: ["g1", "g2"], IFNG: ["g3", "g4"]},
                    target_sum=1e4, ctrl_size=50, n_bins=2, random_state=0)

# B1: normalize_total(1e4) + log1p -> c1 and c2 identical -> c2 score minus c1 score = 0 (float32 tolerance).
expected = [0.0, 0.0, 0.0]
compare(pd.DataFrame({
    "input": "c1 = [1..8], c2 = 10 × c1",
    "score": m3.SCORE_COLUMNS,
    "c1 score": s2.loc["c1"].to_numpy(), "c2 score": s2.loc["c2"].to_numpy(),
    "expected": expected,
    "code gives": (s2.loc["c2"] - s2.loc["c1"]).to_numpy(),
}), "2 library-size normalisation", tol=1e-5)
# The input object must not be modified (docstring: "adata is not modified").
assert ad2.X[1].tolist() == [10, 20, 30, 40, 50, 60, 70, 80]

# %% [markdown]
# ## 3. Score size and sign: an independent oracle (`score_cells`)
#
# **Lines:** `src/m3_signature.py:82` (normalise) and `:89` (scanpy scoring). Brief §5: no test checks a score's actual value.
#
# **Should do:** give a cell whose IFN-α genes are 10× higher than everything else a score of about **+ln 10 ≈ +2.30**, and a cell where everything is equal a score of **0**.
#
# 10 genes: S1, S2 (the IFN-α set), C1…C8. IFN-γ = {C1, C2} here (no overlap with IFN-α, so no-shared = full).
#
# | | S1 | S2 | C1…C8 | total |
# |---|---|---|---|---|
# | c1 | 1 | 1 | 1 each | 10 |
# | c2 | **10** | **10** | 1 each | 28 |
#
# By hand (`n_bins=2`: all genes in one bin; `ctrl_size=50` > 8: controls = every gene not in the set):
# - **c1:** every gene is `log1p(1000)`, so set mean − control mean = **0** for all three scores.
# - **c2, IFN-α:** set genes `s = log1p(10/28 × 10⁴) = log1p(3571.43)`, all 8 controls `c = log1p(1/28 × 10⁴) = log1p(357.14)`.
#   Score = s − c = ln(3572.43 / 358.14) = **+2.3001** (just under ln 10 = 2.3026, because of the +1 in log1p).
# - **c2, IFN-γ {C1, C2}:** set mean = c; controls are S1, S2, C3…C8 → mean (2s + 6c)/8. Score = c − (2s + 6c)/8 =
#   −(s − c)/4 = **−0.5750**. Same value for no-shared.
#
# A sign flip, scoring raw counts, or a log base other than e would each give a different number.

# %%
genes3 = ["S1", "S2"] + [f"C{i}" for i in range(1, 9)]
ad3 = toy([[1] * 10, [10, 10] + [1] * 8], genes3)
s3 = m3.score_cells(ad3, {IFNA: ["S1", "S2"], IFNG: ["C1", "C2"]},
                    target_sum=1e4, ctrl_size=50, n_bins=2, random_state=0)

s = np.log1p(10 / 28 * 1e4)
c = np.log1p(1 / 28 * 1e4)
expected = {("c1", "ifna_score"): 0.0, ("c1", "ifng_score"): 0.0, ("c1", "ifng_noshared_score"): 0.0,
            ("c2", "ifna_score"): s - c,              # +2.3001
            ("c2", "ifng_score"): -(s - c) / 4,       # -0.5750
            ("c2", "ifng_noshared_score"): -(s - c) / 4}
compare(pd.DataFrame({
    "input": [f"{cell}: {ad3.X[int(cell[1]) - 1].astype(int).tolist()}" for cell, _ in expected],
    "score": [col for _, col in expected],
    "expected": [round(v, 4) for v in expected.values()],
    "code gives": [round(float(s3.loc[cell, col]), 4) for cell, col in expected],
}), "3 score oracle", tol=1e-3)

# %% [markdown]
# ## 4. IFN-γ without the shared genes (`score_cells`)
#
# **Line:** `src/m3_signature.py:87` (brief §4).
#
# **Should do:** the no-shared IFN-γ score uses IFN-γ **minus** the genes also in IFN-α (contract step 3, P6).
#
# Toy sets (test plan T4): IFN-α = {A, S}, IFN-γ = {G, S}. Shared = {S}, so no-shared IFN-γ = **{G}** only.
# 10 genes: A, S, G, C1…C7.
#
# | | A | S | G | C1…C7 | total |
# |---|---|---|---|---|---|
# | c1 | 1 | 1 | 1 | 1 each | 10 |
# | c2 | 10 | 10 | **1** | 10 each | 91 |
#
# In c2, G is the only gene that differs. Let `g = log1p(1/91 × 10⁴)`, `c = log1p(10/91 × 10⁴)` (every other gene), so g − c = ln(110.89 / 1099.90) = −2.2944.
# - **No-shared {G}:** g − c = **−2.2944**. If S were wrongly kept, it would equal the full score (−1.1472); if the wrong
#   set were subtracted (IFN-α minus IFN-γ = {A}), it would be 0.
# - **Full IFN-γ {G, S}:** set mean (g + c)/2, controls all = c → (g − c)/2 = **−1.1472**.
# - **IFN-α {A, S}:** set mean c; controls G, C1…C7 → (g + 7c)/8. Score = (c − g)/8 = **+0.2868**.
# - c1: all genes equal → all **0**.

# %%
genes4 = ["A", "S", "G"] + [f"C{i}" for i in range(1, 8)]
ad4 = toy([[1] * 10, [10, 10, 1] + [10] * 7], genes4)
s4 = m3.score_cells(ad4, {IFNA: ["A", "S"], IFNG: ["G", "S"]},
                    target_sum=1e4, ctrl_size=50, n_bins=2, random_state=0)

g = np.log1p(1 / 91 * 1e4)
c = np.log1p(10 / 91 * 1e4)
expected = {("c1", "ifna_score"): 0.0, ("c1", "ifng_score"): 0.0, ("c1", "ifng_noshared_score"): 0.0,
            ("c2", "ifna_score"): (c - g) / 8,          # +0.2868
            ("c2", "ifng_score"): (g - c) / 2,          # -1.1472
            ("c2", "ifng_noshared_score"): g - c}       # -2.2944
compare(pd.DataFrame({
    "input": [f"{cell}: {ad4.X[int(cell[1]) - 1].astype(int).tolist()}" for cell, _ in expected],
    "score": [col for _, col in expected],
    "expected": [round(v, 4) for v in expected.values()],
    "code gives": [round(float(s4.loc[cell, col]), 4) for cell, col in expected],
}), "4 no-shared IFN-gamma", tol=1e-3)

# %% [markdown]
# ## 5. Per-sample score is the mean (`summarise_samples`, `paired_differences`)
#
# **Line:** `src/m3_signature.py:118` (brief §4).
#
# **Should do:** one row per donor × cell type × condition with `n_cells` and the **mean** of the cell scores (B3); the paired difference is stim − ctrl.
#
# | donor | cell type | condition | cell scores | mean | median |
# |---|---|---|---|---|---|
# | A | B cells | ctrl | 1, 2, 9 | **4** | 2 |
# | A | B cells | stim | 3, 5, 7 | **5** | 5 |
#
# Paired stim − ctrl = 5 − 4 = **1** (a median would give 3; a reversed difference would give −1).
# `min_cells=3` is passed so 3 cells per group are enough (the real value, 10, is checked in section 6).

# %%
obs5 = pd.DataFrame({"donor": "A", "cell_type": "B cells", "condition": ["ctrl"] * 3 + ["stim"] * 3},
                    index=[f"c{i}" for i in range(1, 7)])
vals = [1, 2, 9, 3, 5, 7]
cs5 = pd.DataFrame({col: [float(v) for v in vals] for col in m3.SCORE_COLUMNS}, index=obs5.index)

ss5 = m3.summarise_samples(cs5, obs5, min_cells=3).set_index("condition")
pd5 = m3.paired_differences(ss5.reset_index(), control="ctrl", treatment="stim")

# B3: "per-sample score = mean of cell scores"; docstring: "treatment - control".
rows = []
for col in m3.SCORE_COLUMNS:
    rows += [{"input": "ctrl cells 1, 2, 9", "what": f"mean {col}", "expected": 4.0, "code gives": float(ss5.loc["ctrl", col])},
             {"input": "stim cells 3, 5, 7", "what": f"mean {col}", "expected": 5.0, "code gives": float(ss5.loc["stim", col])},
             {"input": "stim − ctrl", "what": col, "expected": 1.0, "code gives": float(pd5[col].iloc[0])}]
rows += [{"input": "ctrl", "what": "n_cells", "expected": 3, "code gives": int(ss5.loc["ctrl", "n_cells"])},
         {"input": "stim", "what": "n_cells", "expected": 3, "code gives": int(ss5.loc["stim", "n_cells"])}]
compare(pd.DataFrame(rows), "5 mean per sample")

# %% [markdown]
# ## 6. Paired 10-cell rule at the boundary (`summarise_samples`)
#
# **Line:** `src/m3_signature.py:123` (brief §4).
#
# **Should do:** keep a donor in a cell type only if **both** conditions have **≥ 10** cells; otherwise drop both of that donor's rows **for that cell type only** (A1; a missing condition counts as 0 cells).
#
# `min_cells=10` passed explicitly.
#
# | cell type | donor | ctrl cells | stim cells | expected |
# |---|---|---|---|---|
# | B cells | A | 9 | 12 | both rows dropped (9 < 10) |
# | B cells | B | 10 | 10 | both kept (10 is not < 10) |
# | B cells | C | 10 | 0 | ctrl row dropped (stim missing = 0) |
# | NK cells | A | 10 | 10 | both kept: A is dropped from B cells only |
#
# A `> 10` instead of `≥ 10` would drop donor B; per-condition (unpaired) filtering would keep A/stim and C/ctrl.

# %%
groups = [("B cells", "A", "ctrl", 9), ("B cells", "A", "stim", 12), ("B cells", "B", "ctrl", 10),
          ("B cells", "B", "stim", 10), ("B cells", "C", "ctrl", 10), ("NK cells", "A", "ctrl", 10),
          ("NK cells", "A", "stim", 10)]
obs6 = pd.DataFrame([(t, d, c) for t, d, c, n in groups for _ in range(n)], columns=["cell_type", "donor", "condition"])
obs6.index = [f"c{i}" for i in range(len(obs6))]
cs6 = pd.DataFrame({col: 0.5 for col in m3.SCORE_COLUMNS}, index=obs6.index)

# A1: < 10 cells in any condition drops the donor from that cell type (both conditions).
expected = {("B cells", "A", "ctrl"): "dropped", ("B cells", "A", "stim"): "dropped",
            ("B cells", "B", "ctrl"): "kept, n=10", ("B cells", "B", "stim"): "kept, n=10",
            ("B cells", "C", "ctrl"): "dropped",
            ("NK cells", "A", "ctrl"): "kept, n=10", ("NK cells", "A", "stim"): "kept, n=10"}

ss6 = m3.summarise_samples(cs6, obs6, min_cells=10)
got = {(r.cell_type, r.donor, r.condition): f"kept, n={r.n_cells}" for r in ss6.itertuples()}
compare(pd.DataFrame({
    "input": [f"{t} / donor {d} / {c}: {n} cells" for t, d, c, n in groups],
    "expected": [expected[(t, d, c)] for t, d, c, _ in groups],
    "code gives": [got.get((t, d, c), "dropped") for t, d, c, _ in groups],
}), "6 paired 10-cell rule")
assert set(got) <= set(expected), f"unexpected extra rows: {set(got) - set(expected)}"

# %% [markdown]
# ## Summary

# %%
print(f"{len(FAILURES)} section(s) with a POSSIBLE BUG: {FAILURES}" if FAILURES else "All examples match.")
assert not FAILURES, FAILURES
