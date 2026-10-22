# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M2a explained by example
#
# Written by a reviewer who did not write `src/m2a_calls.py` (only its signatures and docstrings were read). Each
# section states in one sentence what a key rule should do and runs the **real** code on a tiny input. The `expected`
# value is written first, from `audit/contracts/M2_shared_specific.md` (M2a), `audit/decisions.md` (M2 section,
# Differential expression, Cell-type specificity) and `audit/test_plans/M2a.md`, never from the code. Each table
# shows **input | expected | code gives**.
#
# A mismatch is a **POSSIBLE BUG**: it is recorded, the notebook keeps going, and the last cell fails listing all
# mismatches. Expected values are never edited to match.
#
# Toy cell notation: `↑ 3.0 [2, 4]` = DE up, log2FC 3, 95% CI [2, 4]; `– 0.8 [-0.3, 1.9]` = not DE;
# `n/t` = not tested (baseMean ≤ 1, no q-value).

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from IPython.display import display

from src import m2a_calls as m2a

PROC = ROOT / "data/processed"
CT = m2a.WELL_POWERED + m2a.LOW_POWER
FAILURES = []


def compare(table, name, expected_col="expected", actual_col="code gives"):
    """Show the table; record (not stop on) any row where expected differs from what the code gives."""
    table["match"] = [e == a for e, a in zip(table[expected_col], table[actual_col])]
    display(table)
    if not table["match"].all():
        print(f"POSSIBLE BUG in {name}:")
        display(table[~table["match"]])
        FAILURES.append(name)


# Cell states: (qvalue, log2FC, ci_low, ci_high, tested)
UP = (0.001, 3.0, 2.0, 4.0, True)
DOWN = (0.001, -3.0, -4.0, -2.0, True)
FLAT = (0.5, 0.1, -0.5, 0.7, True)                 # tested, not DE, CI inside (-1, 1)
WIDE_UP = (0.2, 0.8, -0.3, 1.9, True)              # tested, not DE, CI reaches +1
UNTESTED = (np.nan, -0.85, -6.8, 5.1, False)       # like 7SK.2 in B cells: huge CI, never tested


def toy(genes):
    """genes: {gene: {cell_type: state}}; unspecified cell types are FLAT. Returns an m1b_de-shaped frame."""
    rows = []
    for g, cells in genes.items():
        for ct in CT:
            q, lfc, lo, hi, tested = cells.get(ct, FLAT)
            rows.append(dict(run="main", cell_type=ct, gene=g, baseMean=10.0 if tested else 0.07, log2FC=lfc,
                             lfcSE=0.5, ci_low=lo, ci_high=hi, pvalue=q, qvalue=q, tested=tested,
                             cooks_flag=False, log2FC_shrunk=lfc / 2, significant=False, de="none"))
    return pd.DataFrame(rows)


def show_input(genes):
    """8-cell-type view of the toy input: arrow by the DE rule, log2FC and CI."""
    def cell(s):
        q, lfc, lo, hi, tested = s
        if not tested:
            return f"n/t {lfc} [{lo}, {hi}]"
        a = "↑" if q < 0.05 and lfc >= 1 else "↓" if q < 0.05 and lfc <= -1 else "–"
        return f"{a} {lfc} [{lo}, {hi}]"
    display(pd.DataFrame({g: {ct: cell(c.get(ct, FLAT)) for ct in CT} for g, c in genes.items()}).T)


def labels(genes):
    return m2a.label_genes(m2a.derive_calls(toy(genes))).set_index("gene")


# %% [markdown]
# ## 1. One gene per label (`label_genes`)
#
# **Should do:** `shared` = DE all 8 one way; `shared*` = DE in the 6 well-powered, and a low-power cell type rescued by
# a tested CI reaching ±1 that way; `partial` = DE in ≥ 2; `specific` = DE in 1; `specific*` = DE in 1 and a tested
# non-DE cell type's CI reaches ±1; `mixed` = up somewhere and down elsewhere.

# %%
WP_UP = {ct: UP for ct in m2a.WELL_POWERED}
g1 = {
    "G_shared": {ct: UP for ct in CT},
    "G_shared_star": {**WP_UP, "Dendritic cells": UP, "Megakaryocytes": WIDE_UP},
    "G_partial": {"B cells": UP, "NK cells": UP},
    "G_specific": {"CD14+ Monocytes": UP},
    "G_specific_star": {"CD14+ Monocytes": UP, "NK cells": WIDE_UP},
    "G_mixed": {"B cells": UP, "CD14+ Monocytes": DOWN},
}
show_input(g1)
got1 = labels(g1)
display(got1)
# expected: (shared_response, n_cell_types, megakaryocyte_lower_confidence = DE in megakaryocytes)
exp1 = {
    "G_shared": ("shared", 8, True),
    "G_shared_star": ("shared*", 7, False),
    "G_partial": ("partial", 2, False),
    "G_specific": ("specific", 1, False),
    "G_specific_star": ("specific*", 1, False),
    "G_mixed": ("mixed", 2, False),
}
compare(pd.DataFrame({
    "input: gene": list(exp1),
    "expected (label, n_cell_types, Mk flag)": list(exp1.values()),
    "code gives": [(got1.loc[g, "shared_response"], int(got1.loc[g, "n_cell_types"]),
                    bool(got1.loc[g, "megakaryocyte_lower_confidence"])) for g in exp1],
}), "1 labels", expected_col="expected (label, n_cell_types, Mk flag)")

exp_arrows = {"G_shared": "↑" * 8, "G_shared_star": "↑" * 7 + "–", "G_partial": "↑––––↑––",
              "G_specific": "–↑––––––", "G_specific_star": "–↑––––––", "G_mixed": "↑↓––––––"}
exp_star = {"G_specific": False, "G_specific_star": True}
compare(pd.DataFrame({
    "input: gene": list(exp_arrows) + [f"{g} specific_star" for g in exp_star],
    "expected": list(exp_arrows.values()) + list(exp_star.values()),
    "code gives": ["".join(got1.loc[g, CT]) for g in exp_arrows] + [bool(got1.loc[g, "specific_star"])
                                                                    for g in exp_star],
}), "1 arrows and star")

# %% [markdown]
# ## 2. Precedence: up in the 6 well-powered but DOWN in DCs (`label_genes`)
#
# **Should do:** any opposite-direction DE call makes the gene `mixed`, even when it would otherwise be shared.

# %%
g2 = {"G_up6_downDC": {**WP_UP, "Dendritic cells": DOWN, "Megakaryocytes": UP}}
show_input(g2)
got2 = labels(g2)
compare(pd.DataFrame({
    "input: gene": ["G_up6_downDC"],
    "expected (label, n_cell_types)": [("mixed", 8)],
    "code gives": [(got2.loc["G_up6_downDC", "shared_response"], int(got2.loc["G_up6_downDC", "n_cell_types"]))],
}), "2 precedence", expected_col="expected (label, n_cell_types)")

# %% [markdown]
# ## 3. CI rules use only tested cell types (`label_genes`)
#
# **Should do:** a wide CI from a cell type where the gene was not tested (baseMean ≤ 1) neither stars a specific gene
# nor rescues a shared one; the same CI from a tested cell type does.

# %%
WIDE_UP_UNTESTED = (np.nan, 0.8, -0.3, 1.9, False)
g3 = {
    "G_spec_wide_untested": {"CD14+ Monocytes": UP, "Megakaryocytes": UNTESTED},
    "G_spec_wide_tested": {"CD14+ Monocytes": UP, "Megakaryocytes": (0.6, -0.85, -6.8, 5.1, True)},
    "G_rescue_untested": {**WP_UP, "Dendritic cells": UP, "Megakaryocytes": WIDE_UP_UNTESTED},
    "G_rescue_tested": {**WP_UP, "Dendritic cells": UP, "Megakaryocytes": WIDE_UP},
}
show_input(g3)
got3 = labels(g3)
exp3 = {"G_spec_wide_untested": ("specific", False), "G_spec_wide_tested": ("specific*", True),
        "G_rescue_untested": ("partial", False), "G_rescue_tested": ("shared*", False)}
compare(pd.DataFrame({
    "input: gene": list(exp3),
    "expected label": [e[0] for e in exp3.values()],
    "code gives": [got3.loc[g, "shared_response"] for g in exp3],
}), "3 tested-only", expected_col="expected label")
compare(pd.DataFrame({
    "input: gene": ["G_spec_wide_untested", "G_spec_wide_tested"],
    "expected specific_star": [False, True],
    "code gives": [bool(got3.loc[g, "specific_star"]) for g in ["G_spec_wide_untested", "G_spec_wide_tested"]],
}), "3 specific_star", expected_col="expected specific_star")

# %% [markdown]
# ## 4. The DE boundary (`derive_calls`)
#
# **Should do:** DE = q < 0.05 (strict) and |log2FC| ≥ 1 (inclusive, changed by the user 2026-10-02); the input
# `de` column is ignored and re-derived.

# %%
rows4 = [  # cell type, (q, log2FC), expected de
    ("B cells", (0.049, 1.0), "up"),
    ("CD14+ Monocytes", (0.05, 1.0), "none"),       # q exactly 0.05 is not < 0.05
    ("CD4 T cells", (0.049, -1.0), "down"),
    ("CD8 T cells", (0.049, 0.9999), "none"),
    ("FCGR3A+ Monocytes", (0.05, 2.0), "none"),
]
d4 = toy({"G_boundary": {ct: (q, l, l - 0.5, l + 0.5, True) for ct, (q, l), _ in rows4}})
d4["de"] = "up"  # a wrong input `de` everywhere: must not be trusted
got4 = m2a.derive_calls(d4).set_index("cell_type")
compare(pd.DataFrame({
    "input: (q, log2FC)": [r[1] for r in rows4],
    "expected de": [r[2] for r in rows4],
    "code gives": [got4.loc[r[0], "de"] for r in rows4],
}), "4 DE boundary", expected_col="expected de")

# %% [markdown]
# ## 5. LODO stability (`lodo_stability`)
#
# **Should do:** a DE call is `lodo_stable` if every leave-one-donor-out refit keeps the sign (significance not
# needed, 0 counts as a change); a run that is not exactly `lodo_<paired donor>` stops the run.

# %%
donors = {ct: [f"d{i}" for i in range(1, 9)] for ct in m2a.WELL_POWERED}
donors["Dendritic cells"] = [f"d{i}" for i in range(1, 7)]
donors["Megakaryocytes"] = [f"d{i}" for i in range(1, 5)]
calls5 = m2a.derive_calls(toy({"G_tiny_pos": {"B cells": UP}, "G_zero": {"B cells": UP},
                               "G_down": {"NK cells": DOWN}}))
special = {("G_tiny_pos", "B cells", "d3"): 0.001, ("G_zero", "B cells", "d3"): 0.0}
lodo = pd.DataFrame([dict(run=f"lodo_{d}", cell_type=ct, gene=g,
                          log2FC=special.get((g, ct, d), -2.5 if g == "G_down" else 2.5))
                     for ct, ds in donors.items() for d in ds for g in calls5["gene"].unique()])
got5 = m2a.lodo_stability(calls5, lodo, donors).set_index(["gene", "cell_type"])
cases5 = [("G_tiny_pos", "B cells", "up; refits +2.5 ×7, +0.001", True),
          ("G_zero", "B cells", "up; refits +2.5 ×7, 0.0", False),
          ("G_down", "NK cells", "down; refits −2.5 ×8", True),
          ("G_tiny_pos", "CD4 T cells", "not DE", "NA")]
compare(pd.DataFrame({
    "input": [c[2] for c in cases5],
    "expected lodo_stable": [c[3] for c in cases5],
    "code gives": ["NA" if pd.isna(v) else bool(v) for v in
                   (got5.loc[(c[0], c[1]), "lodo_stable"] for c in cases5)],
}), "5 LODO sign", expected_col="expected lodo_stable")

bad = lodo.copy()
bad.loc[(bad["cell_type"] == "B cells") & (bad["run"] == "lodo_d3"), "run"] = "lodo_NOT_A_DONOR"
try:
    m2a.lodo_stability(calls5, bad, donors)
    got_err = "no error"
except ValueError as e:
    got_err = "ValueError"
    print("message:", e)
compare(pd.DataFrame({
    "input": ["B cells: lodo_d3 renamed lodo_NOT_A_DONOR (count still 8)"],
    "expected": ["ValueError"],
    "code gives": [got_err],
}), "5 substituted run")

# %% [markdown]
# ## 6. Real data: labels re-counted independently (`m2_table.csv`)
#
# **Should do:** my own few lines of pandas on `m1b_de.csv`, applying the decided rules, give the same labels as the
# pipeline's table, gene by gene, and the counts recorded in `audit/decisions.md`.

# %%
de = pd.read_csv(PROC / "m1b_de.csv", usecols=["cell_type", "gene", "log2FC", "qvalue", "ci_low", "ci_high",
                                               "tested"])
big = (de["qvalue"] < 0.05) & (de["log2FC"].abs() >= 1)
de["call"] = np.where(big & (de["log2FC"] > 0), "↑", np.where(big & (de["log2FC"] < 0), "↓", "–"))
w = lambda col: de.pivot(index="gene", columns="cell_type", values=col)[CT]
call, lo, hi, tested = w("call"), w("ci_low"), w("ci_high"), w("tested").astype(bool)
up, dn = call == "↑", call == "↓"
keep = (up | dn).any(axis=1)
call, lo, hi, tested, up, dn = (x[keep] for x in (call, lo, hi, tested, up, dn))
n = (up | dn).sum(axis=1)
LP = m2a.LOW_POWER
up_ok = up[m2a.WELL_POWERED].all(axis=1) & (up[LP] | (tested[LP] & (hi[LP] >= 1))).all(axis=1)
dn_ok = dn[m2a.WELL_POWERED].all(axis=1) & (dn[LP] | (tested[LP] & (lo[LP] <= -1))).all(axis=1)
star = (~(up | dn) & tested & ((hi >= 1) | (lo <= -1))).any(axis=1)
mine = pd.Series(np.select(
    [up.any(axis=1) & dn.any(axis=1), n == 8, up_ok | dn_ok, n >= 2, star],
    ["mixed", "shared", "shared*", "partial", "specific*"], "specific"), index=call.index)

table = pd.read_csv(PROC / "m2_table.csv").set_index("gene")
order = ["shared", "shared*", "partial", "specific", "specific*", "mixed"]
decided = dict(zip(order, [21, 87, 856, 241, 1585, 11]))
compare(pd.DataFrame({
    "input: label": order,
    "expected (decisions.md, my count)": [(decided[k], decided[k]) for k in order],
    "code gives": [(decided[k], int((mine == k).sum())) for k in order],
}), "6 my count vs decisions", expected_col="expected (decisions.md, my count)")
compare(pd.DataFrame({
    "input: label": order + ["total genes"],
    "expected (my count)": [int((mine == k).sum()) for k in order] + [len(mine)],
    "code gives": [int((table["shared_response"] == k).sum()) for k in order] + [len(table)],
}), "6 table counts", expected_col="expected (my count)")

same_genes = set(mine.index) == set(table.index)
idx = mine.index.intersection(table.index)
compare(pd.DataFrame({
    "input": ["same gene set", "genes whose label differs", "gene × cell type arrows that differ"],
    "expected": [True, 0, 0],
    "code gives": [same_genes, int((mine[idx] != table.loc[idx, "shared_response"]).sum()),
                   int((call.loc[idx, CT] != table.loc[idx, CT]).to_numpy().sum())],
}), "6 gene-by-gene")

# %%
ex = table[table["shared_response"] == "shared"].head(5)
display(ex[CT + ["n_cell_types", "shared_response"]])
compare(pd.DataFrame({
    "input: shared gene": ex.index,
    "expected arrows (8 same, mine)": ["".join(call.loc[g, CT]) for g in ex.index],
    "code gives": ["".join(ex.loc[g, CT]) for g in ex.index],
}), "6 shared examples", expected_col="expected arrows (8 same, mine)")

# %% [markdown]
# ## Summary

# %%
print("Mismatches (POSSIBLE BUGS):", FAILURES or "none")
assert not FAILURES, FAILURES
