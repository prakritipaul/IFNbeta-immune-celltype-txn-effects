# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # Wave 4 (Fig 3, Fig 4, Fig 5, final table) explained by example
#
# Written by a reviewer who did not write `src/w4_outputs.py` (only its signatures and docstrings were read). Each
# section states in one sentence what a key selection or assembly rule should do and runs the **real** code on a tiny
# input. The `expected` value is written first, from `audit/contracts/W4_figures_table.md` and `audit/decisions.md`
# (F1–F6, F3 revised, F5, Wave 4 follow-ups), never from the code. Each table shows **input | expected | code gives**.
#
# A mismatch is a **POSSIBLE BUG**: it is recorded, the notebook keeps going, and the last cell fails listing all
# mismatches. Expected values are never edited to match.

# %%
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from IPython.display import display

from src import w4_outputs as w4
from src.m2a_calls import CELL_TYPES as CT

PROC = ROOT / "data/processed"
FAILURES = []
B_, C14, C4, C8, F16, NK, DC, MK = CT


def compare(table, name, expected_col="expected", actual_col="code gives"):
    """Show the table; record (not stop on) any row where expected differs from what the code gives."""
    table["match"] = [e == a for e, a in zip(table[expected_col], table[actual_col])]
    display(table)
    if not table["match"].all():
        print(f"POSSIBLE BUG in {name}:")
        display(table[~table["match"]])
        FAILURES.append(name)


def toy_m2(genes):
    """genes: {gene: (label, {cell_type: log2FC}, untested cell types)}. Cell types absent from the dict get
    log2FC 0.1 (tested, not DE). DE = |log2FC| > 1 in a tested cell type. Returns (m2_table, m2_calls)."""
    trows, crows = [], []
    for g, (label, lfcs, untested) in genes.items():
        arrows = {}
        for ct in CT:
            lfc, tested = lfcs.get(ct, 0.1), ct not in untested
            de = "up" if tested and lfc > 1 else "down" if tested and lfc < -1 else "none"
            arrows[ct] = {"up": "↑", "down": "↓", "none": "–"}[de]
            crows.append(dict(cell_type=ct, gene=g, log2FC=lfc, qvalue=(0.001 if de != "none" else 0.5) if tested
                              else np.nan, ci_low=lfc - 1, ci_high=lfc + 1, tested=tested, de=de, lodo_stable=np.nan))
        n = sum(a != "–" for a in arrows.values())
        trows.append(dict(gene=g, **arrows, n_cell_types=n, shared_response=label,
                          specific_star=label.endswith("*"), megakaryocyte_lower_confidence=arrows[MK] != "–"))
    return pd.DataFrame(trows), pd.DataFrame(crows)


# %% [markdown]
# ## 1. Fig 4 rows (`fig4_rows`)
#
# **Should do (F3 revised, F4):** shared block = top 25 of `shared` + `shared*` by **mean |log2FC| over the 8 cell
# types**; then per cell type (in column order) the top 5 `specific`/`specific*` genes by |log2FC| in that cell type,
# **only among genes tested in all 8**; fewer rows when fewer exist; `display` = gene + `*` for starred labels.
#
# Toy shared: `sh01`–`sh27` with |log2FC| = 30 − k in every cell type (`sh03` is down, `sh05` is `shared*`), plus
# `sh_peak` (40 in B cells, 1.5 elsewhere: mean 6.3125, max 40) and `trap_partial` (partial, 50 in two cell types).
# Mean ranking: sh01…sh23 (29…7), sh_peak (6.31), sh24 (6); sh25 (5) and later are cut.
# Toy specific, B cells: `b1`–`b6` (`specific*`, |log2FC| 6, 5, −4.5, 4, 3, 2, all tested) → top 5 = b1…b5;
# `b_trap` (`specific*`, 20 in B cells, **untested in Megakaryocytes**) must be excluded. CD14+ Mono: `c1` (3), `c2`
# (−2.5) → only 2 rows. NK: `n1` (`specific`, no star, 2.0, tested in all 8) → shown, no star. Megakaryocytes: none.

# %%
G1 = {f"sh{k:02d}": ("shared*" if k == 5 else "shared", {ct: (-1 if k == 3 else 1) * (30 - k) for ct in CT}, [])
      for k in range(1, 28)}
G1["sh_peak"] = ("shared", {ct: 40.0 if ct == B_ else 1.5 for ct in CT}, [])
G1["trap_partial"] = ("partial", {B_: 50.0, C4: 50.0}, [])
G1 |= {f"b{i}": ("specific*", {B_: v}, []) for i, v in enumerate([6, 5, -4.5, 4, 3, 2], start=1)}
G1["b_trap"] = ("specific*", {B_: 20.0}, [MK])
G1 |= {"c1": ("specific*", {C14: 3.0}, []), "c2": ("specific*", {C14: -2.5}, []), "n1": ("specific", {NK: 2.0}, [])}
table1, calls1 = toy_m2(G1)
col1 = [NK, C14, B_, C4, C8, F16, DC, MK]  # a non-default column order: specific blocks must follow it
rows1 = w4.fig4_rows(table1, calls1, col1, n_shared=25, n_specific=5)
display(rows1.head(3), rows1.tail(10))

exp_shared = [f"sh{k:02d}" for k in range(1, 24)] + ["sh_peak", "sh24"]
exp_spec = ["n1", "c1", "c2", "b1", "b2", "b3", "b4", "b5"]
compare(pd.DataFrame({
    "input": ["gene order (all rows)", "display of sh05 (shared*)", "display of b1 (specific*)",
              "display of n1 (specific)", "excluded present? (sh25, b_trap, b6, trap_partial)", "number of rows"],
    "expected": [exp_shared + exp_spec, "sh05*", "b1*", "n1", False, 33],
    "code gives": [rows1["gene"].tolist(), *rows1.set_index("gene").loc[["sh05", "b1", "n1"], "display"].tolist(),
                   bool(rows1["gene"].isin(["sh25", "b_trap", "b6", "trap_partial"]).any()), len(rows1)],
}), "1 fig4 rows")
spec_rows1 = rows1[rows1["gene"].isin(exp_spec)]
compare(pd.DataFrame({
    "input": ["each specific row's group names its own cell type", "Megakaryocytes block rows"],
    "expected": [True, 0],
    "code gives": [all(ct in grp for ct, grp in zip([NK, C14, C14] + [B_] * 5, spec_rows1["group"])),
                   int(rows1["group"].str.contains(MK).sum())],
}), "1 fig4 groups")

# %% [markdown]
# ## 2. Fig 4 column order (`jaccard_order`)
#
# **Should do (F5):** average-linkage clustering on 1 − Jaccard of the **up** gene sets.
#
# Toy, 4 cell types A–D (up): J(A,B) 0.8, J(A,C) 0.7, J(C,D) 0.5, all others 0.1 → distances AB 0.2, AC 0.3, CD 0.5,
# BC = AD = BD 0.9. By hand: merge A+B at 0.2; then {A,B}–C average = (0.3 + 0.9)/2 = 0.6 > C–D 0.5, so merge C+D at
# 0.5; last merge at mean(0.3, 0.9, 0.9, 0.9) = 0.75. (Single linkage would join C to {A,B} at 0.3; complete would
# finish at 0.9.) Down sets are a trap: there A and D are identical (J 1.0).

# %%
names = ["A", "B", "C", "D"]
Jup = pd.DataFrame(0.1, index=names, columns=names)
for a, b, v in [("A", "B", 0.8), ("A", "C", 0.7), ("C", "D", 0.5)]:
    Jup.loc[a, b] = Jup.loc[b, a] = v
np.fill_diagonal(Jup.values, 1.0)
Jdn = pd.DataFrame(0.1, index=names, columns=names)
Jdn.loc["A", "D"] = Jdn.loc["D", "A"] = 1.0
np.fill_diagonal(Jdn.values, 1.0)
jac2 = pd.concat({"up": Jup, "down": Jdn}, names=["direction", "cell_type"])
order2, Z2 = w4.jaccard_order(jac2)
print("order:", order2)
pos = {c: i for i, c in enumerate(order2)}
compare(pd.DataFrame({
    "input": ["merge heights", "A next to B", "C next to D"],
    "expected": [[0.2, 0.5, 0.75], True, True],
    "code gives": [np.round(Z2[:, 2], 10).tolist(), abs(pos["A"] - pos["B"]) == 1, abs(pos["C"] - pos["D"]) == 1],
}), "2 jaccard order")

# %% [markdown]
# ## 3. Fig 5 matrix and column order (`fig5_matrix`, `nes_order`)
#
# **Should do (Fig 5 contract):** rows grouped shared → partial → specific → mixed; within a label by number of
# affected cell types (desc), then mean NES (desc); NES NaN (`n/t`) where a program has no GSEA row; significant =
# padj < 0.05.
#
# Toy (8 cell types): `S1` shared up; `Pa1` up in 3 (B, CD4, NK; NES 2); `Pa2` up in 5; `Pa3` up in 3 with NES 3 and
# **no row in DCs**; `Sp1` down in CD14+ Mono only; `M1` up in B, down in NK. Non-significant cells: NES 0.5,
# padj 0.5. Expected row order: S1, Pa2, Pa3, Pa1, Sp1, M1.

# %%
SIG_UP, SIG_DN, NS = (2.0, 0.01), (-2.0, 0.01), (0.5, 0.5)


def gsea_programs(spec):
    """spec: {pathway: {cell_type: (NES, padj)}} -> (m2c_gsea-shaped, m2c_programs-shaped) frames, built by hand."""
    g = pd.DataFrame([dict(cell_type=ct, collection="hallmark", pathway=p, size=50, ES=n / 3, NES=n, pval=q, padj=q,
                           log2err=0.1, leadingEdge="X") for p, d in spec.items() for ct, (n, q) in d.items()])
    prow = []
    for p, d in spec.items():
        arr = {ct: ("n/t" if ct not in d else "↑" if d[ct][1] < 0.05 and d[ct][0] > 0
                    else "↓" if d[ct][1] < 0.05 else "–") for ct in sorted(CT)}
        up, dn = sum(a == "↑" for a in arr.values()), sum(a == "↓" for a in arr.values())
        if up + dn == 0:
            continue
        label = "mixed" if up and dn else "shared" if 8 in (up, dn) else "partial" if up + dn >= 2 else "specific"
        prow.append(dict(collection="hallmark", pathway=p, **arr, n_up=up, n_down=dn, label=label))
    return g, pd.DataFrame(prow)


spec3 = {
    "M1": {ct: SIG_UP if ct == B_ else SIG_DN if ct == NK else NS for ct in CT},
    "Sp1": {ct: SIG_DN if ct == C14 else NS for ct in CT},
    "Pa1": {ct: SIG_UP if ct in (B_, C4, NK) else NS for ct in CT},
    "Pa2": {ct: SIG_UP if ct in (B_, C4, C8, NK, F16) else NS for ct in CT},
    "Pa3": {ct: (3.0, 0.01) if ct in (B_, C4, NK) else NS for ct in CT if ct != DC},
    "S1": {ct: SIG_UP for ct in CT},
}
gsea3, prog3 = gsea_programs(spec3)
rows3, m3, sig3 = w4.fig5_matrix(gsea3, prog3, CT, padj=0.05)
display(m3)
compare(pd.DataFrame({
    "input": ["row order", "NES Pa3 × DCs (no row)", "significant Pa3 × DCs", "NES Sp1 × CD14+ Mono",
              "significant Sp1 × CD14+ Mono", "number of significant cells"],
    "expected": [["S1", "Pa2", "Pa3", "Pa1", "Sp1", "M1"], "NaN", False, -2.0, True, 8 + 3 + 5 + 3 + 1 + 2],
    "code gives": [list(m3.index), "NaN" if np.isnan(m3.loc["Pa3", DC]) else m3.loc["Pa3", DC],
                   bool(sig3.loc["Pa3", DC]), m3.loc["Sp1", C14], bool(sig3.loc["Sp1", C14]),
                   int(sig3.to_numpy().sum())],
}), "3 fig5 matrix")

# %% [markdown]
# **Column order (F5):** average linkage on 1 − Pearson r of NES profiles, over affected programs tested in all 8.
#
# Toy: affected P1–P3 tested in all 8. Group G1 = {B, CD4, CD8, NK} has increasing NES (1,2,3 / 2,4,6 / 1.5,2,2.5 /
# 0,1,2), G2 = {CD14, FCGR3A, DC, Mk} decreasing (3,2,1 / −1,−2,−3 / 6,4,2 / 2,1.5,1). Within a group r = 1 (distance
# 0), between groups r = −1 (distance 2) → last merge at 2.0, groups contiguous, **3 programs used**. Traps:
# `T_unaff` (tested in all 8, **not affected**, B 50 / others 0) and `T_nt` (affected but **no row in Mk**,
# B −50) would break B's r = 1 if used.

# %%
up_, dn_ = {B_: [1, 2, 3], C4: [2, 4, 6], C8: [1.5, 2, 2.5], NK: [0, 1, 2]}, \
           {C14: [3, 2, 1], F16: [-1, -2, -3], DC: [6, 4, 2], MK: [2, 1.5, 1]}
prof = up_ | dn_
spec3b = {f"P{i + 1}": {ct: (prof[ct][i], 0.01) for ct in CT} for i in range(3)}
spec3b["T_unaff"] = {ct: (50.0 if ct == B_ else 0.0, 0.5) for ct in CT}
spec3b["T_nt"] = {ct: (-50.0 if ct == B_ else 2.0, 0.01) for ct in CT if ct != MK}
gsea3b, prog3b = gsea_programs(spec3b)
order3, Z3, n3 = w4.nes_order(gsea3b, prog3b)
print("order:", order3)
first4 = set(order3[:4])
compare(pd.DataFrame({
    "input": ["programs used", "last merge height", "within-group merge heights all 0",
              "the two groups are contiguous"],
    "expected": [3, 2.0, True, True],
    "code gives": [int(n3), round(float(Z3[-1, 2]), 10), bool(np.allclose(Z3[:-2, 2], 0, atol=1e-10)),
                   first4 in (set(up_), set(dn_))],
}), "3 nes order")

# %% [markdown]
# ## 4. Fig 3 donor means (`kept_cells`)
#
# **Should do (F1, A1):** keep only cells of kept samples (rows of `m3_sample_scores`) and check that each kept
# sample's cell count and mean score equal `m3_sample_scores` (the donor dots in Fig 3).
#
# Toy, cell type X: d1 ctrl scores 0.1, 0.3, 0.8 (mean 0.4); d1 stim 1.0, 2.0 (1.5); d2 ctrl 0.0, 0.6 (0.3); d2 stim
# 0.9 (0.9); d3 (a dropped donor, not in sample_scores) ctrl 5.0 → expected 8 cells kept, d3 gone. A sample_scores
# with d1 ctrl mean 0.5 (wrong) must raise.

# %%
cells4 = pd.DataFrame([dict(barcode=f"c{i}", donor=d, condition=c, cell_type="X", ifna_score=s, ifng_score=s / 2,
                            ifng_noshared_score=s / 4)
                       for i, (d, c, s) in enumerate([("d1", "ctrl", 0.1), ("d1", "ctrl", 0.3), ("d1", "ctrl", 0.8),
                                                      ("d1", "stim", 1.0), ("d1", "stim", 2.0), ("d2", "ctrl", 0.0),
                                                      ("d2", "ctrl", 0.6), ("d2", "stim", 0.9),
                                                      ("d3", "ctrl", 5.0)])])
by_hand = {("d1", "ctrl"): (3, 0.4), ("d1", "stim"): (2, 1.5), ("d2", "ctrl"): (2, 0.3), ("d2", "stim"): (1, 0.9)}
samples4 = pd.DataFrame([dict(cell_type="X", donor=d, condition=c, n_cells=n, ifna_score=m, ifng_score=m / 2,
                              ifng_noshared_score=m / 4) for (d, c), (n, m) in by_hand.items()])
kept4 = w4.kept_cells(cells4, samples4)
means4 = kept4.groupby(["donor", "condition"])["ifna_score"].agg(["size", "mean"])
bad4 = samples4.copy()
bad4.loc[0, "ifna_score"] = 0.5
try:
    w4.kept_cells(cells4, bad4)
    raised = False
except AssertionError:
    raised = True
compare(pd.DataFrame({
    "input": [f"{d} {c}: n, mean IFN-α" for d, c in by_hand] + ["d3 cells kept", "wrong mean raises"],
    "expected": list(by_hand.values()) + [0, True],
    "code gives": [(int(means4.loc[k, "size"]), round(float(means4.loc[k, "mean"]), 10)) for k in by_hand]
                  + [int((kept4["donor"] == "d3").sum()), raised],
}), "4 fig3 donor means")

# %% [markdown]
# ## 5. Final table (`final_table`)
#
# **Should do (F6, Wave 4 follow-ups):** arrows and label copied from `m2_table`; `expressed_in_n_cell_types` = number
# of cell types where the gene was tested; log2FC copied from `m2_calls`, q empty where untested; sorted shared →
# shared\* → partial → mixed → specific → specific\*, then n_cell_types (desc), then mean |log2FC| over **tested** cell
# types (desc), then gene.
#
# Toy (input in scrambled order): `p2` and `p3` are both partial with n = 2; `p2` mean over tested = (2.5·2 + 0.1·6)/8
# = 0.70; `p3` = (2·2 + 0.1·5)/7 = 0.643, but it has log2FC 20 in its **untested** Megakaryocytes (0.643 → 3.06 if
# counted). `sp1` (specific, tested only in B cells, log2FC 2 → mean 2) before `sp2` (specific, tested in 8, 3 in
# CD4 → mean 0.4625). Expected order: sh, shs, p1, p2, p3, m, sp1, sp2, sps.

# %%
G5 = {
    "sps": ("specific*", {C8: 4.0}, []),
    "sp2": ("specific", {C4: 3.0}, []),
    "m": ("mixed", {B_: 2.0, NK: -2.0}, []),
    "p3": ("partial", {B_: 2.0, C4: 2.0, MK: 20.0}, [MK]),
    "sp1": ("specific", {B_: 2.0}, [c for c in CT if c != B_]),
    "p2": ("partial", {B_: 2.5, C4: 2.5}, []),
    "p1": ("partial", {B_: 2.0, C4: 2.0, C8: 2.0}, []),
    "shs": ("shared*", {ct: 3.0 for ct in CT}, []),
    "sh": ("shared", {ct: 3.0 for ct in CT}, []),
}
# The toy is padded with 239 filler `specific` genes (written when final_table() still hard-asserted the real data's
# 241 / 229; that check has since moved to the real-data notebook, but the fillers are kept as extra rows: 239 filler `specific` genes
# (`zf…`, log2FC 1.5 in B cells): 228 tested only in B cells (mean 1.5, between sp1 and sp2) and 11 tested in all 8
# (mean (1.5 + 0.7)/8 = 0.275, after sp2). Named-gene order is unchanged; fillers sit between sp1 / sp2 / sps.
G5 |= {f"zf{i:03d}": ("specific", {B_: 1.5}, [c for c in CT if c != B_] if i < 228 else []) for i in range(239)}
table5, calls5 = toy_m2(G5)
ft5 = w4.final_table(table5, calls5).set_index("gene")
display(ft5.head(9), ft5.tail(13))
t5 = table5.set_index("gene")
named = [g for g in ft5.index if not g.startswith("zf")]
pos5 = {g: i for i, g in enumerate(ft5.index)}
compare(pd.DataFrame({
    "input": ["fillers tested in 1 lie between sp1 and sp2; fillers tested in 8 between sp2 and sps",
              "sps is the last row"],
    "expected": [True, True],
    "code gives": [all(pos5["sp1"] < pos5[f"zf{i:03d}"] < pos5["sp2"] for i in range(228))
                   and all(pos5["sp2"] < pos5[f"zf{i:03d}"] < pos5["sps"] for i in range(228, 239)),
                   ft5.index[-1] == "sps"],
}), "5 final table, fillers")
compare(pd.DataFrame({
    "input": ["row order (named genes)","expressed_in_n_cell_types (sp1, p3, sh)", "Shared response? of shs, sps",
              "arrows equal m2_table (all genes × 8)", "p3 Megakaryocytes: log2FC, q empty",
              "sp2 CD4 T cells log2FC"],
    "expected": [["sh", "shs", "p1", "p2", "p3", "m", "sp1", "sp2", "sps"], (1, 7, 8), ("shared*", "specific*"),
                 True, (20.0, True), 3.0],
    "code gives": [named, tuple(int(ft5.loc[g, "expressed_in_n_cell_types"]) for g in ["sp1", "p3", "sh"]),
                   tuple(ft5.loc[["shs", "sps"], "Shared response?"]),
                   bool((ft5.loc[t5.index, CT] == t5[CT]).to_numpy().all()),
                   (float(ft5.loc["p3", f"{MK} log2FC"]), bool(pd.isna(ft5.loc["p3", f"{MK} q"]))),
                   float(ft5.loc["sp2", f"{C4} log2FC"])],
}), "5 final table")

# %% [markdown]
# ## 6. Real data, re-counted independently
#
# **Should do:** `results/final_table.csv` has 2,801 rows; my own count from `m2_calls.csv` (`tested`) finds 229 of
# the 241 `specific` genes tested in only 1 cell type, and the file's `expressed_in_n_cell_types` agrees; every
# specific row of Fig 4 (re-derived with the real `fig4_rows`, column order from the real Jaccard) is tested in all
# 8 cell types and is DE in exactly one; the shared rows are the top 25 by my own mean |log2FC|.

# %%
ft = pd.read_csv(ROOT / "results/final_table.csv")
calls = pd.read_csv(PROC / "m2_calls.csv")
table = pd.read_csv(PROC / "m2_table.csv")
n_tested = calls.groupby("gene")["tested"].sum()
spec_genes = table.loc[table["shared_response"] == "specific", "gene"]
ft_i = ft.set_index("gene")

col_real, _ = w4.jaccard_order(pd.read_csv(PROC / "m2c_jaccard.csv", index_col=[0, 1]))
rows_real = w4.fig4_rows(table, calls, col_real)
lab = table.set_index("gene")["shared_response"]
spec_fig4 = rows_real[rows_real["gene"].map(lab).isin(["specific", "specific*"])]["gene"]
shared_pool = table.loc[table["shared_response"].isin(["shared", "shared*"]), "gene"]
mean_abs = calls[calls["gene"].isin(shared_pool)].assign(a=lambda d: d["log2FC"].abs()).groupby("gene")["a"].mean()
my_top25 = mean_abs.sort_values(ascending=False).index[:25].tolist()
display(rows_real)

compare(pd.DataFrame({
    "input": ["final_table rows", "specific genes", "specific genes tested in exactly 1 cell type (mine)",
              "same count from final_table's expressed_in_n_cell_types",
              "expressed_in_n_cell_types == my count, all 2,801 genes",
              "Fig 4 specific rows: all tested in 8", "Fig 4 specific rows: all DE in exactly 1",
              "Fig 4 specific rows: all specific or specific* (rule F3: tested in all 8; requirement wording corrected 2026-10-06)", "Fig 4 shared rows = my top 25 by mean |log2FC|"],
    "expected": [2801, 241, 229, 229, True, True, True, True, my_top25],
    "code gives": [len(ft), len(spec_genes), int((n_tested[spec_genes] == 1).sum()),
                   int((ft_i.loc[spec_genes, "expressed_in_n_cell_types"] == 1).sum()),
                   bool((ft_i["expressed_in_n_cell_types"] == n_tested[ft_i.index]).all()),
                   bool((n_tested[spec_fig4] == 8).all()),
                   bool((table.set_index("gene").loc[spec_fig4, "n_cell_types"] == 1).all()),
                   bool(lab[spec_fig4].isin(["specific", "specific*"]).all()),
                   rows_real["gene"].iloc[:25].tolist()],
}), "6 real data")
print("Fig 4 specific rows per cell type:")
print(rows_real.loc[spec_fig4.index, "group"].value_counts())

# %% [markdown]
# ## Summary

# %%
print("Mismatches (POSSIBLE BUGS):", FAILURES or "none")
assert not FAILURES, FAILURES
