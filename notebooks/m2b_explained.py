# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M2b explained by example
#
# Written by a reviewer who did not write `src/m2b_validation.py` (only its signatures and docstrings were read).
# Each section states in one sentence what a key calculation should do and runs the **real** code on a tiny input.
# The `expected` value is written first, by hand, from `audit/validation_plan.md`, `audit/decisions.md` (M2 section:
# "M2b checks" and "M2b follow-ups" rows, O5, O7) and `audit/contracts/M2_shared_specific.md` (M2b), never from the
# code. Each table shows **input | expected | code gives**.
#
# A mismatch is a **POSSIBLE BUG**: it is recorded, the notebook keeps going, and the last cell fails listing all
# mismatches. Expected values are never edited to match.

# %%
import sys
from math import comb, log
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from IPython.display import display

from src import m2b_validation as m2b

PROC = ROOT / "data/processed"
FAILURES = []


def compare(table, name, expected_col="expected", actual_col="code gives"):
    """Show the table; record (not stop on) any row where expected differs from what the code gives."""
    table["match"] = [e == a for e, a in zip(table[expected_col], table[actual_col])]
    display(table)
    if not table["match"].all():
        print(f"POSSIBLE BUG in {name}:")
        display(table[~table["match"]])
        FAILURES.append(name)


def de_frame(cell_type, calls, tested=None):
    """Toy main-run DE table: calls = {gene: "up"/"down"/"none"}; log2FC +2 / -2 / 0 to match the call."""
    lfc = {"up": 2.0, "down": -2.0, "none": 0.0}
    df = pd.DataFrame({"run": "main", "cell_type": cell_type, "gene": list(calls),
                       "log2FC": [lfc[c] for c in calls.values()], "de": list(calls.values())})
    df["log2FC_shrunk"] = df["log2FC"]
    df["qvalue"] = np.where(df["de"] == "none", 0.5, 0.001)
    df["tested"] = True if tested is None else [tested[g] for g in df["gene"]]
    df["significant"] = df["de"] != "none"
    return df


# %% [markdown]
# ## 1. P2 / P3 Fisher enrichment and BH (`p2_p3_fisher`, `bh`)
#
# **Should do (P2):** among the genes *tested* in a cell type, a 2×2 of up / not up × IFN-α / not IFN-α; odds ratio
# = ad/bc; one-sided "greater" p = P(X ≥ observed) for the hypergeometric. Untested genes are left out.
#
# Toy: 100 tested genes, 10 of them IFN-α; 10 up-DEGs of which 5 are IFN-α; 5 down-DEGs, none IFN-α; plus 3
# untested genes (one an IFN-α gene, which must not count).
#
# - P2 table [[5, 5], [5, 85]] → OR = 5·85 / (5·5) = **17**; p = Σ_{x=5..10} C(10,x)·C(90,10−x) / C(100,10).
# - P3 table (up vs down) [[5, 5], [0, 5]] → OR = 5·5 / (5·0) = **inf**; p = C(5,5)·C(10,5) / C(15,10) =
#   252 / 3003 = **0.0839**.

# %%
calls = {f"ifn{i}": ("up" if i < 5 else "none") for i in range(10)}
calls |= {f"g{i}": ("up" if i < 5 else "down" if i < 10 else "none") for i in range(90)}
calls |= {"untested_ifn": "none", "untested1": "none", "untested2": "none"}
tested = {g: not g.startswith("untested") for g in calls}
ifna_toy = {f"ifn{i}" for i in range(10)} | {"untested_ifn"}
got1 = m2b.p2_p3_fisher(de_frame("B cells", calls, tested), ifna_toy).iloc[0]

p2_hand = sum(comb(10, x) * comb(90, 10 - x) for x in range(5, 11)) / comb(100, 10)
compare(pd.DataFrame({
    "input": ["n_tested", "n_ifna_tested", "n_up", "up_ifna", "n_down", "down_ifna",
              "P2 odds ratio", "P2 p (greater), 6 s.f.", "P3 odds ratio", "P3 p (greater), 4 d.p."],
    "expected": [100, 10, 10, 5, 5, 0, 17.0, float(f"{p2_hand:.6g}"), np.inf, 0.0839],
    "code gives": [int(got1["n_tested"]), int(got1["n_ifna_tested"]), int(got1["n_up"]), int(got1["up_ifna"]),
                   int(got1["n_down"]), int(got1["down_ifna"]), float(got1["p2_odds_ratio"]),
                   float(f"{got1['p2_p_greater']:.6g}"), float(got1["p3_odds_ratio"]),
                   round(float(got1["p3_p_greater"]), 4)],
}), "1 P2/P3 Fisher")
print(f"hand hypergeometric P2 p = {p2_hand:.4g}")

# %% [markdown]
# **Should do (BH):** p = [0.01, 0.04, 0.03] → sorted 0.01, 0.03, 0.04 → p·3/rank = 0.03, 0.045, 0.04 → running
# minimum from the top = 0.03, 0.04, 0.04 → back in input order **[0.03, 0.04, 0.04]**.

# %%
compare(pd.DataFrame({
    "input": ["[0.01, 0.04, 0.03]"],
    "expected": [[0.03, 0.04, 0.04]],
    "code gives": [[round(float(x), 10) for x in m2b.bh(np.array([0.01, 0.04, 0.03]))]],
}), "1 BH")

# %% [markdown]
# ## 2. P4 expression-matched null (`p4_features`, `p4_null`)
#
# **Should do (features):** control mean = mean DESeq2-normalised count over the cell type's **kept** control
# samples; matching is on rank-percentiles of (control mean, detection).
#
# Toy: 41 genes g01..g41 with counts 1..41 in every kept sample (so size factors are 1 and the control mean of gNN is
# NN) and detection rising in the same order. A **dropped** (kept = False) control sample has the counts reversed:
# if it were wrongly used, every control mean would be 21 and the ranks would tie. Expected: control_mean = 1..41 and
# the two percentile columns both increase with the gene number.

# %%
genes = [f"g{i:02d}" for i in range(1, 42)]
smp = pd.DataFrame({"sample": ["T|d1|ctrl", "T|d1|stim", "T|d2|ctrl", "T|d2|stim", "T|d3|ctrl"],
                    "cell_type": "T", "donor": ["d1", "d1", "d2", "d2", "d3"],
                    "condition": ["ctrl", "stim", "ctrl", "stim", "ctrl"],
                    "n_cells": 100, "total_counts": 861, "kept": [True, True, True, True, False]})
cnt = pd.DataFrame({s: np.arange(1, 42) for s in smp["sample"]}, index=genes)
cnt["T|d3|ctrl"] = np.arange(41, 0, -1)
det = pd.Series(np.linspace(0.01, 0.41, 41), index=genes)
feat = m2b.p4_features(cnt, smp, genes, det)
compare(pd.DataFrame({
    "input": ["control_mean", "pct_mean strictly increasing", "pct_detection strictly increasing"],
    "expected": [list(range(1, 42)), True, True],
    "code gives": [[round(float(x), 6) for x in feat.loc[genes, "control_mean"]],
                   bool(feat.loc[genes, "pct_mean"].is_monotonic_increasing and feat["pct_mean"].is_unique),
                   bool(feat.loc[genes, "pct_detection"].is_monotonic_increasing and feat["pct_detection"].is_unique)],
}), "2 P4 features")

# %% [markdown]
# **Should do (pool):** for each up-DEG, the 30 nearest tested genes (Euclidean on the two percentiles), the focal
# gene excluded; the statistic is the number of IFN-α genes in the matched set; p = (1 + #null ≥ observed) / 1,001.
#
# All 41 genes lie on the diagonal, so the distance from the focal g21 is proportional to |NN − 21|. The 30 nearest
# are exactly g06..g36 minus g21 (|NN − 21| ≤ 15); g05 and g37 (distance 16) are the 31st/32nd.
#
# - **Case A:** IFN-α = {g21} ∪ the 10 far genes (g01–g05, g37–g41). Observed = 1 (g21 is IFN-α). The pool holds no
#   IFN-α gene, and g21 itself is excluded, so every draw gives 0 → p = (1 + 0) / 1001 = **0.000999**; % null genes
#   IFN = 0; % DEGs IFN = 100.
# - **Case B:** IFN-α = the 30 pool genes. Observed = 0; every draw gives 1 → p = (1 + 1000) / 1001 = **1.0**;
#   % null genes IFN = 100; % DEGs IFN = 0.
#
# Either case would break if the pool were 29 or 31 genes, used other features, or included the focal gene.

# %%
far = {f"g{i:02d}" for i in list(range(1, 6)) + list(range(37, 42))}
pool = {f"g{i:02d}" for i in range(6, 37)} - {"g21"}
rows2 = []
for case, ifn, exp in [("A", far | {"g21"}, (1, 0, 0, round(1 / 1001, 9), 0.0, 100.0)),
                       ("B", pool, (0, 1, 1, 1.0, 100.0, 0.0))]:
    res, null = m2b.p4_null(feat, ["g21"], ["g21"], ifn, k=30, n_draws=1000, seed=0)
    rows2.append({"input": f"case {case}: (observed, null min, null max, p, % null IFN, % DEGs IFN)",
                  "expected": exp,
                  "code gives": (int(res["observed"]), int(null.min()), int(null.max()),
                                 round(float(res["p_empirical"]), 9), float(res["pct_null_genes_ifna"]),
                                 float(res["pct_degs_ifna"]))})
compare(pd.DataFrame(rows2), "2 P4 null pool")

# %% [markdown]
# **Should do (p formula on a non-trivial null):** with two up-DEGs g20 and g21 and IFN-α = {g20, g22, g23}, the
# null count varies by draw. Whatever the draws, p must equal (1 + #(null ≥ observed)) / (n_draws + 1), counted here
# from the returned null vector; observed = 1 (only g20 among the up-DEGs is IFN-α). With 9 draws the denominator is
# 10.

# %%
res3, null3 = m2b.p4_null(feat, ["g20", "g21"], ["g20", "g21"], {"g20", "g22", "g23"}, k=30, n_draws=9, seed=0)
print("null counts:", list(null3))
compare(pd.DataFrame({
    "input": ["observed", "number of draws", "p"],
    "expected": [1, 9, (1 + int((null3 >= 1).sum())) / 10],
    "code gives": [int(res3["observed"]), len(null3), float(res3["p_empirical"])],
}), "2 P4 p formula")

# %% [markdown]
# ## 3. R1 / R2 / R3: LODO stability and overlap with other runs (`r1_lodo`, `overlap`)
#
# **Should do (R2/R3 overlap):** genes DE in either direction. Main DE = {A up, B up, C down, D up}; other run DE =
# {A up, B down, E up}; F is DE in neither. n_both = |{A, B}| = 2; Jaccard = 2 / (4 + 3 − 2) = **0.4**; % main
# recovered = 2 / 4 = **50**; sign flips = 1 (B).

# %%
main3 = de_frame("B cells", {"A": "up", "B": "up", "C": "down", "D": "up", "E": "none", "F": "none"})
other3 = de_frame("B cells", {"A": "up", "B": "down", "C": "none", "D": "none", "E": "up", "F": "none"})
got3 = m2b.overlap(main3, other3, "toy").iloc[0]
compare(pd.DataFrame({
    "input": ["n_main", "n_toy", "n_both", "jaccard", "pct_main_recovered", "sign_flips"],
    "expected": [4, 3, 2, 0.4, 50.0, 1],
    "code gives": [int(got3["n_main"]), int(got3["n_toy"]), int(got3["n_both"]), round(float(got3["jaccard"]), 10),
                   round(float(got3["pct_main_recovered"]), 10), int(got3["sign_flips"])],
}), "3 overlap")

# %% [markdown]
# **Should do (R1 LODO):** a DE call is stable if its log2FC keeps the call's sign in **every** one of the cell
# type's 8 leave-one-donor-out refits. Main: A up, B down, C up, D none. A is +1 in all 8; B is −1 in 7 and +0.2 in
# one; C is +0.5 in all 8; D (not DE) flips sign but must not count. Expected: n_de = 3, n_stable = 2 (A, C),
# pct_stable = **66.7**; per call: A True, B False, C True.

# %%
main_r1 = de_frame("B cells", {"A": "up", "B": "down", "C": "up", "D": "none"})
lodo_rows = []
for i in range(1, 9):
    for g, v in {"A": 1.0, "B": 0.2 if i == 5 else -1.0, "C": 0.5, "D": (-1.0) ** i}.items():
        lodo_rows.append({"run": f"lodo_d{i}", "cell_type": "B cells", "gene": g, "log2FC": v})
summ, per_call = m2b.r1_lodo(main_r1, pd.DataFrame(lodo_rows))
s = summ.set_index("cell_type").loc["B cells"]
pc = per_call.set_index("gene")["lodo_stable"]
compare(pd.DataFrame({
    "input": ["n_lodo_runs", "n_de", "n_stable", "pct_stable (1 d.p.)", "per call A, B, C"],
    "expected": [8, 3, 2, 66.7, (True, False, True)],
    "code gives": [int(s["n_lodo_runs"]), int(s["n_de"]), int(s["n_stable"]), round(float(s["pct_stable"]), 1),
                   tuple(bool(pc[g]) for g in "ABC")],
}), "3 R1 LODO")

# %% [markdown]
# ## 4. V9: recovery of the paper's DEGs (`paper_supp1`, `v9`)
#
# **Should do (mapping, real Supp Table 1):** the approved V9 mapping turns Excel-mangled symbols back into names by
# Ensembl ID (ENSG00000184640 → **SEPT9**) and the paper's PSMA2 into our make-unique name **PSMA2.1**; all 3,055
# paper genes then exist in our gene list (decisions, M1b V9 mapping row).

# %%
paper = m2b.paper_supp1()
our_genes = set(pd.read_csv(PROC / "m1b_de.csv", usecols=["gene"])["gene"])
compare(pd.DataFrame({
    "input": ["gene for ENSG00000184640", "paper PSMA2 appears as", "distinct paper genes", "all in our gene list"],
    "expected": ["SEPT9", "PSMA2.1", 3055, True],
    "code gives": [paper.loc[paper["ensembl"] == "ENSG00000184640", "gene"].unique().tolist()[0],
                   sorted(g for g in paper["gene"].unique() if g.startswith("PSMA2"))[0],
                   paper["ensembl"].nunique(), set(paper["gene"]) <= our_genes],
}), "4 V9 mapping")

# %% [markdown]
# **Should do (overlap):** per cell type, % of the paper's DEGs we also call DE, ✓ (pass) if > 70%.
#
# - B cells: paper {X, Y, Z, PSMA2.1}, ours {X, Y, PSMA2.1, W} → in both 3, paper only 1 (Z), ours only 1 (W),
#   recovered 3/4 = **75.0** → pass.
# - CD4 T cells: paper {X, Y, Z}, ours {X, Y} → 2/3 = **66.7** → fail.

# %%
de4 = pd.concat([
    de_frame("B cells", {"X": "up", "Y": "up", "Z": "none", "PSMA2.1": "up", "W": "up"}),
    de_frame("CD4 T cells", {"X": "up", "Y": "up", "Z": "none", "PSMA2.1": "none", "W": "none"}),
], ignore_index=True)
de4["log2FC_shrunk"] = de4["log2FC"] = de4["log2FC"] + 0.1 * np.arange(len(de4))  # non-constant, for the correlation
paper4 = pd.DataFrame({"cell_type": ["B cells"] * 4 + ["CD4 T cells"] * 3,
                       "gene": ["X", "Y", "Z", "PSMA2.1", "X", "Y", "Z"],
                       "ensembl": ["E1", "E2", "E3", "E4", "E1", "E2", "E3"],
                       "paper_log2FC": [2.1, 1.9, 1.5, 2.2, 2.0, 1.8, 1.4]})
ov4 = m2b.v9(de4, paper4)[0].set_index("cell_type")
compare(pd.DataFrame({
    "input": ["B cells", "CD4 T cells"],
    "expected (in_both, paper_only, ours_only, % recovered, V9)": [(3, 1, 1, 75.0, "pass"), (2, 1, 0, 66.7, "fail")],
    "code gives": [(int(ov4.loc[c, "in_both"]), int(ov4.loc[c, "paper_only"]), int(ov4.loc[c, "ours_only"]),
                    round(float(ov4.loc[c, "pct_paper_recovered"]), 1), ov4.loc[c, "V9"])
                   for c in ["B cells", "CD4 T cells"]],
}), "4 V9 overlap", expected_col="expected (in_both, paper_only, ours_only, % recovered, V9)")

# %% [markdown]
# ## 5. Modules: k-means vs the authors' labels, adjusted Rand index (`compare_modules`)
#
# **Should do:** ARI is computed on genes present in both. Toy (a): 6 genes in 3 obvious clusters, authors' labels
# the same grouping under other names → ARI = **1**; an authors' gene missing from ours is ignored.
#
# Toy (b): 4 genes, k = 2, ours {g1, g2} | {g3, g4}, authors {g1, g2, g3} | {g4}. By hand: Σ C(n_ij, 2) = C(2,2) = 1;
# rows Σ C(a_i, 2) = 1 + 1 = 2; columns Σ C(b_j, 2) = C(3,2) = 3; C(4,2) = 6; expected index = 2·3/6 = 1; max =
# (2 + 3)/2 = 2.5 → ARI = (1 − 1)/(2.5 − 1) = **0**.

# %%
za = pd.DataFrame([[5, 0], [5.1, 0], [0, 5], [0, 5.1], [-5, -5], [-5.1, -5]],
                  index=[f"a{i}" for i in range(6)], columns=["s1", "s2"])
auth_a = pd.DataFrame({"gene": [f"a{i}" for i in range(6)] + ["missing"], "module": list("665511") + ["0"]})
zb = pd.DataFrame([[5, 0], [5.1, 0], [-5, 0], [-5.1, 0]], index=["g1", "g2", "g3", "g4"], columns=["s1", "s2"])
auth_b = pd.DataFrame({"gene": ["g1", "g2", "g3", "g4"], "module": ["0", "0", "0", "1"]})
compare(pd.DataFrame({
    "input": ["(a) same grouping, other names (k = 3)", "(b) hand-computed table (k = 2)"],
    "expected": [1.0, 0.0],
    "code gives": [round(float(m2b.compare_modules(za, auth_a, k=3)[2]), 10),
                   round(float(m2b.compare_modules(zb, auth_b, k=2)[2]), 10)],
}), "5 modules ARI")

# %% [markdown]
# ## 6. I6: direction-aligned effect size per donor × cell type (`i6_effects`)
#
# **Should do:** CPM = count / column total × 10⁶; per DE gene log2((stim CPM + 1) / (ctrl CPM + 1)), times the sign of
# the call; median over the cell type's DE genes. x-axes = log of the mean of ctrl and stim n_cells / total_counts
# (natural log assumed; the decisions say "log").
#
# Toy donor, both columns summing to 10⁶ so CPM = count:
#
# | gene | call | ctrl | stim | log2((s+1)/(c+1)) | aligned |
# |---|---|---|---|---|---|
# | G1 | up | 3 | 7 | log2(8/4) = 1 | 1 |
# | G2 | up | 7 | 63 | log2(64/8) = 3 | 3 |
# | G3 | down | 15 | 3 | log2(4/16) = −2 | 2 |
# | G4 | none | 10 | 1000 | (not DE, ignored) | |
#
# Aligned median of [1, 3, 2] = **2**; unaligned median of [1, 3, −2] = **1**; n_de = 3; n_cells 100 and 300 →
# log(200); total_counts 10⁶ both → log(10⁶).

# %%
ctrl = {"G1": 3, "G2": 7, "G3": 15, "G4": 10}
stim = {"G1": 7, "G2": 63, "G3": 3, "G4": 1000}
ctrl["FILL"] = 10**6 - sum(ctrl.values())
stim["FILL"] = 10**6 - sum(stim.values())
cnt6 = pd.DataFrame({"B cells|d1|ctrl": ctrl, "B cells|d1|stim": stim})
smp6 = pd.DataFrame({"sample": list(cnt6.columns), "cell_type": "B cells", "donor": "d1",
                     "condition": ["ctrl", "stim"], "n_cells": [100, 300], "total_counts": [10**6, 10**6],
                     "kept": True})
de6 = de_frame("B cells", {"G1": "up", "G2": "up", "G3": "down", "G4": "none", "FILL": "none"})
g6 = m2b.i6_effects(cnt6, smp6, de6).iloc[0]
compare(pd.DataFrame({
    "input": ["n_de", "median_log2_ratio_aligned", "median_log2_ratio", "log_n_cells", "log_total_counts"],
    "expected": [3, 2.0, 1.0, round(log(200), 9), round(log(10**6), 9)],
    "code gives": [int(g6["n_de"]), round(float(g6["median_log2_ratio_aligned"]), 9),
                   round(float(g6["median_log2_ratio"]), 9), round(float(g6["log_n_cells"]), 9),
                   round(float(g6["log_total_counts"]), 9)],
}), "6 I6")

# %% [markdown]
# ## 7. I5 on the real data: ribosomal / MT- genes among the top 50 shared genes
#
# **Should do:** take genes labelled `shared` or `shared*` in `m2_table.csv`, rank them by mean |log2FC| over the
# 8 cell types (main run, `m1b_de.csv`), keep the top 50, and count RPL/RPS (D4 definition: prefix minus 14
# non-ribosomal genes, none of which is a plausible top IFN gene) and MT- genes. Recounted here with my own few lines;
# expected 0 and 0, and the same numbers as `m2_validation_i5_ribo.csv`.

# %%
table = pd.read_csv(PROC / "m2_table.csv")
shared = table.loc[table["shared_response"].isin(["shared", "shared*"]), "gene"]
de_all = pd.read_csv(PROC / "m1b_de.csv", usecols=["cell_type", "gene", "log2FC"])
score = de_all[de_all["gene"].isin(shared)].groupby("gene")["log2FC"].apply(lambda x: x.abs().mean())
top50 = score.sort_values(ascending=False).head(50).index
my_ribo = [g for g in top50 if g.startswith(("RPL", "RPS"))]
my_mt = [g for g in top50 if g.startswith("MT-")]
saved = pd.read_csv(PROC / "m2_validation_i5_ribo.csv").iloc[0]
print("top 10 shared by mean |log2FC|:", list(top50[:10]))
compare(pd.DataFrame({
    "input": ["n shared + shared* genes", "n ranked", "RPL/RPS in top 50", "MT- in top 50",
              "saved file: (n_shared, n_ranked, n_ribo, n_mt)"],
    "expected": ["(file) same as mine", 50, 0, 0, (len(shared), 50, 0, 0)],
    "code gives": ["(file) same as mine" if int(saved["n_shared"]) == len(shared) else int(saved["n_shared"]),
                   len(top50), len(my_ribo), len(my_mt),
                   (int(saved["n_shared"]), int(saved["n_ranked"]), int(saved["n_ribo"]), int(saved["n_mt"]))],
}), "7 I5")
print(f"shared + shared* = {len(shared)}")

# %% [markdown]
# ## Summary

# %%
print("Mismatches (POSSIBLE BUGS):", FAILURES or "none")
assert not FAILURES, FAILURES
