# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M2c explained by example
#
# Written by a reviewer who did not write `src/m2c_programs.py` or `R/m2c_fgsea.R` (only their signatures and
# docstrings were read). Each section states in one sentence what a key rule should do and runs the **real** code on a
# tiny input. The `expected` value is written first, from `audit/contracts/M2c_programs.md`, `audit/decisions.md`
# (G1–G6, D2, M2c confirmations) and `audit/test_plans/M2c.md`, never from the code. Each table shows
# **input | expected | code gives**.
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

from src import m2c_programs as m2c
from src.m2a_calls import CELL_TYPES as CT

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


def calls_frame(cells):
    """cells: {cell_type: {gene: (qvalue, log2FC, de)}} -> m1b_de-shaped frame on a full gene x cell type grid
    (genes missing from a cell type are tested, not DE: q 0.9, log2FC 0)."""
    genes = sorted({g for d in cells.values() for g in d})
    rows = []
    for ct, d in cells.items():
        for g in genes:
            q, lfc, de = d.get(g, (0.9, 0.0, "none"))
            rows.append(dict(run="main", cell_type=ct, gene=g, baseMean=50.0, log2FC=lfc, lfcSE=0.5,
                             ci_low=lfc - 1, ci_high=lfc + 1, pvalue=q, qvalue=q, tested=True, cooks_flag=False,
                             log2FC_shrunk=lfc, significant=q < 0.05, de=de))
    return pd.DataFrame(rows)


# %% [markdown]
# ## 1. Top 20 up / down (`top_genes`)
#
# **Should do (G5):** per cell type, only DE genes; rank by q ascending, ties by larger |log2FC| (for down genes too);
# at most 20 per direction, fewer rows (no padding) when fewer are DE.
#
# Toy, cell type `A`: 18 fillers `f01`–`f18` (q 1e-39 … 1e-22, in that order), a q tie at 1e-12 between `zz_big`
# (log2FC 5) and `aa_small` (log2FC 1.2), then `cut1` (q 1e-6) and `cut2` (q 0.01), which fall outside the top 20.
# Traps that are **not** DE: `trap_q` (q 1e-50, log2FC 0.8) and `trap_fc` (q 0.2, log2FC 8). Down: `d1` (q 1e-3,
# log2FC −2), `d2` (q 1e-8, −1.1), `d3` (q 1e-3, −4). Cell type `B`: two up genes, nothing down.

# %%
A = {f"f{k:02d}": (10.0 ** -(40 - k), 2.0, "up") for k in range(1, 19)}
A |= {"zz_big": (1e-12, 5.0, "up"), "aa_small": (1e-12, 1.2, "up"), "cut1": (1e-6, 3.0, "up"),
      "cut2": (0.01, 3.0, "up"), "trap_q": (1e-50, 0.8, "none"), "trap_fc": (0.2, 8.0, "none"),
      "d1": (1e-3, -2.0, "down"), "d2": (1e-8, -1.1, "down"), "d3": (1e-3, -4.0, "down")}
B = {"b_late": (1e-4, 1.5, "up"), "b_early": (1e-9, 1.1, "up")}
top = m2c.top_genes(calls_frame({"A": A, "B": B}), n=20)
display(top.head(3), top.tail(6))

exp1 = {
    ("A", "up"): [f"f{k:02d}" for k in range(1, 19)] + ["zz_big", "aa_small"],
    ("A", "down"): ["d2", "d3", "d1"],
    ("B", "up"): ["b_early", "b_late"],
    ("B", "down"): [],
}
got1 = {k: top[(top.cell_type == k[0]) & (top.direction == k[1])].sort_values("rank")["gene"].tolist() for k in exp1}
compare(pd.DataFrame({
    "input: cell type, direction": [f"{a} {d}" for a, d in exp1],
    "expected genes in rank order": [", ".join(v) for v in exp1.values()],
    "code gives": [", ".join(got1[k]) for k in exp1],
}), "1 top genes", expected_col="expected genes in rank order")
compare(pd.DataFrame({
    "input": ["ranks of A up", "total rows", "trap genes present"],
    "expected": [list(range(1, 21)), 25, False],
    "code gives": [top[(top.cell_type == "A") & (top.direction == "up")]["rank"].sort_values().tolist(), len(top),
                   bool(top["gene"].isin(["trap_q", "trap_fc"]).any())],
}), "1 top genes, counts")

# %% [markdown]
# ## 2. Jaccard of DE sets (`jaccard_matrix`)
#
# **Should do (G5, M2c confirmations):** |A ∩ B| / |A ∪ B| on the full DE set of one direction; up and down never mixed;
# two empty sets → NaN (decisions.md), not 1.
#
# Toy: `X` up {a, b, c}, down {e}; `Y` up {b, c, d}, down {e, f, g}; `Z` up {}, down {a, b, c}.
# If directions were mixed, X–Y would be {a,b,c,e} vs {b,c,d,e,f,g} = 3/7 instead of 0.5 (up) and 1/3 (down).

# %%
U, D = (1e-5, 2.0, "up"), (1e-5, -2.0, "down")
jc = calls_frame({"X": {"a": U, "b": U, "c": U, "e": D},
                  "Y": {"b": U, "c": U, "d": U, "e": D, "f": D, "g": D},
                  "Z": {"a": D, "b": D, "c": D}})
ju, jd = m2c.jaccard_matrix(jc, "up"), m2c.jaccard_matrix(jc, "down")
display(ju, jd)
exp2 = {("up", "X", "Y"): 0.5, ("up", "X", "Z"): 0.0, ("up", "Y", "Z"): 0.0, ("up", "X", "X"): 1.0,
        ("up", "Z", "Z"): "NaN", ("down", "X", "Y"): 1 / 3, ("down", "X", "Z"): 0.0, ("down", "Y", "Z"): 0.0,
        ("down", "Y", "X"): 1 / 3, ("down", "Z", "Z"): 1.0}


def jval(m, a, b):
    v = m.loc[a, b]
    return "NaN" if pd.isna(v) else round(float(v), 12)


compare(pd.DataFrame({
    "input: direction, pair": [f"{d} {a}–{b}" for d, a, b in exp2],
    "expected": [v if isinstance(v, str) else round(v, 12) for v in exp2.values()],
    "code gives": [jval(ju if d == "up" else jd, a, b) for d, a, b in exp2],
}), "2 jaccard")

# %% [markdown]
# ## 3. Ranking input (`rank_genes`)
#
# **Should do (G2, T4):** keep only tested genes (baseMean > 1, so 0.99 and 1.00 out, 1.01 in); stat = log2FC / lfcSE;
# sorted by stat descending, ties by gene name (decisions: "tied stats ordered by gene name").
#
# Toy (baseMean, log2FC, lfcSE): `lo099` (0.99, 10, 1), `eq100` (1.00, 9, 1), `hi101` (1.01, 2, 0.5),
# `big_lfc` (50, 6, 4), `tie_b` (50, 1, 1), `tie_a` (50, 0.5, 0.5), `neg` (50, −3, 1.5).
# Ranking by log2FC would put `lo099`/`eq100`/`big_lfc` first.

# %%
toy3 = {"lo099": (0.99, 10, 1), "eq100": (1.00, 9, 1), "hi101": (1.01, 2, 0.5), "big_lfc": (50, 6, 4),
        "tie_b": (50, 1, 1), "tie_a": (50, 0.5, 0.5), "neg": (50, -3, 1.5)}
c3 = calls_frame({"C1": {g: (0.01, v[1], "none") for g, v in toy3.items()}})
c3["baseMean"] = c3["gene"].map(lambda g: toy3[g][0])
c3["lfcSE"] = c3["gene"].map(lambda g: toy3[g][2])
c3["tested"] = c3["baseMean"] > 1
r3 = m2c.rank_genes(c3)
display(r3)
exp3 = [("hi101", 4.0), ("big_lfc", 1.5), ("tie_a", 1.0), ("tie_b", 1.0), ("neg", -2.0)]
compare(pd.DataFrame({
    "input": ["genes in order", "stats in order"],
    "expected": [[g for g, _ in exp3], [s for _, s in exp3]],
    "code gives": [r3["gene"].tolist(), [round(float(s), 12) for s in r3["stat"]]],
}), "3 rank_genes")

# %% [markdown]
# ## 4. GSEA sign, through the real R script (`run_gsea` → `R/m2c_fgsea.R`)
#
# **Should do (G3, T5):** a set whose genes sit near the top of the ranking gets NES > 0, near the bottom NES < 0, both
# padj < 0.05. Toy: 1,000 genes, stat evenly from +3 to −3 (no ties). `TOP` = every 5th gene of ranks 1–150 (30 genes),
# `BOTTOM` = every 5th of ranks 851–1000, `SPREAD` = every 33rd gene over the whole list (no expectation on its
# significance). Cell type `C2` gets the reversed stats, so its signs must flip. Two collections (`toyA`: TOP, BOTTOM,
# SPREAD; `toyB`: five more 30-gene sets) are used again in section 5.

# %%
genes = [f"g{i:04d}" for i in range(1000)]          # g0000 = rank 1 (highest stat)
stat = np.linspace(3, -3, 1000)
ranks = pd.concat([pd.DataFrame({"cell_type": "C1", "gene": genes, "stat": stat}),
                   pd.DataFrame({"cell_type": "C2", "gene": genes[::-1], "stat": stat})], ignore_index=True)
rng = np.random.default_rng(1)
collections = {
    "toyA": {"TOP": genes[0:150:5], "BOTTOM": genes[850:1000:5], "SPREAD": genes[::33][:30]},
    "toyB": {f"B{k}": sorted(rng.choice(genes[k * 100:k * 100 + 400], 30, replace=False).tolist()) for k in range(5)},
}
gsea = m2c.run_gsea(ranks, collections, seed=0)
display(gsea[["cell_type", "collection", "pathway", "size", "NES", "pval", "padj"]])
g4 = gsea.set_index(["cell_type", "pathway"])
exp4 = {("C1", "TOP"): ("NES > 0", True), ("C1", "BOTTOM"): ("NES < 0", True),
        ("C2", "TOP"): ("NES < 0", True), ("C2", "BOTTOM"): ("NES > 0", True)}
compare(pd.DataFrame({
    "input: cell type, set": [f"{c} {p}" for c, p in exp4],
    "expected (sign, padj < 0.05)": list(exp4.values()),
    "code gives": [("NES > 0" if g4.loc[k, "NES"] > 0 else "NES < 0", bool(g4.loc[k, "padj"] < 0.05)) for k in exp4],
}), "4 GSEA sign", expected_col="expected (sign, padj < 0.05)")


# %% [markdown]
# ## 5. BH over both collections within a cell type (`R/m2c_fgsea.R`, the planted bug 4 line)
#
# **Should do (G4):** padj = Benjamini–Hochberg over **all** sets of the cell type (toyA and toyB together, 8 sets),
# separately per cell type. I recompute BH by hand from the script's own `pval` and compare. As a check that this toy
# can tell the right scope from the wrong one, per-collection BH must differ from the script's padj somewhere.

# %%
def bh(p):
    """Benjamini–Hochberg by hand: p_(i) * m / i, cumulative min from the largest p, capped at 1."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adj)
    out[o] = np.minimum(adj, 1)
    return out


mine5 = gsea.groupby("cell_type")["pval"].transform(bh)
per_coll = gsea.groupby(["cell_type", "collection"])["pval"].transform(bh)
display(gsea.assign(my_bh_all=mine5, bh_per_collection=per_coll)[
    ["cell_type", "collection", "pathway", "pval", "padj", "my_bh_all", "bh_per_collection"]])
compare(pd.DataFrame({
    "input": ["sets per cell type", "padj == my BH over both collections (all rows)",
              "per-collection BH differs from padj somewhere (toy discriminates)"],
    "expected": [{"C1": 8, "C2": 8}, True, True],
    "code gives": [gsea.groupby("cell_type").size().to_dict(), bool(np.allclose(gsea["padj"], mine5, rtol=1e-12, atol=0)),
                   bool(not np.allclose(gsea["padj"], per_coll, rtol=1e-12, atol=0))],
}), "5 BH scope")

# %% [markdown]
# ## 6. Program labels (`label_programs`)
#
# **Should do (G4, D2):** affected = padj < 0.05 (0.05 itself is not); direction = sign of NES; `shared` = affected
# in all 8, same sign; `partial` = ≥ 2 same sign; `specific` = 1; `mixed` = up somewhere, down elsewhere; signs of
# unaffected cells ignored; programs affected nowhere are not shown. A cell type with no row for a program is `n/t`,
# counts as not affected, and a cell type with no rows at all still gets an `n/t` column (Codex review fix #1).
#
# Cell order in the strings: B, CD14+ Mono, CD4 T, CD8 T, FCGR3A+ Mono, NK, DC, Mk.

# %%
def gsea_toy(spec):
    """spec: {pathway: {cell_type: (NES, padj)}}; cell types not listed have no row (n/t)."""
    return pd.DataFrame([dict(cell_type=ct, collection="toy", pathway=p, size=30, NES=n, pval=q, padj=q)
                         for p, d in spec.items() for ct, (n, q) in d.items()])


SIG_UP, SIG_DN, NS_DN, NS = (2.0, 0.01), (-2.0, 0.01), (-1.5, 0.3), (0.5, 0.5)
B_, C14, C4, C8, F16, NK, DC, MK = CT
spec6 = {
    "P_shared": {ct: SIG_UP for ct in CT},
    "P_shared_down": {ct: SIG_DN for ct in CT},
    "P_partial": {ct: SIG_UP if ct in (B_, C4, NK) else NS_DN for ct in CT},
    "P_specific": {ct: SIG_DN if ct == C14 else NS for ct in CT},
    "P_mixed": {ct: SIG_UP if ct == B_ else SIG_DN if ct == NK else NS for ct in CT},
    "P_0499": {ct: (2.0, 0.0499) if ct == F16 else NS for ct in CT},
    "P_boundary": {ct: (2.0, 0.05) for ct in CT},
    "P_none": {ct: NS for ct in CT},
    "P_nt": {ct: SIG_UP for ct in CT if ct != DC},
}
got6 = m2c.label_programs(gsea_toy(spec6), CT, padj_threshold=0.05).set_index("pathway")
display(got6)
exp6 = {  # (arrows over CT, n_up, n_down, label)
    "P_shared": ("↑|↑|↑|↑|↑|↑|↑|↑", 8, 0, "shared"),
    "P_shared_down": ("↓|↓|↓|↓|↓|↓|↓|↓", 0, 8, "shared"),
    "P_partial": ("↑|–|↑|–|–|↑|–|–", 3, 0, "partial"),
    "P_specific": ("–|↓|–|–|–|–|–|–", 0, 1, "specific"),
    "P_mixed": ("↑|–|–|–|–|↓|–|–", 1, 1, "mixed"),
    "P_0499": ("–|–|–|–|↑|–|–|–", 1, 0, "specific"),
    "P_nt": ("↑|↑|↑|↑|↑|↑|n/t|↑", 7, 0, "partial"),
}


def row6(df, p):
    if p not in df.index:
        return "not shown"
    r = df.loc[p]
    return ("|".join(r[CT]), int(r["n_up"]), int(r["n_down"]), r["label"])


compare(pd.DataFrame({
    "input: program": list(exp6) + ["P_boundary", "P_none"],
    "expected (arrows, n_up, n_down, label)": list(exp6.values()) + ["not shown", "not shown"],
    "code gives": [row6(got6, p) for p in list(exp6) + ["P_boundary", "P_none"]],
}), "6 labels", expected_col="expected (arrows, n_up, n_down, label)")

# A cell type with no GSEA rows at all (Megakaryocytes): still a column, n/t everywhere.
spec6b = {"P_seven": {ct: SIG_UP for ct in CT if ct != MK}}
got6b = m2c.label_programs(gsea_toy(spec6b), CT, padj_threshold=0.05).set_index("pathway")
display(got6b)
compare(pd.DataFrame({
    "input": ["P_seven (up in 7, Mk has no rows at all)", "columns include all 8 cell types"],
    "expected": [("↑|↑|↑|↑|↑|↑|↑|n/t", 7, 0, "partial"), True],
    "code gives": [row6(got6b, "P_seven"), all(ct in got6b.columns for ct in CT)],
}), "6 n/t cell type")

# %% [markdown]
# ## 7. Real data: labels re-counted independently (`m2c_programs.csv` from `m2c_gsea.csv`)
#
# **Should do:** my own few lines on `m2c_gsea.csv` (padj < 0.05, sign of NES, missing row = n/t = not affected, label
# over the 8 cell types) give shared 5, partial 36, specific 16, mixed 1, the same arrows as `m2c_programs.csv`, and
# padj equals BH over both collections within each real cell type.

# %%
real = pd.read_csv(PROC / "m2c_gsea.csv")
progs = pd.read_csv(PROC / "m2c_programs.csv").set_index(["collection", "pathway"])
real["call"] = np.where(real["padj"] < 0.05, np.where(real["NES"] > 0, "↑", "↓"), "–")
arrows = real.pivot(index=["collection", "pathway"], columns="cell_type", values="call").reindex(columns=CT).fillna("n/t")
n_up, n_dn = (arrows == "↑").sum(axis=1), (arrows == "↓").sum(axis=1)
hit = (n_up + n_dn) > 0
arrows, n_up, n_dn = arrows[hit], n_up[hit], n_dn[hit]
mine7 = pd.Series(np.select([(n_up > 0) & (n_dn > 0), (n_up == 8) | (n_dn == 8), (n_up + n_dn) >= 2],
                            ["mixed", "shared", "partial"], "specific"), index=arrows.index)

order = ["shared", "partial", "specific", "mixed"]
decided = dict(zip(order, [5, 36, 16, 1]))
compare(pd.DataFrame({
    "input: label": order,
    "expected (given, my count)": [(decided[k], decided[k]) for k in order],
    "code gives": [(decided[k], int((mine7 == k).sum())) for k in order],
}), "7 my count vs given", expected_col="expected (given, my count)")
idx = mine7.index.intersection(progs.index)
compare(pd.DataFrame({
    "input": ["same programs", "programs whose label differs", "program × cell type arrows that differ",
              "padj == BH over both collections, per real cell type"],
    "expected": [True, 0, 0, True],
    "code gives": [set(mine7.index) == set(progs.index), int((mine7[idx] != progs.loc[idx, "label"]).sum()),
                   int((arrows.loc[idx, CT] != progs.loc[idx, CT]).to_numpy().sum()),
                   bool(np.allclose(real["padj"], real.groupby("cell_type")["pval"].transform(bh), rtol=1e-9, atol=0))],
}), "7 program-by-program")

# %% [markdown]
# The 5 shared programs, with direction (all 8 arrows the same):

# %%
shared = progs[progs["label"] == "shared"]
display(shared[CT + ["n_up", "n_down", "label"]])
compare(pd.DataFrame({
    "input: shared program": [p for _, p in shared.index],
    "expected arrows (mine)": ["".join(arrows.loc[k, CT]) for k in shared.index],
    "code gives": ["".join(shared.loc[k, CT]) for k in shared.index],
}), "7 shared programs", expected_col="expected arrows (mine)")

# %% [markdown]
# ## Summary

# %%
print("Mismatches (POSSIBLE BUGS):", FAILURES or "none")
assert not FAILURES, FAILURES
