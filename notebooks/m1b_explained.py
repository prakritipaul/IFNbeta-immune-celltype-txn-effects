# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M1b explained by example
#
# Written by a reviewer who did not write `src/m1b_de.py` or `R/m1b_deseq2.R` (only their signatures, docstrings and
# the R header were read). Each section states in one sentence what a key line should do and runs the **real** code on
# a tiny input. The `expected` value is written first, from `audit/contracts/M1b_de.md`, `audit/decisions.md`
# (M1b settings, Differential expression) and `audit/test_plans/M1b.md`, never from the code. Each table shows
# **input | expected | code gives**.
#
# A mismatch is a **POSSIBLE BUG**: it is recorded, the notebook keeps going, and the last cell fails listing all
# mismatches. Expected values are never edited to match.
#
# Toys go through `run_deseq2`, which calls `R/m1b_deseq2.R` as a subprocess, exactly as the pipeline does.

# %%
import subprocess
import sys
import tempfile
from itertools import combinations
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from IPython.display import display

from src import m1b_de

RSCRIPT = "/usr/local/bin/Rscript"
TMP = Path(tempfile.mkdtemp(prefix="m1b_explained_"))  # toy files, outside the repo
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


rng = np.random.default_rng(1)


def nb(mu, disp=0.02):
    """Negative-binomial counts with mean mu and dispersion disp (variance mu + disp * mu^2)."""
    mu = np.asarray(mu, dtype=float)
    return rng.negative_binomial(1 / disp, 1 / (1 + mu * disp))


def toy_samples(donors, cell_type="T1"):
    ids = [f"{cell_type}|{d}|{c}" for d in donors for c in ("ctrl", "stim")]
    return pd.DataFrame({"sample": ids, "cell_type": cell_type, "donor": [i.split("|")[1] for i in ids],
                         "condition": [i.split("|")[2] for i in ids], "kept": True})


# %% [markdown]
# ## 1. The DE rule at its boundaries (`classify`)
#
# **Should do:** significant = q < 0.05 (strict); DE = significant and |log2FC| ≥ 1 (inclusive; the user changed this
# from > 1 on 2026-10-02, so the expected values for |log2FC| = 1 were updated with the requirement), `up` if log2FC > 0,
# `down` if < 0; a NaN q (here: an untested gene) is never significant and is `none`.

# %%
rows1 = [  # (qvalue, log2FC, tested) -> (significant, de), from decisions "Differential expression"
    (0.0499, 1.01, True, (True, "up")),
    (0.05, 2.0, True, (False, "none")),      # q exactly 0.05 is not < 0.05
    (0.01, 1.0, True, (True, "up")),         # |log2FC| exactly 1 counts (>= 1)
    (0.01, 0.9999, True, (True, "none")),    # just below 1 does not
    (0.01, -1.0, True, (True, "down")),
    (0.01, -1.01, True, (True, "down")),     # sign gives the direction; |.| gives the size
    (np.nan, 3.0, False, (False, "none")),   # untested gene: no q-value
]
t1 = pd.DataFrame([r[:3] for r in rows1], columns=["qvalue", "log2FC", "tested"])
t1["cooks_flag"], t1["baseMean"] = False, [10.0] * 6 + [0.5]  # 7 rows: 6 tested + 1 untested
got1 = m1b_de.classify(t1, q_threshold=0.05, lfc_threshold=1.0)
compare(pd.DataFrame({
    "input: (q, log2FC)": [(r[0], r[1]) for r in rows1],
    "expected (significant, de)": [r[3] for r in rows1],
    "code gives": list(zip(got1["significant"].astype(bool), got1["de"])),
}), "1 classify", expected_col="expected (significant, de)")

# %% [markdown]
# ## 2. Label shuffles (`balanced_shuffles`)
#
# **Should do:** list every distinct way to pick half of the donors to have their ctrl/stim labels swapped, counting a
# set and its complement once (swapping {a, b} is the same split as swapping {c, d}): C(8,4)/2 = 35, C(6,3)/2 = 10,
# C(4,2)/2 = 3 (O4). For 4 donors the 3 splits are ab|cd, ac|bd, ad|bc.

# %%
def splits(sets, donors):
    """Turn swap sets into unordered splits {set, complement}, so the choice of representative doesn't matter."""
    return [frozenset([frozenset(s), frozenset(donors) - frozenset(s)]) for s in sets]


rows2 = []
for n, n_expected in [(8, 35), (6, 10), (4, 3)]:
    donors = [f"d{i}" for i in range(1, n + 1)]
    got = m1b_de.balanced_shuffles(donors)
    sp_got = splits(got, donors)
    rows2 += [
        {"input": f"{n} donors: number of shuffles", "expected": n_expected, "code gives": len(got)},
        {"input": f"{n} donors: every swap set is exactly half", "expected": True,
         "code gives": all(len(s) == n // 2 for s in got)},
        {"input": f"{n} donors: no split repeated (incl. as a complement)", "expected": True,
         "code gives": len(set(sp_got)) == len(got)},
    ]
expected_4 = {frozenset([frozenset("ab"), frozenset("cd")]), frozenset([frozenset("ac"), frozenset("bd")]),
              frozenset([frozenset("ad"), frozenset("bc")])}
rows2.append({"input": "a,b,c,d: the splits", "expected": expected_4,
              "code gives": set(splits(m1b_de.balanced_shuffles(list("abcd")), "abcd"))})
compare(pd.DataFrame(rows2), "2 balanced_shuffles")

# %% [markdown]
# ### 2b. One shuffle, on the real samples
#
# **Should do:** in a shuffle run, each sample keeps its donor and its counts; for the donors in the swap set, the
# ctrl sample is relabelled `stim` and the stim sample `ctrl`; the other donors are untouched. Shown for
# Megakaryocytes (4 paired donors → `shuffle_1`..`shuffle_3`) using the run planner (pure; nothing is fitted).

# %%
counts_r = pd.read_csv(PROC / "m1a_counts.csv", index_col=0)
samples_r = pd.read_csv(PROC / "m1a_samples.csv")
counts_c = pd.read_csv(PROC / "m1a_counts_no_outliers.csv", index_col=0)
samples_c = pd.read_csv(PROC / "m1a_samples_no_outliers.csv")
plan = m1b_de._plan_runs(counts_r, samples_r, counts_c, samples_c, runs=("shuffle",))

CT = "Megakaryocytes"
main_mk = samples_r[(samples_r["cell_type"] == CT) & samples_r["kept"]].set_index("sample")
mk_donors = sorted(main_mk["donor"].unique())
swap_sets = {}
for name, c, s, design in plan:
    s_mk = s[(s["cell_type"] == CT) & s["kept"]].set_index("sample")
    if s_mk.empty:
        continue
    flipped = s_mk["condition"] != main_mk.loc[s_mk.index, "condition"]
    swap_sets[name] = (s_mk, flipped, c, design)

name0 = sorted(swap_sets)[0]
s0, f0, c0, design0 = swap_sets[name0]
show = main_mk[["donor", "condition"]].rename(columns={"condition": "main condition"})
show[f"{name0} donor"] = s0.loc[show.index, "donor"]
show[f"{name0} condition"] = s0.loc[show.index, "condition"]
display(show)

rows2b = [
    {"input": f"{name0}: donors unchanged", "expected": True,
     "code gives": bool((show["donor"] == show[f"{name0} donor"]).all())},
    {"input": f"{name0}: counts unchanged", "expected": True,
     "code gives": bool(c0[list(main_mk.index)].equals(counts_r[list(main_mk.index)]))},
    {"input": f"{name0}: design still paired", "expected": "paired", "code gives": design0},
]
got_sets = []
for name, (s_mk, flipped, _, _) in sorted(swap_sets.items()):
    per_donor = flipped.groupby(s_mk["donor"]).agg(["sum", "count"])
    swapped = frozenset(per_donor.index[per_donor["sum"] == 2])
    got_sets.append(swapped)
    rows2b.append({"input": f"{name}: each donor has 0 or both samples flipped; half the donors flipped",
                   "expected": (True, len(mk_donors) // 2),
                   "code gives": (bool(per_donor["sum"].isin([0, 2]).all()), len(swapped))})
rows2b.append({"input": f"{CT}: shuffle runs cover every split of its {len(mk_donors)} donors once",
               "expected": set(splits(m1b_de.balanced_shuffles(mk_donors), mk_donors)),
               "code gives": set(splits(got_sets, mk_donors))})
compare(pd.DataFrame(rows2b), "2b one shuffle")

# %% [markdown]
# ## 3. Contrast direction, through the real R script
#
# **Should do:** control is the reference, so log2FC = log2(stim / ctrl): a gene planted 4× higher in stim comes out
# log2FC ≈ +2 and `up`.
#
# Toy: one cell type, 6 donors × ctrl/stim. 40 null genes (mean 100–1000 per gene, a different level per donor, same
# in ctrl and stim, negative-binomial noise), `up4x` (mean 300 in ctrl, 1200 in stim), `down4x` (1200 → 300),
# and 200 low genes (mean 0.3, so baseMean ≤ 1, untested; used in section 4).

# %%
D6 = [f"d{i}" for i in range(1, 7)]
s3 = toy_samples(D6)
is_stim = (s3["condition"] == "stim").to_numpy()
donor_idx = np.repeat(np.arange(6), 2)
g3 = {}
for k in range(40):
    level = rng.uniform(100, 1000) * rng.lognormal(0, 0.3, 6)
    g3[f"null{k:02d}"] = nb(level[donor_idx])
donor_lvl = rng.lognormal(0, 0.3, 6)[donor_idx]
g3["up4x"] = nb(300 * donor_lvl * np.where(is_stim, 4, 1))
g3["down4x"] = nb(1200 * donor_lvl * np.where(is_stim, 0.25, 1))
for k in range(200):
    g3[f"low{k:03d}"] = nb(np.full(12, 0.3), disp=0.5)
c3 = pd.DataFrame(g3, index=s3["sample"]).T.astype(int)
display(c3.loc[["up4x", "down4x", "null00", "low000"]])

res3, fits3 = m1b_de.run_deseq2(c3, s3, design="paired", run_name="toy", min_basemean=1.0)
r3 = res3.set_index("gene")
display(r3.loc[["up4x", "down4x", "null00"], ["baseMean", "log2FC", "lfcSE", "pvalue", "qvalue", "de"]])
display(fits3)

# Contract: log2FC = treatment vs control (control is the reference level). Tolerance 0.4 = ~3 SE for this toy.
compare(pd.DataFrame({
    "input": ["up4x: log2FC within 0.4 of +2", "up4x: call", "down4x: log2FC within 0.4 of -2", "down4x: call",
              "null genes called DE (of 40)"],
    "expected": [True, "up", True, "down", 0],
    "code gives": [bool(abs(r3.loc["up4x", "log2FC"] - 2) < 0.4), r3.loc["up4x", "de"],
                   bool(abs(r3.loc["down4x", "log2FC"] + 2) < 0.4), r3.loc["down4x", "de"],
                   int((r3.loc[[f"null{k:02d}" for k in range(40)], "de"] != "none").sum())],
}), "3 contrast direction")

# %% [markdown]
# ## 4. The qvalue universe (the planted bug)
#
# **Should do:** within a cell type, qvalue is run on the p-values of **tested** genes only (baseMean > 1);
# untested genes get no q-value (O2, O3). The expected q-values come from my own call of R `qvalue()` on the
# tested genes' p-values, with the decided parameters.
#
# The section 3 toy has 42 tested genes and up to 200 untested low genes, many with a p-value. For contrast, the same
# qvalue call on **all** p-values (the planted bug) is shown, to prove the comparison can tell the two apart.

# %%
def r_qvalue(p, tag):
    pd.DataFrame({"p": p}).to_csv(TMP / f"{tag}_p.csv", index=False)
    r = subprocess.run([RSCRIPT, "-e",
        f"suppressMessages(library(qvalue)); p <- read.csv('{TMP}/{tag}_p.csv')$p; "
        "q <- qvalue(p, lambda = seq(0.05, 0.95, 0.05), pi0.method = 'smoother', smooth.df = 3)$qvalues; "
        f"write.csv(data.frame(q = q), '{TMP}/{tag}_q.csv', row.names = FALSE)"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-1000:]
    return pd.read_csv(TMP / f"{tag}_q.csv")["q"].to_numpy()


tested_genes = r3.index[r3["baseMean"] > 1]
p_tested = r3.loc[tested_genes, "pvalue"]
assert p_tested.notna().all()
q_expected = pd.Series(r_qvalue(p_tested.to_numpy(), "tested"), index=tested_genes)
p_all = r3["pvalue"].dropna()
q_allgenes = pd.Series(r_qvalue(p_all.to_numpy(), "all"), index=p_all.index)

untested_with_p = int((~r3.index.isin(tested_genes) & r3["pvalue"].notna()).sum())
print(f"tested genes: {len(tested_genes)}; untested genes with a p-value: {untested_with_p}")
print(f"toy check: max |q tested-only - q all-genes| on tested genes = "
      f"{(q_expected - q_allgenes[tested_genes]).abs().max():.3g} (must be > 0 for the check to have teeth)")
assert (q_expected - q_allgenes[tested_genes]).abs().max() > 1e-6
display(pd.DataFrame({"pvalue": p_tested, "expected q (tested only)": q_expected,
                      "q if all genes (bug)": q_allgenes[tested_genes], "code gives": r3.loc[tested_genes, "qvalue"]})
        .sort_values("pvalue").head(6))

compare(pd.DataFrame({
    "input": ["tested = baseMean > 1 (same genes)", "tested genes: q equals qvalue on tested p only (to 1e-9)",
              "untested genes: q is NaN"],
    "expected": [True, True, True],
    "code gives": [set(r3.index[r3["tested"].astype(bool)]) == set(tested_genes),
                   bool(np.allclose(r3.loc[tested_genes, "qvalue"], q_expected, rtol=0, atol=1e-9)),
                   bool(r3.loc[~r3.index.isin(tested_genes), "qvalue"].isna().all())],
}), "4 qvalue universe")

# %% [markdown]
# ## 5. Paired vs unpaired design
#
# **Should do:** `paired` fits `~ donor + condition`, `unpaired` fits `~ condition`. When donors differ a lot and the
# stim effect is modest but consistent within each donor, the paired design gives a smaller p-value, because the
# donor term removes the between-donor spread that the unpaired model has to count as noise.
#
# Toy: 6 donors, 200 genes each with its own large donor-to-donor spread (log-normal SD 1, ~e^±2 range), little noise
# within a donor; `target` gets 1.6× in stim in every donor, the others no effect. (200 genes, not 40: R `qvalue`'s smoother fails on too few p-values.)

# %%
g5 = {}
for k in range(200):
    lvl = 300 * rng.lognormal(0, 1, 6)
    g5[f"null{k:03d}"] = nb(lvl[donor_idx], disp=0.01)
tgt = 300 * rng.lognormal(0, 1, 6)
g5["target"] = nb(tgt[donor_idx] * np.where(is_stim, 1.6, 1), disp=0.01)
c5 = pd.DataFrame(g5, index=s3["sample"]).T.astype(int)
display(c5.loc[["target"]])

res_p, _ = m1b_de.run_deseq2(c5, s3, design="paired", run_name="toy_paired", min_basemean=1.0)
res_u, _ = m1b_de.run_deseq2(c5, s3, design="unpaired", run_name="toy_unpaired", min_basemean=1.0)
p_p = float(res_p.set_index("gene").loc["target", "pvalue"])
p_u = float(res_u.set_index("gene").loc["target", "pvalue"])
print(f"target p-value: paired {p_p:.3g}, unpaired {p_u:.3g}")
compare(pd.DataFrame({"input": ["target: paired p < unpaired p"], "expected": [True], "code gives": [p_p < p_u]}),
        "5 paired vs unpaired")

# %% [markdown]
# ## 6. Real data
#
# **Should do:** in the main run, ISG15 (the classic ISG) is `up` in all 8 cell types (CLAUDE.md sanity check).
# The label shuffles keep the counts but break the ctrl/stim link, so they should give at most a handful of DE genes:
# I set the bar, before looking, at ≤ 10 per shuffle run × cell type, and fewer than the main run in every cell type.
# Shuffle runs per cell type: 35, 10 for Dendritic cells, 3 for Megakaryocytes. `m1b_runs.csv` (0.84 GB) is read by
# columns only.

# %%
de_main = pd.read_csv(PROC / "m1b_de.csv")
isg15 = de_main[de_main["gene"] == "ISG15"].set_index("cell_type")
display(isg15[["log2FC", "qvalue", "de"]])
compare(pd.DataFrame({
    "input": ["cell types in m1b_de.csv", "ISG15 `up` in all of them"],
    "expected": [8, 8],
    "code gives": [de_main["cell_type"].nunique(), int((isg15["de"] == "up").sum())],
}), "6 ISG15")

runs = pd.read_csv(PROC / "m1b_runs.csv", usecols=["run", "cell_type", "de"])
runs = runs[(runs["run"] == "main") | runs["run"].str.startswith("shuffle_")]
n_de = (runs.assign(is_de=runs["de"] != "none").groupby(["run", "cell_type"])["is_de"].sum().rename("n_de")
        .reset_index())
main_n = n_de[n_de["run"] == "main"].set_index("cell_type")["n_de"]
sh = n_de[n_de["run"].str.startswith("shuffle_")]
summary = sh.groupby("cell_type")["n_de"].agg(n_shuffles="size", median="median", max="max")
summary["main run DE"] = main_n
display(summary)

expected_n = {ct: 35 for ct in main_n.index}
expected_n.update({"Dendritic cells": 10, "Megakaryocytes": 3})
compare(pd.DataFrame({
    "input: cell type": summary.index,
    "expected (shuffles, max ≤ 10, max < main)": [(expected_n[ct], True, True) for ct in summary.index],
    "code gives": [(int(r["n_shuffles"]), bool(r["max"] <= 10), bool(r["max"] < r["main run DE"]))
                   for _, r in summary.iterrows()],
}), "6 shuffles", expected_col="expected (shuffles, max ≤ 10, max < main)")

# %% [markdown]
# ## Summary

# %%
print("Mismatches (POSSIBLE BUGS):", FAILURES or "none")
assert not FAILURES, FAILURES
