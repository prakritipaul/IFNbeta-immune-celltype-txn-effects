# ---
# jupyter:
#   jupytext:
#     text_representation: {extension: .py, format_name: percent}
#   kernelspec: {display_name: Python 3, language: python, name: python3}
# ---

# %% [markdown]
# # M1a explained by example
#
# Written by a reviewer who did not write `src/m1a_pseudobulk.py` or `R/m1a_varpart.R`. Each section takes one of the
# five lines worth reading from the M1a reviewer brief (`audit/reviews/005_M1a_codex_reviewer_brief.md`, section 4),
# states in one sentence what it should do, and runs the **real** code on a tiny hand-checkable input.
# The `expected` value in each code cell is written first, from `audit/contracts/M1a_pseudobulk.md`,
# `audit/decisions.md` (Wave 1) and `audit/test_plans/M1a.md` (interface), never from the code. Then the code runs, and
# each table shows **input | expected | code gives**.
#
# A mismatch is a **POSSIBLE BUG**: it is recorded, the notebook keeps going so every example runs, and the last
# cell fails listing all mismatches. Expected values are never edited to match.
#
# Thresholds are passed explicitly. The R script runs as a subprocess, exactly as the pipeline calls it.

# %%
import subprocess
import tempfile
import sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from IPython.display import display

from src import m1a_pseudobulk

RSCRIPT = "/usr/local/bin/Rscript"
VARPART = ROOT / "R/m1a_varpart.R"
TMP = Path(tempfile.mkdtemp(prefix="m1a_explained_"))  # toy files, outside the repo

FAILURES = []


def compare(table, name, expected_col="expected", actual_col="code gives"):
    """Show the table; record (not stop on) any row where expected differs from what the code gives."""
    table["match"] = [e == a for e, a in zip(table[expected_col], table[actual_col])]
    display(table)
    if not table["match"].all():
        print(f"POSSIBLE BUG in {name}:")
        display(table[~table["match"]])
        FAILURES.append(name)


def run_varpart(counts, samples, tag, min_basemean=1):
    """Write toy inputs and run R/m1a_varpart.R as the pipeline does. Returns (exit code, out path, stderr tail)."""
    c, s, o = TMP / f"{tag}_counts.csv", TMP / f"{tag}_samples.csv", TMP / f"{tag}_varpart.csv"
    for f in (o, TMP / f"{tag}_varpart_vst.csv", TMP / f"{tag}_varpart_fits.csv"):
        f.unlink(missing_ok=True)
    counts.to_csv(c, index_label="gene")
    samples.to_csv(s, index=False)
    r = subprocess.run([RSCRIPT, str(VARPART), str(c), str(s), str(o), str(min_basemean)],
                       capture_output=True, text=True)
    return r.returncode, o, r.stderr[-1500:]


def toy_adata(rows, genes):
    """rows: (donor, condition, cell_type, counts) per cell."""
    obs = pd.DataFrame([r[:3] for r in rows], columns=["donor", "condition", "cell_type"],
                       index=[f"c{i}" for i in range(1, len(rows) + 1)])
    a = ad.AnnData(X=sp.csr_matrix(np.array([r[3] for r in rows])), obs=obs, var=pd.DataFrame(index=genes))
    a.uns["control"], a.uns["treatment"] = "ctrl", "stim"
    return a


# %% [markdown]
# ## 1. Summing cells into samples (`pseudobulk`)
#
# **Line:** `src/m1a_pseudobulk.py:83` (the sparse matrix product that sums cells into samples).
#
# **Should do:** for each cell type × donor × condition present in the data, add up the raw counts of its cells, gene by
# gene; record how many cells went in (`n_cells`) and their total counts (`total_counts`).
#
# `min_cells=1` here so that the small-group rule (section 2) only bites on the missing condition.
#
# | cell | donor | condition | cell type | g1 | g2 |
# |---|---|---|---|---|---|
# | c1 | A | ctrl | B cells | 1 | 0 |
# | c2 | A | ctrl | B cells | 2 | 3 |
# | c3 | A | stim | B cells | 0 | 5 |
# | c4 | B | ctrl | B cells | 4 | 1 |
# | c5 | B | stim | B cells | 1 | 1 |
# | c6 | B | stim | NK cells | 7 | 0 |
#
# Donor B has no control NK cell, so there is no `NK cells|B|ctrl` column, and `NK cells|B|stim` is unpaired:
# kept = False (interface: "a missing condition counts as 0 cells").

# %%
a1 = toy_adata([
    ("A", "ctrl", "B cells", [1, 0]),
    ("A", "ctrl", "B cells", [2, 3]),
    ("A", "stim", "B cells", [0, 5]),
    ("B", "ctrl", "B cells", [4, 1]),
    ("B", "stim", "B cells", [1, 1]),
    ("B", "stim", "NK cells", [7, 0]),
], genes=["g1", "g2"])

# Contract step 1 + interface: column per sample id "cell_type|donor|condition"; sums by hand from the table above.
# sample: ([g1, g2], n_cells, total_counts, kept)
expected = {
    "B cells|A|ctrl": ([3, 3], 2, 6, True),
    "B cells|A|stim": ([0, 5], 1, 5, True),
    "B cells|B|ctrl": ([4, 1], 1, 5, True),
    "B cells|B|stim": ([1, 1], 1, 2, True),
    "NK cells|B|stim": ([7, 0], 1, 7, False),
}

counts1, samples1 = m1a_pseudobulk.pseudobulk(a1, min_cells=1)
s1 = samples1.set_index("sample")
got = {sid: (counts1[sid].tolist(), int(s1.loc[sid, "n_cells"]), int(s1.loc[sid, "total_counts"]),
             bool(s1.loc[sid, "kept"])) for sid in counts1.columns}
ids = sorted(set(expected) | set(got))
compare(pd.DataFrame({
    "input: sample": ids,
    "expected ([g1, g2], n_cells, total, kept)": [expected.get(i, "absent") for i in ids],
    "code gives": [got.get(i, "absent") for i in ids],
}), "1 summing", expected_col="expected ([g1, g2], n_cells, total, kept)")
# Contract assumption: counts stay integers after summing.
compare(pd.DataFrame({"input": ["counts dtype is integer"], "expected": [True],
                      "code gives": [bool(np.issubdtype(counts1.values.dtype, np.integer))]}), "1 integer counts")

# %% [markdown]
# ## 2. The paired 10-cell rule (A1)
#
# **Line:** `src/m1a_pseudobulk.py:97` (`>= min_cells` in both conditions).
#
# **Should do:** keep a donor's two samples in a cell type only if **both** its control and its stimulated group have
# at least `min_cells` cells; otherwise drop **both** (decision A1).
#
# Real threshold, `min_cells=10`. Each cell has one gene with count 1.
#
# | donor | cell type | ctrl cells | stim cells | expected |
# |---|---|---|---|---|
# | A | B cells | 9 | 12 | both dropped (9 < 10; the 12-cell stim sample goes too) |
# | B | B cells | 10 | 10 | both kept (10 is not < 10: boundary) |

# %%
rows2 = ([("A", "ctrl", "B cells", [1])] * 9 + [("A", "stim", "B cells", [1])] * 12
         + [("B", "ctrl", "B cells", [1])] * 10 + [("B", "stim", "B cells", [1])] * 10)
a2 = toy_adata(rows2, genes=["g1"])

# decisions.md A1: "< 10 cells drops that donor from that cell type (both conditions)".
expected = {"B cells|A|ctrl": (9, False), "B cells|A|stim": (12, False),
            "B cells|B|ctrl": (10, True), "B cells|B|stim": (10, True)}

counts2, samples2 = m1a_pseudobulk.pseudobulk(a2, min_cells=10)
s2 = samples2.set_index("sample")
compare(pd.DataFrame({
    "input: sample": list(expected),
    "expected (n_cells, kept)": list(expected.values()),
    "code gives": [(int(s2.loc[i, "n_cells"]), bool(s2.loc[i, "kept"])) for i in expected],
}), "2 A1 rule", expected_col="expected (n_cells, kept)")

# %% [markdown]
# ## 3. The per-cell-type baseMean > 1 gene filter (R)
#
# **Line:** `R/m1a_varpart.R:67`. This is the **only** check of this boundary: the user dropped test T6.
#
# **Should do:** within each cell type, keep a gene for variance partitioning only if its baseMean (mean of
# size-factor-normalised counts over that cell type's kept samples) is **strictly greater** than 1.
#
# To make baseMean hand-computable, the toy forces every DESeq2 size factor to be exactly 1, so baseMean = plain mean
# of the counts. How: size factors use only genes with no zeros, and each is the median, over those genes, of
# count / gene geometric mean. The toy has 20 mirrored pairs (gene `pNNa` has counts x, `pNNb` has 3600 / x, so in
# every sample one is above its geometric mean and the other equally far below) plus `const1` = 1 everywhere (ratio 1).
# 41 genes, 20 above, 20 below, the middle one is 1: every size factor is 1. The boundary genes contain a zero,
# so they don't affect size factors.
#
# Two cell types (T1, T2), each with 4 donors × 2 conditions = 8 kept samples. Counts of the boundary genes (same in both, except `flip`):
#
# | gene | counts in 8 samples | sum | baseMean | expected |
# |---|---|---|---|---|
# | `below` | 1 1 1 1 1 1 1 0 | 7 | 0.875 | out |
# | `exact` | 2 0 2 0 1 1 1 1 | 8 | **1.000** | out (not > 1) |
# | `const1` | 1 1 1 1 1 1 1 1 | 8 | **1.000** | out (not > 1) |
# | `above` | 2 1 2 0 1 1 1 1 | 9 | 1.125 | in |
# | `flip` in T1 | 2 1 2 0 1 1 1 1 | 9 | 1.125 | in T1 |
# | `flip` in T2 | 1 1 1 1 1 1 1 0 | 7 | 0.875 | out of T2 |
#
# `flip` checks "per cell type": one filter over all cell types would put it in both or neither.
# The 40 pair genes all have baseMean ≥ 30, so they are in.

# %%
DONORS, CONDS = ["d1", "d2", "d3", "d4"], ["ctrl", "stim"]
DIVISORS = [30, 36, 40, 45, 48, 50, 60, 72, 75, 80, 90, 100, 120]  # 3600 / x is an integer between 30 and 120
rng = np.random.default_rng(0)


def pair_genes(n_pairs, n_samples):
    g = {}
    for k in range(1, n_pairs + 1):
        x = rng.choice(DIVISORS, size=n_samples)
        g[f"p{k:02d}a"], g[f"p{k:02d}b"] = x, 3600 // x
    return g


def samples_table(counts, cell_type, n_cells):
    ids = list(counts.columns)
    return pd.DataFrame({
        "sample": ids, "cell_type": cell_type,
        "donor": [i.split("|")[1] for i in ids], "condition": [i.split("|")[2] for i in ids],
        "n_cells": n_cells, "total_counts": counts.sum().to_numpy(), "kept": True,
    })


ids_T1 = [f"T1|{d}|{c}" for c in CONDS for d in DONORS]
ids_T2 = [f"T2|{d}|{c}" for c in CONDS for d in DONORS]
pairs = pair_genes(20, 8)
boundary = {"below": [1, 1, 1, 1, 1, 1, 1, 0], "exact": [2, 0, 2, 0, 1, 1, 1, 1],
            "const1": [1] * 8, "above": [2, 1, 2, 0, 1, 1, 1, 1]}
base = pd.DataFrame({**pairs, **boundary}).T
c_T1 = pd.concat([base, pd.DataFrame({"flip": [2, 1, 2, 0, 1, 1, 1, 1]}).T]).set_axis(ids_T1, axis=1)
c_T2 = pd.concat([base, pd.DataFrame({"flip": [1, 1, 1, 1, 1, 1, 1, 0]}).T]).set_axis(ids_T2, axis=1)
counts3 = pd.concat([c_T1, c_T2], axis=1).astype(int)
n_cells3 = [31, 47, 22, 58, 40, 26, 53, 35]
samples3 = pd.concat([samples_table(c_T1.astype(int), "T1", n_cells3),
                      samples_table(c_T2.astype(int), "T2", n_cells3[::-1])], ignore_index=True)

# Toy check (not the code under test): DESeq2 itself gives size factor exactly 1 for every sample of this toy.
c_T1.astype(int).to_csv(TMP / "sf.csv")
chk = subprocess.run([RSCRIPT, "-e",
    f"suppressMessages(library(DESeq2)); m <- as.matrix(read.csv('{TMP}/sf.csv', row.names = 1)); "
    "cat(all(estimateSizeFactorsForMatrix(m) == 1))"], capture_output=True, text=True)
print("toy size factors all exactly 1:", chk.stdout)
assert chk.stdout.strip().endswith("TRUE"), chk.stderr

# Decision: "keep genes with baseMean > 1", per cell type (contract step 4: "Per cell type, keep genes with baseMean > 1").
pair_names = sorted(pairs)
expected_genes = {"T1": sorted(pair_names + ["above", "flip"]), "T2": sorted(pair_names + ["above"])}
expected_basemean = {("T1", "above"): 1.125, ("T1", "flip"): 1.125, ("T2", "above"): 1.125}

rc3, out3, err3 = run_varpart(counts3, samples3, "basemean", min_basemean=1)
print("R exit status:", rc3, "" if rc3 == 0 else err3)
vp3 = pd.read_csv(out3) if out3.exists() else pd.DataFrame(columns=["cell_type", "gene", "baseMean"])
got_genes = {ct: sorted(vp3.loc[vp3["cell_type"] == ct, "gene"]) for ct in ["T1", "T2"]}
rows = []
for ct in ["T1", "T2"]:
    for g in ["below", "exact", "const1", "above", "flip"]:
        rows.append({"input: cell type": ct, "gene": g, "expected": g in expected_genes[ct],
                     "code gives": g in got_genes[ct]})
compare(pd.DataFrame(rows), "3 baseMean boundary genes")
compare(pd.DataFrame({"input": ["T1: the 40 pair genes all in", "T2: the 40 pair genes all in"],
                      "expected": [True, True],
                      "code gives": [set(pair_names) <= set(got_genes[ct]) for ct in ["T1", "T2"]]}),
        "3 pair genes kept")
v3 = vp3.set_index(["cell_type", "gene"])
compare(pd.DataFrame({
    "input": [f"{ct} {g}" for ct, g in expected_basemean],
    "expected": list(expected_basemean.values()),
    "code gives": [round(float(v3.loc[k, "baseMean"]), 9) if k in v3.index else "absent" for k in expected_basemean],
}), "3 baseMean values")

# %% [markdown]
# ## 4. The transform and the variance model, on a toy with a known answer (R)
#
# **Lines:** `R/m1a_varpart.R:71` (blind VST, parametric fit) and `R/m1a_varpart.R:77` (the model, REML/weights).
#
# **Should do:** per cell type, transform counts with DESeq2 `varianceStabilizingTransformation(blind=TRUE,
# fitType="parametric")`, then split each gene's variance into donor, condition, log(n_cells), log(total_counts) and
# residual with `~ (1|donor) + (1|condition) + log(n_cells) + log(total_counts)`, `REML=FALSE`, `useWeights=TRUE`
# (decisions A2, A3 revised, pinned defaults).
#
# One cell type, 4 donors × 2 conditions, on top of 20 random pair genes:
#
# | gene | d1 ctrl | d2 ctrl | d3 ctrl | d4 ctrl | d1 stim | d2 stim | d3 stim | d4 stim | varies by |
# |---|---|---|---|---|---|---|---|---|---|
# | `donorOnly` | 50 | 200 | 800 | 3000 | 52 | 210 | 760 | 2850 | donor (stim ±5%, up for d1–d2, down for d3–d4) |
# | `condOnly` | 20 | 21 | 20 | 21 | 500 | 520 | 510 | 490 | condition |
#
# The two technical covariates are built to share nothing with donor or condition (so they can't take their variance):
# `n_cells` follows a donor-pair × condition pattern, and a `filler` gene sets each sample's total counts to a
# donor × condition pattern (both patterns are orthogonal to the donor and condition main effects).
#
# Expected: `donorOnly` donor fraction > 0.9; `condOnly` condition fraction > 0.9 (the T7 bar); every row's
# fractions ≥ 0 and sum to 1. The fit type actually used is logged. The code's VST values and fractions should also
# equal an independent run of DESeq2 / variancePartition with exactly the decided settings; if, say, a term were missing or
# `blind=FALSE` were used, the numbers would differ.

# %%
ids4 = [f"T1|{d}|{c}" for c in CONDS for d in DONORS]
g4 = pair_genes(20, 8)
g4["donorOnly"] = np.array([50, 200, 800, 3000, 52, 210, 760, 2850])
g4["condOnly"] = np.array([20, 21, 20, 21, 500, 520, 510, 490])
c4 = pd.DataFrame(g4).T.set_axis(ids4, axis=1)
target_total = np.array([22000, 20000, 22000, 20000, 20000, 22000, 20000, 22000])  # (+,-,+,-) ctrl, reversed stim
c4.loc["filler"] = target_total - c4.sum().to_numpy()
c4 = c4.astype(int)
assert (c4.loc["filler"] > 0).all()
n_cells4 = [30, 30, 50, 50, 50, 50, 30, 30]  # (d1,d2 vs d3,d4) x condition
samples4 = samples_table(c4, "T1", n_cells4)
display(samples4)

rc4, out4, err4 = run_varpart(c4, samples4, "model", min_basemean=1)
print("R exit status:", rc4, "" if rc4 == 0 else err4)
vp4 = pd.read_csv(out4).set_index("gene")
FRACS = ["donor", "condition", "log_n_cells", "log_total_counts", "Residuals"]
display(vp4.loc[["donorOnly", "condOnly"], FRACS].round(3))

# decisions A2 + test plan T7: donor-only gene -> donor > 0.9; condition-only gene -> condition > 0.9.
compare(pd.DataFrame({
    "input": ["donorOnly: donor fraction > 0.9", "condOnly: condition fraction > 0.9",
              "every row: all fractions >= 0", "every row: fractions sum to 1 (to 1e-9)"],
    "expected": [True, True, True, True],
    "code gives": [bool(vp4.loc["donorOnly", "donor"] > 0.9), bool(vp4.loc["condOnly", "condition"] > 0.9),
                   bool((vp4[FRACS] >= 0).all().all()), bool(np.allclose(vp4[FRACS].sum(axis=1), 1, atol=1e-9))],
}), "4 known-answer oracle")

# %% [markdown]
# ### 4b. Same toy, independent run of the decided settings
#
# The expected numbers come from a separate R call written from the decisions (not from the M1a code): DESeq2 size
# factors on all genes, keep baseMean > 1, `varianceStabilizingTransformation(blind=TRUE, fitType="parametric")`,
# then `fitExtractVarPartModel` with the decided formula, `REML=FALSE`, `useWeights=TRUE`. For contrast, the same
# with `blind=FALSE` (design `~ donor + condition`) shows the check can tell the two apart.

# %%
indep = TMP / "independent.R"
indep.write_text(f"""
suppressMessages({{library(DESeq2); library(variancePartition)}})
m <- as.matrix(read.csv("{TMP}/model_counts.csv", row.names = 1, check.names = FALSE))
info <- read.csv("{TMP}/model_samples.csv", check.names = FALSE)
info$donor <- factor(info$donor); info$condition <- factor(info$condition)
rownames(info) <- info$sample
dds <- DESeqDataSetFromMatrix(m[, info$sample], info, design = ~ donor + condition)
dds <- estimateSizeFactors(dds, type = "ratio")
dds <- dds[rowMeans(counts(dds, normalized = TRUE)) > 1, ]
vb <- assay(varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric"))
vu <- assay(varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric"))
d1 <- dds; design(d1) <- ~ 1
fit_type <- attr(dispersionFunction(estimateDispersions(d1, fitType = "parametric", quiet = TRUE)), "fitType")
vp <- fitExtractVarPartModel(vb, ~ (1|donor) + (1|condition) + log(n_cells) + log(total_counts), info,
                             REML = FALSE, useWeights = TRUE)
vp <- as.data.frame(vp)
colnames(vp) <- c("donor", "condition", "log_n_cells", "log_total_counts", "Residuals")[
  match(colnames(vp), c("donor", "condition", "log(n_cells)", "log(total_counts)", "Residuals"))]
write.csv(vp, "{TMP}/indep_varpart.csv")
write.csv(vb, "{TMP}/indep_vst_blind.csv")
write.csv(vu, "{TMP}/indep_vst_unblind.csv")
writeLines(fit_type, "{TMP}/indep_fittype.txt")
""")
r = subprocess.run([RSCRIPT, str(indep)], capture_output=True, text=True)
assert r.returncode == 0, r.stderr[-1500:]

exp_vp = pd.read_csv(TMP / "indep_varpart.csv", index_col=0)
exp_fit = (TMP / "indep_fittype.txt").read_text().strip()

got_vst = pd.read_csv(TMP / "model_varpart_vst.csv").pivot(index="gene", columns="sample", values="vst")
got_fit = pd.read_csv(TMP / "model_varpart_fits.csv").set_index("cell_type").loc["T1", "fit_type_used"]
exp_vst = pd.read_csv(TMP / "indep_vst_blind.csv", index_col=0)[list(got_vst.columns)]
unblind = pd.read_csv(TMP / "indep_vst_unblind.csv", index_col=0)[list(got_vst.columns)]
same_genes = sorted(exp_vst.index) == sorted(got_vst.index) == sorted(vp4.index)
max_vst = float((got_vst.loc[exp_vst.index] - exp_vst).abs().max().max())
max_unblind = float((got_vst.loc[unblind.index] - unblind).abs().max().max())
max_frac = float((vp4.loc[exp_vp.index, FRACS] - exp_vp[FRACS]).abs().max().max())
print(f"max |code VST - blind VST| = {max_vst:.2e};  max |code VST - blind=FALSE VST| = {max_unblind:.2e}  (contrast)")
print(f"max |code fraction - independent fraction| = {max_frac:.2e}")
compare(pd.DataFrame({
    "input": ["same genes analysed", "VST equals blind VST (to 1e-6)", "fractions equal independent fit (to 1e-6)",
              "fit type logged"],
    "expected": [True, True, True, exp_fit],
    "code gives": [same_genes, max_vst < 1e-6, max_frac < 1e-6, got_fit],
}), "4b independent re-run")

# %% [markdown]
# ## 5. The pairing check (R)
#
# **Line:** `R/m1a_varpart.R:60` (one kept row per condition per donor; Codex review 1 fix).
#
# **Should do:** stop with a non-zero exit status, writing nothing, if a kept donor does not have exactly one control
# and one stimulated sample in a cell type, even when the number of rows per donor is still 2.
#
# Same toy as section 4, with two condition labels swapped: `d1`'s stim sample relabelled `ctrl` (d1 now has two
# controls) and `d2`'s ctrl sample relabelled `stim` (two stims). Every donor still has 2 rows.

# %%
bad = samples4.copy()
bad.loc[bad["sample"] == "T1|d1|stim", "condition"] = "ctrl"
bad.loc[bad["sample"] == "T1|d2|ctrl", "condition"] = "stim"
display(bad.groupby("donor")["condition"].agg(list))

# Contract assumption: "Every kept sample has both its donor's conditions (paired design)"; failure = stop.
rc5, out5, err5 = run_varpart(c4, bad, "pairing", min_basemean=1)
print("R stderr (tail):", err5[-400:])
compare(pd.DataFrame({
    "input": ["R exit status is non-zero", "output file written"],
    "expected": [True, False],
    "code gives": [rc5 != 0, out5.exists()],
}), "5 pairing check")

# %% [markdown]
# ## Summary

# %%
print("Mismatches (POSSIBLE BUGS):", FAILURES or "none")
assert not FAILURES, FAILURES
