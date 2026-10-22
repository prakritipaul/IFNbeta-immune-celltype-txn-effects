"""M1b: pseudobulk DE per cell type (R DESeq2 + qvalue), main run plus sensitivity and negative-control reruns.

Contract: audit/contracts/M1b_de.md (decisions O1-O7). Tests: tests/test_m1b.py (written blind by Codex).
The fit itself is R/m1b_deseq2.R; this module validates inputs, builds each run's sample set, calls R,
and validates and collects the results. Invalid input raises ValueError before R runs or anything is written.
"""
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from src.m0_qc import infer_conditions

_ROOT = Path(__file__).resolve().parents[1]
R_SCRIPT = _ROOT / "R/m1b_deseq2.R"
RSCRIPT = "/usr/local/bin/Rscript"
RESULT_COLUMNS = ["run", "cell_type", "gene", "baseMean", "log2FC", "lfcSE", "ci_low", "ci_high", "pvalue", "qvalue",
                  "tested", "cooks_flag", "log2FC_shrunk", "significant", "de"]
FIT_COLUMNS = ["run", "cell_type", "fit_type_used", "n_samples", "n_genes_tested"]
ALL_RUNS = ("main", "no_outliers", "unpaired", "lodo", "shuffle")
Z95 = 1.96  # 95% Wald CI (contract)


def classify(results, q_threshold=0.05, lfc_threshold=1.0):
    """Add the significance and DE calls (decisions: significant = q < 0.05; DE = significant and |log2FC| ≥ 1).

    In:  results with qvalue, log2FC (unshrunk; log2FC_shrunk is never used). Rows with NaN qvalue need a
         documented reason: tested False, cooks_flag True, or baseMean 0 (all-zero gene).
    Out: copy with `significant` (bool) and `de` ("up"/"down"/"none").
         E.g. (q 0.0499, LFC 1.01) → up; (q 0.05, LFC 2) → none; (q 0.01, LFC 1) → significant, none.
    """
    out = results.copy()
    nan_q = out["qvalue"].isna()
    if nan_q.any():
        missing = [c for c in ("tested", "cooks_flag", "baseMean") if c not in out.columns]
        if missing:
            raise ValueError(f"NaN q-values need tested/cooks_flag/baseMean to be explained; missing {missing}")
        explained = ~out["tested"].astype(bool) | out["cooks_flag"].astype(bool) | (out["baseMean"] == 0)
        if (nan_q & ~explained).any():
            raise ValueError(f"{int((nan_q & ~explained).sum())} NaN q-values without a documented reason "
                             "(untested, Cook's outlier or all-zero gene)")
    out["significant"] = (out["qvalue"] < q_threshold).to_numpy()  # NaN < x is False
    big = out["significant"] & (out["log2FC"].abs() >= lfc_threshold)
    out["de"] = np.where(big & (out["log2FC"] > 0), "up", np.where(big & (out["log2FC"] < 0), "down", "none"))
    return out


def balanced_shuffles(donors):
    """Every distinct way to swap ctrl/stim for exactly half of the donors (O4); a set and its complement count once.

    In:  donors (even number). Out: list of frozensets of donors to swap. 8 donors → 35; 6 → 10; 4 → 3.
         E.g. ["a", "b", "c", "d"] → [{a,b}, {a,c}, {a,d}] (the set containing the first donor stands for its complement).
    """
    donors = sorted(donors)
    if len(donors) % 2 or not donors:
        raise ValueError(f"balanced shuffles need an even, non-zero number of donors, got {len(donors)}")
    return [frozenset(s) for s in combinations(donors, len(donors) // 2) if donors[0] in s]


def _validate(counts, samples):
    """Stop unless counts/samples meet the contract. Out: (control, treatment), kept samples.

    Checks: required sample columns; unique ids that match the counts columns exactly; counts numeric,
    finite, non-negative integers; two condition labels (one control); every kept donor × cell type has
    exactly one sample per condition.
    """
    missing = [c for c in ("sample", "cell_type", "donor", "condition", "kept") if c not in samples.columns]
    if missing:
        raise ValueError(f"samples missing columns: {missing}")
    if samples["sample"].duplicated().any() or counts.columns.duplicated().any():
        raise ValueError("duplicated sample ids")
    if set(samples["sample"]) != set(counts.columns):
        raise ValueError("sample ids differ between counts and samples")
    if counts.index.duplicated().any():
        raise ValueError("duplicated gene names in counts")
    if samples["kept"].dtype != bool:
        raise ValueError("kept must be boolean")
    if not all(pd.api.types.is_numeric_dtype(t) and not pd.api.types.is_bool_dtype(t) for t in counts.dtypes):
        raise ValueError("counts must be numeric (raw integer counts)")
    v = counts.to_numpy(dtype=float)
    if not np.isfinite(v).all() or (v < 0).any() or not np.array_equal(v, np.round(v)):
        raise ValueError("counts must be finite non-negative integers (raw counts)")
    kept = samples[samples["kept"]]
    if kept.empty:
        raise ValueError("no kept samples")
    control, treatment = infer_conditions(kept["condition"])
    per = kept.groupby(["cell_type", "donor"])["condition"].agg(lambda c: sorted(c.astype(str)))
    bad = per[per.map(lambda c: c != sorted([control, treatment]))]
    if len(bad):
        raise ValueError(f"pairing violated: each kept donor × cell type needs exactly one {control} and one "
                         f"{treatment} sample; offending: {bad.index.tolist()[:5]}")
    return control, treatment, kept


def run_deseq2(counts, samples, design, run_name, min_basemean=1.0):
    """One R/m1b_deseq2.R call: DESeq2 + qvalue per cell type on the kept samples.

    In:  counts genes × samples (int); samples: sample, cell_type, donor, condition, kept;
         design "paired" (~ donor + condition) or "unpaired" (~ condition); run_name; min_basemean (O2).
    Out: (results, fits). results: RESULT_COLUMNS, one row per cell type × gene (all genes fitted);
         ci = log2FC ∓ 1.96·lfcSE; tested = baseMean > min_basemean; calls from classify().
         fits: FIT_COLUMNS, one row per cell type.
    Routine: 1) validate  2) write kept samples to a temp dir  3) Rscript (non-zero exit → raise, with R's output)
             4) check every kept cell type × gene came back once, sample counts match  5) CI, classify.
    """
    if design not in ("paired", "unpaired"):
        raise ValueError(f"design must be 'paired' or 'unpaired', got {design!r}")
    _, _, kept = _validate(counts, samples)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        counts[kept["sample"]].to_csv(tmp / "counts.csv", index_label="gene")
        kept[["sample", "cell_type", "donor", "condition", "kept"]].to_csv(tmp / "samples.csv", index=False)
        cmd = [RSCRIPT, str(R_SCRIPT), str(tmp / "counts.csv"), str(tmp / "samples.csv"), design,
               str(tmp / "results.csv"), str(tmp / "fits.csv"), repr(float(min_basemean))]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"R/m1b_deseq2.R failed for run {run_name!r}:\n{proc.stdout}\n{proc.stderr}")
        res = pd.read_csv(tmp / "results.csv", keep_default_na=False, na_values=["NA"])
        fits = pd.read_csv(tmp / "fits.csv")

    cell_types = sorted(kept["cell_type"].unique())
    if sorted(fits["cell_type"]) != cell_types or sorted(res["cell_type"].unique()) != cell_types:
        raise ValueError(f"run {run_name}: R returned cell types {sorted(fits['cell_type'])}, expected {cell_types}")
    n_kept = kept.groupby("cell_type").size()
    if not (fits.set_index("cell_type")["n_samples"] == n_kept).all():
        raise ValueError(f"run {run_name}: R fitted a different number of samples than were kept")
    if res.duplicated(["cell_type", "gene"]).any() or (res.groupby("cell_type").size() != len(counts)).any():
        raise ValueError(f"run {run_name}: R did not return every gene exactly once per cell type")
    genes = set(counts.index)
    for ct, got in res.groupby("cell_type")["gene"]:
        missing, extra = genes - set(got), set(got) - genes
        if missing or extra:
            raise ValueError(f"run {run_name}, {ct}: R genes differ from the input; "
                             f"missing {sorted(missing)[:5]}, extra {sorted(extra)[:5]}")
    res["run"] = run_name
    fits["run"] = run_name
    res["ci_low"] = res["log2FC"] - Z95 * res["lfcSE"]
    res["ci_high"] = res["log2FC"] + Z95 * res["lfcSE"]
    return classify(res)[RESULT_COLUMNS], fits[FIT_COLUMNS]


def _read_inputs(in_dir, suffix=""):
    counts = pd.read_csv(Path(in_dir) / f"m1a_counts{suffix}.csv", index_col=0)
    samples = pd.read_csv(Path(in_dir) / f"m1a_samples{suffix}.csv")
    _validate(counts, samples)
    return counts, samples


def _plan_runs(counts, samples, clean_counts, clean_samples, runs):
    """(run_name, counts, samples, design) for every requested run. Pure: nothing is fitted or written.

    lodo_<donor>: main samples minus that donor, in every cell type where it is paired.
    shuffle_<k>:  in each cell type, the k-th balanced split of its paired donors (k = 1..35/10/3) has its
                  two samples' condition labels swapped; donors never change. Cell types with fewer than k splits
                  are absent from shuffle_<k>.
    """
    unknown = set(runs) - set(ALL_RUNS)
    if unknown:
        raise ValueError(f"unknown runs {sorted(unknown)}; allowed {ALL_RUNS}")
    control, treatment, kept = _validate(counts, samples)
    plan = []
    if "main" in runs:
        plan.append(("main", counts, samples, "paired"))
    if "no_outliers" in runs:
        plan.append(("no_outliers", clean_counts, clean_samples, "paired"))
    if "unpaired" in runs:
        plan.append(("unpaired", counts, samples, "unpaired"))
    if "lodo" in runs:
        for donor in sorted(kept["donor"].unique()):
            paired_in = set(kept.loc[kept["donor"] == donor, "cell_type"])
            s = samples.copy()
            s["kept"] = s["kept"] & (s["donor"] != donor) & s["cell_type"].isin(paired_in)
            plan.append((f"lodo_{donor}", counts, s, "paired"))
    if "shuffle" in runs:
        splits = {ct: balanced_shuffles(g["donor"].unique()) for ct, g in kept.groupby("cell_type")}
        for k in range(1, max(len(v) for v in splits.values()) + 1):
            s = samples.copy()
            s["kept"] = s["kept"] & s["cell_type"].map(lambda ct: len(splits[ct]) >= k if ct in splits else False)
            for ct, ct_splits in splits.items():
                if len(ct_splits) >= k:
                    swap = s["kept"] & (s["cell_type"] == ct) & s["donor"].isin(ct_splits[k - 1])
                    s.loc[swap, "condition"] = s.loc[swap, "condition"].map({control: treatment, treatment: control})
            plan.append((f"shuffle_{k}", counts, s, "paired"))
    return plan


def run_m1b(in_dir, out_dir, runs=ALL_RUNS, min_basemean=1.0, max_workers=8):
    """All M1b runs from the M1a files. Inputs are validated before any R call; nothing is written on failure.

    In:  in_dir with m1a_counts.csv, m1a_samples.csv (+ _no_outliers variants); out_dir; runs (subset of ALL_RUNS);
         min_basemean (O2); max_workers: R processes in parallel (each fit is independent and deterministic).
    Out: None. Writes m1b_de.csv (main run), m1b_runs.csv (every run), m1b_fits.csv (run × cell type).
    """
    counts, samples = _read_inputs(in_dir)
    clean_counts, clean_samples = _read_inputs(in_dir, "_no_outliers")
    plan = _plan_runs(counts, samples, clean_counts, clean_samples, runs)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        done = list(pool.map(lambda p: run_deseq2(p[1], p[2], p[3], p[0], min_basemean=min_basemean), plan))

    results = pd.concat([r for r, _ in done], ignore_index=True)[RESULT_COLUMNS]
    fits = pd.concat([f for _, f in done], ignore_index=True)[FIT_COLUMNS]
    if results.duplicated(["run", "cell_type", "gene"]).any():
        raise ValueError("duplicated (run, cell_type, gene) rows")
    for name, c, s, _ in plan:
        expected = sorted(s.loc[s["kept"], "cell_type"].unique())
        got = results[results["run"] == name]
        if sorted(got["cell_type"].unique()) != expected or sorted(fits.loc[fits["run"] == name, "cell_type"]) != expected:
            raise ValueError(f"run {name}: cell types lost or added")
        if (got.groupby("cell_type")["gene"].nunique() != len(c)).any():
            raise ValueError(f"run {name}: genes lost or added")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results.sort_values(["run", "cell_type", "gene"], inplace=True)
    results[results["run"] == "main"].to_csv(out_dir / "m1b_de.csv", index=False)
    results.to_csv(out_dir / "m1b_runs.csv", index=False)
    fits.sort_values(["run", "cell_type"]).to_csv(out_dir / "m1b_fits.csv", index=False)
