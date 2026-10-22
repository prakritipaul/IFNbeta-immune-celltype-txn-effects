"""M1a: pseudobulk raw counts per cell type × donor × condition, plus per-sample covariates.

Contract: audit/contracts/M1a_pseudobulk.md. Tests: tests/test_m1a.py (written blind by Codex).
Variance partitioning is in R/m1a_varpart.R, which reads the files written here.
Every threshold is an argument; invalid input raises ValueError and nothing is written.
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.m0_qc import KANG_CELL_TYPES, _values

M1A_INPUT = "data/processed/m0_cells.h5ad"
REQUIRED_OBS = ["donor", "condition", "cell_type"]
SAMPLE_COLUMNS = ["sample", "cell_type", "donor", "condition", "n_cells", "total_counts", "kept"]
# Contract: the real input's shape and group counts, checked by run_m1a_kang
KANG_SHAPE = (24_562, 15_701)
N_DONORS, N_CELL_TYPES = 8, 8

_ROOT = Path(__file__).resolve().parents[1]


def validate_input(adata):
    """Stop unless adata matches the contract's assumptions.

    In:  adata.
    Out: (control, treatment) labels from adata.uns. Raises ValueError unless: obs has donor,
         condition, cell_type; uns has control and treatment, and every condition label is one of
         them; X is finite, non-negative integers; every donor has cells in both conditions overall.
    """
    missing = [c for c in REQUIRED_OBS if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"obs missing columns: {missing}")
    if "control" not in adata.uns or "treatment" not in adata.uns:
        raise ValueError("uns must hold 'control' and 'treatment' labels")
    control, treatment = str(adata.uns["control"]), str(adata.uns["treatment"])
    extra = set(adata.obs["condition"].astype(str)) - {control, treatment}
    if extra:
        raise ValueError(f"condition labels not in uns control/treatment: {sorted(extra)}")
    v = _values(adata.X)
    if not np.isfinite(v).all():
        raise ValueError("X has non-finite values; raw counts required")
    if (v < 0).any():
        raise ValueError("X has negative values; raw counts required")
    if not np.array_equal(v, np.round(v)):
        raise ValueError("X has non-integer values; raw counts required")
    piped = [c for c in REQUIRED_OBS if adata.obs[c].astype(str).str.contains("|", regex=False).any()]
    if piped:
        raise ValueError(f"labels contain '|', the sample-id delimiter, in obs columns: {piped}")
    obs = adata.obs.astype({"donor": str, "condition": str})
    per_donor = obs.groupby("donor")["condition"].nunique()
    if (per_donor != 2).any():
        raise ValueError(f"donors missing a condition overall (paired design): {per_donor[per_donor != 2].index.tolist()}")
    return control, treatment


def pseudobulk(adata, min_cells=10):
    """Sum raw counts per cell type × donor × condition and apply the small-group rule (A1).

    In:  adata (raw counts in X; obs donor, condition, cell_type; uns control, treatment);
         min_cells: a donor × cell type is kept only if both conditions have >= min_cells cells.
    Out: (counts, samples).
         counts: genes × samples, int64, index "gene", columns = sample ids "cell_type|donor|condition".
         samples: one row per sample (same order as counts' columns): sample, cell_type, donor,
         condition, n_cells, total_counts, kept.
         E.g. donor A with 9 ctrl / 12 stim cells → both A samples kept=False (min_cells=10).
    Routine: 1) validate  2) cells × samples indicator matrix; counts = X.T @ indicator
             3) n_cells, total_counts per sample  4) kept: both conditions >= min_cells (missing = 0)
             5) assert gene totals and cell counts conserved, kept samples paired.
    """
    if not min_cells >= 1:
        raise ValueError(f"min_cells must be >= 1, got {min_cells}")
    control, treatment = validate_input(adata)
    obs = adata.obs[REQUIRED_OBS].astype(str)
    sample_of_cell = obs["cell_type"] + "|" + obs["donor"] + "|" + obs["condition"]
    codes, sample_ids = pd.factorize(sample_of_cell, sort=True)
    indicator = sp.csr_matrix((np.ones(len(codes), dtype=np.int64), (np.arange(len(codes)), codes)),
                              shape=(len(codes), len(sample_ids)))
    X = sp.csr_matrix(adata.X).astype(np.int64)  # integers checked in validate_input
    summed = np.asarray((X.T @ indicator).todense())
    counts = pd.DataFrame(summed, index=pd.Index(adata.var_names, name="gene"), columns=list(sample_ids))

    parts = pd.Series(sample_ids).str.split("|", expand=True)
    samples = pd.DataFrame({
        "sample": list(sample_ids),
        "cell_type": parts[0].to_numpy(),
        "donor": parts[1].to_numpy(),
        "condition": parts[2].to_numpy(),
        "n_cells": np.bincount(codes, minlength=len(sample_ids)),
        "total_counts": summed.sum(axis=0),
    })
    n = samples.pivot_table(index=["cell_type", "donor"], columns="condition", values="n_cells", aggfunc="sum")
    n = n.reindex(columns=[control, treatment]).fillna(0)  # approved by contract: a missing condition counts as 0 cells
    pair_ok = (n >= min_cells).all(axis=1)
    samples["kept"] = pd.MultiIndex.from_frame(samples[["cell_type", "donor"]]).map(pair_ok.to_dict()).astype(bool)

    if not np.array_equal(summed.sum(axis=1), np.asarray(X.sum(axis=0)).ravel()):
        raise ValueError("pseudobulk gene totals differ from the cell-level totals")
    if samples["n_cells"].sum() != adata.n_obs:
        raise ValueError("pseudobulk n_cells does not sum to the number of cells")
    kept_conds = samples[samples["kept"]].groupby(["cell_type", "donor"])["condition"].nunique()
    if (kept_conds != 2).any():
        raise ValueError("a kept donor × cell type lacks one condition (paired design)")
    return counts, samples[SAMPLE_COLUMNS]


def _flow_rows(step, n_in, counts, samples):
    """One cell-accounting row for audit/cell_flow.md."""
    kept = samples[samples["kept"]]
    return (f"| {step} | {n_in} | {int(samples['n_cells'].sum())} | {len(samples)} | {len(kept)} | "
            f"{int(kept['n_cells'].sum())} | {counts.shape[0]} |")


def run_m1a(in_path, out_dir, min_cells=10, audit_dir=None):
    """Pseudobulk with and without the flag_outlier cells. Everything is validated before anything is written.

    In:  in_path (M0 output h5ad, needs obs flag_outlier); out_dir; min_cells (A1);
         audit_dir for cell_flow.md (default out_dir).
    Out: None. Writes m1a_counts.csv, m1a_samples.csv, m1a_counts_no_outliers.csv,
         m1a_samples_no_outliers.csv to out_dir, and replaces the "M1a" section of cell_flow.md.
    Routine: 1) read, validate  2) pseudobulk all cells  3) drop flag_outlier cells, pseudobulk afresh (A1 re-applied)
             4) check counts  5) write.
    """
    adata = ad.read_h5ad(in_path)
    if "flag_outlier" not in adata.obs.columns or adata.obs["flag_outlier"].dtype != bool:
        raise ValueError("obs needs a boolean flag_outlier column")
    counts, samples = pseudobulk(adata, min_cells=min_cells)
    clean = adata[~adata.obs["flag_outlier"].to_numpy()].copy()
    clean_counts, clean_samples = pseudobulk(clean, min_cells=min_cells)
    n_flagged = int(adata.obs["flag_outlier"].sum())
    if clean.n_obs != adata.n_obs - n_flagged:
        raise ValueError("outlier removal did not remove exactly the flagged cells")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "m1a_counts.csv")
    samples.to_csv(out_dir / "m1a_samples.csv", index=False)
    clean_counts.to_csv(out_dir / "m1a_counts_no_outliers.csv")
    clean_samples.to_csv(out_dir / "m1a_samples_no_outliers.csv", index=False)

    section = "\n".join([
        "# M1a cell flow (pseudobulk)",
        "",
        f"min_cells = {min_cells} (A1). Asserted: per-gene sums over samples = cell-level sums; "
        "n_cells sums to cells in; every kept donor × cell type has both conditions.",
        "",
        "| step | cells in | cells in samples | samples | kept samples | cells in kept samples | genes |",
        "|---|---|---|---|---|---|---|",
        _flow_rows("pseudobulk (all cells)", adata.n_obs, counts, samples),
        f"| drop flag_outlier | {adata.n_obs} | {clean.n_obs} | | | | {clean.n_vars} |",
        _flow_rows("pseudobulk (no outliers)", clean.n_obs, clean_counts, clean_samples),
        "",
    ])
    audit_dir = Path(audit_dir) if audit_dir is not None else out_dir
    audit_dir.mkdir(parents=True, exist_ok=True)
    flow = audit_dir / "cell_flow.md"
    old = flow.read_text() if flow.exists() else ""
    old = old.split("# M1a cell flow")[0].rstrip()
    flow.write_text((old + "\n\n" if old else "") + section)


def run_m1a_kang(out_dir=_ROOT / "data/processed", audit_dir=_ROOT / "audit", min_cells=10):
    """run_m1a on the real M0 output, first checking the contract's shape, donors, cell types and conditions.

    In:  out_dir, audit_dir, min_cells.
    Out: as run_m1a. Raises ValueError if the input is not the expected M0 output.
    """
    path = _ROOT / M1A_INPUT
    a = ad.read_h5ad(path, backed="r")
    checks = {
        "shape": (a.shape, KANG_SHAPE),
        "donors": (a.obs["donor"].nunique(), N_DONORS),
        "cell types": (sorted(a.obs["cell_type"].astype(str).unique()), sorted(KANG_CELL_TYPES)),
        "conditions": (a.obs["condition"].nunique(), 2),
    }
    a.file.close()
    bad = {k: v for k, v in checks.items() if v[0] != v[1]}
    if bad:
        raise ValueError(f"input differs from the contract (found, expected): {bad}")
    run_m1a(path, out_dir, min_cells=min_cells, audit_dir=audit_dir)
