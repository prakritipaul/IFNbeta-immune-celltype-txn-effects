"""M0: load Kang 2018, attach the authors' labels, QC metrics, basic filters, QC flags.

Contract: audit/contracts/M0_eda.md. Tests: tests/test_m0.py (written blind by Codex).
Every threshold is an argument; invalid input raises ValueError and nothing is written.
"""
import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

GEO_COLUMNS = ["tsne1", "tsne2", "ind", "stim", "cluster", "cell", "multiplets"]
N_DONORS = 8
N_CELL_TYPES = 8
CONTROL_NAMES = {"ctrl", "control"}  # the one condition label NOT in this set is the treatment
# Genes matching the RPS/RPL prefix that are NOT ribosomal proteins (user decision D4) → 84 ribosomal genes
NON_RIBOSOMAL = {"RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KB1", "RPS6KB2", "RPS6KC1",
                 "RPS6KL1", "RPS10P7", "RPSAP58", "RPL34-AS1", "RPL17-C18orf32", "RPS19BP1"}
TSNE_ATOL = 1e-9

# The real Kang 2018 inputs, checked on every real run (run_m0_kang)
_ROOT = Path(__file__).resolve().parents[1]
KANG_H5AD = _ROOT / "data/kang_2018.h5ad"
KANG_GEO = _ROOT / "data/geo/GSE96583_batch2.total.tsne.df.tsv.gz"
KANG_H5AD_SHA256 = "e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830"
KANG_GEO_SHA256 = "1d57e72e92ca8695250e88cc0f1c3fa8c0be1175d974f8b427c58f1274dc6c09"
KANG_CELL_TYPES = ["B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells",
                   "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells"]  # h5ad vs GEO coordinates differ by ~4e-15 (float parsing)


def _values(X):
    """Stored values of a count matrix.

    In:  X, sparse or dense matrix.
    Out: 1-D array of the stored values (for sparse, the non-zeros only).
    """
    return X.data if sp.issparse(X) else np.asarray(X).ravel()


def _core(barcodes):
    """Barcode sequence without its "-<digits>" suffix (h5ad uses -1/-2, GEO -1/-11).

    In:  barcodes, e.g. ["AAA-1", "AAA-2", "AAA-11"].
    Out: pd.Index of sequences, e.g. ["AAA", "AAA", "AAA"].
    """
    return pd.Index(barcodes).str.replace(r"-\d+$", "", regex=True)


def _check_sha256(path, expected, name):
    """Stop if a file's fingerprint differs from the recorded one.

    In:  path to the file; expected sha256 (None = skip the check); name used in the error.
    Out: None. Raises ValueError on mismatch.
    """
    if expected is not None:
        digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        if digest != expected:
            raise ValueError(f"{name} sha256 {digest} != expected {expected}")


def _donor_id(values):
    """Donor IDs in a common format for comparison.

    In:  donor labels, e.g. ["patient_1016"] (h5ad) or [1016] (GEO).
    Out: array of strings without the "patient_" prefix, e.g. ["1016"].
    """
    return pd.Series(values).astype(str).str.replace(r"^patient_", "", regex=True).to_numpy()


def load_inputs(h5ad_path, geo_path, expected_sha256=None, expected_geo_sha256=None):
    """Load the counts and the authors' metadata, checking both.

    In:  h5ad_path, geo_path; expected sha256 fingerprints (None = skip).
    Out: (adata, geo): AnnData with raw counts in X; DataFrame of GEO metadata indexed by barcode.
    Routine: 1) check fingerprints  2) read h5ad; stop unless X is finite, non-negative integers
             3) read GEO table; stop if any required column is missing.
    """
    _check_sha256(h5ad_path, expected_sha256, "h5ad")
    _check_sha256(geo_path, expected_geo_sha256, "GEO")
    adata = ad.read_h5ad(h5ad_path)
    v = _values(adata.X)
    if not np.isfinite(v).all():
        raise ValueError("X has non-finite values (inf or nan); raw counts required")
    if (v < 0).any():
        raise ValueError("X has negative values; raw counts required")
    if not np.array_equal(v, np.round(v)):
        raise ValueError("X has non-integer values; raw counts required")
    geo = pd.read_csv(geo_path, sep="\t", index_col=0)
    missing = set(GEO_COLUMNS) - set(geo.columns)
    if missing:
        raise ValueError(f"GEO metadata missing columns: {sorted(missing)}")
    return adata, geo


def join_author_metadata(adata, geo):
    """Attach the authors' labels on (barcode sequence, condition). Same cells, same order.

    In:  adata (obs: label, replicate, cluster, tsne1, tsne2; index = barcode);
         geo (index = barcode; stim, ind, cluster, cell, tsne1, tsne2, multiplets).
    Out: copy of adata whose obs is exactly: barcode, condition, donor, cell_type (authors'),
         cluster, tsne1, tsne2.
    Routine: 1) build keys (sequence, condition) for both tables; stop on duplicate GEO keys
             2) look up each cell's GEO row; stop if any cell is unmatched
             3) cross-check donor, cluster, t-SNE (to 1e-9) and singlet status; stop on any mismatch
             4) build the new obs.
    """
    geo_key = pd.MultiIndex.from_arrays([_core(geo.index), geo["stim"].astype(str)])
    if geo_key.duplicated().any():
        raise ValueError(f"GEO has {geo_key.duplicated().sum()} duplicated (barcode, condition) keys")
    cell_key = pd.MultiIndex.from_arrays([_core(adata.obs_names), adata.obs["label"].astype(str)])
    pos = geo_key.get_indexer(cell_key)
    if (pos < 0).any():
        raise ValueError(f"{(pos < 0).sum()} of {adata.n_obs} cells have no GEO row")
    g = geo.iloc[pos]
    obs = adata.obs
    if not (_donor_id(obs["replicate"]) == _donor_id(g["ind"])).all():
        raise ValueError("donor disagrees between h5ad and GEO")
    if not (obs["cluster"].to_numpy() == g["cluster"].to_numpy()).all():
        raise ValueError("cluster disagrees between h5ad and GEO")
    for c in ("tsne1", "tsne2"):
        if not np.allclose(obs[c].to_numpy(float), g[c].to_numpy(float), rtol=0, atol=TSNE_ATOL):
            raise ValueError(f"{c} disagrees between h5ad and GEO")
    if not (g["multiplets"] == "singlet").all():
        raise ValueError("non-singlet cells present")

    out = adata.copy()
    out.obs = pd.DataFrame(
        {
            "barcode": adata.obs_names.to_numpy(),
            "condition": obs["label"].astype(str).to_numpy(),
            "donor": obs["replicate"].astype(str).to_numpy(),
            "cell_type": g["cell"].astype(str).to_numpy(),
            "cluster": obs["cluster"].to_numpy(),
            "tsne1": obs["tsne1"].to_numpy(float),
            "tsne2": obs["tsne2"].to_numpy(float),
        },
        index=adata.obs_names,
    )
    return out


def ribosomal_mask(var_names):
    """Which genes count as ribosomal.

    In:  gene names, e.g. ["RPL13", "RPS6KA1", "MX1"].
    Out: boolean array: RPS/RPL prefix and not in NON_RIBOSOMAL, e.g. [True, False, False].
    """
    names = pd.Index(var_names)
    return np.asarray(names.str.match(r"^RP[SL]") & ~names.isin(NON_RIBOSOMAL))


def compute_qc_metrics(adata):
    """Per-cell and per-gene QC metrics from raw counts.

    In:  adata with raw counts in X.
    Out: copy with obs n_counts (total counts), n_genes (genes > 0), pct_ribo (% counts in ribosomal genes, see ribosomal_mask)
         and var n_cells (cells with count > 0).
    """
    out = adata.copy()
    X = out.X
    out.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
    out.obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
    ribo = ribosomal_mask(out.var_names)
    ribo_counts = np.asarray(X[:, ribo].sum(axis=1)).ravel()
    out.obs["pct_ribo"] = 100 * ribo_counts / out.obs["n_counts"].to_numpy()
    out.var["n_cells"] = np.asarray((X > 0).sum(axis=0)).ravel()
    return out


def basic_filter(adata, min_genes=200, min_cells=3):
    """Cells with >= min_genes genes first, then genes in >= min_cells of the kept cells (all conditions).

    In:  adata; min_genes (per cell); min_cells (per gene).
    Out: filtered copy.
    Routine: 1) keep cells with >= min_genes genes detected
             2) among kept cells, keep genes detected in >= min_cells cells (both conditions pooled).
    """
    n_genes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    out = adata[n_genes >= min_genes].copy()
    n_cells = np.asarray((out.X > 0).sum(axis=0)).ravel()
    return out[:, n_cells >= min_cells].copy()


def mad_outlier(values, n_mads=5.0):
    """True where |v - median| > n_mads * raw (unscaled) MAD, within one group.

    In:  values (one group); n_mads.
    Out: boolean array, same length. E.g. [-1,-1,0,0,1,1,5] → MAD 1, so 5 is not flagged (5 is not > 5).
    """
    v = np.asarray(values, dtype=float)
    med = np.median(v)
    mad = np.median(np.abs(v - med))
    # ponytail: MAD = 0 (all-identical group) flags any deviating value; never seen on real data
    return np.abs(v - med) > n_mads * mad


def flag_qc(adata, n_mads=5.0, min_cells_group=10, min_paired_donors=3):
    """Add boolean flag columns from obs n_counts/n_genes. Never removes cells.

    In:  adata with obs donor, condition, cell_type, n_counts, n_genes; thresholds.
    Out: copy with obs flag_outlier, flag_small_group, flag_underpowered_celltype.
    Routine: 1) flag_outlier: > n_mads MADs above OR below the median on log1p(n_counts) or log1p(n_genes),
                within donor × condition × cell_type; groups < min_cells_group cells are not assessed
             2) flag_small_group: donor × condition × cell_type group has < min_cells_group cells
             3) flag_underpowered_celltype, per cell type: count the "paired donors" (donors with
                >= min_cells_group cells in ctrl AND in stim); flag if that count < min_paired_donors.
                E.g. megakaryocytes: 4 paired donors, 4 >= 3 → not flagged.
    """
    out = adata.copy()
    obs = out.obs
    outlier = np.zeros(out.n_obs, dtype=bool)
    for _, idx in obs.groupby(["donor", "condition", "cell_type"], observed=True).indices.items():
        if len(idx) < min_cells_group:
            continue  # too few cells for a stable median/MAD (user decision D2)
        for metric in ("n_counts", "n_genes"):
            outlier[idx] |= mad_outlier(np.log1p(obs[metric].to_numpy(float)[idx]), n_mads)
    obs["flag_outlier"] = outlier

    group_size = obs.groupby(["donor", "condition", "cell_type"], observed=True)["donor"].transform("size")
    obs["flag_small_group"] = (group_size < min_cells_group).to_numpy()

    # paired donor = both conditions present with >= min_cells_group cells ("left after the rule above")
    big = obs[~obs["flag_small_group"]]
    conds = big.groupby(["cell_type", "donor"], observed=True)["condition"].nunique()
    paired = (conds == 2).groupby(level="cell_type").sum().reindex(obs["cell_type"].unique(), fill_value=0)
    obs["flag_underpowered_celltype"] = (obs["cell_type"].map(paired) < min_paired_donors).to_numpy()
    return out


def infer_conditions(conditions):
    """(control, treatment): exactly 2 labels, exactly one of them in CONTROL_NAMES.

    In:  condition labels, e.g. ["ctrl", "stim", "stim"].
    Out: (control, treatment), e.g. ("ctrl", "stim"). Raises ValueError otherwise (e.g. "a"/"b").
    """
    labels = set(pd.Series(conditions).astype(str))
    control = labels & CONTROL_NAMES
    if len(labels) != 2 or len(control) != 1:
        raise ValueError(f"need exactly one control label {sorted(CONTROL_NAMES)} and one treatment; found {sorted(labels)}")
    (ctrl,) = control
    (treatment,) = labels - control
    return ctrl, treatment


def validate_assumptions(adata, expected_cell_types=None):
    """Stop unless the data matches the contract's assumptions.

    In:  adata with obs condition, donor, cell_type; expected_cell_types (None = count only).
    Out: None. Raises ValueError unless: a valid control/treatment pair; 8 donors, each with both
         conditions; unique gene names; 8 cell types (with the expected names, if given).
    """
    obs = adata.obs
    infer_conditions(obs["condition"])
    if obs["donor"].nunique() != N_DONORS:
        raise ValueError(f"expected {N_DONORS} donors, found {obs['donor'].nunique()}")
    per_donor = obs.groupby("donor", observed=True)["condition"].nunique()
    if (per_donor != 2).any():
        raise ValueError(f"donors missing a condition: {per_donor[per_donor != 2].index.tolist()}")
    if not adata.var_names.is_unique:
        raise ValueError("duplicate gene names")
    if obs["cell_type"].nunique() != N_CELL_TYPES:
        raise ValueError(f"expected {N_CELL_TYPES} cell types, found {obs['cell_type'].nunique()}")
    if expected_cell_types is not None and set(obs["cell_type"].astype(str)) != set(expected_cell_types):
        raise ValueError(f"cell-type names differ from expected: {sorted(set(obs['cell_type']) ^ set(expected_cell_types))}")


def _record(step, before, after, **extra):
    """One cell-accounting row: {step, n_cells_in/out, n_genes_in/out, **extra}."""
    return {"step": step, "n_cells_in": before.n_obs, "n_cells_out": after.n_obs,
            "n_genes_in": before.n_vars, "n_genes_out": after.n_vars, **extra}


def run_m0(h5ad_path, geo_path, out_dir, expected_sha256=None, expected_geo_sha256=None,
           expected_cell_types=None, audit_dir=None, min_genes=200, min_cells=3,
           n_mads=5.0, min_cells_group=10, min_paired_donors=3):
    """All M0 steps. Everything is validated before anything is written.

    In:  input paths; out_dir (for m0_cells.h5ad); audit_dir (for cell_flow.md; default out_dir);
         expected fingerprints and cell-type names; all thresholds.
    Out: list of cell-accounting rows (one per step). Writes m0_cells.h5ad and cell_flow.md.
    Routine: check thresholds → load → join → validate → QC metrics → basic filter → validate → flags
             → check no step gained cells → write outputs.
    """
    thresholds = {"min_genes": min_genes, "min_cells": min_cells, "n_mads": n_mads,
                  "min_cells_group": min_cells_group, "min_paired_donors": min_paired_donors}
    non_finite = sorted(k for k, v in thresholds.items() if not np.isfinite(v))
    if non_finite:
        raise ValueError(f"non-finite thresholds: {non_finite}")
    bad = {k: v for k, v in {"min_genes": min_genes < 0, "min_cells": min_cells < 0, "n_mads": n_mads <= 0,
                              "min_cells_group": min_cells_group < 1, "min_paired_donors": min_paired_donors < 1}.items() if v}
    if bad:
        raise ValueError(f"invalid thresholds: {sorted(bad)}")
    raw, geo = load_inputs(h5ad_path, geo_path, expected_sha256, expected_geo_sha256)
    joined = join_author_metadata(raw, geo)
    validate_assumptions(joined, expected_cell_types)
    qc = compute_qc_metrics(joined)
    filtered = basic_filter(qc, min_genes, min_cells)
    validate_assumptions(filtered, expected_cell_types)
    flagged = flag_qc(filtered, n_mads, min_cells_group, min_paired_donors)
    flagged.uns["control"], flagged.uns["treatment"] = infer_conditions(flagged.obs["condition"])
    flagged.var["n_cells"] = np.asarray((flagged.X > 0).sum(axis=0)).ravel()  # after filtering

    flow = [
        _record("load", raw, raw),
        _record("join", raw, joined, n_matched=joined.n_obs),
        _record("qc_metrics", joined, qc),
        _record("basic_filter", qc, filtered),
        _record("flag_qc", filtered, flagged),
    ]
    for r in flow:
        if r["n_cells_out"] > r["n_cells_in"]:
            raise ValueError(f"step {r['step']} gained cells")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    flagged.write_h5ad(out_dir / "m0_cells.h5ad")
    rows = ["| step | cells in | cells out | genes in | genes out | matched |", "|---|---|---|---|---|---|"]
    rows += [f"| {r['step']} | {r['n_cells_in']} | {r['n_cells_out']} | {r['n_genes_in']} | "
             f"{r['n_genes_out']} | {r.get('n_matched', '')} |" for r in flow]
    audit_dir = Path(audit_dir) if audit_dir is not None else out_dir
    audit_dir.mkdir(parents=True, exist_ok=True)
    (audit_dir / "cell_flow.md").write_text("# M0 cell flow\n\n" + "\n".join(rows) + "\n")
    return flow


def run_m0_kang(out_dir, audit_dir=None, **thresholds):
    """run_m0 on the real Kang 2018 files, always checking both fingerprints and the 8 cell-type names.

    In:  out_dir, audit_dir; optional thresholds (min_genes, min_cells, n_mads, min_cells_group, min_paired_donors).
    Out: as run_m0. The KANG_* constants are read at call time, so they can't be skipped.
    """
    return run_m0(KANG_H5AD, KANG_GEO, out_dir, expected_sha256=KANG_H5AD_SHA256,
                  expected_geo_sha256=KANG_GEO_SHA256, expected_cell_types=KANG_CELL_TYPES,
                  audit_dir=audit_dir, **thresholds)
