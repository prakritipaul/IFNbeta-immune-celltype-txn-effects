"""M3: Hallmark IFN-alpha / IFN-gamma signature scores per cell, summarised per donor x cell type x condition.

Contract: audit/contracts/M3_signature.md. Tests: tests/test_m3.py (written blind by Codex).
Decisions B1-B4 and A1 (audit/decisions.md). Every result-affecting parameter is passed explicitly;
invalid input raises ValueError and nothing is written.
"""
import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from src import m0_qc

IFNA = "HALLMARK_INTERFERON_ALPHA_RESPONSE"
IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
MSIGDB_GMT = "data/msigdb/h.all.v2026.1.Hs.symbols.gmt"
MSIGDB_SHA256 = "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596"
# B4: MSigDB 2026 symbol -> our (2017) var_name
GENE_MAPPING = {"TENT5A": "FAM46A", "WARS1": "WARS", "MARCHF1": "MARCH1", "RIGI": "DDX58",
                "TMT1B": "METTL7B", "PSMA2": "PSMA2.1"}
SCORE_COLUMNS = ["ifna_score", "ifng_score", "ifng_noshared_score"]


def load_gene_sets(gmt_path, expected_sha256, names=(IFNA, IFNG)):
    """The named gene sets from a GMT file whose fingerprint matches.

    In:  gmt_path; expected_sha256 (required); names of the sets to return.
    Out: {name: [genes in file order]}, e.g. {IFNA: ["A", "S"], IFNG: ["G", "S"]}.
    Routine: 1) stop on sha256 mismatch  2) parse "name<TAB>url<TAB>genes..." lines  3) stop if a set is missing.
    """
    digest = hashlib.sha256(Path(gmt_path).read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise ValueError(f"GMT sha256 {digest} != expected {expected_sha256} (wrong MSigDB release?)")
    sets = {}
    for line in Path(gmt_path).read_text().splitlines():
        fields = line.split("\t")
        sets[fields[0]] = fields[2:]
    missing = [n for n in names if n not in sets]
    if missing:
        raise ValueError(f"gene sets missing from GMT: {missing}")
    return {n: sets[n] for n in names}


def map_gene_sets(gene_sets, var_names, mapping=GENE_MAPPING):
    """Rename set genes to our symbols (approved mapping only) and require every gene to be present.

    In:  gene_sets {name: [genes]}; var_names of the data; mapping {set symbol: our symbol}.
    Out: (mapped sets, coverage DataFrame with columns set, n_genes, n_found, pct_found).
         E.g. {IFNA: ["WARS1"]} with var_names ["WARS"] -> {IFNA: ["WARS"]}, coverage 1/1 = 100.0.
    Raises ValueError naming any gene still absent (no fuzzy matching).
    """
    present = set(var_names)
    mapped = {name: [mapping.get(g, g) for g in genes] for name, genes in gene_sets.items()}
    rows, missing = [], {}
    for name, genes in mapped.items():
        absent = [g for g in genes if g not in present]
        if absent:
            missing[name] = absent
        n_found = len(genes) - len(absent)
        rows.append({"set": name, "n_genes": len(genes), "n_found": n_found,
                     "pct_found": 100.0 * n_found / len(genes)})
    if missing:
        raise ValueError(f"set genes absent from var_names after the approved mapping: {missing}")
    return mapped, pd.DataFrame(rows)


def score_cells(adata, gene_sets, target_sum=1e4, ctrl_size=50, n_bins=25, random_state=0):
    """Per-cell IFN-alpha, IFN-gamma and IFN-gamma-without-shared-genes scores (B1). adata is not modified.

    In:  adata with raw counts in X; gene_sets {IFNA: [...], IFNG: [...]} already mapped; scoring parameters.
    Out: DataFrame indexed by obs_names with columns ifna_score, ifng_score, ifng_noshared_score.
    Routine: 1) copy X as float32 into a bare AnnData  2) normalize_total(target_sum) + natural log1p
             3) sc.tl.score_genes three times: IFNA, IFNG, IFNG minus genes also in IFNA (IFNG order kept)
             4) stop on any NaN.
    """
    work = ad.AnnData(X=adata.X.astype(np.float32, copy=True),  # decision: float32 working matrix
                      obs=pd.DataFrame(index=adata.obs_names.copy()),
                      var=pd.DataFrame(index=adata.var_names.copy()))
    sc.pp.normalize_total(work, target_sum=target_sum, exclude_highly_expressed=False, inplace=True)
    sc.pp.log1p(work, base=None)
    ifna = set(gene_sets[IFNA])
    lists = {"ifna_score": list(gene_sets[IFNA]),
             "ifng_score": list(gene_sets[IFNG]),
             "ifng_noshared_score": [g for g in gene_sets[IFNG] if g not in ifna]}
    for score_name, genes in lists.items():
        sc.tl.score_genes(work, genes, ctrl_as_ref=True, ctrl_size=ctrl_size, gene_pool=None,
                          n_bins=n_bins, score_name=score_name, random_state=random_state,
                          copy=False, use_raw=False, layer=None)
    scores = pd.DataFrame({c: work.obs[c].to_numpy(dtype=float) for c in SCORE_COLUMNS}, index=adata.obs_names)
    _check_no_nan(scores)
    return scores


def _check_no_nan(scores):
    """Stop if any score is NaN (contract: one score per cell per set, no NaNs)."""
    n_nan = scores[SCORE_COLUMNS].isna().sum()
    if n_nan.any():
        raise ValueError(f"NaN scores: {n_nan[n_nan > 0].to_dict()}")


def summarise_samples(cell_scores, obs, min_cells=10):
    """Mean score per donor x cell type x condition, keeping only paired donors (A1, B3).

    In:  cell_scores (index = barcodes, SCORE_COLUMNS); obs with donor, condition, cell_type (same index);
         min_cells per donor x cell type x condition.
    Out: one row per kept sample: cell_type, donor, condition, n_cells, mean of each score; sorted.
         A donor x cell type is kept only if EVERY condition in obs has >= min_cells cells (missing = 0),
         e.g. donor A 9 ctrl / 12 stim -> both rows dropped; donor B 10 / 10 -> both kept.
    """
    if not cell_scores.index.equals(obs.index):
        raise ValueError("cell_scores and obs must have the same cells in the same order")
    keys = ["cell_type", "donor", "condition"]
    df = pd.concat([obs[keys].astype(str), cell_scores[SCORE_COLUMNS]], axis=1)
    grouped = df.groupby(keys)
    samples = grouped[SCORE_COLUMNS].mean()
    samples.insert(0, "n_cells", grouped.size())
    # missing condition = 0 cells (A1, interface spec)
    counts = samples["n_cells"].unstack("condition", fill_value=0).reindex(
        columns=sorted(df["condition"].unique()), fill_value=0)
    paired = counts[(counts >= min_cells).all(axis=1)].index
    keep = samples.index.droplevel("condition").isin(paired)
    return samples[keep].reset_index()


def paired_differences(sample_scores, control="ctrl", treatment="stim"):
    """Treatment minus control per donor x cell type, for each score.

    In:  sample_scores from summarise_samples; control/treatment labels.
    Out: one row per donor x cell type: cell_type, donor, then treatment - control of each score
         (e.g. ctrl 2, stim 5 -> 3). Raises if a donor x cell type lacks either condition.
    """
    wide = sample_scores.set_index(["cell_type", "donor", "condition"])[SCORE_COLUMNS].unstack("condition")
    for label in (control, treatment):
        if label not in wide.columns.get_level_values("condition"):
            raise ValueError(f"condition {label!r} absent from sample_scores")
    diff = wide.xs(treatment, axis=1, level="condition") - wide.xs(control, axis=1, level="condition")
    unpaired = diff.index[diff.isna().any(axis=1)].tolist()
    if unpaired:
        raise ValueError(f"unpaired donor x cell type (need both conditions): {unpaired}")
    return diff[SCORE_COLUMNS].reset_index()


def _gene_sets_md(gene_sets, mapped, coverage, mapping):
    """Markdown for audit/m3_gene_sets.md: set sizes, coverage, mappings applied, shared genes."""
    applied = [(name, src, mapping[src]) for name, genes in gene_sets.items() for src in genes if src in mapping]
    shared = [g for g in mapped[IFNG] if g in set(mapped[IFNA])]
    lines = ["# M3 gene sets", "", f"Source: `{MSIGDB_GMT}` (sha256 `{MSIGDB_SHA256}`).", "",
             "## Coverage (after mapping)", "", coverage.to_markdown(index=False), "",
             f"## Approved mappings applied (B4): {len(applied)} set entries, "
             f"{len({s for _, s, _ in applied})} genes", "",
             "| set | MSigDB symbol | our symbol |", "|---|---|---|"]
    lines += [f"| {n} | {s} | {t} |" for n, s, t in applied]
    lines += ["", f"## Genes in both sets: {len(shared)}", "",
              f"IFN-γ no-shared score uses {len(mapped[IFNG]) - len(shared)} of {len(mapped[IFNG])} IFN-γ genes.", "",
              ", ".join(shared), ""]
    return "\n".join(lines)


def run_m3(in_path, gmt_path, out_dir, expected_sha256=MSIGDB_SHA256, audit_dir=None,
           target_sum=1e4, ctrl_size=50, n_bins=25, random_state=0, min_cells=10):
    """All M3 steps. Everything is validated before anything is written.

    In:  in_path (M0 output h5ad, raw counts); gmt_path; out_dir; expected GMT sha256;
         audit_dir (for m3_gene_sets.md; None = not written); all scoring / summary parameters.
    Out: None. Writes out_dir/m3_cell_scores.csv and out_dir/m3_sample_scores.csv.
    Routine: gene sets (sha256, both present) -> read h5ad, raw counts -> map genes (100% or stop)
             -> score -> no NaN -> validate obs (M0 assumptions; uns control/treatment) -> summarise
             (A1, both conditions per kept sample) -> write.
    """
    gene_sets = load_gene_sets(gmt_path, expected_sha256)
    adata = ad.read_h5ad(in_path)
    v = m0_qc._values(adata.X)
    if not (np.isfinite(v).all() and (v >= 0).all() and np.array_equal(v, np.round(v))):
        raise ValueError("X must hold finite non-negative integer raw counts")
    mapped, coverage = map_gene_sets(gene_sets, adata.var_names)
    scores = score_cells(adata, mapped, target_sum=target_sum, ctrl_size=ctrl_size, n_bins=n_bins,
                         random_state=random_state)
    _check_no_nan(scores)  # also guards a replaced score_cells

    m0_qc.validate_assumptions(adata)
    control, treatment = m0_qc.infer_conditions(adata.obs["condition"])
    if "control" in adata.uns and (adata.uns["control"], adata.uns["treatment"]) != (control, treatment):
        raise ValueError(f"uns control/treatment {adata.uns['control']}/{adata.uns['treatment']} "
                         f"disagree with condition labels {control}/{treatment}")
    samples = summarise_samples(scores, adata.obs, min_cells=min_cells)
    paired_differences(samples, control=control, treatment=treatment)  # asserts every kept sample is paired

    cells = pd.concat([adata.obs[["donor", "condition", "cell_type"]].astype(str), scores], axis=1)
    cells.insert(0, "barcode", adata.obs_names.astype(str))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cells.to_csv(out_dir / "m3_cell_scores.csv", index=False)
    samples.to_csv(out_dir / "m3_sample_scores.csv", index=False)
    if audit_dir is not None:
        Path(audit_dir).mkdir(parents=True, exist_ok=True)
        (Path(audit_dir) / "m3_gene_sets.md").write_text(_gene_sets_md(gene_sets, mapped, coverage, GENE_MAPPING))
