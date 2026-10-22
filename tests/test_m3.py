import hashlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import src.m3_signature as m3


def _adata(counts, genes, *, obs=None):
    counts = np.asarray(counts, dtype=float)
    if obs is None:
        obs = pd.DataFrame(index=[f"cell_{i}" for i in range(counts.shape[0])])
    return ad.AnnData(
        X=counts,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )


def _gmt(path, ifna, ifng=None):
    lines = ["\t".join([m3.IFNA, "toy", *ifna])]
    if ifng is not None:
        lines.append("\t".join([m3.IFNG, "toy", *ifng]))
    path.write_text("\n".join(lines) + "\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _scoring_adata():
    genes = ["A", "S", "G"] + [f"B{i:02d}" for i in range(100)]
    row = np.arange(1, len(genes) + 1, dtype=float)
    return _adata(np.vstack([row, row[::-1] + 3]), genes)


# Catches wrong normalisation or failure to subtract matched control-gene expression.
def test_T1_constant_within_cell_scores_zero():
    genes = ["A", "S", "G"] + [f"B{i:02d}" for i in range(100)]
    adata = _adata([[7] * len(genes), [23] * len(genes)], genes)

    result = m3.score_cells(adata, {m3.IFNA: ["A", "S"], m3.IFNG: ["G", "S"]})

    assert list(result.columns) == ["ifna_score", "ifng_score", "ifng_noshared_score"]
    assert result.index.equals(adata.obs_names)
    np.testing.assert_allclose(result.to_numpy(), 0, atol=1e-12)


# Catches scoring raw library-size-dependent counts instead of normalized log counts.
def test_T2_proportional_cells_have_identical_scores():
    genes = ["A", "S", "G"] + [f"B{i:02d}" for i in range(100)]
    row = np.arange(1, len(genes) + 1, dtype=float)
    adata = _adata(np.vstack([row, 10 * row]), genes)

    result = m3.score_cells(adata, {m3.IFNA: ["A", "S"], m3.IFNG: ["G", "S"]})

    np.testing.assert_allclose(result.iloc[0], result.iloc[1], rtol=0, atol=1e-12)


# Catches unseeded scoring or reliance on result-changing score_genes defaults.
def test_T3_deterministic_and_explicit_scoring_parameters(monkeypatch):
    calls = []
    original = sc.tl.score_genes

    def recording_score_genes(*args, **kwargs):
        calls.append(kwargs.copy())
        return original(*args, **kwargs)

    monkeypatch.setattr(sc.tl, "score_genes", recording_score_genes)
    gene_sets = {m3.IFNA: ["A", "S"], m3.IFNG: ["G", "S"]}

    first = m3.score_cells(_scoring_adata(), gene_sets)
    second = m3.score_cells(_scoring_adata(), gene_sets)

    pd.testing.assert_frame_equal(first, second)
    assert len(calls) == 6
    for call in calls:
        assert call["ctrl_size"] == 50
        assert call["n_bins"] == 25
        assert call["random_state"] == 0


# Catches incorrect overlap detection or removal from the IFN-gamma no-shared score.
def test_T4_ifng_shared_gene_subtraction(tmp_path, monkeypatch):
    gmt = tmp_path / "toy.gmt"
    digest = _gmt(gmt, ["A", "S"], ["G", "S"])
    loaded = m3.load_gene_sets(gmt, digest)
    seen_gene_lists = []

    def fake_score_genes(adata, gene_list, *, score_name, **kwargs):
        seen_gene_lists.append((score_name, list(gene_list)))
        adata.obs[score_name] = 0.0

    monkeypatch.setattr(sc.tl, "score_genes", fake_score_genes)
    result = m3.score_cells(_scoring_adata(), loaded)

    assert loaded == {m3.IFNA: ["A", "S"], m3.IFNG: ["G", "S"]}
    assert [genes for _, genes in seen_gene_lists] == [["A", "S"], ["G", "S"], ["G"]]
    assert list(result.columns) == ["ifna_score", "ifng_score", "ifng_noshared_score"]


# Catches omitted approved mappings, fuzzy matching, or silent acceptance of absent genes.
def test_T5_exact_six_mappings_and_unapproved_missing_gene_fails(tmp_path):
    source = list(m3.GENE_MAPPING)
    target = list(m3.GENE_MAPPING.values())
    sets = {m3.IFNA: source[:3], m3.IFNG: source[3:]}

    mapped, coverage = m3.map_gene_sets(sets, target)

    assert mapped == {m3.IFNA: target[:3], m3.IFNG: target[3:]}
    expected_coverage = pd.DataFrame(
        {
            "set": [m3.IFNA, m3.IFNG],
            "n_genes": [3, 3],
            "n_found": [3, 3],
            "pct_found": [100.0, 100.0],
        }
    )
    pd.testing.assert_frame_equal(
        coverage.reset_index(drop=True), expected_coverage, check_dtype=False
    )

    gmt = tmp_path / "missing.gmt"
    digest = _gmt(gmt, source[:3] + ["NOT_APPROVED"], source[3:])
    input_path = tmp_path / "input.h5ad"
    _adata([[1] * len(target)], target).write_h5ad(input_path)
    out_dir = tmp_path / "outputs"

    with pytest.raises(ValueError, match="NOT_APPROVED"):
        m3.run_m3(input_path, gmt, out_dir, expected_sha256=digest)
    assert not out_dir.exists() or list(out_dir.iterdir()) == []


# Catches wrong-release input, absent Hallmark sets, or NaN scores leaking to disk.
def test_T6_invalid_gmt_or_nan_stops_without_outputs(tmp_path, monkeypatch):
    input_path = tmp_path / "input.h5ad"
    obs = pd.DataFrame(
        {"donor": ["D1"], "condition": ["ctrl"], "cell_type": ["T"]},
        index=["bc1"],
    )
    _adata([[1, 2, 3]], ["A", "S", "G"], obs=obs).write_h5ad(input_path)

    valid_gmt = tmp_path / "valid.gmt"
    valid_digest = _gmt(valid_gmt, ["A", "S"], ["G", "S"])
    wrong_hash_out = tmp_path / "wrong_hash"
    with pytest.raises(ValueError, match="sha256|SHA-?256|hash"):
        m3.run_m3(input_path, valid_gmt, wrong_hash_out, expected_sha256="0" * 64)
    assert not wrong_hash_out.exists() or list(wrong_hash_out.iterdir()) == []

    missing_gmt = tmp_path / "missing_set.gmt"
    missing_digest = _gmt(missing_gmt, ["A", "S"])
    missing_set_out = tmp_path / "missing_set"
    with pytest.raises(ValueError, match=m3.IFNG):
        m3.run_m3(input_path, missing_gmt, missing_set_out, expected_sha256=missing_digest)
    assert not missing_set_out.exists() or list(missing_set_out.iterdir()) == []

    def nan_scores(adata, gene_sets, **params):
        return pd.DataFrame(
            {
                "ifna_score": [np.nan],
                "ifng_score": [0.0],
                "ifng_noshared_score": [0.0],
            },
            index=adata.obs_names,
        )

    monkeypatch.setattr(m3, "score_cells", nan_scores)
    nan_out = tmp_path / "nan"
    with pytest.raises(ValueError, match="NaN|nan"):
        m3.run_m3(input_path, valid_gmt, nan_out, expected_sha256=valid_digest)
    assert not nan_out.exists() or list(nan_out.iterdir()) == []


# Catches a <= boundary error or retention of an unpaired donor-condition row.
def test_T7_min_cells_boundary_drops_both_donor_rows():
    rows = []
    for donor, condition, n_cells in [
        ("A", "ctrl", 9),
        ("A", "stim", 12),
        ("B", "ctrl", 10),
        ("B", "stim", 10),
    ]:
        rows.extend((donor, condition, f"{donor}_{condition}_{i}") for i in range(n_cells))
    obs = pd.DataFrame(rows, columns=["donor", "condition", "barcode"]).set_index("barcode")
    obs["cell_type"] = "T"
    scores = pd.DataFrame(
        1.0,
        index=obs.index,
        columns=["ifna_score", "ifng_score", "ifng_noshared_score"],
    )

    result = m3.summarise_samples(scores, obs, min_cells=10)

    assert result[["donor", "condition"]].to_records(index=False).tolist() == [
        ("B", "ctrl"),
        ("B", "stim"),
    ]
    assert result["n_cells"].tolist() == [10, 10]


# Catches median aggregation, an incorrect cell count, or ctrl-minus-stim subtraction.
def test_T8_mean_summary_and_treatment_minus_control_difference():
    obs = pd.DataFrame(
        {
            "donor": ["A"] * 6,
            "condition": ["ctrl"] * 3 + ["stim"] * 3,
            "cell_type": ["T"] * 6,
        },
        index=[f"bc{i}" for i in range(6)],
    )
    values = [1.0, 2.0, 9.0, 4.0, 5.0, 6.0]
    scores = pd.DataFrame(
        {name: values for name in ["ifna_score", "ifng_score", "ifng_noshared_score"]},
        index=obs.index,
    )

    summary = m3.summarise_samples(scores, obs, min_cells=1)
    ctrl = summary.loc[summary["condition"] == "ctrl"].iloc[0]
    assert ctrl["n_cells"] == 3
    assert ctrl["ifna_score"] == 4.0

    paired_input = pd.DataFrame(
        {
            "cell_type": ["T", "T"],
            "donor": ["A", "A"],
            "condition": ["ctrl", "stim"],
            "n_cells": [3, 3],
            "ifna_score": [2.0, 5.0],
            "ifng_score": [2.0, 5.0],
            "ifng_noshared_score": [2.0, 5.0],
        }
    )
    differences = m3.paired_differences(paired_input)
    for column in ["ifna_score", "ifng_score", "ifng_noshared_score"]:
        assert differences.loc[0, column] == 3.0


# Catches lost/duplicated real cells or application of the wrong paired-group filter.
@pytest.mark.realdata
def test_T9_real_data_cell_and_sample_accounting(tmp_path):
    input_path = "data/processed/m0_cells.h5ad"
    gmt_path = "data/msigdb/h.all.v2026.1.Hs.symbols.gmt"
    out_dir = tmp_path / "m3"

    m3.run_m3(input_path, gmt_path, out_dir)

    cell_scores = pd.read_csv(out_dir / "m3_cell_scores.csv")
    sample_scores = pd.read_csv(out_dir / "m3_sample_scores.csv")
    input_barcodes = ad.read_h5ad(input_path, backed="r").obs_names.astype(str)
    assert len(cell_scores) == 24_562
    assert cell_scores["barcode"].is_unique
    assert set(cell_scores["barcode"]) == set(input_barcodes)
    assert len(sample_scores) == 116
    paired_donors = sample_scores.groupby("cell_type")["donor"].nunique()
    assert paired_donors["Dendritic cells"] == 6
    assert paired_donors["Megakaryocytes"] == 4
