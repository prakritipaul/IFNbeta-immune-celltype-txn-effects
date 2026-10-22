import anndata as ad
import numpy as np
import pandas as pd
import pytest

from src import m0_qc


def test_t19_load_inputs_rejects_non_finite_counts(tmp_path):
    for bad_count in (np.inf, np.nan):
        h5ad_path = tmp_path / f"counts-{bad_count}.h5ad"
        ad.AnnData(np.array([[bad_count]])).write_h5ad(h5ad_path)

        with pytest.raises(ValueError):
            m0_qc.load_inputs(h5ad_path, tmp_path / "unused.tsv")


@pytest.mark.parametrize(
    ("threshold", "value"),
    [
        ("min_genes", -1),
        ("min_cells", -1),
        ("n_mads", 0),
        ("min_cells_group", 0),
        ("min_paired_donors", 0),
        ("n_mads", np.nan),
        ("n_mads", np.inf),
        ("min_genes", np.nan),
    ],
)
def test_t20_run_m0_rejects_invalid_thresholds_before_writing(
    tmp_path, monkeypatch, threshold, value
):
    toy = ad.AnnData(
        np.ones((2, 1), dtype=int),
        obs=pd.DataFrame({"condition": ["ctrl", "stim"]}),
    )
    monkeypatch.setattr(m0_qc, "load_inputs", lambda *args, **kwargs: (toy, pd.DataFrame()))
    monkeypatch.setattr(m0_qc, "join_author_metadata", lambda *args: toy)
    monkeypatch.setattr(m0_qc, "validate_assumptions", lambda *args: None)
    monkeypatch.setattr(m0_qc, "compute_qc_metrics", lambda adata: adata)
    monkeypatch.setattr(m0_qc, "basic_filter", lambda adata, *args: adata)
    monkeypatch.setattr(m0_qc, "flag_qc", lambda adata, *args: adata)
    out_dir = tmp_path / "output"

    with pytest.raises(ValueError):
        m0_qc.run_m0("unused.h5ad", "unused.tsv", out_dir, **{threshold: value})

    assert not out_dir.exists()


@pytest.mark.realdata
def test_t21_run_m0_kang_always_checks_fingerprints_and_cell_types(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(m0_qc, "KANG_H5AD_SHA256", "0" * 64)

    with pytest.raises(ValueError):
        m0_qc.run_m0_kang(tmp_path / "output", tmp_path / "audit")


def test_t18_outlier_uses_n_genes_as_well_as_n_counts():
    # Catches outlier detection silently using only n_counts and ignoring n_genes.
    obs = pd.DataFrame(
        {
            "donor": ["d1"] * 10,
            "condition": ["ctrl"] * 10,
            "cell_type": ["T"] * 10,
            "n_counts": [1000] * 10,
            "n_genes": [1, 2, 3, 4, 5, 6, 7, 8, 9, 100_000],
        }
    )
    adata = ad.AnnData(np.ones((10, 1), dtype=int), obs=obs)

    result = m0_qc.flag_qc(
        adata, n_mads=5.0, min_cells_group=10, min_paired_donors=3
    )

    expected = np.array([False] * 9 + [True])
    np.testing.assert_array_equal(result.obs["flag_outlier"].to_numpy(), expected)


def test_t16_outliers_use_cell_type_both_tails_and_skip_small_groups():
    # Catches pooled cell types, one-sided MAD flags, and outlier assessment of groups below 10 cells.
    dc_counts = [900, 950, 975, 1000, 1000, 1000, 1025, 1050, 1075, 1100]
    rows = (
        [("d1", "ctrl", "T", 10)] * 20
        + [("d1", "ctrl", "DC", value) for value in dc_counts]
        + [("d2", "ctrl", "B", value) for value in [100, 101, 99, 102, 98, 100, 101, 99, 100, 1]]
        + [("d3", "ctrl", "Mega", value) for value in [100] * 8 + [10_000]]
    )
    obs = pd.DataFrame(rows, columns=["donor", "condition", "cell_type", "n_counts"])
    obs["n_genes"] = 5
    adata = ad.AnnData(np.ones((len(obs), 1), dtype=int), obs=obs)

    result = m0_qc.flag_qc(
        adata, n_mads=5.0, min_cells_group=10, min_paired_donors=3
    )

    expected = np.zeros(len(obs), dtype=bool)
    expected[39] = True
    np.testing.assert_array_equal(result.obs["flag_outlier"].to_numpy(), expected)


@pytest.mark.realdata
def test_t17_ribosomal_definition_excludes_named_non_ribosomal_genes():
    # Catches kinases, pseudogenes, antisense, and binding proteins counted as ribosomal genes.
    included = ["RPL13", "RPS6", "RPLP0", "RPL22L1"]
    excluded = ["RPS6KA1", "RPSAP58", "RPL34-AS1", "RPS19BP1"]
    var_names = included + excluded
    toy = ad.AnnData(
        np.ones((1, len(var_names)), dtype=int),
        var=pd.DataFrame(index=var_names),
    )

    expected_mask = np.array([True] * 4 + [False] * 4)
    np.testing.assert_array_equal(m0_qc.ribosomal_mask(toy.var_names), expected_mask)
    assert m0_qc.compute_qc_metrics(toy).obs["pct_ribo"].iat[0] == 50.0

    real = ad.read_h5ad("data/kang_2018.h5ad", backed="r")
    try:
        assert int(m0_qc.ribosomal_mask(real.var_names).sum()) == 84
    finally:
        real.file.close()
