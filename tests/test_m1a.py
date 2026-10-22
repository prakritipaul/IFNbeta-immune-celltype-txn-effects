from pathlib import Path
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.m1a_pseudobulk import pseudobulk, run_m1a


REPO_ROOT = Path(__file__).resolve().parents[1]
REAL_DATA = REPO_ROOT / "data/processed/m0_cells.h5ad"
RSCRIPT = "/usr/local/bin/Rscript"
VARPART_SCRIPT = REPO_ROOT / "R/m1a_varpart.R"


def _adata(counts, obs, genes=None):
    counts = np.asarray(counts)
    genes = genes or [f"gene_{i}" for i in range(counts.shape[1])]
    result = ad.AnnData(
        X=sparse.csr_matrix(counts),
        obs=pd.DataFrame(obs, index=[f"cell_{i}" for i in range(len(obs["donor"]))]),
        var=pd.DataFrame(index=genes),
    )
    result.uns["control"] = "ctrl"
    result.uns["treatment"] = "stim"
    return result


def _matrix_sum(matrix, axis):
    return np.asarray(matrix.sum(axis=axis)).ravel()


def _write_r_inputs(tmp_path, counts, samples):
    counts_path = tmp_path / "counts.csv"
    samples_path = tmp_path / "samples.csv"
    out_path = tmp_path / "varpart.csv"
    counts.to_csv(counts_path, index_label="gene")
    samples.to_csv(samples_path, index=False)
    subprocess.run(
        [RSCRIPT, str(VARPART_SCRIPT), str(counts_path), str(samples_path), str(out_path), "1"],
        check=True,
        capture_output=True,
        text=True,
    )
    return pd.read_csv(out_path)


# Catches treating the 10-cell boundary as too small or dropping only one condition.
def test_T1_a1_boundary_and_paired_drop():
    rows = [("A", "ctrl")] * 9 + [("A", "stim")] * 12 + [("B", "ctrl")] * 10 + [("B", "stim")] * 10
    obs = {
        "donor": [donor for donor, _ in rows],
        "condition": [condition for _, condition in rows],
        "cell_type": ["T"] * len(rows),
    }
    _, samples = pseudobulk(_adata(np.ones((len(rows), 1), dtype=int), obs), min_cells=10)

    kept = samples.set_index(["donor", "condition"])["kept"].to_dict()
    assert kept == {
        ("A", "ctrl"): False,
        ("A", "stim"): False,
        ("B", "ctrl"): True,
        ("B", "stim"): True,
    }


# Catches grouping at the wrong level or applying the small-group rule globally.
@pytest.mark.realdata
def test_T2_real_data_kept_samples():
    adata = ad.read_h5ad(REAL_DATA)
    _, samples = pseudobulk(adata, min_cells=10)

    assert len(samples) == 128
    assert int(samples["kept"].sum()) == 116
    dropped_pairs = {
        (cell_type, donor)
        for cell_type, donor in samples.loc[~samples["kept"], ["cell_type", "donor"]].itertuples(index=False, name=None)
    }
    assert dropped_pairs == {
        ("Dendritic cells", "patient_107"),
        ("Dendritic cells", "patient_1039"),
        ("Megakaryocytes", "patient_101"),
        ("Megakaryocytes", "patient_107"),
        ("Megakaryocytes", "patient_1016"),
        ("Megakaryocytes", "patient_1039"),
    }
    assert (~samples["kept"]).sum() == 2 * len(dropped_pairs)


# Catches cells or counts being lost or duplicated during aggregation.
@pytest.mark.realdata
def test_T3_real_data_conservation():
    adata = ad.read_h5ad(REAL_DATA)
    counts, samples = pseudobulk(adata, min_cells=10)

    np.testing.assert_array_equal(counts.sum(axis=1).to_numpy(), _matrix_sum(adata.X, axis=0))
    assert int(samples["n_cells"].sum()) == 24_562


# Catches a reversed outlier filter, mutation of the main run, or stale kept flags.
@pytest.mark.realdata
def test_T4_no_outlier_variant_and_fresh_a1(tmp_path):
    adata = ad.read_h5ad(REAL_DATA)
    expected_counts, expected_samples = pseudobulk(adata, min_cells=10)
    run_m1a(REAL_DATA, tmp_path, min_cells=10)

    main_counts = pd.read_csv(tmp_path / "m1a_counts.csv", index_col=0)
    main_samples = pd.read_csv(tmp_path / "m1a_samples.csv")
    clean_counts = pd.read_csv(tmp_path / "m1a_counts_no_outliers.csv", index_col=0)
    clean_samples = pd.read_csv(tmp_path / "m1a_samples_no_outliers.csv")
    clean = adata[~adata.obs["flag_outlier"].astype(bool)].copy()
    fresh_counts, fresh_samples = pseudobulk(clean, min_cells=10)

    pd.testing.assert_frame_equal(main_counts, expected_counts, check_dtype=False)
    pd.testing.assert_frame_equal(
        main_samples.sort_values("sample").reset_index(drop=True),
        expected_samples.sort_values("sample").reset_index(drop=True),
        check_dtype=False,
    )
    assert int(main_samples["n_cells"].sum()) == 24_562
    assert int(clean_samples["n_cells"].sum()) == 24_048
    np.testing.assert_array_equal(clean_counts.sum(axis=1).to_numpy(), _matrix_sum(clean.X, axis=0))
    pd.testing.assert_frame_equal(clean_counts, fresh_counts, check_dtype=False)
    pd.testing.assert_series_equal(
        clean_samples.set_index("sample")["kept"].sort_index(),
        fresh_samples.set_index("sample")["kept"].sort_index(),
        check_names=False,
    )

    rows = [("A", "ctrl")] * 10 + [("A", "stim")] * 10 + [("B", "ctrl")] * 10 + [("B", "stim")] * 10
    toy_obs = {
        "donor": [donor for donor, _ in rows],
        "condition": [condition for _, condition in rows],
        "cell_type": ["T"] * len(rows),
        "flag_outlier": [True] + [False] * (len(rows) - 1),
    }
    toy_path = tmp_path / "toy.h5ad"
    toy_out = tmp_path / "toy_output"
    _adata(np.ones((len(rows), 1), dtype=int), toy_obs).write_h5ad(toy_path)
    run_m1a(toy_path, toy_out, min_cells=10)

    toy_main = pd.read_csv(toy_out / "m1a_samples.csv").set_index(["donor", "condition"])
    toy_clean = pd.read_csv(toy_out / "m1a_samples_no_outliers.csv").set_index(["donor", "condition"])
    assert toy_main.loc[["A"], "kept"].all()
    assert not toy_clean.loc[["A"], "kept"].any()


# Catches rounding, missing-metadata workarounds, unpaired inputs, and partial writes.
@pytest.mark.parametrize("bad_case", ["negative", "fractional", "missing_obs", "unpaired"])
def test_T5_invalid_input_writes_nothing(tmp_path, bad_case):
    obs = {
        "donor": ["A", "A"],
        "condition": ["ctrl", "stim"],
        "cell_type": ["T", "T"],
        "flag_outlier": [False, False],
    }
    values = [[1], [2]]
    if bad_case == "negative":
        values[0][0] = -1
    elif bad_case == "fractional":
        values[0][0] = 0.5
    elif bad_case == "missing_obs":
        del obs["cell_type"]
    else:
        obs["condition"] = ["ctrl", "ctrl"]
    input_path = tmp_path / "bad.h5ad"
    output_dir = tmp_path / "output"
    _adata(values, obs).write_h5ad(input_path)

    with pytest.raises(ValueError):
        run_m1a(input_path, output_dir, min_cells=1)
    assert not output_dir.exists() or list(output_dir.iterdir()) == []


# Catches >= versus > at baseMean 1 and filtering genes across cell types globally.
@pytest.mark.r
@pytest.mark.skip(reason="Dropped by the user (2026-10-01): toy genes have zero variance, which variancePartition rejects; fixing the toy judged overkill")
def test_T6_basemean_filter_is_strict_and_per_cell_type(tmp_path):
    sample_rows = []
    columns = []
    for cell_type in ["A", "B"]:
        for donor in ["d1", "d2"]:
            for condition in ["ctrl", "stim"]:
                sample = f"{cell_type}|{donor}|{condition}"
                columns.append(sample)
                sample_rows.append([sample, cell_type, donor, condition, 10, 10_000, True])
    samples = pd.DataFrame(
        sample_rows,
        columns=["sample", "cell_type", "donor", "condition", "n_cells", "total_counts", "kept"],
    )
    rng = np.random.default_rng(20261001)
    filler = []
    for _ in range(60):
        values = rng.negative_binomial(12, 0.2, 4) + 1
        filler.extend(np.roll(values, shift) for shift in range(4))
    filler = np.asarray(filler)
    # Every four rows contain all cyclic shifts, so each sample has the same
    # median-of-ratios distribution and DESeq2 normalizes all size factors to 1.
    data = np.vstack(
        [
            [0] * 8,
            [1] * 8,
            [2] * 8,
            [2] * 4 + [1] * 4,
        ]
        + [np.tile(row, 2) for row in filler]
    )
    counts = pd.DataFrame(
        data,
        index=["below", "exact", "above", "cell_specific"] + [f"filler_{i}" for i in range(len(filler))],
        columns=columns,
    )

    result = _write_r_inputs(tmp_path, counts, samples)
    tested = set(result[["cell_type", "gene"]].itertuples(index=False, name=None))
    assert ("A", "below") not in tested and ("B", "below") not in tested
    assert ("A", "exact") not in tested and ("B", "exact") not in tested
    assert ("A", "above") in tested and ("B", "above") in tested
    assert ("A", "cell_specific") in tested and ("B", "cell_specific") not in tested


# Catches wrong model terms, swapped fractions, or variance partitioning raw counts.
@pytest.mark.r
def test_T7_variance_partition_oracle(tmp_path):
    donors = [f"d{i}" for i in range(1, 9)]
    sample_rows = []
    columns = []
    rng = np.random.default_rng(90210)
    n_cells = rng.permutation(np.arange(17, 33))
    total_counts = rng.permutation(np.arange(80_000, 96_000, 1_000))
    for donor_i, donor in enumerate(donors):
        for condition_i, condition in enumerate(["ctrl", "stim"]):
            sample = f"T|{donor}|{condition}"
            columns.append(sample)
            sample_i = 2 * donor_i + condition_i
            sample_rows.append([sample, "T", donor, condition, n_cells[sample_i], total_counts[sample_i], True])
    samples = pd.DataFrame(
        sample_rows,
        columns=["sample", "cell_type", "donor", "condition", "n_cells", "total_counts", "kept"],
    )
    donor_effect = np.repeat([1_000, 8_000, 3_000, 7_000, 2_000, 6_000, 4_000, 5_000], 2)
    condition_effect = np.tile([1_000, 9_000], 8)
    data = np.vstack(
        [donor_effect + rng.integers(0, 3, 16), condition_effect + rng.integers(0, 3, 16)]
        + [rng.integers(4_000, 6_001, 16) for _ in range(1_100)]
    )
    counts = pd.DataFrame(data, index=["gene_X", "gene_Y"] + [f"filler_{i}" for i in range(1_100)], columns=columns)

    result = _write_r_inputs(tmp_path, counts, samples).set_index("gene")
    assert result.loc["gene_X", "donor"] > 0.9
    assert result.loc["gene_Y", "condition"] > 0.9
    fraction_columns = ["donor", "condition", "log_n_cells", "log_total_counts", "Residuals"]
    assert (result[fraction_columns] >= 0).all().all()
    np.testing.assert_allclose(result[fraction_columns].sum(axis=1), 1.0, atol=1e-6)
