import anndata as ad
import numpy as np
import pandas as pd
import pytest

from src import m0_qc


def _assumption_toy():
    donors = [f"patient_{i}" for i in range(8) for _ in range(2)]
    conditions = [condition for _ in range(8) for condition in ("ctrl", "stim")]
    cell_types = [f"type_{i}" for i in range(8) for _ in range(2)]
    obs = pd.DataFrame(
        {"donor": donors, "condition": conditions, "cell_type": cell_types},
        index=[f"cell_{i}" for i in range(16)],
    )
    return ad.AnnData(
        np.ones((16, 1), dtype=int),
        obs=obs,
        var=pd.DataFrame(index=["G1"]),
    )


def _run_toy(tmp_path):
    donors = [f"patient_{i}" for i in range(8) for _ in range(2)]
    conditions = [condition for _ in range(8) for condition in ("ctrl", "stim")]
    cell_types = [f"type_{i}" for i in range(8) for _ in range(2)]
    barcodes = [f"BC{i:02d}-{'1' if condition == 'ctrl' else '2'}" for i, condition in enumerate(conditions)]
    obs = pd.DataFrame(
        {
            "label": conditions,
            "replicate": donors,
            "cluster": np.arange(16),
            "tsne1": np.arange(16, dtype=float),
            "tsne2": np.arange(16, dtype=float) + 0.5,
        },
        index=barcodes,
    )
    adata = ad.AnnData(
        np.ones((16, 200), dtype=int),
        obs=obs,
        var=pd.DataFrame(index=[f"G{i}" for i in range(200)]),
    )
    geo = pd.DataFrame(
        {
            "stim": conditions,
            "ind": donors,
            "cluster": np.arange(16),
            "cell": cell_types,
            "tsne1": obs["tsne1"].to_numpy(),
            "tsne2": obs["tsne2"].to_numpy(),
            "multiplets": ["singlet"] * 16,
        },
        index=[barcode.rsplit("-", 1)[0] + "-1" for barcode in barcodes],
    )
    h5ad_path = tmp_path / "input.h5ad"
    geo_path = tmp_path / "geo.tsv.gz"
    adata.write_h5ad(h5ad_path)
    geo.to_csv(geo_path, sep="\t", compression="gzip")
    return h5ad_path, geo_path, cell_types


# Catches counting paired donors without first excluding donor-condition groups below 10 cells.
def test_T11_underpowered_celltype_uses_viable_groups():
    rows = []
    for donor, size in (("d1", 10), ("d2", 10), ("d3", 9)):
        for condition in ("ctrl", "stim"):
            rows.extend((donor, condition, "target") for _ in range(size))
    obs = pd.DataFrame(rows, columns=["donor", "condition", "cell_type"])
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    obs["n_counts"] = 10
    obs["n_genes"] = 2
    adata = ad.AnnData(
        np.ones((len(obs), 2), dtype=int),
        obs=obs,
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    result = m0_qc.flag_qc(adata, n_mads=5.0, min_cells_group=10, min_paired_donors=3)

    assert result.obs["flag_underpowered_celltype"].all()


# Catches count-only validation accepting misspelled conditions or cell-type names.
def test_T12_validate_assumptions_requires_exact_names():
    expected_cell_types = [f"type_{i}" for i in range(8)]
    for invalid_condition in ("control", "stim "):
        adata = _assumption_toy()
        adata.obs.iloc[0, adata.obs.columns.get_loc("condition")] = invalid_condition
        with pytest.raises(ValueError):
            m0_qc.validate_assumptions(adata, expected_cell_types=expected_cell_types)

    assert m0_qc.infer_conditions(["control", "IFNb"]) == ("control", "IFNb")
    with pytest.raises(ValueError):
        m0_qc.infer_conditions(["a", "b"])
    with pytest.raises(ValueError):
        m0_qc.infer_conditions(["ctrl", "control"])

    adata = _assumption_toy()
    adata.obs.loc[adata.obs["cell_type"] == "type_4", "cell_type"] = "CD4 T cell"
    with pytest.raises(ValueError):
        m0_qc.validate_assumptions(adata, expected_cell_types=expected_cell_types)


# Catches a changed or truncated GEO metadata file passing input validation unnoticed.
def test_T13_load_inputs_checks_geo_sha256(tmp_path):
    h5ad_path, geo_path, _ = _run_toy(tmp_path)

    with pytest.raises(ValueError):
        m0_qc.load_inputs(h5ad_path, geo_path, expected_geo_sha256="0" * 64)


# Catches permissive library-default floating tolerance hiding shifted t-SNE coordinates.
def test_T14_join_rejects_one_millionth_tsne_shift():
    obs = pd.DataFrame(
        {
            "label": ["ctrl"],
            "replicate": ["patient_0"],
            "cluster": [0],
            "tsne1": [1.0],
            "tsne2": [2.0],
        },
        index=["AAA-1"],
    )
    adata = ad.AnnData(
        np.ones((1, 1), dtype=int),
        obs=obs,
        var=pd.DataFrame(index=["G1"]),
    )
    geo = pd.DataFrame(
        {
            "stim": ["ctrl"],
            "ind": ["patient_0"],
            "cluster": [0],
            "cell": ["type_0"],
            "tsne1": [1.0 + 1e-6],
            "tsne2": [2.0],
            "multiplets": ["singlet"],
        },
        index=["AAA-1"],
    )

    with pytest.raises(ValueError):
        m0_qc.join_author_metadata(adata, geo)


# Catches the result H5AD and audit cell-flow report being written to the wrong directories.
def test_T15_run_m0_separates_output_and_audit_dirs(tmp_path):
    h5ad_path, geo_path, cell_types = _run_toy(tmp_path)
    out_dir = tmp_path / "out"
    audit_dir = tmp_path / "audit"

    m0_qc.run_m0(
        h5ad_path,
        geo_path,
        out_dir,
        expected_cell_types=sorted(set(cell_types)),
        audit_dir=audit_dir,
    )

    assert (out_dir / "m0_cells.h5ad").is_file()
    assert (audit_dir / "cell_flow.md").is_file()
    assert not (out_dir / "cell_flow.md").exists()
    assert not (audit_dir / "m0_cells.h5ad").exists()
