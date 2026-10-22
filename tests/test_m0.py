from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from src import m0_qc


REPO_ROOT = Path(__file__).resolve().parents[1]
H5AD_PATH = REPO_ROOT / "data/kang_2018.h5ad"
GEO_PATH = REPO_ROOT / "data/geo/GSE96583_batch2.total.tsne.df.tsv.gz"
EXPECTED_SHA256 = "e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830"
CELL_TYPES = [
    "B cells",
    "CD14+ Monocytes",
    "CD4 T cells",
    "CD8 T cells",
    "Dendritic cells",
    "FCGR3A+ Monocytes",
    "Megakaryocytes",
    "NK cells",
]


def _join_toy():
    obs = pd.DataFrame(
        {
            "label": ["ctrl", "stim", "stim"],
            "replicate": ["patient_1"] * 3,
            "cluster": [1, 2, 3],
            "tsne1": [10.0, 20.0, 30.0],
            "tsne2": [11.0, 21.0, 31.0],
        },
        index=["AAA-1", "AAA-2", "BBB-2"],
    )
    adata = ad.AnnData(np.ones((3, 1), dtype=int), obs=obs, var=pd.DataFrame(index=["G1"]))
    geo = pd.DataFrame(
        {
            "stim": ["ctrl", "stim", "stim"],
            "ind": ["patient_1"] * 3,
            "cluster": [1, 2, 3],
            "cell": ["label_A_ctrl", "label_A_stim", "label_B_stim"],
            "tsne1": [10.0, 20.0, 30.0],
            "tsne2": [11.0, 21.0, 31.0],
            "multiplets": ["singlet"] * 3,
        },
        index=["AAA-1", "AAA-11", "BBB-1"],
    )
    return adata, geo


def _valid_run_inputs(n_cell_types=8, missing_condition=False, count_problem=None):
    n_cells = 16
    barcodes = []
    conditions = []
    donors = []
    for donor in range(8):
        for condition in ("ctrl", "stim"):
            barcodes.append(f"BC{donor:02d}{condition[0].upper()}-{'1' if condition == 'ctrl' else '2'}")
            conditions.append(condition)
            donors.append(f"patient_{donor}")
    if missing_condition:
        conditions[1] = "ctrl"
    clusters = list(range(n_cells))
    obs = pd.DataFrame(
        {
            "label": conditions,
            "replicate": donors,
            "cluster": clusters,
            "tsne1": np.arange(n_cells, dtype=float),
            "tsne2": np.arange(n_cells, dtype=float) + 0.5,
        },
        index=barcodes,
    )
    X = np.ones((n_cells, 200), dtype=float)
    if count_problem == "negative":
        X[0, 0] = -1
    elif count_problem == "fractional":
        X[0, 0] = 0.5
    var_names = [f"G{i}" for i in range(200)]
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=var_names))
    labels = [f"type_{i % n_cell_types}" for i in range(n_cells)]
    geo = pd.DataFrame(
        {
            "stim": conditions,
            "ind": donors,
            "cluster": clusters,
            "cell": labels,
            "tsne1": obs["tsne1"].to_numpy(),
            "tsne2": obs["tsne2"].to_numpy(),
            "multiplets": ["singlet"] * n_cells,
        },
        index=[name.rsplit("-", 1)[0] + "-1" for name in barcodes],
    )
    return adata, geo


def _write_inputs(tmp_path, adata, geo):
    h5ad_path = tmp_path / "input.h5ad"
    geo_path = tmp_path / "geo.tsv.gz"
    adata.write_h5ad(h5ad_path)
    geo.to_csv(geo_path, sep="\t", compression="gzip")
    return h5ad_path, geo_path


# Catches raw-barcode joins and cross-condition author-label swaps.
def test_T1_join_barcode_collision():
    adata, geo = _join_toy()
    result = m0_qc.join_author_metadata(adata, geo)
    assert list(result.obs_names) == ["AAA-1", "AAA-2", "BBB-2"]
    assert result.obs["cell_type"].tolist() == [
        "label_A_ctrl",
        "label_A_stim",
        "label_B_stim",
    ]


# Catches silent unmatched cells and many-to-many metadata joins.
@pytest.mark.parametrize("failure", ["unmatched", "duplicate"])
def test_T2_join_fails_loudly(failure):
    adata, geo = _join_toy()
    if failure == "unmatched":
        geo = geo.drop(index="BBB-1")
    else:
        geo = pd.concat([geo, geo.loc[["AAA-11"]].rename(index={"AAA-11": "AAA-111"})])
    with pytest.raises(ValueError):
        m0_qc.join_author_metadata(adata, geo)


# Catches strict boundary comparisons and filtering genes before cells.
def test_T3_filter_boundaries_and_order():
    X = np.zeros((4, 201), dtype=int)
    X[0, :198] = 1
    X[1:, :200] = 1
    X[0, 200] = 1
    X[1, 200] = 1
    X[2, 200] = 1
    adata = ad.AnnData(X, obs=pd.DataFrame(index=["c199", "c200a", "c200b", "c200c"]),
                       var=pd.DataFrame(index=[f"G{i}" for i in range(200)] + ["THIRD_CELL_REMOVED"]))
    result = m0_qc.basic_filter(adata, min_genes=200, min_cells=3)
    assert result.obs_names.tolist() == ["c200a", "c200b", "c200c"]
    assert result.var_names.tolist() == [f"G{i}" for i in range(200)]


# Catches gene prevalence being calculated separately within each condition.
def test_T4_gene_filter_across_conditions():
    adata = ad.AnnData(
        np.array([[1], [1], [1]], dtype=int),
        obs=pd.DataFrame({"condition": ["ctrl", "ctrl", "stim"]}, index=["a", "b", "c"]),
        var=pd.DataFrame(index=["ACROSS_CONDITIONS"]),
    )
    result = m0_qc.basic_filter(adata, min_genes=1, min_cells=3)
    assert result.var_names.tolist() == ["ACROSS_CONDITIONS"]


# Catches wrong QC axes, denominators, and ribosomal-prefix matching.
def test_T5_qc_metric_arithmetic():
    adata = ad.AnnData(
        np.array([[1, 2, 0], [0, 1, 3]], dtype=int),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["RPL1", "RPS2", "ACTB"]),
    )
    result = m0_qc.compute_qc_metrics(adata)
    np.testing.assert_array_equal(result.obs["n_counts"], [3, 4])
    np.testing.assert_array_equal(result.obs["n_genes"], [2, 2])
    np.testing.assert_allclose(result.obs["pct_ribo"], [100.0, 25.0])
    np.testing.assert_array_equal(result.var["n_cells"], [1, 2, 1])


# Catches inclusive MAD cutoffs, pooled groups, scaled MAD, dropping flags, and off-by-one power rules.
def test_T6_flag_boundaries_and_no_drops():
    boundary = np.array([-1, -1, 0, 0, 1, 1, 5.0])
    assert not m0_qc.mad_outlier(boundary, n_mads=5.0)[-1]
    above = boundary.copy()
    above[-1] = 5.1
    assert m0_qc.mad_outlier(above, n_mads=5.0)[-1]

    grouped_obs = pd.DataFrame(
        {
            "donor": ["d1"] * 7 + ["d2"] * 7,
            "condition": ["ctrl"] * 14,
            "cell_type": ["boundary"] * 14,
            "n_counts": np.expm1(np.r_[above + 2, above + 102]),
            "n_genes": [2] * 14,
        },
        index=[f"g{i}" for i in range(14)],
    )
    grouped = ad.AnnData(np.ones((14, 2), dtype=int), obs=grouped_obs,
                         var=pd.DataFrame(index=["G1", "G2"]))
    grouped_result = m0_qc.flag_qc(
        grouped, n_mads=5.0, min_cells_group=1, min_paired_donors=1
    )
    assert grouped_result.obs.loc["g6", "flag_outlier"]

    rows = []
    for donor in ("d1", "d2", "d3"):
        for condition in ("ctrl", "stim"):
            for cell_type, size in (("two_paired", 9), ("three_paired", 10)):
                if cell_type == "two_paired" and donor == "d3":
                    continue
                rows.extend((donor, condition, cell_type) for _ in range(size))
    obs = pd.DataFrame(rows, columns=["donor", "condition", "cell_type"])
    obs.index = [f"c{i}" for i in range(len(obs))]
    obs["n_counts"] = 10
    obs["n_genes"] = 2
    adata = ad.AnnData(np.ones((len(obs), 2), dtype=int), obs=obs,
                       var=pd.DataFrame(index=["G1", "G2"]))
    result = m0_qc.flag_qc(adata, n_mads=5.0, min_cells_group=10, min_paired_donors=3)
    assert result.n_obs == adata.n_obs
    assert result.obs.loc[result.obs.cell_type == "two_paired", "flag_small_group"].all()
    assert not result.obs.loc[result.obs.cell_type == "three_paired", "flag_small_group"].any()
    assert result.obs.loc[result.obs.cell_type == "two_paired", "flag_underpowered_celltype"].all()
    assert not result.obs.loc[result.obs.cell_type == "three_paired", "flag_underpowered_celltype"].any()


# Catches invalid inputs being worked around or outputs being partially written.
@pytest.mark.parametrize(
    "failure",
    ["negative", "fractional", "missing_condition", "duplicate_genes", "seven_types", "nine_types"],
)
def test_T7_invalid_input_writes_nothing(tmp_path, failure):
    kwargs = {}
    if failure in {"negative", "fractional"}:
        kwargs["count_problem"] = failure
    elif failure == "missing_condition":
        kwargs["missing_condition"] = True
    elif failure == "seven_types":
        kwargs["n_cell_types"] = 7
    elif failure == "nine_types":
        kwargs["n_cell_types"] = 9
    adata, geo = _valid_run_inputs(**kwargs)
    if failure == "duplicate_genes":
        names = adata.var_names.tolist()
        names[-1] = names[0]
        adata.var_names = names
    h5ad_path, geo_path = _write_inputs(tmp_path, adata, geo)
    out_dir = tmp_path / "out"
    with pytest.raises(ValueError):
        m0_qc.run_m0(h5ad_path, geo_path, out_dir)
    assert not out_dir.exists() or list(out_dir.iterdir()) == []


# Catches silent cell losses between pipeline stages and missing join-match accounting.
def test_T8_cell_accounting_chain(tmp_path):
    adata, geo = _valid_run_inputs()
    h5ad_path, geo_path = _write_inputs(tmp_path, adata, geo)
    out_dir = tmp_path / "out"
    flow = m0_qc.run_m0(h5ad_path, geo_path, out_dir)
    for current, following in zip(flow, flow[1:]):
        assert current["n_cells_out"] == following["n_cells_in"]
    written = ad.read_h5ad(out_dir / "m0_cells.h5ad")
    assert flow[-1]["n_cells_out"] == written.n_obs
    join = next(record for record in flow if record["step"] == "join")
    assert join["n_matched"] == adata.n_obs


# Catches wrong real-data files, label drift, failed joins, and output-schema gaps.
@pytest.mark.realdata
def test_T9_real_data_contract(tmp_path):
    out_dir = tmp_path / "out"
    flow = m0_qc.run_m0(H5AD_PATH, GEO_PATH, out_dir, expected_sha256=EXPECTED_SHA256)
    result = ad.read_h5ad(out_dir / "m0_cells.h5ad")
    join = next(record for record in flow if record["step"] == "join")
    assert join["n_matched"] == 24673
    assert result.obs["donor"].nunique() == 8
    assert sorted(result.obs["cell_type"].unique()) == CELL_TYPES
    assert "pct_mito" not in result.obs
    expected_obs = {
        "barcode", "condition", "donor", "cell_type", "cluster", "tsne1", "tsne2",
        "n_counts", "n_genes", "pct_ribo", "flag_outlier", "flag_small_group",
        "flag_underpowered_celltype",
    }
    assert expected_obs <= set(result.obs.columns)
    assert "n_cells" in result.var.columns


# Catches colliding toy barcodes being merged, duplicated, swapped, or relabelled.
def test_T10_toy_cells_have_one_original_condition():
    adata, geo = _join_toy()
    original = adata.obs["label"].copy()
    result = m0_qc.join_author_metadata(adata, geo)
    assert result.n_obs == adata.n_obs
    assert result.obs_names.is_unique
    assert result.obs["condition"].tolist() == original.tolist()
    aaa = result.obs.loc[["AAA-1", "AAA-2"]]
    assert aaa["condition"].tolist() == ["ctrl", "stim"]


# Catches the 234 real collision sequences being merged, duplicated, or assigned the wrong condition.
@pytest.mark.realdata
def test_T10_real_collision_cells_remain_distinct():
    adata, geo = m0_qc.load_inputs(H5AD_PATH, GEO_PATH, expected_sha256=EXPECTED_SHA256)
    original = adata.obs["label"].copy()
    result = m0_qc.join_author_metadata(adata, geo)
    assert result.obs["condition"].tolist() == original.tolist()
    sequences = result.obs_names.str.replace(r"-\d+$", "", regex=True)
    conditions = pd.Series(result.obs["condition"].to_numpy(), index=sequences)
    collision_sequences = conditions.groupby(level=0).nunique().loc[lambda x: x == 2].index
    collision_cells = conditions.loc[conditions.index.isin(collision_sequences)]
    assert len(collision_sequences) == 234
    assert len(collision_cells) == 468
    assert (collision_cells == "ctrl").sum() == 234
    assert (collision_cells == "stim").sum() == 234
