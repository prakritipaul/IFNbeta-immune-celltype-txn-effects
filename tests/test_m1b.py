from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import src.m1b_de as m1b
from src.m1b_de import balanced_shuffles, classify, run_deseq2, run_m1b


REPO_ROOT = Path(__file__).resolve().parents[1]
REAL_INPUT = REPO_ROOT / "data/processed"
GMT = REPO_ROOT / "data/msigdb/h.all.v2026.1.Hs.symbols.gmt"
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


def _samples(cell_types=("T",), donors=tuple(f"d{i}" for i in range(8))):
    rows = []
    for cell_type in cell_types:
        for donor in donors:
            for condition in ("ctrl", "stim"):
                sample = f"{cell_type}|{donor}|{condition}"
                rows.append((sample, cell_type, donor, condition, True))
    return pd.DataFrame(rows, columns=["sample", "cell_type", "donor", "condition", "kept"])


def _nb_counts(samples, seed=20261001, n_null=200, effects=None):
    rng = np.random.default_rng(seed)
    effects = effects or {}
    genes = list(effects) + [f"null_{i}" for i in range(n_null)]
    donor_levels = {d: rng.uniform(0.65, 1.55) for d in samples["donor"].unique()}
    values = np.empty((len(genes), len(samples)), dtype=int)
    for gene_i, gene in enumerate(genes):
        mean = rng.uniform(80, 220)
        for sample_i, row in enumerate(samples.itertuples(index=False)):
            fold = effects.get(gene, 1.0) if row.condition == "stim" else 1.0
            mu = mean * donor_levels[row.donor] * fold
            size = 40.0
            values[gene_i, sample_i] = rng.negative_binomial(size, size / (size + mu))
    return pd.DataFrame(values, index=genes, columns=samples["sample"])


def _write_m1a_inputs(directory, counts, samples, include_clean=True):
    directory.mkdir(parents=True, exist_ok=True)
    counts.to_csv(directory / "m1a_counts.csv", index_label="gene")
    samples.to_csv(directory / "m1a_samples.csv", index=False)
    if include_clean:
        counts.to_csv(directory / "m1a_counts_no_outliers.csv", index_label="gene")
        samples.to_csv(directory / "m1a_samples_no_outliers.csv", index=False)


@pytest.fixture(scope="module")
def real_m1b_output(tmp_path_factory):
    output = tmp_path_factory.mktemp("m1b_real")
    run_m1b(REAL_INPUT, output)
    return output


# Catches a reversed contrast, broken Python/R hand-off, wrong design, or mass false positives.
@pytest.mark.r
def test_T1_simulated_paired_oracle():
    samples = _samples()
    counts = _nb_counts(samples, effects={"up_1": 8, "up_2": 6, "down_1": 1 / 8, "down_2": 1 / 6})

    results, _ = run_deseq2(counts, samples, "paired", "oracle")
    calls = results.set_index("gene")["de"]

    assert (calls[["up_1", "up_2"]] == "up").all()
    assert (calls[["down_1", "down_2"]] == "down").all()
    assert (calls.filter(like="null_") == "none").mean() >= 0.95


# Catches omission of the donor term or samples attached to the wrong donor/condition.
@pytest.mark.r
def test_T2_pairing_removes_donor_baselines():
    samples = _samples()
    rng = np.random.default_rng(314159)
    genes = [f"gene_{i}" for i in range(120)]
    pair_counts = rng.negative_binomial(30, 30 / (30 + rng.uniform(30, 400, (len(genes), 8)))) + 1
    pair_means = np.repeat(pair_counts, 2, axis=1)
    noise_size = 50.0
    counts = rng.negative_binomial(noise_size, noise_size / (noise_size + pair_means)) + 1
    counts = pd.DataFrame(counts, index=genes, columns=samples["sample"])

    results, _ = run_deseq2(counts, samples, "paired", "paired_null")

    np.testing.assert_allclose(results["log2FC"], 0, atol=0.3)
    assert (results["de"] == "none").all()


# Catches an inclusive q boundary, an exclusive LFC boundary, or failure to use absolute log2 fold change.
def test_T3_call_boundaries():
    rows = pd.DataFrame(
        {
            "qvalue": [0.0499, 0.05, 0.01, 0.01, 0.01, 0.01],
            "log2FC": [1.01, 2.0, 1.0, -1.0, 0.9999, -0.9999],
        }
    )

    observed = classify(rows)

    assert observed["significant"].tolist() == [True, False, True, True, True, True]
    assert observed["de"].tolist() == ["up", "none", "up", "down", "none", "none"]


# Catches pre-filtering before DESeq2 or using baseMean >= 1 instead of > 1.
@pytest.mark.r
def test_T4_basemean_filter_after_fit_is_strict():
    samples = _samples(donors=tuple(f"d{i}" for i in range(4)))
    rng = np.random.default_rng(44)
    filler = rng.negative_binomial(20, 0.2, (120, len(samples))) + 1
    data = np.vstack(([0] * 8, [1] * 8, [2] * 8, filler))
    counts = pd.DataFrame(data, index=["below", "exact", "above"] + [f"filler_{i}" for i in range(120)], columns=samples["sample"])

    results, fits = run_deseq2(counts, samples, "paired", "basemean")
    targets = results.set_index("gene").loc[["below", "exact", "above"]]

    assert targets["tested"].tolist() == [False, False, True]
    assert len(results) == len(counts)
    assert int(fits.loc[0, "n_genes_tested"]) == int((results["baseMean"] > 1).sum())


# Catches qvalue correction pooled across cell types.
@pytest.mark.r
def test_T5_qvalue_is_per_cell_type():
    samples = _samples(("A", "B"))
    counts_1 = _nb_counts(samples, seed=51, effects={"signal": 5})
    counts_2 = counts_1.copy()
    b_samples = samples.loc[samples["cell_type"] == "B", "sample"]
    rng = np.random.default_rng(52)
    counts_2.loc[:, b_samples] = rng.negative_binomial(15, 0.1, (len(counts_2), len(b_samples)))

    result_1, _ = run_deseq2(counts_1, samples, "paired", "first")
    result_2, _ = run_deseq2(counts_2, samples, "paired", "second")
    a1 = result_1[result_1["cell_type"] == "A"].set_index("gene")["qvalue"].sort_index()
    a2 = result_2[result_2["cell_type"] == "A"].set_index("gene")["qvalue"].sort_index()

    pd.testing.assert_series_equal(a1, a2, check_names=False, rtol=1e-12, atol=1e-12)


# Catches positional joins, unpaired/duplicate samples, coerced counts, and partial writes.
@pytest.mark.parametrize("bad_case", ["ids", "missing_stim", "two_ctrl", "negative", "fractional", "nonnumeric"])
def test_T6_bad_input_stops_before_writing(tmp_path, monkeypatch, bad_case):
    samples = _samples(donors=("d0", "d1"))
    counts = _nb_counts(samples, n_null=20)
    if bad_case == "ids":
        counts = counts.rename(columns={counts.columns[0]: "unknown_sample"})
    elif bad_case == "missing_stim":
        samples = samples.drop(samples[(samples["donor"] == "d0") & (samples["condition"] == "stim")].index)
        counts = counts[samples["sample"]]
    elif bad_case == "two_ctrl":
        samples.loc[(samples["donor"] == "d0") & (samples["condition"] == "stim"), "condition"] = "ctrl"
    elif bad_case == "negative":
        counts.iloc[0, 0] = -1
    elif bad_case == "fractional":
        counts = counts.astype(float)
        counts.iloc[0, 0] = 0.5
    else:
        counts = counts.astype(object)
        counts.iloc[0, 0] = "not-a-count"
    input_dir, output_dir = tmp_path / "input", tmp_path / "output"
    _write_m1a_inputs(input_dir, counts, samples)
    called = False

    def forbidden(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("R must not run for invalid input")

    monkeypatch.setattr(m1b.subprocess, "run", forbidden)
    with pytest.raises(ValueError):
        run_m1b(input_dir, output_dir, runs=("main",))
    assert not called
    assert not output_dir.exists() or not list(output_dir.iterdir())


# Catches lost cell types, incomplete LODO fits, or duplicate/incomplete balanced shuffles.
@pytest.mark.r
@pytest.mark.realdata
def test_T7_real_data_runs_complete(real_m1b_output):
    fits = pd.read_csv(real_m1b_output / "m1b_fits.csv")
    runs = pd.read_csv(real_m1b_output / "m1b_runs.csv")
    expected_main = {cell_type: 16 for cell_type in CELL_TYPES}
    expected_main["Dendritic cells"] = 12
    expected_main["Megakaryocytes"] = 8
    main = fits[fits["run"] == "main"].set_index("cell_type")["n_samples"].to_dict()
    assert main == expected_main

    lodo = fits[fits["run"].str.startswith("lodo_")].groupby("cell_type")["run"].nunique().to_dict()
    assert lodo == {cell_type: (6 if cell_type == "Dendritic cells" else 4 if cell_type == "Megakaryocytes" else 8) for cell_type in CELL_TYPES}
    shuffled = runs[runs["run"].str.startswith("shuffle_")].groupby("cell_type")["run"].nunique().to_dict()
    expected_shuffles = {cell_type: 35 for cell_type in CELL_TYPES}
    expected_shuffles["Dendritic cells"] = 10
    expected_shuffles["Megakaryocytes"] = 3
    assert shuffled == expected_shuffles

    for n_donors, expected_count in ((8, 35), (6, 10), (4, 3)):
        donors = [f"d{i}" for i in range(n_donors)]
        swaps = balanced_shuffles(donors)
        assert len(swaps) == expected_count
        assert all(len(swap) == n_donors // 2 for swap in swaps)
        assert not any(set(left) == set(donors) - set(right) for i, left in enumerate(swaps) for right in swaps[i + 1:])


# Catches joining R results on gene alone or trusting returned row order.
def test_T8_row_keys_survive_shuffled_r_output(tmp_path, monkeypatch):
    samples = _samples(("A", "B"), donors=("d0", "d1"))
    counts = _nb_counts(samples, n_null=2)
    input_dir, output_dir = tmp_path / "input", tmp_path / "output"
    _write_m1a_inputs(input_dir, counts, samples)
    expected = pd.DataFrame(
        [
            {"run": "main", "cell_type": cell_type, "gene": gene, "baseMean": value, "log2FC": value / 10, "lfcSE": 0.2,
             "ci_low": value / 10 - 0.392, "ci_high": value / 10 + 0.392, "pvalue": 0.2, "qvalue": 0.3,
             "tested": True, "cooks_flag": False, "log2FC_shrunk": value / 20, "significant": False, "de": "none"}
            for value, (cell_type, gene) in enumerate((("A", "null_0"), ("A", "null_1"), ("B", "null_0"), ("B", "null_1")), 11)
        ]
    )

    def fake_run(*args, **kwargs):
        fits = pd.DataFrame({"run": ["main", "main"], "cell_type": ["B", "A"], "fit_type_used": ["parametric", "parametric"], "n_samples": [4, 4], "n_genes_tested": [2, 2]})
        return expected.sample(frac=1, random_state=8).reset_index(drop=True), fits

    monkeypatch.setattr(m1b, "run_deseq2", fake_run)
    run_m1b(input_dir, output_dir, runs=("main",))
    observed = pd.read_csv(output_dir / "m1b_de.csv")

    assert not observed.duplicated(["run", "cell_type", "gene"]).any()
    pd.testing.assert_frame_equal(
        observed.sort_values(["cell_type", "gene"]).reset_index(drop=True),
        expected.sort_values(["cell_type", "gene"]).reset_index(drop=True),
        check_dtype=False,
    )


# Catches disabled Cook's handling or silently dropped, filled, or unexplained NA p-values.
@pytest.mark.r
def test_T9_cooks_outliers_and_unexplained_na():
    s8 = _samples(("eight",), tuple(f"d{i}" for i in range(8)))
    s4 = _samples(("four",), tuple(f"d{i}" for i in range(4)))
    samples = pd.concat([s8, s4], ignore_index=True)
    counts = _nb_counts(samples, seed=99, n_null=150)
    counts.loc["outlier"] = 100
    counts.loc["outlier", "eight|d0|stim"] = 10_000_000
    counts.loc["outlier", "four|d0|stim"] = 10_000_000

    results, _ = run_deseq2(counts, samples, "unpaired", "cooks")
    outliers = results[results["gene"] == "outlier"].set_index("cell_type")

    assert outliers["cooks_flag"].all()
    assert pd.notna(outliers.loc["eight", "pvalue"])
    assert pd.isna(outliers.loc["four", "pvalue"])
    assert outliers.loc["four", "de"] == "none"

    paired, _ = run_deseq2(counts, samples, "paired", "paired_cooks")
    assert not paired["cooks_flag"].any()

    bad = pd.DataFrame({"baseMean": [10.0], "pvalue": [np.nan], "qvalue": [np.nan], "log2FC": [2.0], "cooks_flag": [False], "tested": [True]})
    with pytest.raises(ValueError):
        classify(bad)


# Catches wrong Wald CIs, missing fit logs, or calls made from shrunken rather than MLE LFC.
@pytest.mark.r
def test_T10_reported_statistics_and_unshrunk_calls():
    samples = _samples()
    counts = _nb_counts(samples, seed=1010, effects={"signal": 5})
    results, fits = run_deseq2(counts, samples, "paired", "statistics")
    finite = results[["log2FC", "lfcSE", "ci_low", "ci_high"]].notna().all(axis=1)

    np.testing.assert_allclose(results.loc[finite, "ci_low"], results.loc[finite, "log2FC"] - 1.96 * results.loc[finite, "lfcSE"])
    np.testing.assert_allclose(results.loc[finite, "ci_high"], results.loc[finite, "log2FC"] + 1.96 * results.loc[finite, "lfcSE"])
    assert fits[["run", "cell_type", "fit_type_used"]].notna().all().all()
    assert "log2FC_shrunk" in results
    boundary = pd.DataFrame({"qvalue": [0.01], "log2FC": [1.2], "log2FC_shrunk": [0.8]})
    assert classify(boundary).loc[0, "de"] == "up"


# Catches a reversed contrast, wrong pairing, or a pipeline with no canonical IFN response.
@pytest.mark.r
@pytest.mark.realdata
def test_T11_real_data_ifn_sanity_tripwire(real_m1b_output):
    results = pd.read_csv(real_m1b_output / "m1b_de.csv")
    mappings = {"TENT5A": "FAM46A", "WARS1": "WARS", "MARCHF1": "MARCH1", "RIGI": "DDX58", "TMT1B": "METTL7B", "PSMA2": "PSMA2.1"}
    hallmark = None
    for line in GMT.read_text().splitlines():
        fields = line.split("\t")
        if fields[0] == "HALLMARK_INTERFERON_ALPHA_RESPONSE":
            hallmark = {mappings.get(gene, gene) for gene in fields[2:]}
            break
    assert hallmark is not None and len(hallmark) == 97
    canonical = {"ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "RSAD2", "IRF7", "STAT1"}

    for cell_type in CELL_TYPES:
        calls = results[results["cell_type"] == cell_type].set_index("gene")["de"]
        assert (calls.reindex(sorted(hallmark)) == "up").sum() >= 10
        assert not (calls.reindex(sorted(canonical)) == "down").any()
