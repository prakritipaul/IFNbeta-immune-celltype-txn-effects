import subprocess

import numpy as np
import pandas as pd
import pytest

import src.m1b_de as m1b


# Catches R results that preserve row counts while silently substituting a gene identity.
def test_T12_substituted_gene_name_is_rejected(tmp_path, monkeypatch):
    samples = pd.DataFrame(
        [
            (f"T|{donor}|{condition}", "T", donor, condition, True)
            for donor in ("d0", "d1")
            for condition in ("ctrl", "stim")
        ],
        columns=["sample", "cell_type", "donor", "condition", "kept"],
    )
    counts = pd.DataFrame(
        [[10, 12, 11, 13], [20, 22, 21, 23]],
        index=["gene_a", "gene_b"],
        columns=samples["sample"],
    )

    def fake_r_call(cmd, **kwargs):
        pd.DataFrame(
            {
                "cell_type": ["T", "T"],
                "gene": ["gene_a", "substituted_gene"],
                "baseMean": [11.5, 21.5],
                "log2FC": [0.1, 0.2],
                "lfcSE": [0.1, 0.1],
                "pvalue": [0.5, 0.4],
                "qvalue": [0.5, 0.5],
                "tested": [True, True],
                "cooks_flag": [False, False],
                "log2FC_shrunk": [0.05, 0.1],
            }
        ).to_csv(cmd[5], index=False)
        pd.DataFrame(
            {
                "cell_type": ["T"],
                "fit_type_used": ["parametric"],
                "n_samples": [4],
                "n_genes_tested": [2],
            }
        ).to_csv(cmd[6], index=False)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(m1b.subprocess, "run", fake_r_call)

    with pytest.raises(ValueError):
        m1b.run_deseq2(counts, samples, "paired", "substituted")


# Catches qvalue correction over all genes followed by masking genes with baseMean <= 1.
@pytest.mark.r
def test_T13_qvalue_excludes_untested_genes(tmp_path):
    rng = np.random.default_rng(13013)
    sample_ids = [f"T|d{i}|{condition}" for i in range(8) for condition in ("ctrl", "stim")]
    samples = pd.DataFrame(
        [(f"T|d{i}|{condition}", "T", f"d{i}", condition, True) for i in range(8) for condition in ("ctrl", "stim")],
        columns=["sample", "cell_type", "donor", "condition", "kept"],
    )
    tested = rng.negative_binomial(20, 20 / (20 + rng.uniform(30, 180, (120, 1))), (120, 16)) + 1

    def sparse_genes(n):
        values = np.zeros((n, 16), dtype=int)
        for row in range(n):
            values[row, rng.choice(16, size=6, replace=False)] = 1
        return values

    base = np.vstack([tested, sparse_genes(40)])
    expanded = np.vstack([base, sparse_genes(80)])

    def run(counts, stem):
        counts_path = tmp_path / f"{stem}_counts.csv"
        samples_path = tmp_path / f"{stem}_samples.csv"
        results_path = tmp_path / f"{stem}_results.csv"
        fits_path = tmp_path / f"{stem}_fits.csv"
        genes = [f"tested_{i}" for i in range(120)] + [f"untested_{i}" for i in range(len(counts) - 120)]
        pd.DataFrame(counts, index=genes, columns=sample_ids).rename_axis("gene").to_csv(counts_path)
        samples.to_csv(samples_path, index=False)
        subprocess.run(
            ["/usr/local/bin/Rscript", "R/m1b_deseq2.R", counts_path, samples_path, "paired", results_path, fits_path, "1"],
            check=True,
            capture_output=True,
            text=True,
        )
        return pd.read_csv(results_path).set_index("gene")

    def oracle(results, stem):
        p_path = tmp_path / f"{stem}_p.csv"
        q_path = tmp_path / f"{stem}_q.csv"
        tested_results = results.loc[results["tested"]].dropna(subset=["pvalue"])
        tested_results[["pvalue"]].to_csv(p_path, index_label="gene")
        expression = (
            "d<-read.csv(commandArgs(TRUE)[1]);"
            "d$qvalue<-qvalue::qvalue(d$pvalue,lambda=seq(0.05,0.95,0.05),"
            "pi0.method='smoother',smooth.df=3)$qvalues;"
            "write.csv(d,commandArgs(TRUE)[2],row.names=FALSE)"
        )
        subprocess.run(
            ["/usr/local/bin/Rscript", "-e", expression, p_path, q_path],
            check=True,
            capture_output=True,
            text=True,
        )
        return pd.read_csv(q_path).set_index("gene")["qvalue"]

    base_results = run(base, "base")
    expanded_results = run(expanded, "expanded")
    base_untested = base_results.loc[base_results.index.str.startswith("untested_")]
    expanded_untested = expanded_results.loc[expanded_results.index.str.startswith("untested_")]

    assert len(base_untested) == 40
    assert len(expanded_untested) == 120
    for untested in (base_untested, expanded_untested):
        assert (untested["baseMean"] <= 1).all()
        assert untested["pvalue"].notna().all()
    for stem, results in (("base", base_results), ("expanded", expanded_results)):
        expected = oracle(results, stem)
        observed = results.loc[expected.index, "qvalue"]
        pd.testing.assert_series_equal(observed, expected, check_names=False, rtol=1e-12, atol=1e-12)

    common = [f"tested_{i}" for i in range(120)]
    pd.testing.assert_series_equal(
        base_results.loc[common, "qvalue"],
        expanded_results.loc[common, "qvalue"],
        check_names=False,
        rtol=1e-12,
        atol=1e-12,
    )
