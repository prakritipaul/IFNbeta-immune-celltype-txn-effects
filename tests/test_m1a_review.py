from pathlib import Path
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from src.m1a_pseudobulk import pseudobulk


REPO_ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = "/usr/local/bin/Rscript"


# Catches pairing validation that counts two rows per donor without requiring one row per condition.
@pytest.mark.r
def test_T8_r_scripts_reject_duplicate_condition_rows(tmp_path):
    donors = [f"d{i}" for i in range(1, 9) for _ in range(2)]
    conditions = ["ctrl", "ctrl", "stim", "stim"] + ["ctrl", "stim"] * 6
    sample_ids = [f"s{i}" for i in range(1, 17)]
    samples = pd.DataFrame(
        {
            "sample": sample_ids,
            "cell_type": ["T"] * 16,
            "donor": donors,
            "condition": conditions,
            "n_cells": [41, 57, 46, 68, 52, 73, 61, 49, 77, 55, 64, 83, 58, 71, 89, 66],
            "total_counts": [910, 1370, 1120, 1590, 1280, 1810, 1460, 1190,
                             2030, 1340, 1680, 2210, 1510, 1940, 2390, 1760],
            "kept": [True] * 16,
        }
    )
    samples_path = tmp_path / "samples.csv"
    samples.to_csv(samples_path, index=False)

    counts_path = tmp_path / "counts.csv"
    rng = np.random.default_rng(81)
    counts = rng.negative_binomial(
        n=np.linspace(5, 40, 120)[:, None],
        p=np.linspace(0.15, 0.65, 120)[:, None],
        size=(120, 16),
    )
    counts[:, np.array(conditions) == "stim"] += np.arange(120)[:, None] % 11
    pd.DataFrame(counts, index=[f"G{i}" for i in range(120)], columns=sample_ids).rename_axis("gene").to_csv(counts_path)
    varpart_out = tmp_path / "varpart.csv"
    varpart = subprocess.run(
        [RSCRIPT, str(REPO_ROOT / "R/m1a_varpart.R"), str(counts_path), str(samples_path), str(varpart_out), "1"],
        capture_output=True,
        text=True,
    )

    cancor_out = tmp_path / "cancor.csv"
    cancor = subprocess.run(
        [RSCRIPT, str(REPO_ROOT / "R/m1a_cancor.R"), str(samples_path), str(cancor_out)],
        capture_output=True,
        text=True,
    )

    assert varpart.returncode != 0
    assert cancor.returncode != 0
    assert not varpart_out.exists()
    assert not cancor_out.exists()
    assert any(word in varpart.stderr.lower() for word in ("condition", "pair")), varpart.stderr
    assert any(word in cancor.stderr.lower() for word in ("condition", "pair")), cancor.stderr


# Catches ambiguous pipe-delimited sample IDs caused by pipes inside donor or cell-type labels.
@pytest.mark.parametrize(("donor", "cell_type"), [("d|1", "T"), ("d1", "T|cell")])
def test_T9_pseudobulk_rejects_pipe_in_sample_id_labels(donor, cell_type):
    obs = pd.DataFrame(
        {
            "donor": [donor, donor],
            "condition": ["ctrl", "stim"],
            "cell_type": [cell_type, cell_type],
        },
        index=["cell_1", "cell_2"],
    )
    adata = ad.AnnData(np.ones((2, 1), dtype=int), obs=obs, var=pd.DataFrame(index=["G1"]))
    adata.uns["control"] = "ctrl"
    adata.uns["treatment"] = "stim"

    with pytest.raises(ValueError):
        pseudobulk(adata, min_cells=1)
