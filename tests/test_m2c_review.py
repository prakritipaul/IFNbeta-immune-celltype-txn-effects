import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.m2c_programs import label_programs


REPO_ROOT = Path(__file__).resolve().parents[1]


# A cell type with no testable sets remains part of the program-label universe.
def test_T10_missing_cell_type_is_untested_not_shared():
    cell_types = [f"cell_{i}" for i in range(1, 9)]
    gsea = pd.DataFrame([
        {"cell_type": cell_type, "collection": "H", "pathway": "PROGRAM", "NES": 2, "padj": 0.01}
        for cell_type in cell_types[:7]
    ])

    got = label_programs(
        gsea, cell_types=cell_types, padj_threshold=0.05
    ).set_index("pathway").loc["PROGRAM"]

    assert got["label"] == "partial"
    assert got[cell_types[7]] == "n/t"

    unknown = pd.concat([
        gsea,
        pd.DataFrame([{
            "cell_type": "cell_9", "collection": "H", "pathway": "PROGRAM",
            "NES": 2, "padj": 0.01,
        }]),
    ], ignore_index=True)
    with pytest.raises(ValueError):
        label_programs(unknown, cell_types=cell_types, padj_threshold=0.05)


# A cell type's BH correction pools pathways from every collection.
@pytest.mark.r
def test_T11_fgsea_bh_pools_collections_within_cell_type(tmp_path):
    genes = [f"G{i:04d}" for i in range(1000)]
    ranks_path = tmp_path / "ranks.csv"
    first_gmt = tmp_path / "mostly_small.gmt"
    second_gmt = tmp_path / "mostly_large.gmt"
    out_path = tmp_path / "gsea.csv"
    oracle_path = tmp_path / "oracle.csv"
    pd.DataFrame({
        "cell_type": "toy",
        "gene": genes,
        "stat": np.arange(500, -500, -1),
    }).to_csv(ranks_path, index=False)
    first_gmt.write_text("".join(
        f"SMALL_{i}\tna\t" + "\t".join(genes[i * 20:(i + 1) * 20]) + "\n"
        for i in range(10)
    ))
    second_gmt.write_text("".join(
        f"LARGE_{i}\tna\t" + "\t".join(np.array(genes)[i + 50 * np.arange(20)]) + "\n"
        for i in range(10)
    ))

    subprocess.run(
        ["/usr/local/bin/Rscript", str(REPO_ROOT / "R/m2c_fgsea.R"),
         str(ranks_path), str(first_gmt), str(second_gmt), str(out_path), "11011"],
        check=True, capture_output=True, text=True,
    )
    got = pd.read_csv(out_path)
    got[["collection", "pval"]].to_csv(oracle_path, index=False)
    expression = (
        "d<-read.csv(commandArgs(TRUE)[1]);"
        "d$pooled<-p.adjust(d$pval,method='BH');"
        "d$per_collection<-unsplit(lapply(split(d$pval,d$collection),p.adjust,method='BH'),d$collection);"
        "write.csv(d,commandArgs(TRUE)[2],row.names=FALSE)"
    )
    oracle_out = tmp_path / "oracle_adjusted.csv"
    subprocess.run(
        ["/usr/local/bin/Rscript", "-e", expression, str(oracle_path), str(oracle_out)],
        check=True, capture_output=True, text=True,
    )
    expected = pd.read_csv(oracle_out)

    assert np.allclose(got["padj"], expected["pooled"], rtol=1e-12, atol=1e-15)
    assert not np.allclose(expected["pooled"], expected["per_collection"], rtol=1e-12, atol=1e-15)
