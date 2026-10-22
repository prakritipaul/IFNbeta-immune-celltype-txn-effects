import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

import src.m3_signature as m3


# Catches float32 scoring that only occurs accidentally when the input matrix is already float32.
def test_T10_score_cells_uses_float32_working_matrix(monkeypatch):
    genes = ["A", "S", "G"] + [f"B{i:02d}" for i in range(100)]
    counts = np.arange(1, 2 * len(genes) + 1, dtype=np.float64).reshape(2, len(genes))
    adata = ad.AnnData(X=counts, var=pd.DataFrame(index=genes))
    observed_dtypes = []
    original = sc.tl.score_genes

    def recording_score_genes(scoring_adata, *args, **kwargs):
        observed_dtypes.append(scoring_adata.X.dtype)
        return original(scoring_adata, *args, **kwargs)

    monkeypatch.setattr(sc.tl, "score_genes", recording_score_genes)

    m3.score_cells(adata, {m3.IFNA: ["A", "S"], m3.IFNG: ["G", "S"]})

    assert observed_dtypes == [np.dtype("float32")] * 3
