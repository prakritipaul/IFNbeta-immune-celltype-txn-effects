import pandas as pd
import pytest

import src.m2a_calls as m2a
from tests.test_m2a import _gene_rows, _lodo_rows


# Catches traceability checks that reuse the label function they are checking.
def test_T11_traceability_independently_rederives_labels(monkeypatch):
    calls = pd.DataFrame(_gene_rows("g", {"B cells": "up"}))
    original = m2a.label_genes

    def wrong_label_genes(input_calls):
        table = original(input_calls)
        table["shared_response"] = "shared"
        return table

    monkeypatch.setattr(m2a, "label_genes", wrong_label_genes)
    with pytest.raises(ValueError):
        m2a.check_traceability(m2a.label_genes(calls), calls)


# Catches LODO validation that checks run count but not exact run identity.
def test_T12_lodo_rejects_missing_and_unexpected_run_names():
    calls = pd.DataFrame(_gene_rows("g", {"B cells": "up"}))
    paired_donors = {
        cell_type: [f"d{i}" for i in range(6 if cell_type == "Dendritic cells" else 4 if cell_type == "Megakaryocytes" else 8)]
        for cell_type in m2a.WELL_POWERED + m2a.LOW_POWER
    }
    unexpected = _lodo_rows(calls)
    replace = (unexpected["cell_type"] == "B cells") & (unexpected["run"] == "lodo_d0")
    unexpected.loc[replace, "run"] = "lodo_NOT_A_DONOR"

    with pytest.raises(ValueError, match=r"(?s)(?=.*lodo_d0)(?=.*lodo_NOT_A_DONOR)"):
        m2a.lodo_stability(calls, unexpected, paired_donors)

    substituted = _lodo_rows(calls)
    replace = (substituted["cell_type"] == "Dendritic cells") & (substituted["run"] == "lodo_d0")
    substituted.loc[replace, "run"] = "lodo_d6"

    with pytest.raises(ValueError):
        m2a.lodo_stability(calls, substituted, paired_donors)
