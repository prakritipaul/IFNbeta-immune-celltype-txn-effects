from pathlib import Path

import pandas as pd
import pytest

from src.m2a_calls import (
    LOW_POWER,
    WELL_POWERED,
    check_traceability,
    derive_calls,
    label_genes,
    lodo_stability,
    run_m2a,
)


CELL_TYPES = WELL_POWERED + LOW_POWER
REPO_ROOT = Path(__file__).resolve().parents[1]


def _gene_rows(gene, de_by_cell=None, ci_by_cell=None):
    de_by_cell = de_by_cell or {}
    ci_by_cell = ci_by_cell or {}
    rows = []
    for cell_type in CELL_TYPES:
        de = de_by_cell.get(cell_type, "none")
        log2fc = {"up": 2.0, "down": -2.0, "none": 0.0}[de]
        ci_low, ci_high = ci_by_cell.get(cell_type, (log2fc - 0.2, log2fc + 0.2))
        rows.append(
            {
                "cell_type": cell_type,
                "gene": gene,
                "log2FC": log2fc,
                "qvalue": 0.01 if de != "none" else 0.5,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "de": de,
                "tested": True,
            }
        )
    return rows


def _lodo_rows(calls, values=None):
    values = values or {}
    rows = []
    for row in calls.itertuples(index=False):
        n_runs = 6 if row.cell_type == "Dendritic cells" else 4 if row.cell_type == "Megakaryocytes" else 8
        default = 0.2 if row.de != "down" else -0.2
        for run_i in range(n_runs):
            rows.append(
                {
                    "run": f"lodo_d{run_i}",
                    "cell_type": row.cell_type,
                    "gene": row.gene,
                    "log2FC": values.get((row.cell_type, run_i), default),
                }
            )
    return pd.DataFrame(rows)


# Catches wrong shared/specific/mixed label logic, precedence, arrows, or DE counts.
def test_T1_label_oracle():
    rows = []
    rows += _gene_rows("shared", {cell_type: "up" for cell_type in CELL_TYPES})
    rows += _gene_rows(
        "shared_star",
        {cell_type: "up" for cell_type in WELL_POWERED},
        {cell_type: (0.2, 1.0) for cell_type in LOW_POWER},
    )
    rows += _gene_rows("partial", {"B cells": "up", "CD4 T cells": "up"})
    rows += _gene_rows("specific", {"B cells": "up"})
    rows += _gene_rows(
        "specific_star",
        {"B cells": "up"},
        {"CD4 T cells": (-1.0, 0.2)},
    )
    rows += _gene_rows("mixed", {"B cells": "up", "CD4 T cells": "down"})
    rows += _gene_rows(
        "mixed_precedence",
        {**{cell_type: "up" for cell_type in WELL_POWERED}, "Dendritic cells": "down"},
        {"Megakaryocytes": (0.2, 1.2)},
    )

    table = label_genes(pd.DataFrame(rows)).set_index("gene")

    assert table["shared_response"].to_dict() == {
        "shared": "shared",
        "shared_star": "shared*",
        "partial": "partial",
        "specific": "specific",
        "specific_star": "specific*",
        "mixed": "mixed",
        "mixed_precedence": "mixed",
    }
    assert table["n_cell_types"].to_dict() == {
        "shared": 8,
        "shared_star": 6,
        "partial": 2,
        "specific": 1,
        "specific_star": 1,
        "mixed": 2,
        "mixed_precedence": 7,
    }
    assert table.loc["mixed", ["B cells", "CD4 T cells", "NK cells"]].tolist() == ["↑", "↓", "–"]


# Catches wrong q/LFC inclusivity or trusting M1b's supplied DE column.
def test_T2_de_boundaries_rederived():
    cases = {
        "up_exact": (0.0499, 1.0, "down", "up"),
        "down_exact": (0.0499, -1.0, "up", "down"),
        "q_equal": (0.05, 2.0, "up", "none"),
        "lfc_inside": (0.01, 0.9999, "up", "none"),
        "both_fail": (0.05, -0.9999, "down", "none"),
    }
    rows = []
    for gene, (qvalue, log2fc, supplied_de, _) in cases.items():
        gene_rows = _gene_rows(gene)
        for row in gene_rows:
            row.update(qvalue=qvalue, log2FC=log2fc, de=supplied_de)
        rows += gene_rows

    calls = derive_calls(pd.DataFrame(rows)).set_index(["gene", "cell_type"])

    for gene, (_, _, _, expected) in cases.items():
        assert (calls.loc[gene, "de"] == expected).all()


# Catches strict CI boundaries, absolute-value direction errors, or wrong specific-star scope.
def test_T3_ci_rules_boundary_and_direction():
    rows = []
    rows += _gene_rows("up_boundary", {cell_type: "up" for cell_type in WELL_POWERED}, {cell_type: (0.1, 1.0) for cell_type in LOW_POWER})
    rows += _gene_rows("down_boundary", {cell_type: "down" for cell_type in WELL_POWERED}, {cell_type: (-1.0, -0.1) for cell_type in LOW_POWER})
    rows += _gene_rows("up_inside", {cell_type: "up" for cell_type in WELL_POWERED}, {cell_type: (0.1, 0.999) for cell_type in LOW_POWER})
    rows += _gene_rows("down_inside", {cell_type: "down" for cell_type in WELL_POWERED}, {cell_type: (-0.999, -0.1) for cell_type in LOW_POWER})
    rows += _gene_rows("opposite_only", {cell_type: "up" for cell_type in WELL_POWERED}, {cell_type: (-1.2, 0.5) for cell_type in LOW_POWER})
    rows += _gene_rows("specific_other", {"B cells": "up"}, {"CD4 T cells": (0.1, 1.0)})
    rows += _gene_rows("own_ci_only", {"B cells": "up"}, {"B cells": (-1.2, 2.2)})
    rows += _gene_rows("specific_untested_ci", {"B cells": "up"}, {"CD4 T cells": (0.1, 1.0)})
    rows += _gene_rows("specific_nan_ci", {"B cells": "up"}, {"CD4 T cells": (float("nan"), float("nan"))})
    rows += _gene_rows("shared_untested_ci", {cell_type: "up" for cell_type in WELL_POWERED}, {cell_type: (0.1, 1.0) for cell_type in LOW_POWER})
    rows += _gene_rows("shared_nan_ci", {cell_type: "up" for cell_type in WELL_POWERED}, {cell_type: (0.1, 1.0) for cell_type in LOW_POWER})
    for row in rows:
        row["tested"] = True
        if (row["gene"], row["cell_type"]) in {
            ("specific_untested_ci", "CD4 T cells"),
            ("shared_untested_ci", "Dendritic cells"),
        }:
            row["tested"] = False
        if row["gene"] == "shared_nan_ci" and row["cell_type"] == "Dendritic cells":
            row["ci_low"] = row["ci_high"] = float("nan")

    labels = label_genes(pd.DataFrame(rows)).set_index("gene")["shared_response"]

    assert labels.to_dict() == {
        "up_boundary": "shared*",
        "down_boundary": "shared*",
        "up_inside": "partial",
        "down_inside": "partial",
        "opposite_only": "partial",
        "specific_other": "specific*",
        "own_ci_only": "specific",
        "specific_untested_ci": "specific",
        "specific_nan_ci": "specific",
        "shared_untested_ci": "partial",
        "shared_nan_ci": "partial",
    }


# Catches requiring LODO significance or treating zero as retaining a call's sign.
def test_T4_lodo_sign():
    calls = pd.DataFrame(
        _gene_rows("g", {"B cells": "up", "CD14+ Monocytes": "up", "CD4 T cells": "up", "CD8 T cells": "down"})
    )
    values = {
        ("B cells", 7): 0.001,
        ("CD14+ Monocytes", 7): 0.0,
        ("CD4 T cells", 7): -0.001,
        ("CD8 T cells", 7): -0.001,
    }

    paired_donors = {
        cell_type: [f"d{i}" for i in range(6 if cell_type == "Dendritic cells" else 4 if cell_type == "Megakaryocytes" else 8)]
        for cell_type in CELL_TYPES
    }
    observed = lodo_stability(calls, _lodo_rows(calls, values), paired_donors).set_index("cell_type")

    assert bool(observed.loc["B cells", "lodo_stable"])
    assert not bool(observed.loc["CD14+ Monocytes", "lodo_stable"])
    assert not bool(observed.loc["CD4 T cells", "lodo_stable"])
    assert bool(observed.loc["CD8 T cells", "lodo_stable"])
    assert pd.isna(observed.loc["NK cells", "lodo_stable"])


# Catches stability computed from missing, duplicated, or wrong-count LODO refits.
def test_T5_incomplete_lodo_stops():
    calls = pd.DataFrame(_gene_rows("g", {"B cells": "up"}))
    valid = _lodo_rows(calls)
    paired_donors = {
        cell_type: [f"d{i}" for i in range(6 if cell_type == "Dendritic cells" else 4 if cell_type == "Megakaryocytes" else 8)]
        for cell_type in CELL_TYPES
    }
    b_first = valid.index[valid["cell_type"] == "B cells"][0]
    cases = [
        valid.drop(b_first),
        pd.concat([valid, valid.loc[[b_first]]], ignore_index=True),
        valid[~((valid["cell_type"] == "B cells") & (valid["run"] == "lodo_d7"))],
    ]

    for bad_lodo in cases:
        with pytest.raises(ValueError, match=r"B cells.*lodo|lodo.*B cells"):
            lodo_stability(calls, bad_lodo, paired_donors)


# Catches silently deduplicating or filling a missing gene-by-cell-type DE row.
def test_T6_join_integrity():
    valid = pd.DataFrame(_gene_rows("g"))
    bad_inputs = [
        pd.concat([valid, valid.iloc[[0]]], ignore_index=True),
        valid.iloc[1:].copy(),
    ]

    for bad in bad_inputs:
        with pytest.raises(ValueError, match=r"gene|cell.type|duplicate|missing"):
            derive_calls(bad)


# Catches traceability checks that miss stale arrows, counts, or response labels.
def test_T7_traceability_tamper_test():
    calls = pd.DataFrame(_gene_rows("g", {"B cells": "up", "CD4 T cells": "up"}))
    table = label_genes(calls)
    changes = [("B cells", "↓"), ("n_cell_types", 7), ("shared_response", "shared")]

    for column, value in changes:
        tampered = table.copy()
        tampered.loc[0, column] = value
        with pytest.raises(ValueError, match=column.replace("_", ".")):
            check_traceability(tampered, calls)


# Catches early gene filtering, duplicate table rows, or missing arrow columns.
def test_T8_output_shape():
    calls = derive_calls(
        pd.DataFrame(_gene_rows("somewhere", {"B cells": "up"}) + _gene_rows("never"))
    )
    table = label_genes(calls)

    assert calls.groupby("gene").size().to_dict() == {"never": 8, "somewhere": 8}
    assert table["gene"].tolist() == ["somewhere"]
    assert table.loc[0, CELL_TYPES].tolist() == ["↑", "–", "–", "–", "–", "–", "–", "–"]


# Catches dropping megakaryocyte calls or applying its lower-confidence flag elsewhere.
def test_T9_megakaryocyte_flag():
    calls = pd.DataFrame(
        _gene_rows("mega_only", {"Megakaryocytes": "up"})
        + _gene_rows("b_only", {"B cells": "up"})
    )

    table = label_genes(calls).set_index("gene")

    assert table.loc["mega_only", "Megakaryocytes"] == "↑"
    assert bool(table.loc["mega_only", "megakaryocyte_lower_confidence"])
    assert not bool(table.loc["b_only", "megakaryocyte_lower_confidence"])


# Catches large-scale filtering or shared/mixed/single-cell-type classification errors.
@pytest.mark.realdata
def test_T10_real_data(tmp_path):
    run_m2a(REPO_ROOT / "data/processed", tmp_path)
    calls = pd.read_csv(tmp_path / "m2_calls.csv")
    table = pd.read_csv(tmp_path / "m2_table.csv")

    assert len(table) == 2801
    assert table["shared_response"].value_counts().to_dict() == {
        "shared": 21,
        "shared*": 87,
        "partial": 856,
        "specific": 241,
        "specific*": 1585,
        "mixed": 11,
    }
    assert (table["n_cell_types"] == 1).sum() == 1826
    assert (calls["de"] != "none").groupby(calls["gene"]).sum().eq(8).sum() == 21
