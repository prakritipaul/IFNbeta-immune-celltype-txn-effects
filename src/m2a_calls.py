"""M2a: shared vs cell-type-specific calls per gene, LODO stability, and the final table.

Contract: audit/contracts/M2_shared_specific.md (part M2a; decisions O1-O4 and the specific* rule).
Tests: tests/test_m2a.py (written blind by Codex). Invalid input raises ValueError before anything is written.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from src.m1b_de import classify

WELL_POWERED = ["B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells"]
LOW_POWER = ["Dendritic cells", "Megakaryocytes"]
CELL_TYPES = WELL_POWERED + LOW_POWER
EXPECTED_LODO = {**{ct: 8 for ct in WELL_POWERED}, "Dendritic cells": 6, "Megakaryocytes": 4}  # paired donors (A1)
CALL_COLUMNS = ["cell_type", "gene", "log2FC", "qvalue", "ci_low", "ci_high", "tested", "de"]
ARROW = {"up": "↑", "down": "↓", "none": "–"}


def derive_calls(m1b_de, q_threshold=0.05, lfc_threshold=1.0):
    """Long calls per gene × cell type, `de` re-derived from qvalue and unshrunk log2FC (the input `de` is ignored).

    In:  m1b_de with exactly one row per gene for each of the 8 cell types.
    Out: CALL_COLUMNS. E.g. (q 0.0499, log2FC -1.0) → down; (q 0.05, log2FC 2) → none.
    """
    if set(m1b_de["cell_type"]) != set(CELL_TYPES):
        raise ValueError(f"cell_type values must be exactly the 8 cell types, got {sorted(set(m1b_de['cell_type']))}")
    dup = m1b_de.duplicated(["gene", "cell_type"], keep=False)
    if dup.any():
        raise ValueError(f"duplicated gene × cell_type rows, e.g. {m1b_de.loc[dup, ['gene', 'cell_type']].iloc[0].tolist()}")
    per_gene = m1b_de.groupby("gene").size()
    if (per_gene != len(CELL_TYPES)).any():
        raise ValueError(f"missing gene × cell_type rows for {int((per_gene != len(CELL_TYPES)).sum())} genes, "
                         f"e.g. {per_gene[per_gene != len(CELL_TYPES)].index[0]}")
    return classify(m1b_de.drop(columns="de"), q_threshold, lfc_threshold)[CALL_COLUMNS]


def lodo_stability(calls, lodo, paired_donors):
    """Add `lodo_stable`: a DE call keeps its sign (log2FC > 0 for up, < 0 for down) in every LODO refit (O4).

    In:  lodo rows run ("lodo_<donor>"), cell_type, gene, log2FC. Each cell type's runs must be exactly
         {"lodo_<donor>"} over its paired donors, each run holding every gene of `calls` once.
         paired_donors: {cell_type: donors} (required; run_m2a reads it from m1a_samples.csv). Expected runs
         never come from `lodo` itself (Codex re-check: that would validate the input against itself).
         Significance in the refit is not required.
    Out: calls + lodo_stable (boolean; NA where de == "none"). Zero counts as a sign change.
    """
    if set(lodo["cell_type"]) != set(calls["cell_type"]):
        raise ValueError(f"lodo cell types {sorted(set(lodo['cell_type']))} differ from calls")
    runs_by_ct = lodo.groupby("cell_type")["run"].agg(set)
    if set(paired_donors) != set(runs_by_ct.index):
        raise ValueError(f"paired_donors cell types {sorted(paired_donors)} differ from lodo {sorted(runs_by_ct.index)}")
    expected = {ct: {f"lodo_{d}" for d in paired_donors[ct]} for ct in runs_by_ct.index}
    for ct, sub in lodo.groupby("cell_type"):
        dup = sub.duplicated(["run", "gene"])
        if dup.any():
            raise ValueError(f"{ct}: duplicated lodo rows in run {sub.loc[dup, 'run'].iloc[0]}")
        runs = runs_by_ct[ct]
        if runs != expected[ct]:
            raise ValueError(f"{ct}: lodo runs missing {sorted(expected[ct] - runs)}, unexpected {sorted(runs - expected[ct])}")
        if len(runs) != EXPECTED_LODO[ct]:
            raise ValueError(f"{ct}: expected {EXPECTED_LODO[ct]} lodo runs, got {len(runs)} ({sorted(runs)})")
        genes = set(calls.loc[calls["cell_type"] == ct, "gene"])
        for run, rs in sub.groupby("run"):
            if set(rs["gene"]) != genes:
                raise ValueError(f"{ct}: lodo run {run} has {len(set(rs['gene']) ^ genes)} missing or extra genes")
    is_de = calls["de"] != "none"
    m = calls.loc[is_de, ["cell_type", "gene", "de"]].merge(lodo, on=["cell_type", "gene"], validate="one_to_many")
    if m["log2FC"].isna().any():
        bad = m[m["log2FC"].isna()].iloc[0]
        raise ValueError(f"{bad.cell_type}: lodo run {bad.run} has NaN log2FC for DE gene {bad.gene}")
    m["same_sign"] = np.where(m["de"] == "up", m["log2FC"] > 0, m["log2FC"] < 0)
    stable = m.groupby(["cell_type", "gene"])["same_sign"].all()
    out = calls.copy()
    out["lodo_stable"] = pd.array([pd.NA] * len(out), dtype="boolean")
    out.loc[is_de, "lodo_stable"] = stable.reindex(pd.MultiIndex.from_frame(out.loc[is_de, ["cell_type", "gene"]])).to_numpy()
    return out


def label_genes(calls):
    """One row per gene DE in >= 1 cell type: arrows per cell type, n_cell_types, shared_response, flags.

    Labels, in precedence order (O1, O2, specific* rule):
      mixed      up in some cell types and down in others
      shared     DE in all 8 in one direction
      shared*    DE in all 6 well-powered in one direction; each low-power cell type DE that way or its 95% CI
                 reaching |log2FC| >= 1 that way (ci_high >= 1 for up, ci_low <= -1 for down)
      partial    DE in >= 2
      specific*  DE in 1, and some non-DE cell type's CI reaches |log2FC| >= 1 (either direction)
      specific   DE in 1
    CI rules use only cell types where the gene was tested (`tested`, baseMean > 1); a NaN CI never counts.
    """
    if "tested" not in calls.columns:
        raise ValueError("calls need a `tested` column: CI rules use only cell types where the gene was tested")
    def wide(col):
        return calls.pivot(index="gene", columns="cell_type", values=col)[CELL_TYPES]

    de, tested = wide("de"), wide("tested").astype(bool)
    ci_lo, ci_hi = wide("ci_low").where(tested), wide("ci_high").where(tested)  # untested CI → NaN → never counts
    up, dn = de == "up", de == "down"
    n_up, n_dn = up.sum(axis=1), dn.sum(axis=1)
    n = n_up + n_dn
    mixed = (n_up > 0) & (n_dn > 0)
    shared_up = (n_dn == 0) & up[WELL_POWERED].all(axis=1) & (up[LOW_POWER] | (ci_hi[LOW_POWER] >= 1)).all(axis=1)
    shared_dn = (n_up == 0) & dn[WELL_POWERED].all(axis=1) & (dn[LOW_POWER] | (ci_lo[LOW_POWER] <= -1)).all(axis=1)
    shared = shared_up | shared_dn
    reach = ((ci_hi >= 1) | (ci_lo <= -1)) & (de == "none")
    label = np.select([mixed, shared & (n == 8), shared, n >= 2, reach.any(axis=1)],
                      ["mixed", "shared", "shared*", "partial", "specific*"], "specific")
    table = de.replace(ARROW)
    table["n_cell_types"] = n
    table["shared_response"] = label
    table["specific_star"] = table["shared_response"] == "specific*"
    table["megakaryocyte_lower_confidence"] = de["Megakaryocytes"] != "none"  # O3
    table = table[n > 0].reset_index()
    table.columns.name = None
    return table


def check_traceability(table, calls):
    """I2: re-derive every table cell from calls (`de` from qvalue/log2FC by the DE rule) and require 100% agreement.

    Raises ValueError naming the first disagreeing gene and column.
    """
    big = (calls["qvalue"] < 0.05) & (calls["log2FC"].abs() >= 1.0)
    rederived = np.where(big & (calls["log2FC"] > 0), "up", np.where(big & (calls["log2FC"] < 0), "down", "none"))
    if (rederived != calls["de"].to_numpy()).any():
        bad = calls[rederived != calls["de"].to_numpy()].iloc[0]
        raise ValueError(f"calls de disagrees with the DE rule: {bad.gene} in {bad.cell_type}")
    # Independent re-derivation (not label_genes): one gene at a time, straight from the rules.
    rows = []
    for gene, g in calls.groupby("gene", sort=True):
        g = g.set_index("cell_type").loc[CELL_TYPES]
        de_types = [ct for ct in CELL_TYPES if g.at[ct, "de"] != "none"]
        if not de_types:
            continue
        directions = {g.at[ct, "de"] for ct in de_types}

        def ci_reaches(ct, direction):  # untested or NaN CI never counts
            if not g.at[ct, "tested"]:
                return False
            return g.at[ct, "ci_high"] >= 1 if direction == "up" else g.at[ct, "ci_low"] <= -1

        if len(directions) == 2:
            label = "mixed"
        else:
            d = directions.pop()
            all_wp = all(g.at[ct, "de"] == d for ct in WELL_POWERED)
            lp_ok = all(g.at[ct, "de"] == d or ci_reaches(ct, d) for ct in LOW_POWER)
            if all_wp and lp_ok:
                label = "shared" if len(de_types) == 8 else "shared*"
            elif len(de_types) >= 2:
                label = "partial"
            elif any(ci_reaches(ct, "up") or ci_reaches(ct, "down") for ct in CELL_TYPES if ct not in de_types):
                label = "specific*"
            else:
                label = "specific"
        row = {"gene": gene, **{ct: ARROW[g.at[ct, "de"]] for ct in CELL_TYPES}}
        row.update(n_cell_types=len(de_types), shared_response=label, specific_star=label == "specific*",
                   megakaryocyte_lower_confidence="Megakaryocytes" in de_types)
        rows.append(row)
    expected = pd.DataFrame(rows).set_index("gene")
    got = table.set_index("gene")
    if list(got.index) != list(expected.index):
        raise ValueError(f"table genes differ from re-derived genes ({len(got)} vs {len(expected)})")
    if list(got.columns) != list(expected.columns):
        raise ValueError(f"table columns {list(got.columns)} differ from {list(expected.columns)}")
    for col in expected.columns:
        diff = got[col].astype(str) != expected[col].astype(str)
        if diff.any():
            g = diff[diff].index[0]
            raise ValueError(f"{col} disagrees for gene {g}: table {got.at[g, col]!r}, re-derived {expected.at[g, col]!r}")


def run_m2a(in_dir, out_dir):
    """Read m1b_de.csv, the lodo rows of m1b_runs.csv and paired donors from m1a_samples.csv; write m2_calls.csv and m2_table.csv (nothing on failure)."""
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    m1b = pd.read_csv(in_dir / "m1b_de.csv")
    if set(m1b["run"]) != {"main"}:
        raise ValueError(f"m1b_de.csv must hold only the main run, got {sorted(set(m1b['run']))}")
    runs = pd.read_csv(in_dir / "m1b_runs.csv", usecols=["run", "cell_type", "gene", "log2FC"])
    lodo = runs[runs["run"].str.startswith("lodo_")]
    del runs
    samples = pd.read_csv(in_dir / "m1a_samples.csv")
    kept = samples[samples["kept"]]
    n_cond = kept.groupby(["cell_type", "donor"])["condition"].nunique()
    if (n_cond != 2).any():
        raise ValueError(f"kept samples without both conditions: {n_cond[n_cond != 2].index.tolist()}")
    paired = {ct: set(d.index.get_level_values("donor")) for ct, d in n_cond.groupby(level="cell_type")}
    calls = lodo_stability(derive_calls(m1b), lodo, paired)
    table = label_genes(calls)
    check_traceability(table, calls)
    out_dir.mkdir(parents=True, exist_ok=True)
    calls.to_csv(out_dir / "m2_calls.csv", index=False)
    table.to_csv(out_dir / "m2_table.csv", index=False)
