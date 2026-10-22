# M2a test plan

**Status:** APPROVED by the user (2026-10-02), with the DE rule changed to |log2FC| ≥ 1. Codex proposed 11 blind from the contract; Claude merged the two CI-rule tests (T3) and added the label-precedence case to T1. Codex writes exactly the approved tests, in `tests/test_m2a.py`, before the implementation exists.

| # | Test | Bug it catches | Expected answer from | From |
|---|---|---|---|---|
| T1 | **Label oracle:** toy genes across all 8 cell types: all 8 up → `shared`; 6 well-powered up, DCs and megakaryocytes not DE with `ci_high ≥ 1` → `shared*`; 2 DE same direction → `partial`; 1 DE, other CIs inside (−1, 1) → `specific`; 1 DE and another cell type's CI reaching ±1 → `specific*`; up in one, down in another → `mixed`. *Precedence:* 6 well-powered up but DCs DE **down** → `mixed`, not shared. Arrows and `n_cell_types` checked | Wrong label logic or precedence | O1, O2, specific\* rule | Codex + Claude |
| T2 | **DE boundaries re-derived:** (q 0.0499, log2FC ±1.0) → up / down (\|log2FC\| = 1 counts); q = 0.05, \|log2FC\| = 0.9999, or a failing threshold → none; the input `de` column is not trusted | `≤` vs `<` on q; `>` vs `≥` on log2FC; shrunken LFC; trusting M1b's column | DE rule | Codex |
| T3 | **CI rules at the boundary and by direction:** `ci_high = 1.000` counts towards up-shared, `ci_low = −1.000` towards down-shared, 0.999 / −0.999 don't; a CI reaching only the opposite direction doesn't count. For specific\*: any non-DE cell type can trigger it (not only DCs / megakaryocytes), the DE cell type's own CI can't. *Added (user, 2026-10-02):* a CI from a cell type where the gene was not tested (baseMean ≤ 1) or a NaN CI never triggers a star nor rescues shared | Strict vs inclusive; absolute-value logic; wrong scope | O1; specific\* rule | Codex (2 merged) |
| T4 | **LODO sign:** up call with all LODO log2FC > 0 (incl. +0.001, not significant) → `lodo_stable = True`; one refit ≤ 0 → False; mirrored for down | Requiring significance; zero counted as same sign | O4 | Codex |
| T5 | **Incomplete LODO stops:** one expected refit missing, one duplicated, or the wrong count per cell type (8 / 6 / 4) → error naming the cell type and run, no outputs | Stability from incomplete runs | Contract assumption | Codex |
| T6 | **Join integrity:** a duplicated and a missing gene × cell type row in `m1b_de` → error, nothing deduplicated or filled | Silent duplication / loss | Contract assumption | Codex |
| T7 | **Traceability tamper test (I2):** alter one arrow, one `n_cell_types` or one label after building the table → the traceability check fails and names it | Partial checking; stale outputs | I2 (100% agreement) | Codex |
| T8 | **Output shape:** a gene DE somewhere and a gene never DE → calls file has 8 rows for both; table has only the first, once, 8 cell-type columns of ↑ / ↓ / – | Filtering too early; duplicates; missing columns | Output schema | Codex |
| T9 | **Megakaryocyte flag:** a gene DE only in megakaryocytes keeps its call and is flagged lower confidence; the same call in B cells isn't flagged | Dropping low-power calls; flag applied everywhere | O3 | Codex |
| T10 | **Real data:** table has 2,801 genes; `shared` 21, `shared*` 87, `partial` 856, `specific` 241, `specific*` 1,585, `mixed` 11; 1,826 with `n_cell_types` = 1 | Large-scale label or filter errors | Counts computed independently by Claude from `m1b_de.csv` with the tested-only CI rule (2026-10-02); the agent's implementation agrees | Codex |

## Interface the tests call

Naming/interface only. Invalid input raises `ValueError` naming the violated assumption, before anything is written.

```python
# src/m2a_calls.py
WELL_POWERED = ["B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells"]
LOW_POWER = ["Dendritic cells", "Megakaryocytes"]

derive_calls(m1b_de, q_threshold=0.05, lfc_threshold=1.0) -> pd.DataFrame
    # Long: cell_type, gene, log2FC, qvalue, ci_low, ci_high, de ("up"/"down"/"none") re-derived from qvalue and log2FC.
    # Raises on duplicated or missing gene × cell type rows (all 8 cell types × every gene).
lodo_stability(calls, lodo) -> pd.DataFrame
    # lodo: rows run ("lodo_<donor>"), cell_type, gene, log2FC. Adds lodo_stable (bool; NaN where de == "none").
    # Raises if any cell type lacks exactly its expected LODO runs (8 / 6 / 4) or has duplicates.
label_genes(calls) -> pd.DataFrame
    # One row per gene DE in >= 1 cell type: one column per cell type ("↑"/"↓"/"–"), n_cell_types,
    #   shared_response ("shared"/"shared*"/"partial"/"specific"/"specific*"/"mixed"), specific_star (bool),
    #   megakaryocyte_lower_confidence (bool: DE in Megakaryocytes).
check_traceability(table, calls) -> None
    # Re-derives every table cell and label from calls; raises ValueError naming the first disagreement.
run_m2a(in_dir, out_dir) -> None
    # Reads m1b_de.csv and the lodo rows of m1b_runs.csv (usecols only); writes m2_calls.csv and m2_table.csv. Nothing on failure.
```

## Added after Codex review 1 (user-approved), in `tests/test_m2a_review.py`

| # | Test | Bug it catches | Expected answer from |
|---|---|---|---|
| T11 | `check_traceability` rejects a table produced by a consistently wrong `label_genes` (monkeypatched to call a one-cell-type gene `shared`) | Traceability that re-uses the function it checks (review #1) | I2: independent re-derivation |
| T12 | `lodo_stability` with one expected `lodo_<donor>` replaced by `lodo_NOT_A_DONOR` (count unchanged) → `ValueError` naming the missing and unexpected runs | LODO checked by count, not identity (review #2) | Contract: every expected LODO run exists |
