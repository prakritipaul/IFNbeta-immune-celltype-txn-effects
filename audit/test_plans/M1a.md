# M1a test plan

**Status:** APPROVED by the user (2026-10-01). Codex proposed 12 blind from the contract; Claude merged overlaps and added T4b and T7; the user dropped Codex's hand-summed toy and shuffled-metadata tests. Codex writes exactly the approved tests, in `tests/test_m1a.py`, before the implementation exists.

| # | Test | Bug it catches | Expected answer from | From |
|---|---|---|---|---|
| T1 | A1 boundary: donor A 9 ctrl / 12 stim → both A samples `kept=False`; donor B 10 / 10 → both kept | `<=` vs `<`; dropping only the small sample | A1 | Codex |
| T2 | Real data: 128 samples, 116 kept; dropped pairs exactly DC × {107, 1039}, megakaryocytes × {101, 107, 1016, 1039} | Wrong grouping; rule applied globally | A1 + `data_facts.md` table | Codex + Claude |
| T3 | Real data: every gene's sum over samples = its sum in `m0_cells.h5ad`; `n_cells` sums to 24,562 | Lost or duplicated cells | I1 identity | Codex |
| T4 | No-outlier variant: `n_cells` sums to 24,048 (24,562 − 514), gene totals match the non-flagged cells; main run unchanged. **(b)** A1 is re-applied after removing outliers (a group can fall below 10) | Reversed filter; stale `kept` copied from the main run | Contract; 514 flagged cells | Codex + Claude |
| T5 | Bad input stops with an error and writes no outputs: negative or fractional count; missing `obs` column; a kept donor lacking one condition | Rounding, silent workarounds, unpaired samples | Contract assumptions | Codex |
| ~~T6~~ **Dropped by the user 2026-10-01** (toy genes had zero variance; fixing it judged overkill; the baseMean filter is now untested) ~~baseMean filter, per cell type: genes at baseMean just below 1, exactly 1, just above 1 → only the last enters `m1a_varpart.csv`; a gene can pass in one cell type and fail in another~~ | `>=` vs `>`; one filter for all cell types | Decision (baseMean > 1) | Codex + Claude |
| T7 | Variance-partitioning oracle: simulated pseudobulk where gene X varies only by donor and gene Y only by condition (small noise) → donor fraction of X > 0.9, condition fraction of Y > 0.9; every row's fractions ≥ 0 and sum to 1 | Wrong model terms; swapped columns; untransformed input | Construction; variance arithmetic | Claude + Codex |

## Interface the tests call

Naming/interface only; no result-changing choices. Invalid input raises `ValueError` naming the violated assumption, before any file is written.

```python
# src/m1a_pseudobulk.py
pseudobulk(adata, min_cells=10) -> tuple[pd.DataFrame, pd.DataFrame]
    # adata: raw counts in X; obs has donor, condition, cell_type; uns['control'], uns['treatment'].
    # Returns (counts, samples). counts: genes × samples, int, columns = sample ids "cell_type|donor|condition",
    #   one column per donor × cell type × condition present in the data.
    # samples: one row per sample: sample, cell_type, donor, condition, n_cells, total_counts, kept.
    # kept = False for both samples of a donor × cell type when either condition has < min_cells cells (A1);
    #   a missing condition counts as 0 cells.
    # Raises: negative / fractional / non-finite counts; missing obs column; a donor with no cells in one condition overall.
run_m1a(in_path, out_dir, min_cells=10) -> None
    # Writes m1a_counts.csv, m1a_samples.csv, m1a_counts_no_outliers.csv, m1a_samples_no_outliers.csv
    #   (the last two from cells with flag_outlier == False; A1 applied afresh). Writes nothing if validation fails.
M1A_INPUT = "data/processed/m0_cells.h5ad"
```

```
# R/m1a_varpart.R (called as a subprocess)
Rscript R/m1a_varpart.R <counts.csv> <samples.csv> <out.csv> <min_basemean>
    # Uses kept samples only; per cell type: DESeq2 size factors → keep genes with baseMean > min_basemean
    #   → varianceStabilizingTransformation(blind=TRUE) → variancePartition ~ (1|donor) + (1|condition) + log(n_cells) + log(total_counts).
    # out.csv columns: cell_type, gene, baseMean, donor, condition, log_n_cells, log_total_counts, Residuals.
    # Non-zero exit status on invalid input.
```

## Added after Codex review 1 (user-approved 2026-10-01), in `tests/test_m1a_review.py`

| # | Test | Bug it catches | Expected answer from |
|---|---|---|---|
| T8 | `R/m1a_varpart.R` and `R/m1a_cancor.R` given kept samples where one donor has two control rows and another two treatment rows (row counts still 2 per donor) → non-zero exit, no output file | Pairing checked by row count, not one row per condition (review #1) | Contract assumption: every kept sample has both its donor's conditions |
| T9 | `pseudobulk` with `\|` inside a donor or cell-type label → `ValueError` before aggregation | Sample ID `cell_type\|donor\|condition` split wrongly (review #3) | Sample-ID schema |
