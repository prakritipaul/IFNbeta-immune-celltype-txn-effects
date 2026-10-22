# M1a reviewer brief (Codex, commit 02a32b3)

## 1. What M1a does, in 3 lines

Input: 24,562 cells × 15,701 genes of raw counts, with donor, condition, cell type and the 514 precomputed outlier flags.
Pseudobulk: sums those cells into 128 donor × cell-type × condition samples; the approved 10-cell paired rule keeps 116 samples containing 24,451 cells (the no-outlier sensitivity input keeps 116 containing 23,941 of 24,048 cells).
Outputs: 15,701-gene count matrices and sample tables, plus variance fractions for 51,581 cell-type × gene combinations across all 8 cell types, their VST values, fit types and 128 term-pair correlations.

## 2. Assumptions coded in

- A sample is usable only when both conditions for that donor × cell type have at least 10 cells ([src/m1a_pseudobulk.py:97](../../src/m1a_pseudobulk.py#L97)). For example, 9 control/12 stimulated drops both samples; this preserves pairing but excludes 111 cells from the main analysis.
- The input is raw, finite, non-negative integer counts ([src/m1a_pseudobulk.py:43](../../src/m1a_pseudobulk.py#L43)). This matters because summing normalized expression would not produce valid DESeq2 counts.
- Every kept donor has exactly one sample in each of two conditions within each cell type ([R/m1a_varpart.R:60](../../R/m1a_varpart.R#L60)). This prevents a duplicated control row from masquerading as a pair.
- Genes enter variance partitioning separately by cell type only when DESeq2 normalized baseMean is strictly greater than 1 ([R/m1a_varpart.R:66](../../R/m1a_varpart.R#L66)). Thus a gene can be analysed in B cells but omitted in NK cells.
- Variance is fitted on blind DESeq2 VST values with donor and condition as random effects and log cell count and log total count as fixed effects ([R/m1a_varpart.R:71](../../R/m1a_varpart.R#L71), [R/m1a_varpart.R:77](../../R/m1a_varpart.R#L77)). These choices determine what each reported fraction means.

## 3. Decisions the user didn't make

None found. The result-affecting choices in the implementation trace to the approved contract or `audit/decisions.md`, including the 10-cell rule, all-cells main analysis, no-outlier sensitivity input, transform, model, gene filter and pinned library parameters.

## 4. The lines worth reading yourself

- [src/m1a_pseudobulk.py:83](../../src/m1a_pseudobulk.py#L83): the sparse matrix multiplication that actually sums cells into samples.
- [src/m1a_pseudobulk.py:97](../../src/m1a_pseudobulk.py#L97): the `>= 10 in both conditions` rule; reversing this comparison would silently change included samples.
- [R/m1a_varpart.R:67](../../R/m1a_varpart.R#L67): the strict, per-cell-type `baseMean > 1` gene filter.
- [R/m1a_varpart.R:71](../../R/m1a_varpart.R#L71): the blind VST and explicit parametric fit request.
- [R/m1a_varpart.R:77](../../R/m1a_varpart.R#L77): the variance model and its explicit REML/weight settings; wrong terms could still produce plausible fractions.

## 5. What nothing checks

- The strict baseMean boundary has no running test: the user dropped T6. A future `> 1` → `>= 1` change could admit genes exactly at 1 without a failure.
- No numerical oracle checks `canCorPairs`; T8 checks only that malformed pairs are rejected. For example, swapping two term labels in the 128-row correlation output could remain unnoticed.
- The notebook's PCA coordinates and centroid-distance calculation are not asserted. For example, measuring distance after scaling genes instead of in the approved unscaled VST space could yield a convincing but different outlier plot.
- Regression sensitivity was not demonstrated because the user explicitly skipped Wave 1's planted-bug check. Passing tests therefore have not been shown to catch a controlled realistic defect.

## 6. Verdict

The targeted suite passes: **12 passed, 1 skipped** (T6), including real-data count conservation, the variance-partition oracle and the two review regressions. I independently confirmed 128/116 total/kept samples, 24,562 input cells, variance fractions summing to 1 within `2.3 × 10⁻15`, and parametric fits logged for all 8 cell types. No unapproved result-affecting choice or confirmed defect was found.

Before approving M1a, the user must inspect the per-cell-type PCA/distances and the variance-fraction plus collinearity plots: those scientific judgments are intentionally not automated. Approval should retain the explicit evidence limits above, especially the untested baseMean boundary and skipped planted-bug check.
