# M0 reviewer brief (Codex, commit 142cec1)

## 1. What M0 does, in 3 lines

Reads the Kang count matrix (24,673 cells × 15,706 genes) and GEO metadata (29,065 droplets), then joins all 24,673 cells to the authors' labels using barcode sequence + condition.
Computes counts, detected genes and % ribosomal; the approved filters remove 111 cells and 5 genes, leaving 24,562 × 15,701.
Writes `m0_cells.h5ad`, cell accounting and EDA/QC reports; 514 outliers and 7 small donor–condition–cell-type groups are flagged, not removed.

## 2. Assumptions coded in

- [src/m0_qc.py:109](../../src/m0_qc.py#L109): identity is barcode sequence + condition. This keeps the 234 barcode sequences present in both runs as 468 distinct cells.
- [src/m0_qc.py:118](../../src/m0_qc.py#L118): donor, cluster and t-SNE must agree between files, and every joined cell must be a singlet. A mismatched source stops rather than silently relabelling cells.
- [src/m0_qc.py:180](../../src/m0_qc.py#L180): cells are filtered first at ≥200 detected genes, then genes at ≥3 cells across the survivors. Reversing that order can change the retained gene list.
- [src/m0_qc.py:214](../../src/m0_qc.py#L214): outliers are assessed within donor × condition × cell type, on log counts **or** log genes, with groups under 10 skipped. For example, a high-RNA dendritic cell is compared with dendritic cells, not pooled T cells.
- [src/m0_qc.py:255](../../src/m0_qc.py#L255): the result must have one recognized control and one treatment, 8 donors with both conditions, unique gene symbols, and the exact 8 author cell-type names.

## 3. Decisions she didn't make

None found. Result-affecting thresholds, grouping, ribosomal-gene definition, retained cells/donors/cell types, and condition handling are recorded in the contract or `audit/decisions.md`.

## 4. The lines worth reading yourself

- [src/m0_qc.py:113](../../src/m0_qc.py#L113): constructs the cell-side join key; a plausible wrong key would silently swap or lose stimulated collision barcodes.
- [src/m0_qc.py:167](../../src/m0_qc.py#L167): defines % ribosomal as ribosomal counts / total counts; a wrong axis or denominator would still yield plausible percentages.
- [src/m0_qc.py:218](../../src/m0_qc.py#L218): unions the two MAD tests; the planted bug removed `n_genes`, passed the original tests, and was caught in review.
- [notebooks/00_eda.py:183](../../notebooks/00_eda.py#L183): checks the notebook's filtering agrees in shape with the independently rerun pipeline.
- [notebooks/00_eda.py:346](../../notebooks/00_eda.py#L346): computes stimulated-minus-control composition shifts; reversing this subtraction would reverse the biological interpretation without crashing.

## 5. What nothing checks

- No test executes or visually judges `00_eda.ipynb`. For example, swapped colours or obscured points in the t-SNE could survive while every numerical test passes.
- No assertion pins the real-data final counts and flags to 24,562 cells, 15,701 genes, 514 outliers and 7 small groups. I independently confirmed the saved shape and 514 flags, but a future intentional threshold change could alter them without a focused regression failure.
- The notebook recomputes tables using the same functions as the pipeline. For example, its per-metric outlier assertion would not detect a shared error in `mad_outlier`; the toy oracle tests, not that notebook assertion, provide the independent check.
- The tests establish code behavior, not scientific adequacy: they cannot decide whether the 514 retained outliers, the 4 paired megakaryocyte donors, or the unmodelled condition/10x-run confounding make later conclusions credible.

## 6. Verdict

**Ready for owner review, not automatic approval.** All 35 tests pass (including real-data tests); the planted one-metric outlier bug was caught by review and T18 now fails on that bug. Before approving M0, inspect the executed notebook's authors' t-SNE, the two donor × condition × cell-type outlier panels, the composition plot (especially patient_107), and covariate plots (% ribosomal by condition/donor). Confirm that these visuals support the already-recorded choices to keep all cells, donors and cell types and to treat 10x-run confounding as a caveat.
