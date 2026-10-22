# Planted-bug checks

A small, realistic bug is planted on a throwaway branch without telling the reviewers, to test whether the tests and reviews catch it. The bug never reaches `planning`/`main`. The bug was written down before any check ran.

| # | Module | Bug planted | Tests | Codex review (not told) | Gap closed by |
|---|---|---|---|---|---|
| 1 | M0 | `flag_qc` outlier loop checks only `n_counts`; `n_genes` dropped ("tidy flag_qc loop") | **Missed**: 24/24 passed | **Caught** as finding #1 (high): "a cell with ordinary total counts but log1p(n_genes) more than 5 MADs … receives flag_outlier=False" | T18: passes on the real code, **fails on the bugged code** |
| 2 | M1b | `R/m1b_deseq2.R`: qvalue run on **all** genes with a p-value (incl. baseMean ≤ 1), untested genes masked afterwards ("qvalue universe") | **Missed**: 15/15 fast tests passed (real-data T7/T11 not run) | **Caught** as finding #1 (high): toy p = 0.0001 gives q 0.019 tested-only vs 0.040 with all genes | T13: passes on the real code, **fails on the bugged code** |
| 3 | M2a | `label_genes`: upward shared\* rescue uses `ci_high > 1` instead of `≥ 1` (boundary) | **Caught**: T1 and T3 failed (2 of 12) | **Caught** as finding #1 (medium), with a reproducer; also noted the independent traceability check would abort the real run | Not needed: tests already catch it |
| 4 | M2c | `R/m2c_fgsea.R`: BH run **per collection** (Hallmark and KEGG separately) instead of over both within the cell type | **Missed**: 10/10 passed (T7 supplied precomputed padj) | **Caught** as finding #1 (high) | T11: passes on the real code, **fails on the bugged code** |

**Why the tests missed it:** every outlier test (T6, T16) holds `n_genes` constant, so dropping the genes metric changed no expected value.

**M1b, why the tests missed it:** T4 checks that untested genes get no q-value, and T5 that cell types don't share a correction, but nothing checked *which* p-values qvalue was run on within a cell type.

**M2c, why the tests missed it:** T7 checked the significance boundary on padj values handed in directly, so nothing exercised the R script's BH step; Codex had flagged exactly this weakness in its review.
