# Codex review 1 and re-check: M2c (Wave 3)

## Review 1

1. pytest result: **6 passed, 3 failed** in 42.25s. T5, T6, and T9 failed because R could not open a port in the sandbox; this is not counted as a code failure.

2. Findings

| # | Severity | file:line | Finding | Concrete failing scenario or reproducer |
|---|---|---|---|---|
| 1 | medium | `src/m2c_programs.py:86-95` | `label_programs` derives the cell-type universe only from returned GSEA rows. If one cell type has no testable sets, it disappears and a program affected in the remaining seven is incorrectly labeled `shared`, contrary to labeling across all 8. | Seven significant rows for `cell_1`–`cell_7`, with `cell_8` absent, produce `n_up=7, label="shared"`. |
| 2 | low | `notebooks/m2c_programs.py:206` | The notebook identifies a specific program’s cell type using `!= "–"`, so it mistakes `n/t` for an affected cell. This already misreports two real programs and risks incorrect FINDINGS text. | `HALLMARK_ANGIOGENESIS` is printed as B cells, but B cells are `n/t`; the affected cell is CD14+ Monocytes. `KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC` is similarly printed as B cells instead of CD4 T cells. |

3. Proposed NEW tests

- Supply eight ranked cell types where one has no gene set passing the 15–500 matched-gene limit. Assert it remains represented as `n/t` and seven same-direction hits are `partial`, not `shared`. This fails currently.
- Exercise the notebook’s specific-cell reporting with an earlier `n/t` cell and a later `↑`/`↓` cell. Assert the reported cell is the arrow-bearing cell. This fails currently.

Test-suite weakness noted: T7 does not test BH calculation or scope—it supplies precomputed `padj` values directly. The current R implementation nevertheless correctly applies BH within each cell type across both collections.

4. Unapproved decisions: **none**. The notebook’s “open question” comment at line 29 is stale; the implemented `average` linkage matches the approved decision.
## Re-check

| Finding | Status | Evidence |
|---|---|---|
| Missing cell type falsely yields `shared` | Fixed | Reproducer returns `partial`, `n_up=7`, and `cell_8="n/t"`; unknown cell types raise `ValueError`. |
| Notebook mistakes `n/t` for affected cell | Fixed | Arrow-only lookup confirmed. Real output reports `HALLMARK_ANGIOGENESIS → CD14+ Monocytes` and KEGG ARVC → `CD4 T cells`. |
| Regression checks | Pass | `n/t` marker retained; T6 uses its revised distributed 501-gene toy; R implementation unchanged. |

Pytest: **7 passed, 3 sandbox-only failures**. T5, T6, and T9 failed solely because R could not open a port; not counted as code failures.

New issues: **none**.

**Verdict: PASS**