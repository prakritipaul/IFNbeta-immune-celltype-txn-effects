# M3 reviewer brief (Codex, commit 02a32b3)

## 1. What M3 does, in 3 lines

- Input: 24,562 filtered cells × 15,701 genes plus MSigDB 2026.1 Hallmark IFN-α (97 genes) and IFN-γ (200 genes), both at 100% coverage after six approved symbol mappings.
- Processing: scores every cell for IFN-α, full IFN-γ, and IFN-γ with the 73 shared genes removed (127 genes remain).
- Outputs: `m3_cell_scores.csv` has 24,562 unique cells × 7 columns; `m3_sample_scores.csv` has 116 donor × cell-type × condition rows × 7 columns, retaining 24,451 cells after the approved paired 10-cell rule.

## 2. Assumptions coded in

- The GMT must have the approved SHA-256 and both named Hallmark sets ([src/m3_signature.py:34](../../src/m3_signature.py#L34), [src/m3_signature.py:41](../../src/m3_signature.py#L41)); this prevents a different release from silently changing the scores.
- Only the six approved symbol mappings are allowed, and every mapped gene must exist ([src/m3_signature.py:55](../../src/m3_signature.py#L55), [src/m3_signature.py:65](../../src/m3_signature.py#L65)); for example, `PSMA2` must map to `PSMA2.1`, not be dropped.
- Counts are cast to float32, normalized to 10,000 per cell, natural-log transformed, then scored with the approved fixed Scanpy parameters ([src/m3_signature.py:79](../../src/m3_signature.py#L79), [src/m3_signature.py:82](../../src/m3_signature.py#L82), [src/m3_signature.py:89](../../src/m3_signature.py#L89)); these choices set the numerical score scale.
- The IFN-γ specificity control removes every gene also in IFN-α ([src/m3_signature.py:84](../../src/m3_signature.py#L84)); this changes the full 200-gene set to a 127-gene control.
- Sample scores are cell means, and a donor is retained for a cell type only when every condition has at least 10 cells ([src/m3_signature.py:117](../../src/m3_signature.py#L117), [src/m3_signature.py:121](../../src/m3_signature.py#L121)); for example, 9 control/12 stimulated drops both rows.

## 3. Decisions the user didn't make

None found. Every result-affecting choice I found is in the approved contract or `audit/decisions.md`.

## 4. The lines worth reading yourself

- [src/m3_signature.py:56](../../src/m3_signature.py#L56): applies the exact gene-symbol substitutions; a plausible wrong mapping changes a signature without crashing.
- [src/m3_signature.py:82](../../src/m3_signature.py#L82): performs per-cell library-size normalization; scoring raw counts would turn depth into apparent IFN response.
- [src/m3_signature.py:87](../../src/m3_signature.py#L87): constructs IFN-γ without shared genes; subtracting the wrong set would undermine the specificity comparison.
- [src/m3_signature.py:118](../../src/m3_signature.py#L118): takes the mean cell score per sample; a median here would quietly change every donor-level result.
- [src/m3_signature.py:123](../../src/m3_signature.py#L123): applies the paired 10-cell threshold; `> 10` instead of `>= 10` would silently discard boundary samples.

## 5. What nothing checks

- No test supplies a heterogeneous toy matrix with an independently calculated expected Scanpy score. The current tests check invariants and passed parameters, so code that negated all three computed scores after scoring could still pass: zero remains zero, proportional cells remain equal, and real-data row counts remain correct.
- The notebook reports P5, P6, and I6 but deliberately assigns no pass/fail criteria. For example, it displays score-versus-library-size correlations but nothing asserts how large a within-cell-type/condition correlation is acceptable; that scientific judgment remains with the user.
- The planted-bug check was skipped for Wave 1 by the user's recorded decision, so M3 has no evidence that its suite catches an unknown realistic mutation.

## 6. Verdict

The targeted suite passes: **10/10 tests** (`tests/test_m3.py` and `tests/test_m3_review.py`, including the real-data rerun). I independently confirmed the saved row counts, unique barcodes, 116-sample pairing pattern, 100% gene-set coverage, and absence of missing/non-finite scores. The paired median IFN-α increase is positive in all eight cell types (0.452 to 0.898), which is directionally sensible but not an independent biological validation.

Before approving M3, the user should inspect the notebook's three scientific diagnostics: donor-paired increases (P5), whether IFN-α still exceeds IFN-γ after shared genes are removed (P6), and within-cell-type/condition dependence on library size (I6). Approval should acknowledge the missing independent numerical-score oracle and the skipped planted-bug check.
