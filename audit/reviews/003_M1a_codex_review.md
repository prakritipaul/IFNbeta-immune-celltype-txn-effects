# Codex review 1 and re-check: M1a (Wave 1)

## Review 1

1. Pytest: **9 passed, 1 skipped** in 63.21s. The skip is the user-approved dropped T6. Snapshot: `23b171382c39442d5f9950536ee150046819d1df`.

2. Findings

| # | Severity | file:line | Finding | Concrete failing scenario or reproducer |
|---|---|---|---|---|
| 1 | High | [R/m1a_varpart.R:59](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-acfad106b7a1d751e/R/m1a_varpart.R:59), [R/m1a_cancor.R:22](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-acfad106b7a1d751e/R/m1a_cancor.R:22) | Pairing validation checks only two rows per donor, not one row for each condition. Malformed unpaired metadata is accepted and produces apparently valid variance fractions. `m1a_cancor.R` has no pairing check. | Constructed 16 samples where `d0` had two controls, `d1` two treatments, and the other donors were paired. `m1a_varpart.R` exited successfully and wrote all three outputs. |
| 2 | Medium | [notebooks/m1a_pseudobulk.py:32](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-acfad106b7a1d751e/notebooks/m1a_pseudobulk.py:32), [notebooks/m1a_pseudobulk.py:79](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-acfad106b7a1d751e/notebooks/m1a_pseudobulk.py:79), [notebooks/m1a_pseudobulk.py:106](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-acfad106b7a1d751e/notebooks/m1a_pseudobulk.py:106) | The notebook hardcodes `ctrl`/`stim` instead of obtaining the declared control/treatment labels. This contradicts the requested condition-label robustness. | Valid outputs using labels `control`/`treated` fail at the pivot selection or colour lookup. |
| 3 | Low | [src/m1a_pseudobulk.py:75](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-acfad106b7a1d751e/src/m1a_pseudobulk.py:75), [src/m1a_pseudobulk.py:83](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-acfad106b7a1d751e/src/m1a_pseudobulk.py:83) | Sample metadata is reconstructed by splitting an unescaped `\|`-delimited ID. A delimiter inside any label silently corrupts the output schema. | With cell type `A\|B`, the returned metadata becomes `cell_type=A`, `donor=B`, `condition=d`; both samples are incorrectly marked dropped. |

3. Proposed new tests

- Pass duplicated same-condition rows for one donor to both R scripts; require a non-zero exit and no output files.
- Execute the notebook against otherwise-valid fixtures whose condition labels are `control`/`treated`; require successful completion and both labels in the plots/tables.
- Call `pseudobulk` with `|` in a donor or cell-type label; require an explicit `ValueError` before aggregation.

4. Unapproved decisions

None. The hardcoded condition names are an implementation violation, not a documented decision.

## Re-check after fixes

| Finding | Status | Evidence |
|---|---|---|
| M1a #1 — R pairing validation | FIXED | Both scripts rejected duplicated same-condition rows, exited non-zero, and wrote no outputs. |
| M1a #2 — notebook condition labels | FIXED | Notebook ran end-to-end with `control`/`treated`; both labels appeared correctly in outputs and plots. |
| M1a #3 — pipe-delimited labels | FIXED | Donor and cell-type labels containing `\|` raised `ValueError` before aggregation. |
| M3 #2 — explicit float32 cast | NOT FIXED | `src/m3_signature.py` is absent from this branch; the reproducer fails during collection with `ModuleNotFoundError`. No M3 source fix appears in the diff. |

Pytest: **12 passed, 1 skipped in 62.49s**. The skip is the user-approved dropped T6.

No new M1a issue was introduced by the fixes.

**Overall: PARTIAL** — M1a passes step 6; the approved M3 fix is missing from this branch.