# Codex review 1 and re-checks: M2a (Wave 3)

## Review 1

1. pytest: **10 passed in 9.96s**. Command exited 0; only benign Arrow CPU-cache warnings.

2. Findings

| # | Severity | file:line | Finding | Concrete failing scenario or reproducer |
|---|---|---|---|---|
| 1 | medium | [src/m2a_calls.py:122](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-ad1f077ceace0031e/src/m2a_calls.py:122) | I2 traceability is not independent: it calls `label_genes`, the function it is supposed to check. A shared defect in label construction and validation passes unnoticed. T7 only tampers with the finished table, so it misses this weakness. | Monkeypatching `label_genes` to return a table incorrectly labeling a one-cell-type gene as `shared` makes `check_traceability` return successfully. |
| 2 | medium | [src/m2a_calls.py:48](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-ad1f077ceace0031e/src/m2a_calls.py:48) | LODO validation checks only the number of distinct run names, not that the expected donor-specific refits exist. Thus a missing donor can be replaced by an arbitrary `lodo_*` run while satisfying the contract’s count check. T5 tests deletion and duplication, but not identity substitution. | Replacing B-cell `lodo_d7` with `lodo_NOT_A_DONOR` preserves eight runs; `lodo_stability` accepts it and reports the call stable. |

3. Proposed NEW tests

- **Independent traceability oracle:** replace or monkeypatch `label_genes` with a consistently wrong implementation, then require `check_traceability` to reject its table by independently reconstructing arrows, counts, flags, and label precedence from primitive call fields. This fails currently.
- **Exact LODO identity:** construct the expected donor-run set, replace one expected run with `lodo_NOT_A_DONOR` while preserving the required count, and require a `ValueError` naming the missing/unexpected run. This fails currently.

4. Unapproved decisions

None found. The DE boundary, tested-only directional CI rules, label precedence, megakaryocyte flag, schemas, arrows, and validation-before-writing behavior otherwise match the approved contract and decisions.
## Re-check

| Finding | Status | Evidence |
|---|---|---|
| Independent traceability | FIXED | Wrong `label_genes` output was rejected. Independent logic correctly covers precedence, directional/tested-only CI rules, stars, arrows, counts, and megakaryocyte flag. |
| Exact LODO identity | PARTIAL | Explicit `paired_donors` rejects missing/unexpected runs, and production passes it. But the optional fallback accepted a DC substitution `lodo_d0 → lodo_d6`; it infers expectations from the input being validated. |

Pytest: **12 passed in 25.21s**.

New issue: **medium** — `lodo_stability(..., paired_donors=None)` remains a validation loophole for low-power cell types. Remove the fallback or require `paired_donors`.

Verdict: **PARTIAL**

Audit followed the [continuous-independent-audit skill](/Users/prakritipaul/.agents/skills/continuous-independent-audit/SKILL.md).
Follow-up: the fallback was removed (paired_donors required) and Codex's DC substitution case added to T12; 12/12 pass.
