# Codex review of the planted-bug snapshot: M2a (not told about the bug)

| # | Severity | file:line | Finding | Evidence |
|---|---|---|---|---|
| 1 | Medium | `src/m2a_calls.py:105` | Upward `shared*` rescue uses strict `ci_high > 1`; the approved rule is inclusive (`≥ 1`). | The boundary reproducer expects `shared*` at `ci_high == 1.0` but receives `partial`; `test_T3_ci_rules_boundary_and_direction` fails. The independent traceability path correctly uses `>= 1` at line 143, so such input would ultimately abort rather than write output. Downward and `specific*` boundaries are correctly inclusive. |

Verdict: M2a is not fully correct as written; all other requested rules checked conform, but the upward shared-rescue boundary must be inclusive.