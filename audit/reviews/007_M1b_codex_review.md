# Codex review 1 and re-check: M1b (Wave 2)

## Review 1

1. Pytest: `14 passed, 2 deselected in 159.04s` using the requested command.

2. Findings

| # | Severity | file:line | Finding | Concrete failing scenario or reproducer |
|---|---|---|---|---|
| 1 | medium | [src/m1b_de.py:131](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-ad831d5818f7dad3c/src/m1b_de.py:131) | R-output validation checks gene counts and uniqueness, but not that returned gene identities equal the input genes. M2 could therefore receive silently substituted/missing genes. | A fake R result replacing `real_a, real_b` with `wrong_a, wrong_b` was accepted; `run_deseq2()` returned `ACCEPTED_GENES ['wrong_a', 'wrong_b']`. |
| 2 | low | [src/m1b_de.py:213](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-ad831d5818f7dad3c/src/m1b_de.py:213) | `m1b_de.csv` contains `run`, although the approved contract says `m1b_runs.csv` has the same schema **plus** `run`. The real file has 15 columns, including constant `run="main"`. | An exact-schema assertion derived from the contract fails. Existing T8 instead includes `run` in its expected main-file schema, so it codifies the mismatch. |
| 3 | low | [notebooks/m1b_de.py:44](/Users/prakritipaul/Git/IFNbeta-immune-celltype-txn-effects/.claude/worktrees/agent-ad831d5818f7dad3c/notebooks/m1b_de.py:44) | The notebook reloads all columns of `m1b_runs.csv` solely to display its shape/head/tail, after already loading selected columns. The real file is 785 MB, creating avoidable multi-GB memory pressure and risking failure of the required return-package notebook. | Instrument `pd.read_csv` to reject unrestricted reads of `m1b_runs.csv`; execution fails at line 46. |

3. Proposed NEW tests

- R-boundary gene-set test: fake an R response with the correct row count but one substituted gene; require `run_deseq2()` to raise `ValueError`. This fails currently.
- Exact output-schema test: require `m1b_de.csv` to omit `run`, while `m1b_runs.csv` includes it; also assert exact column order and key uniqueness. This fails currently.
- Notebook large-file access test: execute/instrument the notebook and forbid a full unrestricted `read_csv` of `m1b_runs.csv`; require chunked/head-tail access. This fails currently.

4. Unapproved decisions

None beyond finding #2’s contract-schema deviation. Contrast direction, reference level, designs, shuffle/LODO construction, qvalue scope, post-fit baseMean filtering, CI, classification, Cook’s handling, and `lfcShrink(type="normal")` matched the approved decisions in the reviewed paths.

Audit method followed: [continuous-independent-audit skill](file:///Users/prakritipaul/.agents/skills/continuous-independent-audit/SKILL.md).
## Re-check after fixes

| Finding | Status | Evidence |
|---|---|---|
| #1 Gene identity validation | FIXED | Substituted-gene reproducer/T12 raises `ValueError`; all 8 real cell types exactly match the 15,701 input genes. |
| #2 `run` column schema | RESOLVED-BY-CONTRACT | Contract now specifies constant `run="main"`; both output files have the matching 15-column schema. |
| #3 Notebook memory pressure | FIXED | Instrumented execution passed; all `m1b_runs.csv` reads use `usecols`, `nrows`, or `skiprows`. |

Pytest: `15 passed, 2 deselected in 138.00s (0:02:17)`.

Diff reviewed at `9f06e21e`; R code unchanged, no new issue found.

**Verdict: PASS**