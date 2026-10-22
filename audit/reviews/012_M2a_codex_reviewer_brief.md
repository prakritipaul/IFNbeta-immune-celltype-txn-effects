# M2a reviewer brief (Codex, commit 827bea4)

## 1. What M2a does, in 3 lines

Input: 125,608 main-DE rows (15,701 genes × 8 cell types), 910,658 LODO rows from `m1b_runs.csv`, and the 116 kept pseudobulk samples.
Output: `m2_calls.csv`, 125,608 rows with re-derived up/down/none calls and LODO stability; 4,947 calls are DE and all 4,947 retain direction in every LODO refit.
Output: `m2_table.csv`, 2,801 genes × 13 columns: 21 shared, 87 shared*, 856 partial, 241 specific, 1,585 specific*, and 11 mixed.

## 2. Assumptions coded in

- Every gene must occur exactly once in each of the eight named cell types ([src/m2a_calls.py:27](../../src/m2a_calls.py#L27), [src/m2a_calls.py:32](../../src/m2a_calls.py#L32)). This prevents a missing cell type from making a gene look specific.
- DE means `qvalue < 0.05` and `|log2FC| >= 1`, re-derived rather than trusting M1b's `de` column ([src/m2a_calls.py:21](../../src/m2a_calls.py#L21)). For example, q = 0.05 is not DE, while log2FC = 1 exactly can be.
- “Shared” requires one direction in all six well-powered types; DC and megakaryocytes may be rescued only by a tested CI reaching the same-direction threshold ([src/m2a_calls.py:99](../../src/m2a_calls.py#L99), [src/m2a_calls.py:105](../../src/m2a_calls.py#L105)). Thus OAS2, DE in seven types, is `shared*`, not `shared`.
- A one-type call becomes `specific*` when a tested, non-DE cell type has a CI reaching either ±1 ([src/m2a_calls.py:108](../../src/m2a_calls.py#L108)). Untested or missing CIs do not count.
- LODO stability means sign only—not significance—and requires the exact omitted-donor run set: 8, 6, or 4 refits by cell type ([src/m2a_calls.py:54](../../src/m2a_calls.py#L54), [src/m2a_calls.py:62](../../src/m2a_calls.py#L62), [src/m2a_calls.py:73](../../src/m2a_calls.py#L73)). A zero refit effect is unstable.
- Any megakaryocyte DE call is retained but flagged lower confidence ([src/m2a_calls.py:115](../../src/m2a_calls.py#L115)); 27 genes carry that flag.

## 3. Decisions the user didn't make

None found. The thresholds, tested-only CI rule, label precedence, low-power handling, LODO definition, and megakaryocyte flag all trace to the approved contract or `audit/decisions.md`.

## 4. The lines worth reading yourself

- [src/m2a_calls.py:36](../../src/m2a_calls.py#L36): rebuilds DE calls from q-value and unshrunk log2FC; a wrong boundary silently changes arrows and gene inclusion.
- [src/m2a_calls.py:73](../../src/m2a_calls.py#L73): defines “same sign” for every LODO refit; `>= 0` here would make zero effects look stable.
- [src/m2a_calls.py:105](../../src/m2a_calls.py#L105): implements same-direction `shared*` CI rescue; the planted `> 1` bug looked plausible but lost exact-boundary genes.
- [src/m2a_calls.py:109](../../src/m2a_calls.py#L109): sets label precedence (`mixed` before shared/partial/specific); reordering can plausibly relabel conflicting responses.
- [src/m2a_calls.py:145](../../src/m2a_calls.py#L145): independently rebuilds every final label for traceability rather than calling `label_genes` again.

## 5. What nothing checks

- Numeric input validity is not asserted. For example, a row with `ci_low > ci_high`, an infinite log2FC, or a non-boolean `tested` value could yield a plausible star/DE label instead of stopping; current upstream files do not show an observed mismatch, but no M2a test guards this.
- The two output writes are not atomic ([src/m2a_calls.py:194](../../src/m2a_calls.py#L194)). For example, a disk error after `m2_calls.csv` is replaced but before `m2_table.csv` is written could leave outputs from different runs; tests cover validation failures, not I/O interruption.
- The real-data test freezes aggregate label counts, not individual gene identities. A swap between two genes that preserves all six totals would escape T10, although the internal traceability assertion still checks the table against the same input on each run.

## 6. Verdict

No correctness defect found. `pytest tests/test_m2a.py tests/test_m2a_review.py -q` passes **12/12** on commit `827bea4`; the planted inclusive-CI bug was caught by T1 and T3. Independent checks found zero DE-rule mismatches and reproduced all output dimensions and label totals; classic ISGs are shared/shared* (for example, ISG15 and MX1 are shared, OAS2 and STAT1 shared*).

Before approving M2a, the user should inspect the explained examples/output display and decide whether the scientifically dominant result is credible: **1,585 of 2,801 genes are `specific*`**, versus 241 unstarred specific genes, because uncertainty in another tested cell type reaches ±1. They should also notice that all 4,947 DE calls are LODO direction-stable. These are consequences of approved rules, not implementation failures, but they drive the interpretation.
