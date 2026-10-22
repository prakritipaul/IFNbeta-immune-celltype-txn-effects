# Codex review of the planted-bug snapshot: M1b (not told about the bug)

| # | Severity | file:line | Finding | Evidence / reproducer |
|---|---|---|---|---|
| 1 | High | `R/m1b_deseq2.R:86–88` | `qvalue` is run on every non-NA p-value, including genes with `baseMean <= 1`, then those genes are masked afterward. O2/O3 require running `qvalue` only on tested genes (`baseMean > 1`) with non-NA p-values. This can change q-values and DE calls. | Synthetic R reproducer changed q-values substantially—for `p=0.0001`, correct tested-only q was `0.0193`, versus `0.0400` with the current universe. The present main real-data output happens to be numerically unaffected (all eight cell types had identical tested-gene q-values within floating-point tolerance), but that does not guarantee reruns or other inputs are unaffected. |

Verdict: M1b’s model, contrast, thresholds, CI, shrinkage, Cook’s handling, LODO, and shuffles conform, but the q-value testing universe is implemented incorrectly and should be fixed before approval.