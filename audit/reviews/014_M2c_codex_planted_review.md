# Codex review of the planted-bug snapshot: M2c (not told about the bug)

| # | Severity | file:line | Finding | Evidence |
|---|---|---|---|---|
| 1 | High | `R/m2c_fgsea.R:55` | BH correction is performed separately per collection, not across all Hallmark and KEGG sets within each cell type. | `ave(res$pval, collection_of[res$pathway], ...)` partitions p-values by collection. Decision D2 requires one BH adjustment over both collections per cell type. This can change `padj`, affected calls, and program labels. |

Verdict: M2c is not correct as written; all reviewed rules align except the multiple-testing scope.