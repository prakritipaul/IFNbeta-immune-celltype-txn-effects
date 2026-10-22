# Contract M2c: shared programs (top genes and GSEA)

**Status:** APPROVED 2026-10-02 (G1–G6; test plan T1–T9 incl. the T9 tripwire). **Wave 3**, alongside M2a/M2b. Tier: results, with a small test set (G6).

## Goal
Go beyond single genes: find **groups of cell types that share a response program**, at two levels:
1. **Top genes:** which cell types share their strongest up- and down-regulated genes.
2. **Programs (GSEA):** which gene sets (Hallmark, KEGG) go up, down, or both across cell types, and which cell types run the same programs.

## Inputs
| File | Used for |
|---|---|
| `data/processed/m1b_de.csv` (main run: `log2FC`, `qvalue`, `de`, Wald `stat` or `log2FC/lfcSE`) | Ranking and top genes |
| `data/msigdb/h.all.v2026.1.Hs.symbols.gmt` (sha256 `eecaf6da…f596`, 6 mappings B4) | Hallmark sets |
| `data/msigdb/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt` (sha256 `4a87c8d1…fa5e`) | KEGG sets |

## Part 1: top genes shared across cell types (G5)
- Per cell type: **top 20 up** and **top 20 down** DE genes, ranked by q, then |log2FC| (fewer if fewer are DE, e.g. megakaryocytes have 4 down).
- For every gene in any top list: the cell types where it is in the top 20, in which direction.
- **Cell-type similarity:** Jaccard index between the full DE gene sets of every pair of cell types (up and down separately), hierarchically clustered → groups of cell types.

## Part 2: programs (GSEA, G1–G4)
- Per cell type: rank **all tested genes** (baseMean > 1) by the DESeq2 Wald statistic (= log2FC / lfcSE); run R `fgsea` on Hallmark (50) + KEGG legacy (186), set size 15–500 (fgsea defaults, passed explicitly), fixed seed.
- A program is **affected** in a cell type if BH-adjusted p < 0.05 (BH over all sets within the cell type); direction = sign of NES.
- **Program labels** across the 8 cell types, same scheme as genes but without the CI rule: `shared` (affected in all 8, same direction) / `partial` (≥ 2, same direction) / `specific` (1) / `mixed` (up in some, down in others).
- **Cell-type groups:** hierarchical clustering of cell types on their NES profiles (correlation distance).

## Outputs
| File | Contents |
|---|---|
| `src/m2c_programs.py`, `R/m2c_fgsea.R` | Logic (tested) |
| `data/processed/m2c_top_genes.csv` | cell_type, direction, rank, gene, log2FC, qvalue; plus a gene × cell type sharing table |
| `data/processed/m2c_jaccard.csv` | Pairwise Jaccard of DE sets (up / down) |
| `data/processed/m2c_gsea.csv` | cell_type, collection, pathway, size, NES, pval, padj, leading edge |
| `data/processed/m2c_programs.csv` | One row per program affected anywhere: per-cell-type ↑ / ↓ / –, label |
| `notebooks/m2c_programs.py` | Top-gene sharing table, Jaccard heatmap with dendrogram, NES heatmap (programs × cell types) with dendrogram, labelled programs |

## Assumptions (asserted)
- Gene-set files match their sha256; every gene-set gene is matched by exact symbol (after B4's mappings) or reported as absent, never fuzzily.
- Every cell type's ranking has no ties broken arbitrarily without a fixed seed; no NaN statistics among tested genes.

## Open decisions
None beyond G1–G6 (decided). Low-stakes details fixed here: set size 15–500, `eps = 0` (exact p-values), seed 0.

## Return package
Explained-by-example notebook, outputs head/tail, Codex tests, `FINDINGS.md` section. No separate Codex brief unless the user asks (G6: small module).
