# M2c test plan

**Status:** APPROVED by the user (2026-10-02), T9 kept. Claude drafted 5; Codex proposed 12 blind from the contract and critiqued the draft (drop the real-data test; don't require a random set to be non-significant). Merged to 9; T9 kept as a softened tripwire for the user to decide.

| # | Test | Bug it catches | Expected answer from | From |
|---|---|---|---|---|
| T1 | **Top-20 ordering:** 25 up DE genes with known q (ties broken by larger \|log2FC\|) and 4 down → exactly 20 up rows ranked 1–20 and 4 down rows, no padding | Wrong sort key or direction; off-by-one; padding | G5 | Both |
| T2 | **Only DE genes, and sharing joins:** non-DE genes with extreme q / log2FC never enter a top list; toy lists (A up in cell types 1 and 2, down in 3; B up in 1) → the sharing table records exactly those memberships and directions | Ranking all tested genes; join inflation; lost direction | G5; construction | Codex |
| T3 | **Jaccard:** {a,b,c} vs {b,c,d} → 0.5; identical → 1; disjoint → 0; symmetric, diagonal 1, up and down separate; 8 cell types → 28 pairs | Wrong formula; mixing directions; dropped cell type | Set arithmetic | Both |
| T4 | **Ranking input:** baseMean 0.99 / 1.00 / 1.01 → only 1.01 ranked; Wald stat = log2FC / lfcSE (2 / 0.5 → 4); NaN stat, zero lfcSE or a duplicated gene → error naming cell type and gene | Wrong filter; ranking by log2FC or q; silent drops | Contract; arithmetic | Codex |
| T5 | **GSEA sign oracle (through R):** set A's genes all at the top, set B's all at the bottom → NES(A) > 0, NES(B) < 0, both padj < 0.05; same seed → identical output | Reversed NES; wrong ranking column; seed ignored | Construction | Both |
| T6 | **Gene sets:** sets with 14 / 15 / 500 / 501 matched genes → only 15 and 500 tested; near-miss symbols (case, punctuation) reported absent, never matched; one changed byte in a GMT → checksum error before GSEA | Size-boundary error; fuzzy matching; wrong file | Contract (15–500, exact match, sha256) | Codex |
| T7 | **BH scope:** hand-computable p-values in two cell types → BH within each cell type over all its sets; padj 0.0499 affected, 0.05 not | Global or per-collection BH; `≤` vs `<` | Contract G4; arithmetic | Codex |
| T8 | **Program labels** *(+ user D2: a program untested in a cell type shows `n/t`, not `–`, and counts as not affected)*: 8 × same sign significant → `shared`; 2–7 → `partial`; 1 → `specific`; significant + and − → `mixed`; none significant → not in the output | Label logic; counting non-significant signs | G4 | Both |
| T9 | **Real-data tripwire:** Hallmark IFN-α response has positive NES in all 8 cell types and padj < 0.05 in all 6 well-powered ones. *(The top-20-contains-a-canonical-ISG half failed in CD14+ Mono and CD4 T and was removed by the user, 2026-10-02; see decisions.)* | Pipeline broken on real data | M3 (IFN-α score up in all 8), P1 (all 10 ISGs up in 7 cell types, 8 in megakaryocytes) | Claude; Codex would drop it |

## Interface the tests call

```python
# src/m2c_programs.py
top_genes(calls, n=20) -> pd.DataFrame          # cell_type, direction, rank, gene, log2FC, qvalue
jaccard_matrix(calls, direction) -> pd.DataFrame  # cell type × cell type
label_programs(gsea, padj_threshold=0.05) -> pd.DataFrame  # one row per program affected anywhere: 8 columns ↑/↓/–, label
run_m2c(in_dir, out_dir, seed=0) -> None
```

```
# R/m2c_fgsea.R
Rscript R/m2c_fgsea.R <ranks.csv: cell_type,gene,stat> <gmt> [<gmt> ...] <out.csv> <seed>
    # fgsea(pathways, stats, minSize=15, maxSize=500, eps=0) per cell type; BH within cell type.
```

## Added after Codex review 1 (user-approved)

| # | Test | Bug it catches | Expected answer from |
|---|---|---|---|
| T10 | `label_programs` with GSEA rows for 7 cell types only (the 8th has no testable sets): a program significant in those 7 in one direction → `partial`, and the 8th shows `n/t` | Cell-type universe taken from the GSEA rows (review #1) | G4: labels across all 8 cell types; D2 |

## Added after the planted-bug check (M2c)

| # | Test | Bug it catches | Expected answer from |
|---|---|---|---|
| T11 | Through `R/m2c_fgsea.R`, one cell type, two toy GMT collections: `padj` equals R's `p.adjust(pval, "BH")` computed by the test over **all** returned sets of that cell type together, and differs from per-collection BH on this toy | BH run per collection (planted bug 4, missed by T1–T10; T7 supplied precomputed padj) | Contract: BH over all sets within the cell type |
