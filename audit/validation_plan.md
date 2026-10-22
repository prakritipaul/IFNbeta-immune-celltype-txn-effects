# Validation plan

Agreed with the user on 2026-09-29, before any results exist. Source of every decision: `audit/decisions.md`. **Never change a criterion after seeing results without asking the user.** Items marked **OPEN** still need a user decision.

## Shared definitions

- **Cells:** the 24,673 cells in `data/kang_2018.h5ad` (all are authors' singlets), minus whatever QC removes after the EDA.
- **Cell types:** the authors' per-cell labels from GEO (`cell`), 8 types: CD4 T, CD8 T, B, NK, CD14+ monocytes, FCGR3A+ monocytes, dendritic cells, megakaryocytes.
- **Main vs exploratory cell types:** **OPEN (D3)**, decided after the EDA from donor × condition × cell type counts.
- **DE (Kang 2018 rule):** per cell type, pseudobulk raw counts per donor × condition, R DESeq2 `~ donor + condition`, genes with baseMean > 1, unshrunk log2FC, FDR by R `qvalue`. **Significant** = q < 0.05. **DE** = q < 0.05 and |log2FC| ≥ 1, up or down.
- **IFN gene set:** MSigDB `HALLMARK_INTERFERON_ALPHA_RESPONSE` (version frozen at first use and recorded).
- **Canonical ISGs:** ISG15, IFI6, IFIT1, IFIT2, IFIT3, MX1, OAS1, RSAD2, IRF7, STAT1.
- **Specificity (Kang 2018 Fig 3b):** a gene's pattern = the cell types where it is DE. A "specific" call is flagged **specific\*** if, in any cell type where the gene is not DE, the 95% CI for log2FC reaches |log2FC| ≥ 1. The final table uses ↑ / ↓ / –.

**Reporting style:** numbers and plots, not verdicts. The only pass mark is V9 (✓ if > 70%). The user may revise any criterion later; log changes in `audit/decisions.md`.

## 1. Positive controls: do we find the IFN response?

| ID | Check | How | Reported as | Expected answer from |
|---|---|---|---|---|
| P1 | Canonical ISGs DE per cell type | Count of the 10 ISGs that are DE | Bar chart (ISGs DE per cell type) + histogram (x = ISGs DE 0–10, y = proportion of cell types). No pass mark | IFN biology; Kang 2018 |
| P2 | IFN enrichment among up-DEGs (Q1) | Fisher test: up-DEGs vs all genes tested in that cell type, for the Hallmark IFN-α set | q and odds ratio per cell type | IFN biology |
| P3 | Up vs down | Fisher test: are up-DEGs more IFN-enriched than down-DEGs? | Odds ratio and p, with up/down list sizes | IFN genes are induced, not repressed |
| P4 | Expression-matched null (Q2) | Per DEG, ~30 nearest genes on control mean + detection rate, drawn from **all** tested genes, focal gene excluded; one per DEG, without replacement; 1,000 draws; same enrichment on each | Empirical p per cell type (number, no marker), plus diagnostics: % null genes that are DEGs, % null genes that are IFN, % DEGs that are IFN | Construction (expression bias). If it fails: caveat "enrichment not beyond expression-matched null"; nothing else changes |
| P5 | IFN-α signature rises (Fig 3) | Score per cell, summarised per donor × cell type × condition; gene-set coverage reported | Coverage %; median donor-paired stim − ctrl per cell type. Numbers, no pass marks | Type I IFN biology |
| P6 | IFN-γ specificity control | Same score with Hallmark IFN-γ; repeat with genes shared by both sets removed | IFN-α vs IFN-γ increase per cell type. Numbers, no pass marks | Type I vs II biology |

## 2. Negative control: do we find things that aren't there?

| ID | Check | How | Reported as | Expected answer from |
|---|---|---|---|---|
| N1 | Balanced label shuffles (D5) | All 35 distinct ways of swapping ctrl/stim labels for exactly 4 of 8 donors; rerun DE on each. Same rule for the smaller cell types: half their donors swapped, 10 shuffles for Dendritic cells (6 donors), 3 for Megakaryocytes (4) (O4) | Per cell type: 35 shuffle DE counts plotted next to the real count. No pass mark. Hundreds of shuffle DEGs = stop and investigate | Null by construction |

## 3. Robustness

| ID | Check | How | Reported as | Expected answer from |
|---|---|---|---|---|
| R1 | Leave-one-donor-out | 8 DE refits, each dropping one donor | Whether the genes in the final table keep their direction in all 8 | Resampling |
| R3 | Outlier sensitivity | DE rerun without the 514 flagged cells | Overlap of DE calls with the main run; genes whose call changes | Tests whether possible same-donor doublets or low-quality cells drive any call |
| R2 | Unpaired sensitivity | DE with `~ condition` | Overlap of DE calls with the paired run | Design sensitivity |

## 4. Reproducing Kang 2018

| ID | Check | How | Reported as | Expected answer from |
|---|---|---|---|---|
| V9 | DE genes vs Supp Table 1 | Per cell type: % of the paper's DEGs that we also call DE; correlation of log2FC | Overlap % **with ✓ if > 70%**; correlation as a number | `data/kang2018_supp/MOESM4.xlsx` (3,055 genes) |
| V9b | Unshrunk LFC confirmed | Paper's log2FC vs our unshrunk and shrunken estimates | Which matches better | Supp Table 1 |
| F3a | Recreate Fig 3a | Authors' t-SNE coordinates, ctrl blue / stim red, and cell types | Side by side with published | Kang 2018 Fig 3a |
| F3b | Recreate Fig 3b | Pseudobulk log2 CPM per donor × cell type × condition (9 groups incl. all PBMCs), genes DE in ≥ 1 cell type, rows z-scored, ordered by the authors' modules | Side by side with published | Kang 2018 Fig 3b; `data/kang2018_code/final.gene.gct.txt` |
| F3d | Recreate Fig 3d | Cell-type proportions per donor, ctrl vs stim | Side by side with published | Kang 2018 Fig 3d |
| M1 | Modules reproduce | Authors' 7-cluster labels vs our k-means (k = 7) on our DEGs' per-donor pseudobulk expression. Labelled "method inferred" | Overlap table + adjusted Rand index | Authors' GCT file |
| M2 | Module biology | KEGG enrichment (Fisher, KEGG gene sets) of our modules | Compared with Supp Table 2, at least Cluster III (401 genes, Influenza A / Measles) | `data/kang2018_supp/MOESM5.xlsx` |

## 4b. Shared programs (M2c, added 2026-10-02)

| ID | Check | How | Reported as |
|---|---|---|---|
| G-top | Top genes shared | Top 20 up and top 20 down per cell type; which cell types share each; Jaccard clustering of full DE lists | Table + heatmap |
| G-sea | Shared programs | fgsea per cell type (Hallmark + KEGG legacy, Wald-stat ranking), BH < 0.05; programs labelled shared / partial / specific / mixed; cell types clustered by NES | Table + NES heatmap |

## 5. Data and code integrity

| ID | Check | How | Reported as |
|---|---|---|---|
| I1 | Cell accounting (incl. old V1) | Every filter, join and aggregation logs n in → n out, and asserts what must hold (e.g. joins count **matches**, not rows; pseudobulk counts sum to the raw total). Hand-summed pseudobulk spot check is a Codex test | `audit/cell_flow.md` |
| I2 | Heatmap/table traceability (V11) | Every table cell and heatmap tile matches its DESeq2 coefficient and the one DE rule | 100% agreement (automatic) |
| I3 | Reading the data correctly (V13) | Both UMAP panels use the same cells and coordinates; our plot of the authors' t-SNE matches their Fig 3a | Automatic check + visual |
| I4 | Simulated-data oracle (V14) | Fake counts with planted effects (known genes, known cell types) through the full Python → R pipeline; planted genes found, unplanted not | Codex test plan for the DE module |
| I5 | Ribosomal/mito sanity (V10) | Number of RPL/RPS/MT- genes among the top 50 shared DEGs | A number, no pass mark |
| I6 | Covariates (V15) | Signature scores and effect sizes vs library size and % mito; with vs without covariates **only if M1a adds one** | Correlation plots, no pass mark |

## Dropped (and why)

Fig 3f (test not specified in the paper) · edgeR (little added value for strong IFN signal) · shrinkage on/off (paper used unshrunk) · Q4, enrichment on shuffles (other enrichment checks suffice) · control-vs-control splits (a misunderstanding; overlaps N1) · split-donor replication (4 donors per half too few) · V12 annotation check (authors' labels used) · ΔLFC specificity test (paper's rule + specific\* flag instead).

## Still open

1. **D3:** main vs exploratory cell types, decided after the EDA from donor × condition × cell type counts.
