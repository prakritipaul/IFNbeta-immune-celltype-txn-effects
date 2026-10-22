# Figure captions

Methods-style captions: what is plotted, from which input, by which method. No interpretation (results are in `FINDINGS.md`). Data: Kang et al. 2018 PBMCs, 8 donors, control (ctrl) vs IFN-β-stimulated (stim); 24,562 cells × 15,701 genes after filtering (cells with ≥ 200 genes detected; genes detected in ≥ 3 cells); cell-type labels from the authors (GEO GSE96583). **DE** throughout = pseudobulk DESeq2 (`~ donor + condition`, per cell type), Wald test, unshrunk log2 fold change (stim vs ctrl), q-values by `qvalue` within each cell type, called at q < 0.05 and |log2FC| ≥ 1. Pseudobulk samples with < 10 cells in either condition were excluded with their pair, leaving 8 paired donors per cell type except dendritic cells (6) and megakaryocytes (4).

---

**`fig1_composition`. Cell numbers per donor, condition and cell type.**
Stacked bars of cell counts (top) and cell-type proportions (bottom) for each of the 8 donors, ctrl and stim side by side, coloured by cell type. Input: 24,562 cells after filtering. Donor order and colours follow Kang 2018 Fig 3d.

**`fig2_umap`. UMAP of all cells, coloured by condition (left) and cell type (right).**
Counts were normalised to 10,000 per cell and log-transformed (natural log, pseudocount 1); 2,000 highly variable genes (Seurat method) were scaled (maximum 10) and reduced by PCA (50 components, ARPACK). A k-nearest-neighbour graph (k = 15, 30 PCs, Euclidean distance) was embedded with UMAP (min_dist 0.5, spread 1.0, spectral initialisation, seed 0). No batch or condition integration was applied. Both panels show the same 24,562 cells at identical coordinates.

**`m1b_concordance`. DE genes per cell type compared with Kang 2018 Supplementary Table 1.**
For each cell type, genes DE in our analysis only (blue), in both lists (purple) and in the paper's list only (red), as counts (left) and as percentages of the union (right). Overlap is by gene symbol regardless of direction; 10 paper symbols were mapped to our names (5 Excel-converted dates restored by Ensembl ID; 5 de-duplicated `.1` names). Paper list: 3,055 unique genes across 8 cell-type sheets.

**`m2b_p1_isgs`. Canonical interferon-stimulated genes called DE per cell type (check P1).**
Left: number of the 10 canonical ISGs (ISG15, IFI6, IFIT1, IFIT2, IFIT3, MX1, OAS1, RSAD2, IRF7, STAT1) called DE in each cell type. Right: distribution of that number across the 8 cell types (proportion of cell types).

**`m2b_p4_null`. Hallmark IFN-α genes among up-regulated genes vs an expression-matched null (check P4).**
For each cell type, grey histogram: number of MSigDB Hallmark IFN-α response genes (2026.1; 6 symbols mapped to the dataset's names) in each of 1,000 random gene sets, each built by replacing every up-regulated DE gene with one of its 30 nearest genes on rank-percentiles of log mean normalised control expression and control detection rate (fraction of the cell type's control cells with count > 0), drawn from all tested genes (baseMean > 1), excluding the gene itself, without reuse within a set (seed 0). Coloured line: observed number of Hallmark IFN-α genes among the up-regulated DE genes. Empirical p = (1 + number of null sets ≥ observed) / 1,001.

**`m2b_n1_shuffles`. DE genes with shuffled vs real condition labels (check N1).**
Number of DE genes + 1 (log scale) per cell type. Grey points: label shuffles, in which the ctrl and stim labels were swapped within exactly half of a cell type's paired donors, every distinct split once (35 shuffles for 8 donors, 10 for dendritic cells, 3 for megakaryocytes), each refitted with the same DESeq2 and qvalue pipeline; points are spread horizontally. Diamonds: real labels.

**`m2b_r2_r3`. Agreement of DE calls between the main analysis and two sensitivity analyses (checks R2, R3).**
Left: Jaccard index between the sets of DE genes in the main analysis and in (blue) the unpaired design (`~ condition`) or (red) pseudobulk recomputed after removing 514 outlier cells (> 5 median absolute deviations from their donor × condition × cell type median in log total counts or log genes detected). Right: percentage of the main analysis' DE calls also called DE in each sensitivity analysis.

**`m2b_v9b_scatter`. log2 fold changes reported by Kang 2018 vs ours (check V9b).**
For each cell type, every gene DE in Supplementary Table 1 that was tested in our analysis (n = 105–2,135 per cell type): the paper's log2FC (x) against ours (y), unshrunk maximum-likelihood estimates (top row) or shrunken by DESeq2 `lfcShrink(type = "normal")` (bottom row). Coloured: DE in both analyses; grey: DE in the paper only. Dashed line: identity.

**`m2b_f3b_heatmap`. Recreation of Kang 2018 Fig 3b.**
Rows: the 2,801 genes DE in at least one cell type. Columns: 144 pseudobulk samples, grouped as all PBMCs (sum over all cells per donor × condition) and the 8 cell types, each group with 8 ctrl then 8 stim donors (the paper shows stim first), including the 12 samples excluded from DE by the 10-cell rule. Values: log2(CPM + 1), z-scored across samples within each row (−3 to 3). Rows are ordered by the authors' 7 gene modules (`final.gene.gct.txt`, column `7 Clusters`; left colour bar and legend), labelled as in that file: 0 (6 genes), 1 (487), 2 (378), 3 (357), 4 (274), 5 (256), 6 (605); label 3 corresponds to Supplementary Table 2 Cluster III (matched on size), other labels have no established correspondence to the paper's numbering. The 438 genes absent from the authors' file are placed at the bottom (no module).

**`m2b_i6`. Per-donor effect size vs pseudobulk sample size (check I6).**
Each point is one paired donor × cell type (58 points). y: median over the cell type's DE genes of the direction-aligned log2((stim CPM + 1) / (ctrl CPM + 1)) (sign multiplied by the direction of the DE call). x: log of the mean number of cells (left) or mean total counts (right) of the donor's ctrl and stim samples. Colour: cell type.

**`fig3_ifn_signature`. IFN signature scores in control and stimulated cells, per cell type.**
Per-cell scores from `scanpy.tl.score_genes` (counts normalised to 10,000 per cell, log1p; 50 control genes per scored gene from 25 expression bins; seed 0) for (A) the MSigDB 2026.1 Hallmark IFN-α response set (97 genes), (B1) the full Hallmark IFN-γ response set (200 genes, 73 of them shared with IFN-α) and (B2) the IFN-γ set without the 73 shared genes (127 genes; specificity control). Violins: kernel density (Scott bandwidth) of the scores of ctrl (blue) and stim (red) cells of the kept samples (24,451 cells; the 111 cells of donors excluded from a cell type by the 10-cell rule are not shown). Dots: each donor's mean score, joined ctrl→stim; 8 donors per cell type, dendritic cells 6, megakaryocytes 4. y-axis shared within a row.

**`m0_paper_fig3a`. Recreation of Kang 2018 Fig 3a: the authors' t-SNE.**
Each dot is one cell (24,562 cells after filtering), placed at the t-SNE coordinates from the authors' GEO metadata (`tsne1`, `tsne2`; not recomputed). Left: coloured by condition (ctrl blue, stim red). Right: coloured by the authors' cell-type label. Cells drawn in a random order (seed 0) so neither condition hides the other. Exported unchanged from the executed `notebooks/00_eda.ipynb` (section 10).

**`m0_paper_fig3d`. Recreation of Kang 2018 Fig 3d: cell-type proportions per donor.**
Stacked bars of the fraction of each donor's cells in each authors' cell type, for ctrl (left) and stim (right); 8 donors in the paper's order, 24,562 cells after filtering. Exported unchanged from the executed `notebooks/00_eda.ipynb` (section 9).

**`fig4_de_heatmap`. Top shared and cell-type-specific IFN-β-responsive genes.**
Rows (61 genes): of the 108 genes labelled shared or shared\* (labels assigned per cell type from the DE calls; see M2a), the 25 with the largest mean |log2FC| over the 8 cell types (a ranking for display only), then, per cell type, the top 5 single-cell-type genes (`specific` or `specific*`) by |log2FC| in that cell type among genes tested (baseMean > 1) in all 8 cell types (fewer where fewer exist: CD8 T cells 4, megakaryocytes 2). \* after a gene name = shared\* or specific\* (label relies on a 95% CI). Colour: unshrunk DESeq2 log2FC (stim vs ctrl), diverging, centred at 0, on two scales with separate colour bars: shared block ±11 (symmetric about the maximum |log2FC|), specific block ±3, where values beyond ±3 saturate (1 cell: SOCS3 in CD8 T cells, +3.3). Dots: DE (q < 0.05, |log2FC| ≥ 1). All shown genes are tested (baseMean > 1) in every cell type. Right bar: `Shared response?` label. Columns ordered by average-linkage hierarchical clustering of the Jaccard index between the cell types' up-regulated DE gene sets (distance 1 − J; dendrogram above).

**`fig4_readme_shared` and `fig4_readme_specific`. Fig 4 split into its two blocks, for the README.**
Pixel crops of `fig4_de_heatmap.png`, unchanged otherwise. `fig4_readme_shared`: the title, dendrogram, shared block (top 25 shared/shared\* genes, colour bar ±11) and legend, with the column labels strip of the full figure appended underneath. `fig4_readme_specific`: the specific block (top single-cell-type genes per cell type, colour bar ±3) with its column labels. Data, rows, colours and column order as in `fig4_de_heatmap`.

**`fig5_programs_heatmap`. Response programs across cell types (GSEA).**
Rows: the 58 Hallmark (H) and KEGG legacy (K) gene sets with BH-adjusted p < 0.05 in at least one cell type, grouped by label (shared: affected in all 8, same sign; partial: 2–7, same sign; specific: 1; mixed: opposite signs), within a group by number of affected cell types, then mean NES. Colour: normalised enrichment score (fgsea multilevel, tested genes ranked by DESeq2 Wald statistic), diverging, centred at 0, limits ±4; dots: BH-adjusted p < 0.05 within cell type; grey hatched n/t: not tested (set size outside 15–500 in that cell type). Columns ordered by average-linkage hierarchical clustering of the cell types' NES profiles (distance 1 − Pearson r) over the 46 affected sets tested in all 8 cell types (dendrogram above).
