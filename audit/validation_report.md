# Validation report (M2b)

- **Sources:** notebook `notebooks/m2_validation.ipynb`, logic in `src/m2b_validation.py`, outputs `data/processed/m2_validation_*.csv`, plots `figures/m2b_*.png`.
- **Definitions:** `audit/validation_plan.md`, plus the "M2b checks" row in `audit/decisions.md`.
- **DE rule:** calls re-derived with q < 0.05 and |log2FC| ≥ 1. This changes 0 calls in every run.
- **Pass marks:** only V9 has one.

| ID | Result | Plot / table |
|---|---|---|
| P1 | All 10 canonical ISGs are up in 7 cell types. Megakaryocytes have 8/10 (IRF7 and STAT1 not DE). No ISG goes down | `m2b_p1_isgs.png` |
| P2 | Hallmark IFN-α among up-DEGs: odds ratio 44 (FCGR3A+ Mono) to 174 (CD4 T). One-sided p, BH across 8: q ≤ 1.8e-17 in every cell type (Megakaryocytes 14 of 23 up genes, the least significant) | `p2_p3_fisher` |
| P3 | Up vs down: IFN-α genes among down-DEGs are 0 in 6 cell types and 1 in CD14+ and FCGR3A+ Mono (odds ratio 83 and 45; inf elsewhere). One-sided p ranges from 2e-23 (CD14+) to 0.041 (Megakaryocytes, 4 down genes) | `p2_p3_fisher` |
| P4 | Observed IFN-α genes among up-DEGs exceed all 1,000 expression-matched draws in every cell type: p = 0.001, the minimum. Observed vs null mean: CD4 T 78 vs 4.6; CD14+ 84 vs 15.1; Megakaryocytes 14 vs 3.9. Diagnostics: 1.2–4.0% of null genes are IFN-α (Megakaryocytes 17%); 4.6–32% of null genes are DEGs; 6.9–61% of up-DEGs are IFN-α. CD14+ needed the next-nearest fallback 112 times in 1.2 M picks | `m2b_p4_null.png`, `p4_matched_null` |
| N1 | Shuffles give 0–2 DE genes (at most 2 in B, 1 in CD14+, 0 elsewhere), against 27–2,335 with the real labels | `m2b_n1_shuffles.png` |
| R1 | 100% of DE calls keep their sign in every LODO refit, in all 8 cell types (4,947 calls) | `r1_lodo` |
| R2 | The unpaired run recovers 89% (Megakaryocytes) to 99% (B) of main calls. Jaccard: 0.50 (Megakaryocytes), 0.67 (DC), 0.79–0.93 in the others. 0 sign flips | `m2b_r2_r3.png` |
| R3 | The no-outlier run recovers 85% (Megakaryocytes) to 99.7% (CD14+) of main calls. Jaccard 0.85–0.99. 0 sign flips | `m2b_r2_r3.png` |
| V9 | **Pass in 2 of 8 cell types** (CD14+ 84.8%, CD4 T 79.9%); fail in 6 (24.3–69.5%). Unchanged from M1b; no opposite signs | `v9_overlap`, `m1b_concordance.png` |
| V9 correlation | Every paper DEG we tested (n = 105–2,135): Pearson 0.983–0.994 in 7 cell types, Megakaryocytes 0.92. DE in both lists: 0.970–0.995, Megakaryocytes 0.93 | `v9_correlation` |
| V9b | Correlation is similar for unshrunk and shrunken (r ≥ 0.97 for both, except Megakaryocytes and shrunken CD8 T). The median \|paper − ours\| is smaller for shrunken values in most cell types (CD4 T 0.06 vs 0.14 on DE-in-both genes; Megakaryocytes 0.45 vs 1.28); unshrunk is closer only in CD8 T. **Recorded as favouring shrunken values; the calls stay unshrunk (decision)** | `m2b_v9b_scatter.png` |
| F3b | 2,801 DE genes × 144 columns (9 groups × 8 donors × 2 conditions). 2,363 genes are in the authors' file; 438 (16%) are "no module" | `m2b_f3b_heatmap.png`, `f3b_zscores` |
| Modules | Our k-means vs the authors' clusters: **ARI 0.33** (2,363 genes; "method inferred"). Authors' label 1 maps mostly to our module 0 (356/487); label 2 to our 2 (225/378); label 6 to our 1 (313/605). Label 3 (Cluster III) splits between our module 3 (177) and module 6 (160; 72% of module 6) | `modules_overlap` |
| Module biology | KEGG genes get the 6 B4 renames (3 more genes reach the background). KEGG legacy has no Measles, Influenza A, Herpes simplex or Hepatitis C sets, the top 4 Cluster III terms. Our module 6 (mostly authors' Cluster III) ranks RIG-I-like receptor signalling 1st (q 3e-5) and cytosolic DNA sensing 2nd (q 7e-4); both are Cluster III terms. Apoptosis, also a Cluster III term, is 4th (q 0.033), after proteasome (q 0.010). Our module 3 has no term with q < 0.05 | `kegg_modules` |
| I5 | **0 RPL/RPS and 0 MT- genes** among the top 50 of the 108 shared + shared\* genes (`m2_table.csv`), ranked by mean \|log2FC\| over 8 cell types (3.2–7.9). The top 10 are CXCL10, IFIT1, CCL8, IFIT3, RSAD2, IFIT2, CMPK2, ISG15, TNFSF10, MX1. MT2A, at 25th, is a metallothionein, not a mitochondrial gene | `i5_ribo`, notebook I5 cell |
| I6 | Effect size is the median direction-aligned log2((stim CPM + 1)/(ctrl CPM + 1)) over DE genes: 1.5–2.9 per donor, Megakaryocytes 5.3–7.5. Spearman against log n_cells: −0.96 (B), −0.93 (CD8 T), −0.86 (CD4 T), −0.83 (NK), −0.67 (FCGR3A+), −0.52 (CD14+), −0.26 (DC), +0.20 (Megakaryocytes). Against log total_counts: −0.98 to −0.40. Each uses 4–8 donors. **Smaller samples show larger effects, and this holds with the pseudocount** | `m2b_i6.png`, `i6_spearman` |
