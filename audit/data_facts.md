# Data facts (written by notebooks/00_eda.ipynb)

Read this instead of guessing. Source files and fingerprints are in `CLAUDE.md`.

- h5ad: 24673 cells × 15706 genes; `X` raw integer counts (sparse, float32 storage); no `.raw`, no layers.
- Gene names: symbols in `var_names`, unique. **No `MT-` genes** (so no % mito). 98 genes with the RPS/RPL prefix; 84 counted as ribosomal after excluding 14 non-ribosomal genes (decision D4).
- Key ISGs present: ISG15, IFI6, IFIT1, IFIT2, IFIT3, MX1, OAS1, OAS2, OAS3, RSAD2, IRF7, STAT1. Missing: none.
- Conditions: `ctrl` (control) / `stim` (treatment); stored as `uns['control']` / `uns['treatment']` in the M0 output.
- Donors: patient_1039, patient_107, patient_101, patient_1016, patient_1244, patient_1256, patient_1488, patient_1015 (column `donor`, from h5ad `replicate`); GEO `ind` holds the bare number.
- Cell types (authors', GEO `cell`): B cells, CD14+ Monocytes, CD4 T cells, CD8 T cells, Dendritic cells, FCGR3A+ Monocytes, Megakaryocytes, NK cells.
- Join: key (barcode without suffix, condition); 24673/24673 matched; 234 sequences occur in both conditions among these singlets (313 across all GEO droplets).
- The h5ad's own `cell_type` (re-analysis) agrees with the authors' on 91.0% of cells; not used.
- Basic filters: cells 24673 → 24562; genes 15706 → 15701.
- M0 output: `data/processed/m0_cells.h5ad`, 24562 cells × 15701 genes.

## Cells per donor × condition

| donor        |   ctrl |   stim |   total |
|:-------------|-------:|-------:|--------:|
| patient_1039 |    408 |    584 |     992 |
| patient_107  |    517 |    525 |    1042 |
| patient_101  |    858 |   1166 |    2024 |
| patient_1016 |   1723 |   1635 |    3358 |
| patient_1244 |   1879 |   1464 |    3343 |
| patient_1256 |   2108 |   2026 |    4134 |
| patient_1488 |   2030 |   2549 |    4579 |
| patient_1015 |   2738 |   2352 |    5090 |

## Cells per cell type × donor × condition

| cell_type         |   ('patient_101', 'ctrl') |   ('patient_101', 'stim') |   ('patient_107', 'ctrl') |   ('patient_107', 'stim') |   ('patient_1015', 'ctrl') |   ('patient_1015', 'stim') |   ('patient_1016', 'ctrl') |   ('patient_1016', 'stim') |   ('patient_1039', 'ctrl') |   ('patient_1039', 'stim') |   ('patient_1244', 'ctrl') |   ('patient_1244', 'stim') |   ('patient_1256', 'ctrl') |   ('patient_1256', 'stim') |   ('patient_1488', 'ctrl') |   ('patient_1488', 'stim') |
|:------------------|--------------------------:|--------------------------:|--------------------------:|--------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|
| B cells           |                       107 |                       135 |                        44 |                        54 |                        453 |                        341 |                        126 |                        116 |                         26 |                         39 |                        116 |                         93 |                        216 |                        208 |                        218 |                        279 |
| CD14+ Monocytes   |                       202 |                       253 |                       209 |                       169 |                        779 |                        671 |                        369 |                        329 |                        112 |                        150 |                        415 |                        290 |                        390 |                        375 |                        294 |                        348 |
| CD4 T cells       |                       293 |                       382 |                       158 |                       186 |                        848 |                        734 |                        422 |                        353 |                        182 |                        272 |                       1068 |                        854 |                       1001 |                        954 |                       1182 |                       1484 |
| CD8 T cells       |                        74 |                       105 |                        20 |                        15 |                        179 |                        133 |                        536 |                        489 |                         26 |                         29 |                         58 |                         32 |                        131 |                        112 |                         58 |                         43 |
| Dendritic cells   |                        21 |                        34 |                         8 |                        14 |                         30 |                         30 |                         19 |                         15 |                          9 |                         17 |                         47 |                         33 |                         26 |                         32 |                         43 |                         54 |
| FCGR3A+ Monocytes |                        80 |                       144 |                        32 |                        37 |                        252 |                        228 |                        115 |                        128 |                         31 |                         43 |                         56 |                         35 |                         69 |                         97 |                        101 |                        144 |
| Megakaryocytes    |                        11 |                         9 |                         6 |                         4 |                         20 |                         20 |                          8 |                         10 |                          4 |                         11 |                         14 |                         13 |                         14 |                         17 |                         18 |                         29 |
| NK cells          |                        70 |                       104 |                        40 |                        46 |                        177 |                        195 |                        128 |                        195 |                         18 |                         23 |                        105 |                        114 |                        261 |                        231 |                        116 |                        168 |
