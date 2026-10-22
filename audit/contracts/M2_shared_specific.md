# Contract M2: shared vs cell-type-specific calls, and the validation checks

**Status:** APPROVED 2026-10-02 (O1–O7 in `audit/decisions.md`, M2 section). **Wave 3.** Two parts (decision O6): **M2a** produces the calls and the final table (result tier: all 8 steps + planted-bug check); **M2b** computes the validation checks as numbers and plots (report tier: steps 1, 4, 7, 8).

## Goal
Turn M1b's per-cell-type DE results into the answer to the project question: for every gene DE in at least one cell type, which cell types respond, in which direction, and is the response **shared** or **cell-type specific**. Then report every agreed validation check (`audit/validation_plan.md`) on those calls.

## Inputs
| File | Used for |
|---|---|
| `data/processed/m1b_de.csv` (main run) | Calls, CIs, the table |
| `data/processed/m1b_runs.csv`, `m1b_fits.csv` | R1 (LODO), R2 (unpaired), R3 (no outliers), N1 (shuffles) |
| `data/processed/m1a_counts.csv`, `m1a_samples.csv` | F3b heatmap (log2 CPM), I6 (covariates), P4 (control mean) |
| `data/processed/m0_cells.h5ad` | P4 detection rate (fraction of cells expressing) |
| `data/msigdb/h.all.v2026.1.Hs.symbols.gmt` (+ 6 mappings, B4) | P2–P4 (Hallmark IFN-α) |
| `data/kang2018_supp/MOESM4.xlsx` (Supp Table 1), `MOESM5.xlsx` (Supp Table 2), `data/kang2018_code/final.gene.gct.txt` | V9, V9b, F3b, module comparison (M1/M2 checks) |
| `data/msigdb/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt` (O5) | Module biology |

## M2a: calls and the final table
1. **Per gene × cell type:** `up` / `down` / `none` from M1b (DE = q < 0.05 and |log2FC| ≥ 1; unshrunk).
2. **Pattern** = the set of cell types where the gene is DE (Kang 2018 Fig 3b).
3. **Shared response?** label per gene: `shared` / `shared*` = DE in all 6 well-powered cell types in one direction, and in DCs and megakaryocytes DE in that direction or with the 95% CI reaching |log2FC| ≥ 1 in it (`*` = CI rule needed); `partial` = DE in ≥ 2, not shared; `specific` / `specific*` = DE in 1; `mixed` = up in some, down in others (O1, O2). Megakaryocyte calls flagged lower confidence (O3).
4. **specific\* flag** (decided): a gene called specific gets an asterisk if, in any cell type where it is not DE, the 95% CI for log2FC reaches |log2FC| ≥ 1.
5. **LODO stability (R1):** for each DE call, does log2FC keep its sign in every leave-one-donor-out refit of that cell type? Column `lodo_stable`.
6. **I2 traceability:** every table cell is re-derived from `m1b_de.csv` and the rule, and must agree 100% (assert).

**Outputs:** `data/processed/m2_calls.csv` (gene × cell type, long: `de`, `log2FC`, `qvalue`, `ci_low/high`, `lodo_stable`); `data/processed/m2_table.csv` (one row per gene DE anywhere: 8 cell-type columns ↑ / ↓ / –, `shared_response`, `specific_star`, `n_cell_types`, flags); `notebooks/m2_calls.py`.

## M2b: validation checks (numbers and plots; only V9 has a pass mark)
| ID | What |
|---|---|
| P1 | 10 canonical ISGs DE per cell type: bar chart + histogram |
| P2 / P3 | Fisher enrichment of Hallmark IFN-α among up-DEGs (q, odds ratio); up vs down |
| P4 | Expression-matched null: 1,000 draws of ~30 genes matched on control mean + detection rate, focal gene excluded; empirical p per cell type + diagnostics |
| N1 | Shuffle DE counts (35 / 10 / 3) plotted next to the real count |
| R1 / R2 / R3 | LODO direction stability; overlap of calls with the unpaired and no-outlier runs |
| V9 / V9b | Recovery of Supp Table 1 (already computed in M1b; restated) + log2FC correlation; paper LFC vs our unshrunk and shrunken |
| F3b | Recreate Fig 3b: pseudobulk log2 CPM, 9 groups incl. all PBMCs, genes DE in ≥ 1 cell type, z-scored rows, ordered by the authors' modules |
| Modules | Authors' 7 clusters vs our k-means (k = 7, seed fixed): overlap table + adjusted Rand index ("method inferred") |
| Module biology | KEGG enrichment (Fisher) of our modules vs Supp Table 2, at least Cluster III |
| I5 | RPL/RPS/MT genes among the top 50 shared DE genes |
| I6 | DE effect sizes vs `n_cells` / `total_counts` (correlation plots) |

**Outputs:** `data/processed/m2_validation_*.csv`, `audit/validation_report.md` (every check: result, plot, pass/fail only for V9), `notebooks/m2_validation.py`.

## Assumptions (asserted; failure = stop and report)
- Every gene × cell type in the table exists in `m1b_de.csv` exactly once; calls re-derived from it agree 100%.
- Every expected LODO and shuffle run exists (8/6/4 and 35/10/3).

## Open decisions
None. O1–O7 decided.

## Return package
M2a: Codex reviewer brief (+ HTML), explained-by-example, outputs head/tail, Codex tests, planted bug. M2b: explained-by-example, `validation_report.md`. Both: `FINDINGS.md` section.
