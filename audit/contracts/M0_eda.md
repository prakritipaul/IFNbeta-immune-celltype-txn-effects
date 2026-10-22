# Contract M0: EDA, join, basic filtering, QC flags

**Status:** COMPLETE: approved by the user 2026-10-01 (contract approved 2026-09-29, O1–O4 as recommended). Test plan: `audit/test_plans/M0.md` (approved). **Wave 0**, one agent (Claude writes, Codex tests). Nothing downstream starts until the user has reviewed M0's outputs.

## Goal
Learn the data before any analysis: see it (head/tail, distributions, outliers), attach the authors' labels, apply the agreed basic filters, and **flag** anything concerning for the user to decide on. M0 drops nothing beyond the agreed basic filters.

## Inputs
| File | What | Check on load |
|---|---|---|
| `data/kang_2018.h5ad` | 24,673 cells × 15,706 genes, raw counts in `X` | sha256 `e6a5adac…9830`; `X` all non-negative integers |
| `data/geo/GSE96583_batch2.total.tsne.df.tsv.gz` | Authors' per-cell metadata, 29,065 rows | 29,065 rows; columns `tsne1 tsne2 ind stim cluster cell multiplets` |

## Outputs
| File | Contents |
|---|---|
| `src/m0_qc.py` | Functions: load, join, QC metrics, basic filters, flags. No hidden defaults: every threshold is an argument |
| `notebooks/00_eda.ipynb` | Executed notebook: everything in the `CLAUDE.md` EDA checklist, with plots |
| `data/processed/m0_cells.h5ad` | Raw counts after basic filters. `obs`: `barcode, condition, donor, cell_type` (authors'), `cluster, tsne1, tsne2, n_counts, n_genes, pct_ribo`, one boolean column per QC flag. `var`: gene symbol, `n_cells` |
| `audit/data_facts.md` | Facts other agents read instead of guessing: column names, label spellings, counts, gene naming, what's missing |
| `audit/qc_flags.md` | Table of flags: what, evidence (number + plot), proposed action, why |
| `audit/cell_flow.md` | n cells / n genes in → out for every step, including join matches |

## Steps (each logged in `cell_flow.md` with an assertion)
1. Load both files; verify the input checks above.
2. **Join** on (barcode without suffix, condition). Assert: all 24,673 cells match; donor, cluster and t-SNE agree 100%; all are singlets.
3. Compute QC metrics from raw counts: total counts, genes detected, % ribosomal (`RPL`/`RPS`). Ignore the re-analysis's `nCount_SCT`, `seurat_clusters` etc.
4. **Basic filters (agreed):** keep cells with ≥ 200 genes, then genes detected in ≥ 3 cells.
5. **Flags only** (no drops), thresholds as proposed in `CLAUDE.md`: outlier cells (> 5 MADs, per donor × condition), groups < 10 cells, cell types with < 3 paired donors, weak donors, composition shifts, covariate trends.
6. Write outputs.

## Assumptions (each asserted in code; failure = stop and report, never work around)
- Both conditions present for all 8 donors.
- Gene names are unique symbols (`var_names`), no duplicates.
- The authors' labels have exactly 8 cell types, matching the names in `CLAUDE.md`.

## Validation checks owned by M0
I1 cell accounting (join counts matches, not rows) · I3 reading the data correctly (authors' t-SNE plotted vs their Fig 3a) · Fig 3d recreation (proportions per donor) · EDA inputs for D3.

## Open decisions (user)
| # | Decision | Options | Recommendation |
|---|---|---|---|
| O1 | **No `MT-` genes in the file**, so % mito can't be computed | (a) Drop % mito from QC and covariates, and say so. (b) Rebuild counts from GEO's raw matrices (`GSE96583_RAW.tar`), which may include `MT-` genes: a bigger job and a new join | **(a)**. Cells were already QC'd by the authors (singlets only). % ribo and library size remain as covariates. Note it as a caveat |
| O2 | Gene filter "≥ 3 cells" counted across | (a) all cells. (b) each condition separately | **(a)** all cells: the standard, and it doesn't bias towards genes present in one condition |
| O3 | QC flags computed | (a) after the basic filters. (b) before | **(a)** after: flags then describe the cells we'll actually use |
| O4 | Code layout | (a) logic in `src/m0_qc.py`, notebook calls it. (b) everything in the notebook | **(a)**: Codex can test functions; the notebook stays readable |

## Return package (at the review gate)
Reviewer brief (non-author agent) · `data_facts.md`, `qc_flags.md`, `cell_flow.md` · shape + head/tail of `m0_cells.h5ad` `obs`/`var` · the notebook · Codex test results.
