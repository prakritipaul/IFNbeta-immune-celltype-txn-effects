# Contract Figs 1–2: composition and UMAP

**Status:** APPROVED 2026-10-01 (decisions in `audit/decisions.md`, Wave 1 section). **Wave 1**, parallel with M1a and M3. Tier: displays results (steps 1, 4, 7 explained-by-example only, 8).

## Goal
Fig 1: cells by condition, donor and cell type. Fig 2: one UMAP, coloured by condition and by cell type.

## Inputs
`data/processed/m0_cells.h5ad` (all 24,562 cells, outliers included per the keep-all decision).

## Outputs
| File | Contents |
|---|---|
| `src/figures.py` | Plot functions + the shared colours/orders (`CT_COLORS`, `COND_COLORS`, `DONOR_ORDER`, `STACK_ORDER`, moved from `00_eda.py`, which then imports them) |
| `notebooks/fig1_fig2.py` | Executed: builds both figures |
| `figures/fig1_composition.{png,pdf}` | Per donor, ctrl and stim bars side by side, stacked by cell type; counts and proportions panels. Kang colours and donor order |
| `figures/fig2_umap.{png,pdf}` | Two panels, same embedding: condition; cell type |
| `data/processed/f2_umap.csv` | `barcode, condition, umap1, umap2` (lets I3 check both panels use the same cells and coordinates) |

## Steps
1. Fig 1 from `obs` counts. Assert: bar totals = cells per donor × condition in `data_facts.md`.
2. Fig 2: normalise → HVGs → PCA → neighbours → UMAP with the C2 recipe, integration per C1. Assert: both panels plot the same 24,562 cells and coordinates.

## Validation checks owned
I3 (both UMAP panels: same cells, same coordinates).

## Open decisions
None. C1–C2 decided (`audit/decisions.md`).

## Return package
Explained-by-example notebook (Claude subagent) · the figures · shape + head/tail of `f2_umap.csv`.
