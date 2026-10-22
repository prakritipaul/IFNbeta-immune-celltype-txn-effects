# Contract Wave 4: Fig 3, Fig 4, Fig 5 and the final table

**Status:** APPROVED 2026-10-05 (F1–F6 in `audit/decisions.md`; Fig 5 added by the user). **Wave 4**, one agent. Tier: displays results (steps 1, 4, 7 explained-by-example only, 8). Nothing here computes new results: every number comes from an approved module's output.

## Inputs
| File | Used for |
|---|---|
| `data/processed/m3_cell_scores.csv`, `m3_sample_scores.csv` | Fig 3 |
| `data/processed/m2_table.csv`, `m2_calls.csv` | Fig 4, final table |
| `data/processed/m2c_jaccard.csv` (up-gene Jaccard) | Column order of Fig 4 |
| `data/processed/m2c_gsea.csv`, `m2c_programs.csv` | Fig 5 and its column order |
| `src/figures.py` | Shared colours and orders (CT_COLORS, COND_COLORS) |

## Fig 3: IFN signature, control vs stimulated, per cell type (F1, F2)
Per cell type: violins of per-cell scores, ctrl and stim side by side (COND_COLORS), with each donor's mean overlaid as a dot and joined ctrl→stim by a line (paired). Row A: Hallmark IFN-α score. Row B1: full Hallmark IFN-γ; row B2: IFN-γ without the shared genes (specificity control) (F2 revised 2026-10-22). Megakaryocytes and DCs show their 4 / 6 kept donors.
**Output:** `figures/fig3_ifn_signature.{png,pdf}`.

## Fig 4: heatmap of the top shared and cell-type-specific genes (F3–F5)
Rows: the top 25 shared genes (shared + shared\*, by mean |log2FC|) and the top 5 specific genes per cell type **among genes tested in all 8 cell types** (`specific*`, marked \*, plus the clean `specific` FTL in megakaryocytes), by |log2FC| (F3, revised by the user after seeing the first version); no gene-set filter (F4); title "Top shared and cell-type-specific IFN-β-responsive genes". Columns: 8 cell types, ordered by average-linkage clustering of the up-regulated DE gene sets (Jaccard, from M2c), with the dendrogram drawn above (F5, user). Colour: unshrunk log2FC, diverging blue–white–red centred at 0, symmetric limits; cells that are DE marked with a dot; untested cells grey. Row groups labelled (shared / specific to <cell type>) with a side bar for the `shared_response` label.
**Output:** `figures/fig4_de_heatmap.{png,pdf}`.

## Fig 5: response programs across cell types (added by the user after this contract was drafted, 2026-10-05)
Rows: the 58 programs affected in at least one cell type (M2c), grouped by label (shared / partial / specific / mixed). Columns: 8 cell types, ordered by average-linkage clustering of their NES profiles (1 − Pearson correlation, programs tested in all 8; from M2c), dendrogram above. Colour: NES, diverging, centred at 0; significant cells (BH < 0.05) marked with a dot; untested (`n/t`) cells grey and hatched.
**Output:** `figures/fig5_programs_heatmap.{png,pdf}`.

## Final table (F6)
`results/final_table.csv`: one row per gene DE in ≥ 1 cell type (2,801): `gene`, 8 cell-type columns (↑ / ↓ / –), `Shared response?` (M2a labels: shared, shared*, partial, specific, specific*, mixed), `n_cell_types`, `megakaryocyte_lower_confidence`, and per cell type the log2FC and q (so every arrow can be traced). Sorted: shared → shared* → partial → mixed → specific → specific*, then by n_cell_types and mean |log2FC|. A short `results/README.md` explains the columns and footnotes (shared\*, specific\*, megakaryocyte flag, 10x-run confound).

## Assertions
- Every arrow in the final table equals M2a's `m2_table.csv` (100%); every Fig 4 cell equals the log2FC in `m2_calls.csv`; every Fig 5 cell equals the NES in `m2c_gsea.csv`, and its labels equal `m2c_programs.csv`.
- Fig 3 donor means equal `m3_sample_scores.csv`.

## Open decisions
None. F1–F6 decided.

## Return package
Explained-by-example notebook (Claude subagent), the figures with captions in `figures/CAPTIONS.md` and takeaways in `FINDINGS.md`, head/tail of `final_table.csv`, `FINDINGS.md` section.
