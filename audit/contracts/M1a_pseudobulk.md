# Contract M1a: pseudobulk + variance partitioning

**Status:** APPROVED 2026-10-01 (decisions in `audit/decisions.md`, Wave 1 section). **Wave 1**, parallel with M3 and Figs 1–2. Tier: produces results (all 8 steps + planted-bug check).

## Goal
Sum raw counts into pseudobulk samples for M1b's DE, and measure how much variance donor, condition and the technical covariates explain, so the user can decide whether M1b's `~ donor + condition` needs another term.

## Inputs
| File | Check on load |
|---|---|
| `data/processed/m0_cells.h5ad` (24,562 cells × 15,701 genes, raw counts; `obs`: `donor, condition, cell_type, flag_outlier, …`; `uns['control']`, `uns['treatment']`) | Shape as stated; `X` non-negative integers; 8 donors, 8 cell types, 2 conditions |

## Outputs
| File | Contents |
|---|---|
| `src/m1a_pseudobulk.py` | Pseudobulk + covariate table + writers. Every threshold is an argument |
| `R/m1a_varpart.R` | Variance partitioning (reads the files below, writes results) |
| `data/processed/m1a_counts.csv` | Genes × samples, integer counts. Sample id = `cell_type|donor|condition` |
| `data/processed/m1a_samples.csv` | One row per sample: `sample, cell_type, donor, condition, n_cells, total_counts, kept` (`kept` per A1) |
| `data/processed/m1a_counts_no_outliers.csv` + `…_samples_no_outliers.csv` | Same, without the 514 `flag_outlier` cells (input for M1b's R3) |
| `data/processed/m1a_varpart.csv` | Per cell type × gene: variance fraction for each term + residual |
| `notebooks/m1a_pseudobulk.py` | Executed: sample table, PCA per cell type (outlier samples, plotted with distances, no automatic flag, as in M0), variance-fraction violins, covariate vs condition plots |

## Steps (each logged in `audit/cell_flow.md` with an assertion)
1. Sum raw counts per donor × cell type × condition (up to 128 samples). Assert: sample counts sum to the cell-level total for every gene; `n_cells` sums to 24,562.
2. Apply the small-group rule (A1); record which samples are dropped and why.
3. Covariates per sample: `n_cells`, `total_counts`. (10x run = condition, so it can't be a model term: caveat only. % ribo stays out, per F6.)
4. Per cell type, keep genes with baseMean > 1 (the DE gene filter), transform (A3), fit the A2 model per gene.
5. Report. No numeric rule (A4): the user judges from the plots at the gate whether M1b needs a covariate.

## Assumptions (asserted; failure = stop and report)
- Every kept sample has both its donor's conditions (paired design).
- Counts are integers after summing (no normalisation before DESeq2).

## Validation checks owned
I1 (pseudobulk sums to the raw total) · outlier-sample flag (QC list in `CLAUDE.md`) · input to I6.

## Open decisions
None. A1–A4 decided (`audit/decisions.md`).

## Return package
Codex reviewer brief (+ HTML) · `m1a_explained` notebook · shape + head/tail of every output · Codex test results · planted-bug result.
