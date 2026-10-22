# Contract M1b: pseudobulk differential expression

**Status:** APPROVED 2026-10-01 (decisions O1–O7 in `audit/decisions.md`, M1b section). **Wave 2**, one agent. Tier: produces results (all 8 steps + planted-bug check, O7).

## Goal
Per cell type, test every gene for an IFN-β effect with the paper's method, and produce the reruns that the robustness and negative-control checks need. M1b makes **no** shared/specific calls (that is M2).

## Inputs
| File | Check on load |
|---|---|
| `data/processed/m1a_counts.csv` + `m1a_samples.csv` (128 samples, 116 `kept`) | Integer counts; sample ids match between files; kept samples paired (one per condition per donor × cell type) |
| `data/processed/m1a_counts_no_outliers.csv` + `…_samples_no_outliers.csv` | Same checks (input for R3) |

## Method (all decided: `audit/decisions.md`)
Per cell type, kept samples only: R DESeq2 1.40.2, `~ donor + condition`, size factors `type="ratio"`, `fitType="parametric"` with the fallback logged, Wald test, `betaPrior=FALSE` (unshrunk MLE log2FC), Cook's handling at defaults (O1). Fit all genes, then keep baseMean > 1 (O2); FDR by R `qvalue` per cell type on their p-values (O3). Shrunken `lfcShrink(type="normal")` log2FC also written, for V9b only (O5). **Significant** = q < 0.05; **DE** = q < 0.05 and |log2FC| ≥ 1.

## Outputs
| File | Contents |
|---|---|
| `src/m1b_de.py` | Writes the inputs for each run, calls R, validates and collects results |
| `R/m1b_deseq2.R` | One DESeq2 + qvalue fit per cell type for a given sample set and design |
| `data/processed/m1b_de.csv` | Main run. Per cell type × gene: `run` (always `main`), `baseMean, log2FC, lfcSE, ci_low, ci_high` (95% Wald CI, for M2's specific\* flag), `pvalue, qvalue, significant, de` (`up` / `down` / `none`), `cooks_flag`, `log2FC_shrunk` |
| `data/processed/m1b_runs.csv` | Same columns, for every run (below) |
| `data/processed/m1b_fits.csv` | Per run × cell type: dispersion fit actually used, genes tested, samples |
| `notebooks/m1b_de.py` | DE counts per cell type (up/down), volcano plots, top genes (CLAUDE.md sanity check: ISGs should dominate), outputs head/tail |

## Runs (each is the same R script on a different sample set or design)
| Run | What | For check |
|---|---|---|
| `main` | All kept samples, `~ donor + condition` | Everything downstream |
| `no_outliers` | Kept samples rebuilt without the 514 flagged cells | R3 |
| `unpaired` | Main samples, `~ condition` | R2 |
| `lodo_<donor>` | Main samples minus one donor (one run per paired donor: 8, 6 or 4) | R1 |
| `shuffle_<k>` | Within-donor ctrl/stim labels swapped for half of each cell type's paired donors, every distinct split: 35 (8 donors), 10 (Dendritic cells), 3 (Megakaryocytes) (O4) | N1 |

## Assumptions (asserted; failure = stop and report)
- Every gene × cell type with any non-zero count in the main run has a log2FC; a NaN p-value appears only for a documented DESeq2 reason (Cook's outlier or all-zero gene), and is counted in the notebook.
- No run silently drops a cell type: if a fit fails, the script stops and reports.

## Validation checks owned
I4 simulated-data oracle (planted DE genes found, unplanted not; Codex test) · inputs for P1–P4, N1, R1–R3, V9, V9b (computed in M2) · DESeq2 fit types logged.

## Open decisions
None. O1–O7 decided.

## Return package
Codex reviewer brief (+ HTML) · `m1b_explained` notebook · shape + head/tail of every output · Codex test results.
