# notebooks/

Each notebook exists twice: the `.py` file is the source, in jupytext percent format and readable as a diff, and the `.ipynb` file is the executed copy with its outputs. The notebooks hold no analysis logic of their own. They call the module in `src/` (and through it the scripts in `R/`), show the results, and assert that the numbers agree with their sources.

## Analysis notebooks (one per module, in wave order)

| Notebook | Module | Calls | What it shows |
|---|---|---|---|
| `00_eda` | M0 · Wave 0 | `src/m0_qc.py` | Visual EDA: structure, counts, filters, QC flags, donor × condition × cell type sizes, covariates. Writes `audit/data_facts.md` |
| `m1a_pseudobulk` | M1a · Wave 1 | `src/m1a_pseudobulk.py`, `R/m1a_varpart.R`, `R/m1a_cancor.R` | Pseudobulk samples and variance partitioning; why DE uses no technical covariate |
| `m3_signature` | M3 · Wave 1 | `src/m3_signature.py` | Hallmark IFN-α and IFN-γ signature scores per cell and per donor |
| `fig1_fig2` | Figs 1–2 · Wave 1 | `src/figures.py` | Fig 1 (composition) and Fig 2 (UMAP) |
| `m1b_de` | M1b · Wave 2 | `src/m1b_de.py`, `R/m1b_deseq2.R` | Paired pseudobulk DESeq2 per cell type, sensitivity runs, comparison with the paper |
| `m2_calls` | M2a · Wave 3 | `src/m2a_calls.py` | Shared vs cell-type-specific labels per gene; leave-one-donor-out stability |
| `m2_validation` | M2b · Wave 3 | `src/m2b_validation.py` | Validation checks from `audit/validation_plan.md` (ISGs, nulls, shuffles, robustness, agreement with the paper) |
| `m2c_programs` | M2c · Wave 3 | `src/m2c_programs.py`, `R/m2c_fgsea.R` | Shared top genes and GSEA programs (Hallmark + KEGG) across cell types |
| `w4_figures_table` | Wave 4 | `src/w4_outputs.py` | Figs 3–5 and `results/final_table.csv` |

## Explained-by-example notebooks (`*_explained`)

One per module: `m0_`, `m1a_`, `m3_`, `fig12_`, `m1b_`, `m2a_`, `m2b_`, `m2c_`, `w4_explained`. Each takes the lines of code worth checking (those named in the module's Codex reviewer brief) and, for each one:
1. builds a tiny, hand-checkable input.
2. shows the **expected** output, written by a Claude subagent from the contract and the decisions *before* the code was run. That subagent did not write the code.
3. shows what the real function returns, with an `assert` that the two match.

A mismatch is reported as a possible bug, and the expected value is never edited to match. These notebooks are checks, not findings, so they are not summarised in `FINDINGS.md`.
