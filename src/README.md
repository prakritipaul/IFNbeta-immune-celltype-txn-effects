# src/

The Python logic, one file per module, written by Claude against a contract the user approved first and tests Codex wrote blind. Notebooks in `notebooks/` call these files.

Every parameter that can change a result is an explicit argument, with no hidden defaults. Invalid input raises `ValueError` before anything is written.

| File | Module | What it does | Contract | Tests |
|---|---|---|---|---|
| `m0_qc.py` | M0 | Loads the h5ad, joins the authors' GEO labels on (barcode, condition), computes QC metrics, applies the basic filters and sets QC flags. Writes `data/processed/m0_cells.h5ad` | `M0_eda.md` | `test_m0*.py` |
| `m1a_pseudobulk.py` | M1a | Sums raw counts per cell type × donor × condition, applies the 10-cell paired rule and writes the inputs for `R/m1a_varpart.R` | `M1a_pseudobulk.md` | `test_m1a*.py` |
| `m3_signature.py` | M3 | Scores Hallmark IFN-α and IFN-γ per cell, then summarises per donor × cell type × condition | `M3_signature.md` | `test_m3*.py` |
| `figures.py` | Figs 1–2 | Composition bars and the UMAP, plus the shared colours and orders | `F12_figures.md` | notebook asserts |
| `m1b_de.py` | M1b | Validates inputs, builds each DE run (main, sensitivity, shuffles, leave-one-donor-out), calls `R/m1b_deseq2.R` and collects the results | `M1b_de.md` | `test_m1b*.py` |
| `m2a_calls.py` | M2a | Labels each gene shared / shared\* / partial / specific / specific\* / mixed, and computes leave-one-donor-out sign stability | `M2_shared_specific.md` | `test_m2a*.py` |
| `m2b_validation.py` | M2b | The validation checks defined in `audit/validation_plan.md` | `M2_shared_specific.md` | report tier: numbers and plots |
| `m2c_programs.py` | M2c | Finds the top genes shared across cell types, ranks genes for GSEA, calls `R/m2c_fgsea.R` and labels programs | `M2c_programs.md` | `test_m2c*.py` |
| `w4_outputs.py` | Wave 4 | Builds Figs 3–5 and the final table. It computes nothing new: every value shown is asserted equal to its source file | `W4_figures_table.md` | `w4_explained` notebook |

Contracts are in `audit/contracts/`, and tests in `tests/`.
