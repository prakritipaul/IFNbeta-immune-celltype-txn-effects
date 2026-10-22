# M0 Codex review of the planted-bug snapshot (not told about the bug)

Finding #1 is the planted bug (caught). Findings #2–#5 are real observations about the actual code (the planted change touched only the outlier loop), to be decided by the user.

| # | Severity | file:line | Finding | Concrete failing scenario |
|---:|---|---|---|---|
| 1 | High | `src/m0_qc.py:206` | Outlier detection evaluates only `n_counts`; the approved decision requires the union of outliers from log counts **and log genes detected**. | A cell with ordinary total counts but `log1p(n_genes)` more than 5 raw MADs from its donor × condition × cell-type median receives `flag_outlier=False`. |
| 2 | Medium | `src/m0_qc.py:63` | Required input dimensions are never checked: 24,673 × 15,706 for the h5ad and 29,065 GEO rows. | A truncated GEO table containing every h5ad cell, or a structurally valid subset h5ad retaining eight donors and cell types, can pass and produce results despite violating the contract. |
| 3 | Medium | `src/m0_qc.py:74` | The raw-count validation accepts positive infinity because `inf >= 0` and `round(inf) == inf`. | An `X` entry of `np.inf` passes `load_inputs`, then yields infinite library size and invalid downstream QC metrics instead of raising `ValueError`. |
| 4 | Medium | `src/m0_qc.py:264` | Threshold arguments are not validated, contradicting the module guarantee that invalid input raises `ValueError`. | `run_m0(..., min_genes=-1, min_cells=-1, n_mads=-1)` runs and writes output; the negative MAD threshold flags every assessed cell. |
| 5 | Medium | `src/m0_qc.py:264` | Production identity checks are optional: both hashes and the expected cell-type names default to `None`, so the approved GEO fingerprint and eight exact author labels can be bypassed. | Calling `run_m0` with an altered GEO file containing eight misspelled but unique cell-type labels, while omitting the optional arguments, can write a contract-invalid output. |