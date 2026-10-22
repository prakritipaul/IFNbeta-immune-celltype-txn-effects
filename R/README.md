# R/

The statistics that follow Kang 2018 run in R, as they did in the paper (R 4.3, Bioconductor). Python (`src/`) writes the inputs to CSV, calls one of these scripts with `Rscript`, and reads its output CSV back. Each script's header gives its exact usage, inputs and outputs.

| Script | Called by | What it does |
|---|---|---|
| `m1a_varpart.R` | `src/m1a_pseudobulk.py` | Variance partitioning per gene (`variancePartition`): `~ (1\|donor) + (1\|condition) + log(n_cells) + log(total_counts)` on variance-stabilised pseudobulk |
| `m1a_cancor.R` | `src/m1a_pseudobulk.py` | Canonical correlations between those terms, to check whether the technical covariates are just donor in disguise |
| `m1b_deseq2.R` | `src/m1b_de.py` | Pseudobulk DE: DESeq2 `~ donor + condition` (or unpaired for a sensitivity run), Wald test, then `qvalue` per cell type on genes with baseMean > 1 |
| `m2c_fgsea.R` | `src/m2c_programs.py` | Preranked GSEA (`fgsea`) per cell type on Hallmark + KEGG, with BH correction across both collections within a cell type |

The packages and their versions are listed in the main README under Setup.
