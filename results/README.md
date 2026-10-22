# Final table: `final_table.csv`

One row per gene DE in at least one cell type (**2,801 genes**). DE = pseudobulk DESeq2 (`~ donor + condition`), q < 0.05 (`qvalue`, within cell type) and |log2FC| ≥ 1, stim vs ctrl. Built by `src/w4_outputs.py::final_table` from `data/processed/m2_table.csv` (arrows, labels) and `m2_calls.csv` (log2FC, q); every arrow is asserted equal to both.

| Column | Meaning |
|---|---|
| `gene` | Gene symbol as in the dataset |
| 8 cell-type columns | ↑ DE up, ↓ DE down, – not DE (or not tested) |
| `Shared response?` | M2a label, see below |
| `n_cell_types` | Number of cell types where the gene is DE |
| `expressed_in_n_cell_types` | Number of cell types where the gene was tested (baseMean > 1, from `m2_calls.csv` `tested`) |
| `megakaryocyte_lower_confidence` | True if the gene is DE in megakaryocytes (footnote 3) |
| `<cell type> log2FC` | Unshrunk DESeq2 log2 fold change, stim vs ctrl, verbatim from `m2_calls.csv` |
| `<cell type> q` | q-value; **empty = not tested** there (baseMean ≤ 1); the log2FC of an untested gene is shown but not used for any call |

**Labels in `Shared response?`**
- `shared`: DE in all 8 cell types, same direction (21 genes).
- `shared*`: DE in the 6 well-powered cell types (B, CD14+ Mono, CD4 T, CD8 T, FCGR3A+ Mono, NK), and in dendritic cells / megakaryocytes either DE or with a tested 95% CI reaching |log2FC| ≥ 1 the same way (87).
- `partial`: DE in ≥ 2 cell types, same direction, not shared (856).
- `mixed`: up in some cell types, down in others (11).
- `specific`: DE in 1 cell type, and no other tested cell type's 95% CI reaches a ≥ 2-fold change (241).
  A `specific` gene with expressed_in_n_cell_types = 1 is specific because only that cell type expresses it (229 of 241 clean specific genes); see FINDINGS.md.
- `specific*`: DE in 1 cell type, but another tested cell type's CI reaches ≥ 2-fold, so "specific" may be a power gap (1,585).

**Sort order:** shared → shared\* → partial → mixed → specific → specific\*, then `n_cell_types` (descending), then mean |log2FC| over the cell types where the gene was tested (descending), then gene.

**Footnotes**
1. \* (shared\*, specific\*): the label relies on 95% confidence intervals, taken only from cell types where the gene was tested.
2. Dendritic cells have 6 paired donors and megakaryocytes 4 (pseudobulk samples with < 10 cells dropped with their pair); the other six cell types have 8.
3. Megakaryocyte calls are lower confidence: 4 donors, 24% recovery of the paper's DE genes (V9), technical covariates partly track condition.
4. Condition is confounded with the 10x run (one pool and run per condition): a run effect cannot be separated from IFN-β in any result.
