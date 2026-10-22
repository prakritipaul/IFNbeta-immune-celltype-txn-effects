# M0 cell flow

| step | cells in | cells out | genes in | genes out | matched |
|---|---|---|---|---|---|
| load | 24673 | 24673 | 15706 | 15706 |  |
| join | 24673 | 24673 | 15706 | 15706 | 24673 |
| qc_metrics | 24673 | 24673 | 15706 | 15706 |  |
| basic_filter | 24673 | 24562 | 15706 | 15701 |  |
| flag_qc | 24562 | 24562 | 15701 | 15701 |  |

# M1a cell flow (pseudobulk)

min_cells = 10 (A1). Asserted: per-gene sums over samples = cell-level sums; n_cells sums to cells in; every kept donor × cell type has both conditions.

| step | cells in | cells in samples | samples | kept samples | cells in kept samples | genes |
|---|---|---|---|---|---|---|
| pseudobulk (all cells) | 24562 | 24562 | 128 | 116 | 24451 | 15701 |
| drop flag_outlier | 24562 | 24048 | | | | 15701 |
| pseudobulk (no outliers) | 24048 | 24048 | 128 | 116 | 23941 | 15701 |
