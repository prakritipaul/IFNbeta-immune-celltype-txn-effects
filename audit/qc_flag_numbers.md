# QC flag numbers (written by notebooks/00_eda.ipynb; proposed actions are added at review)

## Outlier cells (> 5 raw MADs above or below, on log counts or log genes, within donor × condition × cell type; groups < 10 cells not assessed)

| donor        |   ctrl |   stim |
|:-------------|-------:|-------:|
| patient_101  |     16 |     24 |
| patient_107  |      7 |     10 |
| patient_1015 |     65 |     49 |
| patient_1016 |     38 |     37 |
| patient_1039 |      8 |     11 |
| patient_1244 |     45 |     18 |
| patient_1256 |     44 |     32 |
| patient_1488 |     48 |     62 |

By cell type:

| cell_type         |   flagged |   cells |   pct |
|:------------------|----------:|--------:|------:|
| B cells           |        85 |    2571 |   3.3 |
| CD14+ Monocytes   |        14 |    5355 |   0.3 |
| CD4 T cells       |       172 |   10373 |   1.7 |
| CD8 T cells       |        92 |    2040 |   4.5 |
| Dendritic cells   |        33 |     432 |   7.6 |
| FCGR3A+ Monocytes |        59 |    1592 |   3.7 |
| Megakaryocytes    |         9 |     208 |   4.3 |
| NK cells          |        50 |    1991 |   2.5 |

## Groups with < 10 cells (donor × condition × cell type)

|                                             |   cells |
|:--------------------------------------------|--------:|
| ('Dendritic cells', 'patient_107', 'ctrl')  |       8 |
| ('Dendritic cells', 'patient_1039', 'ctrl') |       9 |
| ('Megakaryocytes', 'patient_101', 'stim')   |       9 |
| ('Megakaryocytes', 'patient_107', 'ctrl')   |       6 |
| ('Megakaryocytes', 'patient_107', 'stim')   |       4 |
| ('Megakaryocytes', 'patient_1016', 'ctrl')  |       8 |
| ('Megakaryocytes', 'patient_1039', 'ctrl')  |       4 |

## Cell types by donors with both conditions ≥ 10 cells

|                   |   donors with both conditions ≥ 10 cells | flag_underpowered_celltype   |
|:------------------|-----------------------------------------:|:-----------------------------|
| B cells           |                                        8 | False                        |
| CD14+ Monocytes   |                                        8 | False                        |
| CD4 T cells       |                                        8 | False                        |
| CD8 T cells       |                                        8 | False                        |
| Dendritic cells   |                                        6 | False                        |
| FCGR3A+ Monocytes |                                        8 | False                        |
| Megakaryocytes    |                                        4 | False                        |
| NK cells          |                                        8 | False                        |

## Donor sizes

| donor        |   cells |   fraction of median donor |
|:-------------|--------:|---------------------------:|
| patient_101  |    2024 |                       0.6  |
| patient_107  |    1042 |                       0.31 |
| patient_1015 |    5090 |                       1.52 |
| patient_1016 |    3358 |                       1    |
| patient_1039 |     992 |                       0.3  |
| patient_1244 |    3343 |                       1    |
| patient_1256 |    4134 |                       1.23 |
| patient_1488 |    4579 |                       1.37 |

## Composition shift, stim − ctrl (percentage points)

| donor        |   B cells |   CD14+ Monocytes |   CD4 T cells |   CD8 T cells |   Dendritic cells |   FCGR3A+ Monocytes |   Megakaryocytes |   NK cells |
|:-------------|----------:|------------------:|--------------:|--------------:|------------------:|--------------------:|-----------------:|-----------:|
| patient_1039 |       0.3 |              -1.8 |           2   |          -1.4 |               0.7 |                -0.2 |              0.9 |       -0.5 |
| patient_107  |       1.8 |              -8.2 |           4.9 |          -1   |               1.1 |                 0.9 |             -0.4 |        1   |
| patient_101  |      -0.9 |              -1.8 |          -1.4 |           0.4 |               0.5 |                 3   |             -0.5 |        0.8 |
| patient_1016 |      -0.2 |              -1.3 |          -2.9 |          -1.2 |              -0.2 |                 1.2 |              0.1 |        4.5 |
| patient_1244 |       0.2 |              -2.3 |           1.5 |          -0.9 |              -0.2 |                -0.6 |              0.1 |        2.2 |
| patient_1256 |       0   |               0   |          -0.4 |          -0.7 |               0.3 |                 1.5 |              0.2 |       -1   |
| patient_1488 |       0.2 |              -0.8 |          -0   |          -1.2 |               0   |                 0.7 |              0.3 |        0.9 |
| patient_1015 |      -2   |               0.1 |           0.2 |          -0.9 |               0.2 |                 0.5 |              0.1 |        1.8 |

## Technical covariates: median stim / ctrl ratio per cell type

| cell_type         |   n_counts |   n_genes |   pct_ribo |
|:------------------|-----------:|----------:|-----------:|
| B cells           |       1.06 |      1.07 |       0.84 |
| CD14+ Monocytes   |       1.13 |      1    |       0.53 |
| CD4 T cells       |       1.07 |      1.05 |       0.91 |
| CD8 T cells       |       1.07 |      1.03 |       0.85 |
| Dendritic cells   |       0.9  |      0.93 |       0.69 |
| FCGR3A+ Monocytes |       1.06 |      0.98 |       0.72 |
| Megakaryocytes    |       0.97 |      1.06 |       1.03 |
| NK cells          |       1.09 |      1.02 |       0.79 |

## Technical covariates by donor: each donor's median / overall median

| donor        |   n_counts |   n_genes |   pct_ribo |
|:-------------|-----------:|----------:|-----------:|
| patient_101  |       1.03 |      1.04 |       0.84 |
| patient_107  |       1.14 |      1.04 |       0.63 |
| patient_1015 |       1.04 |      1.04 |       0.8  |
| patient_1016 |       0.85 |      0.92 |       0.89 |
| patient_1039 |       1.13 |      1.09 |       0.81 |
| patient_1244 |       1.06 |      1.01 |       1.26 |
| patient_1256 |       0.99 |      1.01 |       1.04 |
| patient_1488 |       0.98 |      0.95 |       1.32 |
