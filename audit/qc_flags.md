# M0 QC flags: for the user to decide

Numbers: `audit/qc_flag_numbers.md` (written by the notebook). Plots: `notebooks/00_eda.ipynb`. Nothing below has been dropped; each row is a proposal.

| # | Flag | Evidence | Proposed action | Why |
|---|---|---|---|---|
| F1 | **Outlier cells** (after D2: per donor × condition × cell type, both tails, groups < 10 cells not assessed) | **514 cells (2.1%)**: 379 high, 135 low. By cell type: CD4 T 172 (1.7%), CD8 T 92 (4.5%), B 85 (3.3%), FCGR3A+ mono 59 (3.7%), NK 50 (2.5%), DC 33 (7.6%), CD14+ mono 14 (0.3%), megakaryocytes 9 (4.3%). Notebook §7 | **Decided: keep all**; sensitivity check without them in M1b | Pseudobulk sums cells, so one abnormal cell can dominate a small group; dropping 2% barely changes large groups |
| F2 | **Small groups** (< 10 cells) | 7 groups: DC `patient_107`/`patient_1039` ctrl (8, 9 cells); megakaryocytes in 5 donor × condition groups (4–9 cells). Notebook §8 | Input to **D3** | Too few cells for a reliable pseudobulk sample |
| F3 | **Paired donors per cell type** (both conditions ≥ 10 cells) | 8 for six cell types; **DC 6; megakaryocytes 4**. Notebook §8 | **Decided (D3): all 8 main** | Matches the paper; "shared" requires DE in all 8 |
| F4 | **Weak donors** | `patient_107` (1,042 cells) and `patient_1039` (992) are about 0.3× the median donor | **Keep** | Both have both conditions and most cell types; pseudobulk sums cells, so smaller donors just contribute noisier samples |
| F5 | **Composition shifts** | Most donors move < 3 percentage points; `patient_107` CD14+ mono −8.2 pp (stim − ctrl). Notebook §9 | **Keep**, and note as a caveat | Could be IFN-driven cell loss or annotation; affects proportions (Fig 3d), not per-cell-type DE |
| F6 | **% ribosomal tracks condition** (84 ribosomal genes, D4) | Median stim/ctrl: CD14+ mono 0.53, DC 0.69, FCGR3A+ mono 0.72; counts and genes ≈ 1.0. By donor: `patient_1488` 1.32×, `patient_1244` 1.26×, `patient_107` 0.63× the overall median | **Don't use % ribo as a DE covariate**; report it | Likely biology (strong ISG transcription dilutes the ribosomal share). As a covariate it would absorb part of the IFN effect |
| F7 | **Basic filters** | Removed 111 cells (< 200 genes) and 5 genes (< 3 cells) | Done (agreed) | — |
| F8 | **No % mito** | The h5ad has no `MT-` genes | Caveat (agreed, O1) | — |

## Decisions

Q1 → D2 (per cell type, both tails), D3 (all 8 main), D4 (84 ribosomal genes) are decided; see `audit/decisions.md`.

F1: keep all cells (sensitivity check in M1b). F4: keep both small donors. F5: keep `patient_107`, shift noted as a caveat. F6: % ribo stays out of the DE model. **All M0 QC flags decided.**
