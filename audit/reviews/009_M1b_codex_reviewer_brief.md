# M1b reviewer brief (Codex, commit b254e07)

## 1. What M1b does, in 3 lines

- Input: raw pseudobulk counts for 15,701 genes in 116 kept donor × condition × cell-type samples (16 samples for six cell types, 12 dendritic, 8 megakaryocyte).
- It fits treatment versus control separately in all 8 cell types with paired DESeq2 (`~ donor + condition`), then repeats no-outlier, unpaired, leave-one-donor-out and balanced-label-shuffle runs.
- Output: 125,608 main rows (8 × 15,701), 4,788,805 all-run rows across 305 fits/46 named runs, and 305 fit records; main calls range from 27 DE genes in megakaryocytes to 2,335 in CD14+ monocytes.

## 2. Assumptions coded in

- Kept samples must form exact donor pairs, one control and one treatment per cell type ([src/m1b_de.py:91](../../src/m1b_de.py#L91)); this prevents an apparently paired fit from using incomplete or duplicated pairs.
- Treatment is explicitly contrasted against `ctrl`/`control` ([R/m1b_deseq2.R:73](../../R/m1b_deseq2.R#L73)); reversing this would turn every up-response into a down-response.
- Size-factor, Wald-fit, parametric-trend, no-beta-prior and Cook's settings are fixed explicitly ([R/m1b_deseq2.R:66](../../R/m1b_deseq2.R#L66), [R/m1b_deseq2.R:70](../../R/m1b_deseq2.R#L70), [R/m1b_deseq2.R:74](../../R/m1b_deseq2.R#L74)); these determine fitted effects and p-values.
- Genes are fitted first, but only `baseMean > 1` genes enter per-cell-type q-value estimation ([R/m1b_deseq2.R:84](../../R/m1b_deseq2.R#L84)); a gene exactly at 1 is not tested.
- Calls use strict cutoffs: q < 0.05 and unshrunk |log2FC| > 1 ([src/m1b_de.py:46](../../src/m1b_de.py#L46)); for example q = 0.05 or log2FC = 1 produces no DE call.

## 3. Decisions the user didn't make

None found. The result-affecting settings and rerun definitions trace to approved O1–O7 and the contract.

## 4. The lines worth reading yourself

- [R/m1b_deseq2.R:64](../../R/m1b_deseq2.R#L64): selects the paired model for the main analysis.
- [R/m1b_deseq2.R:73](../../R/m1b_deseq2.R#L73): fixes the biological direction as treatment minus control.
- [R/m1b_deseq2.R:84](../../R/m1b_deseq2.R#L84): defines which fitted genes enter multiple testing.
- [R/m1b_deseq2.R:87](../../R/m1b_deseq2.R#L87): performs q-value correction on that tested universe, separately inside each cell type.
- [src/m1b_de.py:186](../../src/m1b_de.py#L186): swaps condition labels within donors for the negative-control runs.

## 5. What nothing checks

- No independent oracle checks the numerical log2FC or p-value of a real gene. The simulation checks direction and calls, but not that, for example, the reported ISG15 log2FC agrees with an independently calculated paired estimate.
- No fingerprint ties the shipped CSVs to the current inputs/code. Tests can regenerate results, but a stale copied `m1b_de.csv` could still look structurally valid and pass the notebook's displays.
- M1b deliberately does not judge robustness. For example, whether a CD14+ monocyte call survives all leave-one-donor-out fits or appears in shuffled-label fits is deferred to M2.

## 6. Verdict

Ready for the user's scientific gate, with no process failure found. Fresh non-real-data execution: **16 passed, 2 real-data tests deselected**; the planted q-value-universe bug was caught by independent review and its T13 regression test now passes. Cheap checks on the shipped outputs found all 15,701 genes in every main cell type, no duplicate main keys, all 305 fits parametric, no unexplained main q-value NAs, and 14–84 Hallmark IFN-α genes called up per cell type (canonical ISGs were never down). Before approval, the user should inspect the notebook's DE counts/top-gene tables and volcano plots, especially the small megakaryocyte result (23 up, 4 down), and remember that robustness/negative-control interpretation remains an M2 gate.
