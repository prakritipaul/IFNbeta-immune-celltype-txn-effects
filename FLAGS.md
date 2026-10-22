# Flags

Things the user needs to know or decide. **Open** = waiting on the user; **To fix** = agreed, not yet done; **Done**.
Details: `audit/qc_flags.md` (QC numbers), `audit/reviews/002_M0_codex_reviewer_brief.html` (code review).

## M0 (EDA)

| # | Flag | Example | Status |
|---|---|---|---|
| 1 | **Outliers are mostly normal myeloid cells** | 357 cells flagged; 353 for having *too much* RNA. 95 of 432 DCs (22%) flagged, because DCs normally carry more RNA than T cells and the cutoff pools all cell types | **Done (D2):** per donor × condition × cell type, both tails, groups < 10 skipped → 514 flagged (379 high, 135 low); DCs 7.6%. **Done:** keep all cells; DE rerun without flagged cells as a sensitivity check (M1b) |
| 2 | **Two cell types have few paired donors** | Donors with both conditions ≥ 10 cells: 8 for six cell types, DC 6, megakaryocytes 4 | **Done (D3):** all 8 main, as in the paper. "Shared" requires DE in all 8 |
| 3 | **"Ribosomal" includes non-ribosomal genes** | `RPS6KA1` is a kinase, not a ribosomal protein. Effect: median 0.11 pp of % ribo | **Done (D4):** 84 genes (14-gene exclusion list) |
| 4 | **% ribo drops after stimulation** | Median stim/ctrl: CD14+ mono 0.53, DC 0.70 | **Done:** % ribo stays out of the DE model. With 84 genes: CD14+ mono 0.53, DC 0.69 |
| 5 | **Two small donors** | `patient_107` (1,042 cells) and `patient_1039` (992): about 0.3× the median donor | **Done:** keep both |
| 6 | **One donor's composition shifts** | `patient_107`: CD14+ monocytes −8.2 percentage points in stim; others < 3 | **Done:** keep; caveat in the write-up |
| 7 | **No % mito** | The h5ad has no `MT-` genes | Done: caveat (O1) |
| 8 | **All `FutureWarning`s silenced in the notebook** | A pandas behaviour change would go unannounced. Claude added this without approval (rule 0) | Done |
| 9 | **Thresholds rely on code defaults** | Editing a default in `m0_qc.py` would change the output, but the plots would still draw the old cutoff | Done |
| 10 | **Some numbers aren't produced by the notebook** | "95 of 432 DCs" was computed outside the notebook | Done |
| 11 | **Two filtering runs aren't checked against each other** | The notebook and `run_m0` both give 24,562 cells today, but nothing asserts it | Done |
| 12 | **Condition is confounded with 10x run** | One pool per condition, one run each (paper: "Two pools, IFN-β-treated and control…"), so a run-to-run technical difference looks exactly like an IFN-β effect | Caveat: can't be modelled out; stated in the write-up |
| 13 | **`fillna(0)` in the outlier table** | Blank would mean "none flagged" | **Done (D1):** approved |
| 14 | **EDA checklist gaps** (Codex brief) | GEO columns not plotted; covariates not compared by donor | **Done (D5):** both added. New: `patient_1488` and `patient_1244` have ~1.3× the median % ribo; `patient_107` 0.63× |
