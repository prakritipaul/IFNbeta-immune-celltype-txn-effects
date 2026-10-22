# M1b test plan

**Status:** APPROVED by the user (2026-10-01). Codex proposed 12 blind from the contract; Claude merged three input-failure tests into T6 and added the shuffle-structure check to T7; the user added T11. Codex writes exactly the approved tests, in `tests/test_m1b.py`, before the implementation exists.

| # | Test | Bug it catches | Expected answer from | From |
|---|---|---|---|---|
| T1 | **Oracle (I4):** simulated paired counts, 8 donors, some genes planted strongly up, some strongly down, many null → full Python → R run: planted genes `up`/`down` with the right sign; nearly all null genes `none` | Reversed contrast, broken Python↔R hand-off, wrong design, mass false positives | Simulation parameters; validation plan I4 | Codex |
| T2 | **Pairing matters:** 8 donors with very different baselines; ctrl and stim differ only by small random within-pair noise (no condition effect) → log2FC ≈ 0, no DE genes. *(Revised: identical pairs left DESeq2 unable to fit dispersion.)* | Donor term missing; samples attached to the wrong donor or condition | No treatment effect by construction | Codex |
| T3 | **Call boundaries:** (q 0.0499, LFC 1.01) → significant, `up`; (q 0.05, LFC 2) → not significant, `none`; (q 0.01, LFC 1) → significant, `up` *(rule changed to \|log2FC\| ≥ 1 by the user, 2026-10-02)*; (q 0.01, LFC −1.01) → `down` | `<=` vs `<`; `>=` vs `>`; missing absolute value | DE rule in decisions | Codex |
| T4 | **baseMean after fitting:** genes with fitted baseMean just below 1, exactly 1, just above 1 → only the last is tested by qvalue and output as tested; all genes were fitted by DESeq2 | Pre-filtering; `>= 1` | O2 | Codex |
| T5 | **qvalue per cell type:** same p-values for cell type A in two runs, very different p-values for cell type B → A's q-values identical | One correction across cell types | O3 | Codex |
| T6 | **Bad input stops before R runs, nothing written:** sample ids that don't match between counts and metadata; a donor with ctrl but no stim, or two ctrl rows; negative, fractional or non-numeric counts | Positional joins, unpaired leftovers, coerced counts | Contract input checks | Codex (3 merged) |
| T7 | **Runs complete (real data):** main fits use 16 samples in six cell types, 12 in DCs, 8 in megakaryocytes; LODO 8 / 6 / 4 refits; 35 shuffles for each 8-donor cell type, 10 for DCs, 3 for megakaryocytes; each shuffle swaps exactly half the donors and no two shuffles are the same split (a swap set and its complement count as the same) | Lost cell types, incomplete or duplicated shuffles | Contract; O4; C(8,4)/2 = 35, C(6,3)/2 = 10, C(4,2)/2 = 3 | Codex + Claude |
| T8 | **Row keys:** two cell types share gene names and R returns rows in shuffled order → every `(run, cell_type, gene)` row carries its own values, no duplicates or losses | Joining on gene alone; trusting row order | Constructed unique values | Codex |
| T9 | **Cook's outliers, unpaired design (`~ condition`):** an extreme single-sample gene in an 8-donor and in a 4-donor cell type → DESeq2's default handling kept (replaced in the 8-donor type, p-value NA in the 4-donor type), `cooks_flag` set, NA reported, never dropped or filled; an unexplained NA stops the run. Paired fits: `cooks_flag` never set. *(Revised: DESeq2 computes Cook's only in design cells with ≥ 3 samples, so never under `~ donor + condition`.)* | Cook's handling turned off; NAs silently removed | O1 (corrected) | Codex |
| T10 | **Reported statistics:** `ci_low/high = log2FC ∓ 1.96 × lfcSE`; the fit type used is logged for every run × cell type; `log2FC_shrunk` exists but calls use the unshrunk `log2FC` | Wrong CI; unlogged fallback; calls on shrunken LFC | Wald CI identity; O5 | Codex |
| T11 | **IFN sanity tripwire (real data, main run):** in every cell type, ≥ 10 of the 97 Hallmark IFN-α genes (after the 6 approved symbol mappings) are `up`, and none of the 10 canonical ISGs (ISG15, IFI6, IFIT1, IFIT2, IFIT3, MX1, OAS1, RSAD2, IRF7, STAT1) is `down` | Reversed contrast, wrong pairing, broken pipeline producing no IFN signal | Kang 2018 Supp Table 1: 45–84 of 97 IFN-α genes and 10/10 canonical ISGs up in every cell type; 10 is a deliberately low bar | User |

## Interface the tests call

Naming/interface only; no result-changing choices. Invalid input raises `ValueError` naming the violated assumption, before R runs or anything is written.

```python
# src/m1b_de.py
classify(results, q_threshold=0.05, lfc_threshold=1.0) -> pd.DataFrame
    # Adds `significant` (qvalue < q_threshold) and `de` ("up"/"down"/"none": significant and |log2FC| >= lfc_threshold, by sign of log2FC).
    # Rows with NaN qvalue: significant=False, de="none" (and must carry cooks_flag or an all-zero reason).
balanced_shuffles(donors) -> list[frozenset]
    # All distinct ways to swap ctrl/stim for exactly half of the donors; a set and its complement count once (8 donors → 35).
run_deseq2(counts, samples, design, run_name, min_basemean=1.0) -> tuple[pd.DataFrame, pd.DataFrame]
    # counts: genes × samples (int); samples: sample, cell_type, donor, condition, kept; design "paired" (~ donor + condition) or "unpaired" (~ condition).
    # Calls R/m1b_deseq2.R on kept samples; returns (results, fits). results rows: run, cell_type, gene, baseMean, log2FC, lfcSE,
    #   ci_low, ci_high, pvalue, qvalue, tested (baseMean > min_basemean), cooks_flag, log2FC_shrunk, significant, de.
    # fits rows: run, cell_type, fit_type_used, n_samples, n_genes_tested.
run_m1b(in_dir, out_dir, runs=("main", "no_outliers", "unpaired", "lodo", "shuffle")) -> None
    # Reads m1a_* files from in_dir; writes m1b_de.csv (main), m1b_runs.csv (all runs), m1b_fits.csv. Nothing on failure.
```

```
# R/m1b_deseq2.R (called by run_deseq2; tests may call it directly)
Rscript R/m1b_deseq2.R <counts.csv> <samples.csv> <paired|unpaired> <results_out.csv> <fits_out.csv> <min_basemean>
    # Per cell type: DESeq2 (ratio size factors, parametric fit with logged fallback, Wald, betaPrior=FALSE, Cook's defaults),
    #   fit all genes, then qvalue (lambda 0.05–0.95 by 0.05, smoother, smooth.df 3) on genes with baseMean > min_basemean.
    # Non-zero exit on invalid input.
```

## Added after Codex review 1 (user-approved), in `tests/test_m1b_review.py`

| # | Test | Bug it catches | Expected answer from |
|---|---|---|---|
| T12 | `run_deseq2` given an R result with the right number of rows but a substituted gene name (per cell type) → `ValueError`, nothing returned | Gene identities not checked at the R boundary (review #1) | Contract: results are per cell type × input gene |

## Added after the planted-bug check (M1b)

| # | Test | Bug it catches | Expected answer from |
|---|---|---|---|
| T13 | R script on a toy cell type with tested genes (baseMean > 1) plus many untested genes (baseMean ≤ 1) with non-NA p-values: the tested genes' q-values equal `qvalue()` run on the tested genes' p-values alone (same parameters), and adding or removing untested genes doesn't change them | qvalue run on all genes and untested ones masked afterwards (planted bug 1, missed by T1–T12) | O2/O3: qvalue only on tested genes |
