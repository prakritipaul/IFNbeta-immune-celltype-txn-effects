# Audit ledger

One row per verified behaviour. ✅ verified · ⏳ incomplete · ❌ confirmed defect · 🔄 stale · 🚧 blocked

## M0 (EDA, join, filters, QC flags) — approved by the user 2026-10-01

| Module / file | Behaviour | Tests | Review | Status |
|---|---|---|---|---|
| `m0_qc.py` load | Fingerprints (h5ad + GEO); finite, non-negative integer counts | T7, T9, T13, T19, T21 | Codex re-check 3 | ✅ |
| `m0_qc.py` join | (sequence, condition) key; all 24,673 matched; donor/cluster/t-SNE (1e-9) agree; 234 collision cells kept distinct | T1, T2, T10, T14 | Codex review 1 + re-check | ✅ |
| `m0_qc.py` QC metrics | n_counts, n_genes, pct_ribo (84 ribosomal genes) | T5, T17 | Codex re-check 2 | ✅ |
| `m0_qc.py` basic filter | Cells ≥ 200 genes, then genes ≥ 3 cells (all conditions) → 24,562 × 15,701 | T3, T4 | Codex review 1 | ✅ |
| `m0_qc.py` flags | Outliers > 5 raw MADs both tails per donor × condition × cell type (groups < 10 skipped), on counts **or** genes; small groups; underpowered cell types; nothing dropped | T6, T11, T16, T18 | Codex re-check 2; planted bug 1 caught | ✅ |
| `m0_qc.py` conditions / assumptions | Control = `ctrl`/`control`, other = treatment; 8 donors × 2 conditions; 8 named cell types | T12, T7 | Codex re-check | ✅ |
| `m0_qc.py` run | Thresholds validated (incl. non-finite) before writing; outputs and cell accounting | T8, T15, T20 | Codex re-check 3 | ✅ |
| `00_eda.py` | Visual EDA; per-group outlier plots reproduce `flag_outlier` (assert); Fig 3a/3d recreations match the paper | notebook asserts | Codex brief 002; explained-by-example (9/9) | ✅ |

**Suite:** 35 passed (`pytest tests/ -q`), commit at M0 approval.

## Wave 1 (M1a, M3, Figs 1–2)

| Module / file | Behavior | Tests | Review | Status |
|---|---|---|---|---|
| `m1a_pseudobulk.py` pseudobulk | Raw counts summed per donor × cell type × condition (128 samples); totals conserved; `\|` in labels rejected | T1, T3, T9 | Codex review 1 + re-check | ✅ |
| `m1a_pseudobulk.py` A1 rule | < 10 cells in either condition drops the pair → 116 kept; re-applied in the no-outlier variant | T1, T2, T4 | Codex review 1 + re-check | ✅ |
| `m1a_varpart.R` filter | baseMean > 1 per cell type | T6 **dropped by user** | explained-by-example only | ⏳ untested by suite |
| `m1a_varpart.R` / `m1a_cancor.R` model | varianceStabilizingTransformation(blind), parametric fit logged; variancePartition oracle; pairing one row per condition | T7, T8 | Codex review 1 + re-check | ✅ |
| `m3_signature.py` | Gene sets (sha256, 6 mappings, 100% coverage); float32 normalise + log1p + score_genes; no-shared IFN-γ; mean per sample; A1 rule | T1–T10 | Codex review 1 + re-check | ✅ |
| `figures.py` / `fig1_fig2.py` | Fig 1 totals = data_facts; UMAP C2 recipe, no integration; I3 same cells/coords | notebook asserts | explained-by-example (4/4) | ✅ |
| All result modules | Planted-bug check | — | **skipped by user** | — |

Explained by example: M1a 9/9, M3 6/6 (incl. an independent hand-computed score), Figs 4/4. Briefs: `audit/reviews/005`, `006`.

**Suite:** 57 passed, 1 skipped (T6), Wave 1 merged into `planning`.

## Wave 2 (M1b)

| Module / file | Behavior | Tests | Review | Status |
|---|---|---|---|---|
| `m1b_deseq2.R` model | Paired `~ donor + condition` (unpaired for R2), control reference, ratio size factors, parametric fit logged, Wald, betaPrior FALSE, Cook's defaults | T1, T2, T9, T10 | Codex review 1 + re-check | ✅ |
| `m1b_deseq2.R` testing | Fit all genes, then baseMean > 1; qvalue per cell type on tested genes only | T4, T5, T13 | Planted bug 2 (qvalue universe): tests missed, Codex caught; T13 added | ✅ |
| `m1b_de.py` calls / IO | DE rule q < 0.05 and \|log2FC\| ≥ 1 (changed 2026-10-02; T3 being updated); input validation; R output must return exactly the input genes | T3, T6, T8, T12 | Codex review 1 + re-check | ✅ |
| `m1b_de.py` runs | main, no_outliers, unpaired, LODO (8/6/4), shuffles (35/10/3); 305 fits | T7 | Codex review 1 | ✅ |
| Real data | IFN tripwire: ≥ 10 IFN-α genes up per cell type, no canonical ISG down | T11 | — | ✅ |
| V9 (paper overlap) | > 70% of paper's DE genes recovered | notebook | user-accepted as reported | ❌ 2/8 pass (diagnostics in FINDINGS) |

Explained by example: M1b 6/6. Brief: `audit/reviews/009`.

**Suite:** 65 passed, 1 skipped (M1a T6), 10 real-data tests deselected in the merge check (T7/T11 passed on the real run). M1b merged into `planning` with the user's approval.

## Wave 3 (M2a, M2b)

| Module / file | Behavior | Tests | Review | Status |
|---|---|---|---|---|
| `m2a_calls.py` labels | DE re-derived (≥ 1); shared / shared* / partial / specific / specific* / mixed with precedence; tested-only CI rules; megakaryocyte flag | T1–T3, T8–T10 | Codex review 1 + re-check (2 fixes + 1 follow-up) | ✅ |
| `m2a_calls.py` LODO / integrity | Sign stability per call; exact `lodo_<donor>` runs (paired donors required); independent traceability | T4–T7, T11, T12 | Codex re-check; planted bug 3 caught by tests and Codex | ✅ |
| `m2b_validation.py` | P1–P4, N1, R1–R3, V9, V9b, F3b, modules, KEGG, I5, I6 | report tier (explained-by-example, all match) | — | ✅ (V9 ❌ 2/8, accepted) |

| `m2c_programs.py` / `m2c_fgsea.R` | Top 20 up/down + sharing; Jaccard; Wald-stat ranking; fgsea (Hallmark + KEGG), BH within cell type over both collections; program labels with n/t; cell-type universe explicit | T1–T11 | Codex review 1 + re-check; planted bug 4 missed by tests, caught by Codex → T11 | ✅ |

Explained by example: M2a 23/23, M2b all match, M2c all match. Briefs: `audit/reviews/012` (M2a). Wave 3 merged.

## Wave 4 (Figs 3–5, final table)

| Module / file | Behavior | Tests | Review | Status |
|---|---|---|---|---|
| `w4_outputs.py` Fig 3 | Per-cell violins + paired donor means (IFN-α; IFN-γ without shared genes); donor means equal `m3_sample_scores.csv` | notebook assertions | explained-by-example (all match) | ✅ |
| `w4_outputs.py` Fig 4 | Top 25 shared + top 5 specific per cell type among genes tested in all 8; columns by up-gene Jaccard; two colour scales; every cell equals `m2_calls` log2FC | notebook assertions | explained-by-example | ✅ |
| `w4_outputs.py` Fig 5 | 58 programs by label; n/t; columns by NES profiles; every cell equals `m2c_gsea` NES | notebook assertions | explained-by-example | ✅ |
| `results/final_table.csv` | 2,801 genes × 29 columns; arrows equal `m2_table` 100%; `expressed_in_n_cell_types` (229 of 241 `specific` = 1) | notebook assertions | explained-by-example | ✅ |

Display tier (no Codex tests). **All waves merged into `planning`.**
