# Decisions

Current, user-approved decisions only. Agents read this and never re-ask. When a decision changes, **replace** its row; don't append a contradicting one. History lives in git.

## Setup and process

| Decision | Reason |
|---|---|
| Environment: `uv` venv (Python) + R 4.3 with Bioconductor DESeq2 and qvalue | Installed and working; R matches the paper's tools |
| Rules live only in project `CLAUDE.md` (`AGENTS.md` links to it for Codex) | One source of truth |
| Module plan: M1 split into M1a (pseudobulk + variance partitioning) and M1b (DE) | Covariate choice for DE needs variance partitioning first |

## Data and cells

| Decision | Reason |
|---|---|
| Cells = all 24,673 cells in `data/kang_2018.h5ad` (all are authors' singlets, all found in the GEO metadata), before QC filtering. GEO-only rows (doublets, ambiguous, 6 unlabelled) are not used | h5ad has no counts for them |
| **Cell-type labels = the authors'** (GEO GSE96583 `batch2.total.tsne.df.tsv`, column `cell`), not the h5ad's `cell_type` | h5ad's `cell_type` is a later re-analysis; agrees with the authors on only 91% of cells |
| Join on (barcode without suffix, condition) | Barcode sequences occur in both runs (313 across all droplets, 234 among our singlets), renamed differently in the two tables |
| CD14+ and FCGR3A+ monocytes are separate cell types | They respond differently to IFN |
| M0 (EDA): % mito dropped from QC and covariates (the h5ad has no `MT-` genes; stated as a caveat); gene filter "≥ 3 cells" counted across all cells; QC flags computed after the basic filters; outlier flags use raw (unscaled) MAD, > 5 MADs **above or below** the median on log counts / log genes, computed per **donor × condition × cell type** (groups < 10 cells not assessed); "ribosomal" = `RPS`/`RPL` prefix minus 14 non-ribosomal genes (`RPS6KA1`–`KA5`, `RPS6KB1`, `RPS6KB2`, `RPS6KC1`, `RPS6KL1`, `RPS10P7`, `RPSAP58`, `RPL34-AS1`, `RPL17-C18orf32`, `RPS19BP1`) → 84 genes; `fillna(0)` in the notebook's outlier table approved (blank means none flagged); EDA also plots every GEO column and summarises each covariate by donor; logic in `src/m0_qc.py`, called by the notebook | Authors already QC'd cells (singlets only); library size and % ribo remain. Testable code, readable notebook |
| M0 review fixes: weak donors, composition shifts and covariate trends are reported as numbers and plots for the user to judge (no automatic flags); GEO file verified by sha256; conditions: exactly 2 labels, the control is whichever is `ctrl` or `control` (defined once in the code), the other is the treatment; stored in the output as `uns['control']`/`uns['treatment']` and cell-type names must match the 8 authors' names on real data; t-SNE agreement tolerance explicit (absolute 1e-9); `cell_flow.md` written to `audit/` | User accepted Codex review 1 recommendations |
| **Keep all 8 donors**, including the two small ones (`patient_1039`, `patient_107`, ~0.3× the median donor); keep `patient_107` despite its −8.2 pp CD14+ mono composition shift (noted as a caveat); **% ribo stays out of the DE model** | Both conditions present for every donor; the shift affects proportions (Fig 3d), not per-cell-type DE; % ribo tracks condition, so as a covariate it would absorb part of the IFN effect |
| M0 input guards (Codex planted-bug review #3–5): reject non-finite counts; reject invalid thresholds; the real-data run always checks both fingerprints and the 8 cell-type names | Cheap protections against silent bad input |
| **Keep all cells**, including the 514 flagged outliers; `flag_outlier` stays as a column. Sensitivity check in M1b: rerun DE without flagged cells and compare calls | Matches the paper (all singlets kept); low-side flags can't be confirmed as damaged without % mito; high side may hold ~390 same-donor doublets demuxlet can't see (estimate), so the sensitivity run tests whether they matter |
| D3: **all 8 cell types are main results** (as in the paper), including DC (6 paired donors) and megakaryocytes (4). The definition of "shared" is in the M2 section (O1) | User chose (b); matches Kang 2018, which ran DE on all 8 |

## Wave 1 (M1a, M3, Figs 1–2)

| Decision | Reason |
|---|---|
| A1: a donor × condition × cell type group with < 10 cells drops **that donor from that cell type** (both conditions): 128 → 116 samples; DC keeps 6 donors, megakaryocytes 4. Same rule for M3's per-sample summaries | Matches the paired-donor definition D3 was decided on; an unpaired leftover adds nothing under `~ donor + condition` |
| A2–A3: variance partitioning with R `variancePartition`, per cell type, `~ (1\|donor) + (1\|condition) + log(n_cells) + log(total_counts)`, on DESeq2 `vst(blind=TRUE)`, genes with baseMean > 1 | Standard tool and input after DESeq2 |
| A4 (decided at the gate): **M1b keeps `~ donor + condition`, no technical covariate.** Variance partitioning shows `log n_cells` and `log total_counts` correlate r = 0.92–0.997 with each other and are 72–98% explained by donor, so their fractions re-label donor variation; they barely track condition (\|r\| ≤ 0.19; megakaryocytes ~0.31). `canCorPairs` plot added as evidence | No effect separate from donor, which is already in the model; matches Kang 2018 |
| Codex review 1 (Wave 1): fix R pairing check (one row per condition), condition names read from data in the M1a notebook, reject the `\|` character in labels, explicit float32 cast in M3; **skip** atomic writes in M3 (validation already precedes writing; only I/O errors would leave partial files) | User accepted Claude's recommendations |
| Wave 1: **planted-bug check skipped** for M1a and M3; Wave 1 merged into `planning` after the reviewer briefs and explained-by-example notebooks | User trusts the checks so far (blind tests + Codex review + re-check) |
| M1a test T6 (baseMean boundary) dropped; the filter is checked only by review | User: fixing the zero-variance toy is overkill |
| Outlier pseudobulk samples (M1a PCA): plots + distances, no automatic flag | Same as M0's weak-donor precedent |
| B1: `scanpy.tl.score_genes(ctrl_size=50, n_bins=25, random_state=0)` after `normalize_total(target_sum=1e4)` + `log1p` | Standard (Seurat AddModuleScore); I6 checks library-size dependence |
| B2: MSigDB **2026.1.Hs** Hallmark (`h.all.v2026.1.Hs.symbols.gmt`, sha256 `eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596`) | Latest release |
| B3: per-sample score = mean of cell scores | Simple; matches pseudobulk summing |
| B4: map 6 set genes to our names (5 renamed after 2017, plus PSMA2; WARS1 is in both sets, so 7 set entries): TENT5A→FAM46A, WARS1→WARS, MARCHF1→MARCH1, RIGI→DDX58, TMT1B→METTL7B, PSMA2→PSMA2.1 (confirm it is PSMA2) → 100% coverage of both sets | Same genes under post-2017 names |
| A3 (revised): transform = DESeq2 **`varianceStabilizingTransformation(blind=TRUE)`** for all 8 cell types (not `vst`) | `vst` needs ≥ 1,000 genes with baseMean > 5; megakaryocytes have 612. One method everywhere; DESeq2's own suggestion |
| DESeq2 dispersion trend: `fitType="parametric"` with DESeq2's automatic fallback to local allowed, **logged per cell type** (also in M1b) | Default DESeq2 behaviour, as in the paper; the switch is recorded, not silent |
| Pinned defaults (approved): size factors `type="ratio"`; variancePartition `REML=FALSE`, `useWeights=TRUE`; M3 `score_genes(ctrl_as_ref=True)`; `normalize_total(exclude_highly_expressed=False)`, natural-log `log1p`, scoring on a float32 working matrix (cast explicitly); Figs: HVGs on all cells (no `batch_key`), `n_bins=20`, PCA `svd_solver="arpack"`, neighbors `metric="euclidean"`, UMAP `spread=1.0`, `init_pos="spectral"` | Standard defaults, now explicit |
| M1a notebook PCA / outlier-sample distances on the same variance-stabilised data as variance partitioning | Consistency |
| PSMA2.1 confirmed as PSMA2 (B4) | No plain `PSMA2` in the data; `.1` is the make-unique suffix; expression and IFN response (log2FC 1–2 in the same cell types) match the paper's PSMA2 |
| **V9 result accepted as reported:** passes (> 70% recovery) in 2 of 8 cell types (CD14+ Mono, CD4 T); fails in 6 (24–70%). Criterion unchanged; method unchanged. Diagnostics (in `FINDINGS.md`): pooled qvalue and all-8-donors for DCs/megakaryocytes each explain part of the gap, neither all | User: our calls are a near-subset of the paper's with no sign flips; method stays as decided |
| V9 symbol mapping (paper Supp Table 1 → our names): 5 Excel date-mangled symbols restored by Ensembl ID (SEPT9, MARCH2, SEPT6, SEPT4, MARCH1); 5 make-unique names CLN3, MRPS17, NDUFB8, PSMA2, TAP2 → `.1` → all 3,055 paper genes matched | User approved; otherwise these genes could never count as overlap. Only PSMA2 checked individually (B4) |
| C1: UMAP **without integration** | The stim/ctrl separation within cell types is the IFN response Fig 2 shows |
| C2: `normalize_total(1e4)`, `log1p`, 2,000 HVGs (`flavor="seurat"`), `scale(max_value=10)`, PCA 50 comps, neighbors k = 15 on 30 PCs, UMAP `min_dist=0.5`, `random_state=0` | Standard scanpy recipe, every parameter explicit |

## M1b (DE) settings

| Decision | Reason |
|---|---|
| O1: DESeq2 Cook's-distance handling at **defaults** (`cooksCutoff` on, `minReplicatesForReplace=7`). DESeq2 only computes Cook's distances for samples in design cells with ≥ 3 samples, so under `~ donor + condition` (every sample alone in its cell) **it never applies to paired fits**, here or in the paper; it acts only in the unpaired run (R2: 61 genes flagged on real data). Affected genes reported in `cooks_flag`. Single-donor-driven calls are caught by LODO (R1) instead | Exactly the paper's DESeq2 behaviour (best for V9); premise corrected after implementation (`DESeq2:::recordMaxCooks`) |
| O2: fit DESeq2 on all genes, **then** keep baseMean > 1 for qvalue (no pre-filter); `padj` not used | DESeq2's recommended practice; baseMean from the fit, as in Supp Table 1 |
| O3: qvalue **per cell type**, `lambda = seq(0.05, 0.95, 0.05)`, `pi0.method = "smoother"`, `smooth.df = 3` (defaults, explicit) | Matches the paper's per-cell-type tables |
| O4: N1 label shuffles in **all 8 cell types**: within-donor ctrl/stim swaps for exactly half of each cell type's paired donors, every distinct split once (a swap set and its complement count once): 35 with 8 donors, 10 for Dendritic cells (6), 3 for Megakaryocytes (4). LODO (R1) refits once per paired donor (8, 6 or 4) | User: include the small cell types too |
| O5: shrunken log2FC for V9b via `lfcShrink(type="normal")` | Built in; DESeq2's default shrinkage in the paper's era |
| O6: Wald test, `betaPrior=FALSE`, passed explicitly | Defaults; unshrunk MLE as decided |
| O7: planted-bug check **done** for M1b | DE is the core result |
| M1b tests T2 and T9 rewritten (user-approved): T2's toy gets small within-pair noise (identical pairs make DESeq2's dispersion fit impossible); T9 tests Cook's handling on the unpaired design, where DESeq2 applies it. `m1b_runs.csv` keeps every gene of every run (0.84 GB, gitignored) as the contract specifies | Wrong test premises, not wrong code |
| Codex review 1 (M1b): fix #1 (check R returns exactly the input genes, new test T12), #3 (notebook reads `m1b_runs.csv` by columns only); #2 resolved by updating the contract (`m1b_de.csv` keeps a constant `run` column) | User accepted Claude's recommendations |
| M1b test T11 (user): real-data tripwire, ≥ 10 of 97 IFN-α genes `up` in every cell type and no canonical ISG `down`. The only pass mark besides V9; P1/P2 stay numbers | Paper has 45–84/97 and 10/10 in every cell type; catches a broken pipeline |

## M2 (shared vs specific) settings

| Decision | Reason |
|---|---|
| O1: **shared** = DE in all 6 well-powered cell types (B, CD14+ Mono, CD4 T, CD8 T, FCGR3A+ Mono, NK) in the same direction, **and** in Dendritic cells and Megakaryocytes either DE in that direction or with the 95% CI reaching \|log2FC\| ≥ 1 in that direction; genes that needed the CI rule are flagged (`shared*`). **108 genes** on the M1b results with the tested-only CI rule below (120 under the first, any-CI reading; 21 if all 8 had to be DE). Replaces "shared requires DE in all 8" (D3), which still holds that all 8 cell types are main results | New evidence: megakaryocytes (4 donors) have 27 DE genes, which would cap "shared" at 27; same logic as specific\* |
| O2: `shared_response` labels: `shared` / `shared*` (O1), `partial` (DE in ≥ 2 cell types, not shared), `specific` / `specific*` (DE in 1), `mixed` (up in some, down in others) | Informative; mixed-direction genes (11) get their own label |
| O3: megakaryocyte calls marked **lower confidence** in the final table (footnote: 4 donors, V9 24%, covariates partly track condition) | Flag, not a drop |
| O4: R1 `lodo_stable` = log2FC keeps its sign in every LODO refit of that cell type (significance not required) | As in the validation plan |
| O5: KEGG gene sets = MSigDB 2026.1 C2 KEGG legacy (`c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt`, sha256 `4a87c8d1260e637ec174f5960885519d46fa29da650e70d3a3641a4d41c4fa5e`) | Closest to the KEGG the paper's DAVID used |
| O6: M2 split into **M2a** (calls + table: all 8 steps + planted bug) and **M2b** (validation checks: report tier, steps 1, 4, 7, 8) | Checks scaled to risk |
| O7: k-means for the module check: k = 7, `random_state=0`, `n_init=25` | Paper's 7 modules; reproducible; robust to an unlucky start |
| M2b checks (user accepted Claude's recommendations): **P4** control mean = mean DESeq2-normalised count over the cell type's kept control samples; detection = fraction of that cell type's control cells (all cells of kept donors) with count > 0; matching on rank-percentiles of log control mean and detection (Euclidean), exactly 30 nearest per up-DEG from all tested genes, focal excluded, one draw per DEG without reuse within a draw (next-nearest unused if a pool runs out); statistic = number of IFN-α genes in the matched set vs among up-DEGs; p = (1 + #null ≥ observed)/1,001; seed 0. **F3b** log2(CPM + 1); includes the 12 samples dropped by A1 (paper's 8 donors per group); all-PBMC = sum over all cells; our DE genes, those missing from the authors' file shown as "no module"; z-score ddof 0. **KEGG** background = genes tested in ≥ 1 cell type; one-sided Fisher + BH (as DAVID). **P2** BH across the 8 cell types; one-sided "greater" headline. **I5** shared genes (incl. `shared*`) ranked by mean \|log2FC\|. **I6** per donor × cell type: median log2(stim CPM / ctrl CPM) over that cell type's DE genes vs log n_cells and log total_counts (mean of ctrl and stim), Spearman per cell type. **V9 correlation** reported on both gene sets, leading with every paper DEG we tested | Definitions the validation plan didn't fix |
| **V9b:** correlation can't separate unshrunk from shrunken log2FC (r ≥ 0.97 both), but the median \|paper − ours\| favours shrunken in most cell types (e.g. CD4 T 0.06 vs 0.14; megakaryocytes 0.45 vs 1.28). **Calls stay on unshrunk log2FC**; V9b is recorded as favouring shrunken *values*, which contradicts the earlier note "confirmed later by V9b" | Paper's DESeq2 version probably shrank by default; switching the DE rule to shrunken LFC would call fewer genes and widen the V9 gap |
| **CI rules (shared\*, specific\*) use only cell types where the gene was tested** (baseMean > 1); a gene with no CI (untested or all-zero) never triggers or rescues anything. Real data: shared 21 + shared\* 87 = **108**, partial 856, specific 241, specific\* 1,585, mixed 11. A missing LODO log2FC for a DE call stops the run | Untested CIs are uninformative (e.g. −6.8 to 5.1) and starred 1,823 of 1,826 specific genes |
| M2b follow-ups: **I6** uses log2(CPM + 1) (as F3b) and leads with the direction-aligned median (sign × median); x-axis = log of the ctrl/stim mean; **KEGG** genes get the 6 B4 renames; **P4** tie-breaking (stable sort, random DE order per draw, average ranks) confirmed | Zeros made the unpseudocounted ratio infinite |
| M2c contract and test plan (T1–T9) approved, **T9 real-data tripwire kept** | User |
| **M2c D2:** a program not tested in a cell type (set size outside 15–500 there) is shown as **`n/t`**, distinct from `–` (tested, not affected); it counts as not affected for labels. **Confirmed:** fgsea 1.26 defaults passed explicitly (`sampleSize=101`, `scoreType="std"`, `gseaParam=1`, `nPermSimple=1000`), `nproc=1` for reproducibility, `fgseaMultilevel` once per cell type on both collections, BH within cell type; average linkage for Jaccard and NES trees; cell-type clustering on programs tested in all 8 (count printed; display only); Jaccard of two empty sets = NaN; tied stats ordered by gene name | User accepted Claude's recommendations |
| **M2c T6 toy changed** (user): its 500-gene set is placed among ordinary ranks instead of all at the top, so fgsea (eps = 0) doesn't compute an extreme exact p-value (the run had stalled for > 2 h). The size-limit check (14 / 15 / 500 / 501) is unchanged; eps = 0 stays for real data | Test-only change; no effect on results |
| **Wave 4 F1–F4, F6 (user accepted recommendations, 2026-10-05):** F1 Fig 3 = per-cell violins (ctrl, stim) + donor means joined ctrl→stim; F2 Fig 3 row B = IFN-γ shown twice: B1 full Hallmark set, B2 without the 73 shared genes (revised by the user 2026-10-22: show both; Codex also recommended the full set as published); F3 (revised after seeing the figure) Fig 4 rows = top 25 shared (by mean \|log2FC\|) + top 5 specific per cell type **among genes tested in all 8 cell types** (`specific*`, marked \*, except FTL in megakaryocytes, a clean `specific`), by \|log2FC\|; fewer where fewer exist (megakaryocytes 2); F4 no gene-set filter on Fig 4 rows, title "top shared and cell-type-specific IFN-β-responsive genes"; F6 final table with M2a labels in "Shared response?", log2FC and q per cell type, megakaryocyte flag, `results/README.md` | Every cell type appears in Fig 4; uncertainty visible; every arrow traceable. F3 revised because 229 of 241 clean `specific` genes are expressed only in their own cell type (grey everywhere else): the user chose rows showing a cell-type-specific *response* in genes every cell type expresses |
| **Wave 4 approved and merged** (2026-10-06): Figs 3–5 and `results/final_table.csv`; all project modules complete | User |
| **Wave 4 follow-ups (user, 2026-10-06):** Fig 4 uses **two colour scales** (shared block ±11, specific block ±3, each with its own colour bar). The final table gets a column **`expressed_in_n_cell_types`** (cell types where the gene was tested, baseMean > 1) plus a README note that a `specific` gene with value 1 is specific because only that cell type expresses it; **M2a's labels are not changed** | Low-touch: makes specific-expression vs specific-response visible without reopening an approved module |
| **Wave 4, F5 (user, 2026-10-05): both groupings of cell types.** Fig 4 (gene heatmap) columns ordered by Jaccard clustering of up-regulated DE gene sets; **new Fig 5** (program heatmap, GSEA NES) columns ordered by NES-profile clustering. Fig 5 was **added after the original contract** (outputs now 5 figures + table; CLAUDE.md updated) | The two groupings answer different questions (shared genes vs shared programs) and disagree on DCs, NK and CD4 T; the user wants both shown |
| **M2c approved and merged** (2026-10-02); Wave 3 complete | User |
| **M2b approved and merged** (2026-10-02): every check in `audit/validation_report.md` reported; V9 fails in 6 of 8 (accepted, causes in FINDINGS); I6 shows larger per-donor effects in smaller samples (low-count caveat) | User |
| **M2a approved and merged** (2026-10-02). Codex review 1 (M2c): fix #1 `label_programs` always labels across all 8 cell types (a cell type with no testable sets shows `n/t`, never disappears), new test T10; fix #2 the notebook names a specific program's cell type from its arrow, not from `!= "–"` | User |
| **Codex review 1 (M2a):** fix #1 traceability rebuilds arrows, counts and labels with its own code (independent of `label_genes`), new test T11; fix #2 LODO runs must be exactly `lodo_<donor>` for each paired donor, new test T12 | User accepted |
| **M2c T9 reduced to its IFN-α half** (Hallmark IFN-α response: positive NES in all 8, padj < 0.05 in the 6 well-powered; passes, NES 2.8–3.0). The second half ("each top-20 up list contains one of the 10 canonical ISGs") **failed** in CD14+ monocytes and CD4 T cells and was removed, not widened: ranking by q favours consistent genes (CD4 T: ISG15 #22), and those lists are still IFN genes (ISG20, OASL, PLSCR1, IFITM2, BST2, HERC5) | User chose (a): the program-level check is the real tripwire; widening the gene list after seeing the result would move the bar |
| **M2c (shared programs), added by the user:** (G5) per cell type the **top 20 up and top 20 down** genes (by q, then \|log2FC\|) and which cell types share them, plus Jaccard clustering of the full DE lists; (G1) GSEA on MSigDB Hallmark + KEGG legacy; (G2) genes ranked by DESeq2 Wald statistic; (G3) R `fgsea`; (G4) a program is affected if BH-adjusted p < 0.05 within the cell type, labelled across all 8 cell types like genes (shared / partial / specific / mixed), and cell types clustered by NES; (G6) result tier with a small test set (~5). Reverses "threshold-free GSEA dropped" | User: common programs, not just common genes |

## Differential expression (follows Kang 2018)

| Decision | Reason |
|---|---|
| Pseudobulk: raw counts summed per donor × cell type × condition; R DESeq2, `~ donor + condition`, per cell type. Unpaired `~ condition` as a sensitivity check | Paper: "individuals as replicates"; paired removes donor variation |
| Keep genes with baseMean > 1 | Paper's "filtered out" is a typo: Supp Table 1's minimum baseMean is 1.00 |
| FDR by R `qvalue`. **Significant** = q < 0.05. **DE** = q < 0.05 and \|log2FC\| ≥ 1, up or down | Exactly the paper's rule: Supp Table 1 has 3,055 unique genes, all \|log2FC\| ≥ 1 ("logFC > 2" meant fold change > 2); its BH padj reaches 0.135, so selection used qvalue. **Changed by the user (2026-10-02) from \|log2FC\| > 1 to ≥ 1**: \|log2FC\| = 1 counts as DE; changes 0 calls on our data (no significant gene at exactly 1.000) |
| **Ribosomal genes stay in the DE test** like any other gene (no gene-class exclusions); check I5 only counts them among the top DE genes | User decision; % ribo is still excluded as a *covariate* (F6) |
| Unshrunk (MLE) log2FC | Supp Table 1: DE genes at baseMean 1–2 have median \|log2FC\| 1.95, which shrinkage would prevent. V9b later found shrunken *values* closer to the paper's (see M2 section); calls stay unshrunk |

## Cell-type specificity

| Decision | Reason |
|---|---|
| A gene's pattern = the cell types where it is DE (Kang 2018 Fig 3b). No formal between-cell-type test | Consistency with the paper |
| A "specific" call is flagged **specific\*** if, in any cell type where the gene is not DE, the 95% CI for log2FC reaches \|log2FC\| ≥ 1. Decided per gene. Final table uses ↑ / ↓ / – | Low power can fake specificity; errs towards extra asterisks |

## IFN signature (Fig 3)

| Decision | Reason |
|---|---|
| Signature = MSigDB Hallmark IFN-α response; Hallmark IFN-γ scored as a relative specificity control | IFN-β is type I (shares IFN-α's receptor and genes); IFN-γ's MHC-II genes track cell identity |

## Reproducing Kang 2018

| Decision | Reason |
|---|---|
| Recreate Fig 3a (authors' t-SNE coordinates, which the h5ad carries), 3b and 3d in a validation report. Fig 3f dropped | 3f's test isn't specified well enough |
| Fig 3b: pseudobulk log2 CPM per donor × cell type × condition, 9 groups incl. all PBMCs; genes DE in ≥ 1 cell type; rows z-scored; ordered by the authors' 7 modules | Faithful to the authors' heatmap file |
| Modules: authors' 7-cluster labels (`data/kang2018_code/final.gene.gct.txt`), plus our k-means (k = 7), compared by overlap table + adjusted Rand index, labelled "method inferred" | Method never stated; file format suggests k-means in Morpheus |
| Module biology: KEGG enrichment (Fisher) of our modules vs Supp Table 2, at least Cluster III | DAVID is web-only |
| V9: per cell type, % of the paper's DEGs (Supp Table 1) we also call DE, ✓ if > 70%; plus log2FC correlation | Gene-level reproduction; extra genes don't count against us |

## Validation (details in `validation_plan.md`)

| Decision | Reason |
|---|---|
| Positive controls: 10 canonical ISGs DE (bar chart per cell type + histogram, x = ISGs DE, y = proportion of cell types, no pass mark); Fisher enrichment of up-DEGs (q, odds ratio); up vs down enrichment; expression-matched null from all tested genes (1,000 draws, empirical p as a number, with diagnostic columns); IFN-α/γ scores as numbers | User prefers numbers and plots to pass marks |
| Negative control: all 35 balanced within-donor label shuffles, DE counts plotted next to the real count, no pass mark | The only false-positive check |
| Robustness: leave-one-donor-out; genes in the final table should keep direction in all 8 refits | Single-donor effects |
| Integrity: cell accounting (counts matches, not rows); heatmap/table traceability; same cells in both UMAP panels + authors' t-SNE vs their Fig 3a; simulated-data oracle test (Codex, DE module); number of ribosomal genes in the top 50 shared DEGs; covariate correlation plots, with/without-covariate comparison only if M1a adds one | Catch silent errors |
| Dropped: Fig 3f, edgeR, shrinkage on/off, enrichment on shuffles, control-vs-control splits, split-donor replication, annotation check (authors' labels used), GSEA, ΔLFC specificity test | Out of scope (reasons in git history) |
