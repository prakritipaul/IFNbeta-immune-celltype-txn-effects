# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**Question:** Which transcriptional responses to IFN-β are shared across PBMC cell types, and which are cell-type specific?

**Data:** Kang et al. 2018, PBMCs from 8 donors, control vs IFN-β. Load `data/kang_2018.h5ad` with `anndata.read_h5ad` (the file `pertpy.data.kang_2018()` downloads, from `https://scverse-exampledata.s3.eu-west-1.amazonaws.com/pertpy/kang_2018.h5ad`; `import pertpy` 1.0.3 fails with statsmodels 0.15). **Cell-type labels come from the authors' GEO metadata** (`data/geo/GSE96583_batch2.total.tsne.df.tsv.gz`, column `cell`), not the h5ad's `cell_type` (a later re-analysis).

**Join trap:** 313 barcode sequences occur in both the control and stimulated runs (all 29,065 GEO droplets; 234 among the 24,673 singlets we use). The h5ad suffixes all stimulated barcodes `-2`; GEO keeps `-1` and renames the 313 colliding stimulated ones `-11`. Join on (barcode without suffix, condition), never on the raw barcode: a raw-barcode join silently matches only control cells. Expected result: all 24,673 cells match, donor/cluster/t-SNE agree 100%. Before writing code that touches the data, read `audit/data_facts.md` (written by the EDA module) for column names, label spellings and group sizes. Don't guess them.

**Sanity check:** the response must be dominated by classic ISGs (ISG15, IFIT family, MX1, OAS family, STAT1 signalling). If ribosomal or mitochondrial genes dominate, something is wrong, so investigate before continuing.

**Outputs:** only these.
1. Fig 1: cells by condition / donor / cell type.
2. Fig 2: UMAP coloured by condition and by cell type.
3. Fig 3: IFN signature score, control vs stimulated, per cell type.
4. Fig 4: DE heatmap of the top shared and cell-type-specific IFN-response genes; cell types grouped by their shared DE genes (Jaccard clustering).
5. Fig 5: heatmap of response programs (GSEA enrichment scores, Hallmark + KEGG) across cell types; cell types grouped by their program profiles. *Added by the user after the original contract (2026-10-05), together with the M2c "shared programs" analysis (added 2026-10-02): the user wanted to see common programs, not only common genes, and both ways of grouping cell types.*
6. Table: gene × cell type (↑/↓/–) plus a "Shared response?" column.

**Decisions already made:**
- DE follows Kang 2018: pseudobulk (raw counts summed per donor × cell type × condition), R DESeq2 with `~ donor + condition`, R `qvalue` for FDR (both in R, as in the paper). Significant = q < 0.05; DE = q < 0.05 and |log2FC| ≥ 1. Cells are not replicates.
- CD14+ and FCGR3A+ monocytes stay separate cell types.
- The IFN signature in Fig 3 uses a fixed published gene set (MSigDB Hallmark IFN-α response), not our own DE genes, to avoid circularity. Hallmark IFN-γ is scored as a specificity control. Full decision log: `audit/decisions.md`.

## Environment

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv --only-binary llvmlite,numba -r requirements.txt
.venv/bin/python -m pytest tests/                       # all tests
.venv/bin/python -m pytest tests/test_x.py::test_name   # one test
```

DE runs in R 4.3 (`/usr/local/bin/Rscript`) with Bioconductor DESeq2 1.40.2, qvalue 2.32.0 and variancePartition 1.30.2 and fgsea 1.26.0:
`Rscript -e 'BiocManager::install(c("DESeq2", "qvalue", "variancePartition", "fgsea"))'`. Python writes pseudobulk counts and sample metadata to files, R reads them and writes results to files.

`--only-binary` is required: this is an Intel Mac, and recent `llvmlite` has no x86_64 macOS wheel. Without the flag it tries to build from source and fails for lack of CMake.

## How we work

The user's rules, for Claude and Codex alike. Keep every message to the user short: the user's attention is the scarcest resource, so spend it on decisions, not on reading.

0. **The user decides; agents propose.** Ask *before* coding any choice that could change a number, a gene list, what a figure shows, or what is included. That covers thresholds, filters, methods, model design, small or missing groups, and multiple-testing scope. It also covers:
   - **Error workarounds:** `dropna`, `fillna`, `errors="coerce"`, `try/except`, silenced warnings, skipped groups. When code fails, stop and report the cause. Don't make it run.
   - **Library defaults:** pass every parameter that affects results explicitly. Never rely on a default without having asked.

   Subagents return questions to the coordinator; they never pick an answer. There are no defaults "for now". Only cosmetic choices (names, colours, file layout) are exempt.
1. **Validation before building.** Controls and pass criteria are agreed in `audit/validation_plan.md` before results exist. They are never changed afterwards without asking.
2. **EDA first, and visual.** `notebooks/00_eda.ipynb` covers the EDA checklist below. The facts go in `audit/data_facts.md`, which agents read instead of guessing.
3. **Cell accounting.** Every filter, join or aggregation logs n in → n out and asserts what must hold. For example: a join keeps the cell count, and pseudobulk cell counts sum to the total. The running table goes in `audit/cell_flow.md`.
4. **Contracts before parallel work.** Each module has a contract of one page or less: inputs, output file + schema, assumptions, validation checks, open decisions. The user approves it before dispatch. One worktree per agent, and modules connect only through the agreed output files. At most about 3 agents per wave.
5. **Return package per module**, reviewed by the user at a gate after each wave before merging:
   - **Reviews folder:** `audit/reviews/` keeps only the latest Codex reviews, numbered in the order generated (`001_<file>`, `002_<file>`, …); older ones live in git history.
   - **Reviewer brief** written by **Codex** (a different model from the code's author), always also rendered to HTML for the user (`.venv/bin/python scripts/md2html.py <brief>.md`, which writes `<brief>.html` with the same headers and file links that open in VS Code at the line): assumptions coded in (with `file:line`), any unapproved decision (should be none; any entry is a process failure), the 3–5 lines worth reading, and what nothing checks.
   - **Explained by example** (`notebooks/<module>_explained.py`, executed): for each line worth reading, a tiny input, the **expected** output written by a **Claude subagent** (not the author) from the contract/decisions *before* running, and what the **real function gives**, with an `assert` that they match. A mismatch is reported as a possible bug; the expected value is never edited to match. This lets the user verify behaviour without reading the syntax.
   - **Outputs shown:** shape plus head and tail of every output file.
   - **Figure captions:** every file saved in `figures/` gets a caption in `figures/CAPTIONS.md`, written like a paper's figure legend: what is plotted, the input, the method and its parameters, n. No interpretation; the figure is embedded in `FINDINGS.md` with a short caption, and its main takeaway goes in a row of the bucketed "Figures at a glance" tables at the top of `FINDINGS.md`.
   - **Findings:** the key takeaways of the module's analysis notebook go in `FINDINGS.md` under the matching section (EDA, QC + Analysis Sanity Checks · Concordance with published paper · Validation Checks · Biological Findings), in bullets: why, method, results, interpretation, with ⚠️ **Analysis caveat** blocks for confounders and limits (numbers, not process; explained-by-example notebooks excluded).
   - **Checks:** validation checks marked pass/fail, 1–2 diagnostic plots, and Codex audit results.
6. **Codex audit, user-designed tests.** Follow the `continuous-independent-audit` skill. The user approves 5–10 tests in `audit/test_plans/<module>.md` before any are written. Status goes in `audit/ledger.md`.

**How to ask** (so the user doesn't have to rubber-stamp):
- **Batch questions** into one message per gate, not a stream of pings.
- **High-stakes first.** Mark each decision high- or low-stakes, and lead with the ones that change results.
- **One compact card per decision:**

  | Option | What it does | Who uses it, and why | Effect here |
  |---|---|---|---|

  End each card with a recommendation and one line of reasoning. For a library parameter, the rows are its default and the common alternatives.
- **Record each answer** in `audit/decisions.md` (decision, reason) and never re-ask. Keep only current decisions: when one changes, replace its row rather than appending a contradicting one. History lives in git.

**Module workflow** (every module, in this order). **Each check is done by a different model, or at least a different context, from the step it checks.**
1. **The user approves the contract.**
2. **Test plan:** Codex proposes blind, then Claude, then the user; the user approves 5–10 tests.
3. **Codex writes the tests**, before any code exists.
4. **Claude writes the code** and runs the tests. Diagnose each failure (wrong code, wrong test or wrong requirement?); decisions go to the user.
5. **Codex reviews the code.** Each fix gets a new test that fails on the old code and passes on the new; decisions go to the user.
6. **Codex re-checks** the fixed code: re-runs its reproducers and the full suite.
7. **On the final code and notebooks:** **Codex** writes the reviewer brief; a **Claude subagent** writes the explained-by-example notebook.
8. **The user** reads the brief, examples and outputs, makes the open decisions, and approves the module.

**Scale to risk:**
- **Produces results** (e.g. pseudobulk, DE, shared/specific calls): all 8 steps, plus a planted-bug check.
- **Displays results** (figures): steps 1, 4, 7 (explained-by-example only) and 8.

**Planted-bug check** (result-producing modules, after step 6): Claude introduces one small, realistic bug on a throwaway branch (e.g. `>= min_cells` → `> min_cells`) without telling the reviewers, then runs the tests and asks Codex to review that snapshot.
- Caught → the checks have teeth.
- Missed → a gap to close, with a new test.

The bug never reaches `main`; the outcome is recorded in `audit/planted_bugs.md`.

**Module plan:**

```
Wave 0  M0 EDA (incl. covariate vs condition/donor/cell type confounding) → user review
Wave 1  M1a pseudobulk + variance partitioning | M3 IFN signature | Fig 1 + Fig 2   (all need only M0)
Wave 2  M1b pseudobulk DE (covariates chosen from M1a)
Wave 3  M2a shared/specific calls | M2b validation checks | M2c shared programs (added by the user mid-project)
Wave 4  Fig 3, Fig 4, Fig 5 (added by the user after the contract), final table
```

## EDA checklist (M0)

**Structure**
- Head and tail of `obs` (cell metadata) and `var` (gene metadata). dtypes, NA counts and duplicates per column.
- Number of unique cell barcodes, donors, genes and cell types. Are barcodes unique across control and stimulated? The same barcode can appear in both lanes.
- A list of every metadata column: what it holds, its unique values or range, and a plot suited to it (bar plot for categorical, histogram for numeric).

**Counts matrix**
- Does `X` hold raw integer counts, or is it normalized? Is there a `.raw` or a counts layer? Pseudobulk needs raw counts.
- Gene IDs: symbols or Ensembl? Any duplicated symbols? Are ISG15, MX1, IFIT1/2/3, OAS1/2/3 and STAT1 present under those names? Mito (`MT-`) and ribosomal (`RPS`/`RPL`) genes present?

**Filtering** (in one place, reporting how many cells and genes each step removes)
- Valid cells: expressing ≥ 200 genes. Valid genes: expressed in ≥ 3 cells.
- Show the distributions before and after, with the cutoffs drawn on the plots.

**Distributions and outliers**, by donor and by condition:
- Cells per donor, as a stacked bar plot split by condition and by cell type.
- Genes detected per donor: for each donor, two side-by-side violins (control vs stimulated).
- Per-cell distributions (total counts, genes detected, % ribo) as violins, with outliers flagged. (No % mito: the h5ad has no `MT-` genes.)
- Library size, control vs stimulated. A systematic difference matters for normalization.

**Design checks**
- Does every donor have both conditions? This must hold for the paired design.
- Donor × condition × cell type cell-count table, flagging groups too small to pseudobulk (e.g. < 10 cells).
- Doublet or ambiguous-donor cells still present? Kang 2018 assigned donors with demuxlet.
- Were cell types annotated on the pooled data? The IFN response can shift clustering, so check that each label's marker genes hold in both conditions.

## QC flags (end of M0)

The EDA ends with `audit/qc_flags.md`: one table of everything concerning, each row with **what, evidence (number + plot), proposed action (keep / drop / caveat), and why**. Agents flag and recommend; **only the user drops cells, samples or cell types**. Decisions go in `audit/decisions.md`.

Flag at least the following. The thresholds are starting proposals for the user to confirm:
- **Outlier cells:** more than 5 raw (unscaled) MADs from the median on log total counts or log genes detected, computed per donor × condition, since stimulation shifts these.
- **Underpowered groups:** donor × condition × cell type with < 10 cells. They make an unreliable pseudobulk sample.
- **Underpowered cell types:** fewer than 3 donors with both conditions left after the rule above. Too few pairs for paired DE, so the cell type is dropped from the table or clearly marked.
- **Weak donors:** a donor with far fewer cells than the others overall or in most cell types, or one missing a condition.
- **Outlier samples:** a pseudobulk sample far from the rest of its cell type × condition in PCA.
- **Covariate trends:** a technical covariate that differs by condition, donor or cell type (plot it, and give an effect size or correlation, not just a p-value), and any covariate perfectly confounded with condition.
- **Composition shifts:** cell-type proportions that change sharply between control and stimulated within donors. This may be biology or an annotation artifact.
- **Annotation doubts:** a cell type whose marker genes don't hold in both conditions.

## Technical covariates

The EDA lists every technical covariate that exists in the data, and records which of the following are missing:
- **Sequencing lane / batch / 10x run.** If control and stimulated were run in separate lanes, condition is confounded with lane.
- **Library size** (total counts per cell) and **genes detected per cell.**
- **% ribosomal** reads. (% mitochondrial is unavailable: the h5ad has no `MT-` genes.)
- **Doublet / demuxlet status.**
- At pseudobulk level: **cells per sample** and **total counts per sample.**

How they are used:
- **Check each covariate against condition, donor and cell type first.** A covariate perfectly confounded with condition cannot go in the same model: it would absorb the effect, or make the model unidentifiable. Report it as a caveat instead.
- **Variance partitioning** at pseudobulk level: for each gene, the fraction of variance explained by condition, cell type, donor and the technical covariates.
- **Correlation of results with covariates:** check whether DE effect sizes or signature scores track a technical covariate (e.g. score vs library size or % ribo). If they do, the "biology" may be technical.
- Include non-confounded covariates in the DE model only if they explain meaningful variance. Record the choice and its reason in the module's contract.
