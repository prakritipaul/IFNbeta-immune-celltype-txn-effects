# audit/

The paper trail of the workflow: what the data contains, what the user decided, what each module promised, how it was tested and reviewed, and what the checks found. Agents read the first three files before writing code.

**Read before coding**

| File | What it holds | Written by |
|---|---|---|
| `data_facts.md` | Exact column names, label spellings, group sizes, which genes exist | M0 EDA notebook |
| `decisions.md` | Every current user decision, with its reason; replaced, never contradicted | Updated at each gate |
| `qc_flags.md` (+ `qc_flag_numbers.md`) | Everything concerning in the data, each with evidence and a proposed action; only the user drops anything | M0 |

**Per module**

| Folder / file | What it holds |
|---|---|
| `contracts/` | One page per module, approved by the user before coding: inputs, output files and schema, assumptions, checks, open decisions |
| `test_plans/` | The 5–10 tests the user approved. Codex proposed blind, then Claude added, then the user chose |
| `reviews/` | Codex output, numbered in the order generated (`0NN_<module>_codex_review`, `_codex_planted_review`, `_codex_reviewer_brief`), each as `.md` and rendered `.html`. Only the latest are kept; older ones are in git history |
| `m3_gene_sets.md` | The MSigDB gene sets M3 used and their symbol mappings |

**Checks and status**

| File | What it holds |
|---|---|
| `validation_plan.md` | Controls and pass criteria, agreed before any results existed |
| `validation_report.md` | Each check's result (M2b) |
| `cell_flow.md` | Cells in → cells out at every filter, join and aggregation |
| `ledger.md` | Audit status of every behaviour: which tests and reviews cover it |
| `planted_bugs.md` | Bugs planted on purpose, whether tests and Codex caught them, and the test added for each miss |
