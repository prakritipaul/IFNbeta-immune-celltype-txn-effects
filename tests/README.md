# tests/

99 pytest tests. Run them with `.venv/bin/python -m pytest tests/`.

| Files | Written by | When |
|---|---|---|
| `test_<module>.py` | **Codex, blind**, from the contract and the test plan the user approved (`audit/test_plans/<module>.md`) | Before any code for that module existed |
| `test_<module>_review.py` | Added after Codex's code review, or after a planted bug the suite missed | Each test fails on the old (or bugged) code and passes on the fixed code |
| `test_m0_decisions.py` | Codex: tests for the user's M0 QC decisions, input guards, and the planted-bug gap (T18) | After the M0 gate and planted-bug check |
| `conftest.py` | Puts the repo on the import path; defines the `realdata` and `r` markers | |

Modules: `m0` (EDA, join, QC), `m1a` (pseudobulk), `m3` (signature scores), `m1b` (DE), `m2a` (shared/specific calls), `m2c` (programs).

Display-only modules (figures, final table) have no tests here. They are checked by their notebook asserts and their `*_explained` notebooks.

Tests marked `realdata` read the downloaded Kang 2018 files, and tests marked `r` run an R script. To run only the fast synthetic tests, use `pytest -m "not realdata and not r"`. Which behaviour each test covers is listed in `audit/ledger.md`.
