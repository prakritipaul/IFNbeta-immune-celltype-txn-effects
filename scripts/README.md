# scripts/

| Script | What it does |
|---|---|
| `md2html.py` | Renders a Markdown file (e.g. a Codex reviewer brief) as HTML next to it. Relative `file:line` links become `vscode://` links that open the file at that line. Usage: `.venv/bin/python scripts/md2html.py audit/reviews/<file>.md` |
