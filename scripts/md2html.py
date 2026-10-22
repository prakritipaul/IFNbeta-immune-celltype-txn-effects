"""Render a review markdown file as a readable HTML page next to it.

Usage: .venv/bin/python scripts/md2html.py audit/reviews/<file>.md
Relative links like [src/m0_qc.py:95](../../src/m0_qc.py#L95) become vscode://file/...:95,
so clicking one opens the file at that line in VS Code.
"""
import re
import sys
from pathlib import Path

import markdown

CSS = """
:root { --bg:#fcfcfb; --ink:#2b2b28; --muted:#5c5c57; --rule:#e6e6e3; --accent:#2a78d6; --code:#f3f3f1; }
@media (prefers-color-scheme: dark) { :root { --bg:#1a1a19; --ink:#f2f2ef; --muted:#c3c2b7; --rule:#34342f; --accent:#6da7ec; --code:#262624; } }
body { background:var(--bg); color:var(--ink); font:16px/1.6 -apple-system, "Segoe UI", sans-serif; max-width:52rem; margin:2rem auto; padding:0 16px; }
h1 { font-size:1.6rem; border-bottom:1px solid var(--rule); padding-bottom:.4rem; }
h2 { font-size:1.2rem; margin-top:2rem; color:var(--accent); }
a { color:var(--accent); } li { margin:.35rem 0; }
code { background:var(--code); padding:.1em .35em; border-radius:4px; font-size:.9em; }
table { border-collapse:collapse; margin:1rem 0; } th, td { border:1px solid var(--rule); padding:.35rem .6rem; text-align:left; vertical-align:top; }
nav { color:var(--muted); font-size:.9rem; border:1px solid var(--rule); border-radius:6px; padding:.5rem 1rem; }
del { color:var(--muted); }
"""


def vscode_links(text, md_path):
    def repl(m):
        target, line = m.group(1), m.group(2)
        abs_path = (md_path.parent / target).resolve()
        return f"(vscode://file{abs_path}{':' + line if line else ''})"
    return re.sub(r"(?<=\])\((?!https?:|vscode:)([^()\[\]#\s]+\.(?:py|md|ipynb|R))(?:#L(\d+))?\)", repl, text)


def main(md_file):
    md_path = Path(md_file).resolve()
    text = vscode_links(md_path.read_text(), md_path)
    text = re.sub(r"~~(.+?)~~", r"<del>\1</del>", text)
    md = markdown.Markdown(extensions=["tables", "toc", "fenced_code"])
    body = md.convert(text)
    title = next((l[2:].strip() for l in text.splitlines() if l.startswith("# ")), md_path.stem)
    html = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width, initial-scale=1'><title>{title}</title>"
            f"<style>{CSS}</style></head><body><nav><strong>Contents</strong>{md.toc}</nav>{body}</body></html>")
    out = md_path.with_suffix(".html")
    out.write_text(html)
    print(out)


if __name__ == "__main__":
    for f in sys.argv[1:]:
        main(f)
