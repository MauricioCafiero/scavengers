#!/usr/bin/env python
"""Render paper_info/METHODS_OUTLINE.pdf from paper_info/METHODS_OUTLINE.md.

Markdown -> HTML (python-markdown, tables + fenced code) -> paginated PDF via PyMuPDF's
Story API. Run with the repo venv (.venv/bin/python); pymupdf and markdown are the only
requirements beyond the markdown source, which the venv reads nothing else from.

    .venv/bin/python code/render_methods_pdf.py [md_path [pdf_path]]

Defaults to the METHODS_OUTLINE pair in paper_info/. The PDF is a rendering of the md —
regenerate it after every edit to the md, or it rots silently.
"""
import sys
from pathlib import Path

import fitz  # PyMuPDF
import markdown

REPO = Path(__file__).resolve().parent.parent

CSS = """
body { font-family: sans-serif; font-size: 9.5px; line-height: 1.45; }
h1 { font-size: 17px; margin-bottom: 4px; }
h2 { font-size: 13px; margin-top: 14px; border-bottom: 0.5px solid #999; }
h3 { font-size: 10.5px; margin-top: 8px; }
code, pre { font-family: monospace; font-size: 8.5px; background-color: #f2f2f2; }
pre { padding: 4px; }
table { border-collapse: collapse; margin: 6px 0; }
th, td { border: 0.5px solid #888; padding: 2px 5px; font-size: 8.5px; text-align: left; }
li { margin-bottom: 3px; }
"""


def main() -> None:
    md_path = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "paper_info" / "METHODS_OUTLINE.md"
    pdf_path = Path(sys.argv[2]) if len(sys.argv) > 2 else md_path.with_suffix(".pdf")

    text = md_path.read_text()
    html = markdown.markdown(text, extensions=["tables", "fenced_code"])

    story = fitz.Story(html=html, user_css=CSS)
    writer = fitz.DocumentWriter(str(pdf_path))
    page_rect = fitz.paper_rect("a4")
    body = page_rect + (40, 42, -40, -40)
    more = 1
    while more:
        dev = writer.begin_page(page_rect)
        more, _ = story.place(body)
        story.draw(dev)
        writer.end_page()
    writer.close()

    doc = fitz.open(str(pdf_path))
    doc.set_metadata({
        "title": "Methods outline — fragment-based peptide design",
        "author": "M. Cafiero",
        "producer": "PyMuPDF (code/render_methods_pdf.py)",
    })
    doc.saveIncr()
    print(f"{pdf_path}: {doc.page_count} pages, {pdf_path.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()