"""Render runs/folds_md/cofold_binding_tables.md to a PDF, one table per page.

The PDF that sat beside that markdown until 2026-10-19 was built by matplotlib (its metadata says
so) by a script that was never committed, so it went four days stale across a restructure and a
whole campaign with no way to rebuild it. This is that script, committed: the markdown is the
source of truth and the PDF is derived from it, so regenerating is `python code/render_cofold_pdf.py`
rather than an afternoon of remembering.

Cells keep their unicode (−, •, ↑, ↓) and lose their markdown: `**bold**` becomes a bold cell and
`` `code` `` loses its backticks, because a PDF table has no use for the markup.

Needs matplotlib, which the bare `python3` on this machine does not have -- use the repo venv (or the
conda openmm-md env, which also carries it). The plain `python code/render_cofold_pdf.py` in an
earlier version of this docstring fails with ModuleNotFoundError.

    .venv/bin/python code/render_cofold_pdf.py       # -> runs/folds_md/cofold_binding_tables.pdf
    .venv/bin/python code/render_cofold_pdf.py --out other.pdf
"""
import argparse
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "runs", "folds_md", "cofold_binding_tables.md")


def tables(path):
    """Yield (caption, header, rows, legend) per markdown table.

    The legend is the prose that FOLLOWS the table, which is where this document defines what the
    numbers are -- "Each cell is enclosed fraction / ligand-centroid separation" under Table 3,
    the ΔG/residence and •/↑/↓ key under Tables 1 and 2. A first version of this script took only
    the caption above each table, so every legend was dropped and the PDF showed bare numbers.
    Two paragraphs are enough to carry the key and the caveat under it.
    """
    lines = open(path).read().split("\n")
    caption, out = None, []
    i = 0
    while i < len(lines):
        if lines[i].startswith("**Table "):
            # Only the bolded span is the caption; the same line often runs on into prose
            # ("...treatments.** Added 2026-10-18. Tables 1 and 2 set..."), which would otherwise
            # be chopped mid-sentence at whatever the line break happened to be.
            m = re.match(r"\*\*(.+?)\*\*", lines[i])
            caption = m.group(1).strip() if m else re.sub(r"\*\*", "", lines[i]).strip()
        if lines[i].startswith("|") and i + 1 < len(lines) and set(lines[i + 1].replace("|", "").replace(" ", "")) <= {"-", ":"}:
            header = [c.strip() for c in lines[i].strip().strip("|").split("|")]
            i += 2
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            # prose paragraphs after the table, up to the next table or heading
            paras, cur, j = [], [], i
            while j < len(lines) and len(paras) < 2:
                ln = lines[j]
                if ln.startswith(("|", "#")) or ln.startswith("**Table "):
                    break
                if ln.strip():
                    cur.append(ln.strip())
                elif cur:
                    paras.append(" ".join(cur)); cur = []
                j += 1
            if cur and len(paras) < 2:
                paras.append(" ".join(cur))
            out.append((caption, header, rows, "\n".join(paras)))
            caption = None
            continue
        i += 1
    return out


def clean(cell):
    """(text, bold) -- markdown stripped, unicode kept."""
    bold = "**" in cell
    text = cell.replace("**", "").replace("`", "").replace("\\*", "*")
    return text, bold


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", default=os.path.join(REPO, "runs", "folds_md", "cofold_binding_tables.pdf"))
    p.add_argument("--src", default=SRC)
    a = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    import textwrap

    tabs = tables(a.src)
    print(f"{len(tabs)} tables in {os.path.relpath(a.src, REPO)}")
    with PdfPages(a.out) as pdf:
        for caption, header, rows, legend in tabs:
            legend_lines = []
            for para in legend.split("\n"):
                legend_lines += textwrap.wrap(clean(para)[0], width=168) + [""]
            while legend_lines and legend_lines[-1] == "":
                legend_lines.pop()
            # One landscape page per table; height follows the row count plus the legend beneath it.
            fig_h = max(4.0, 0.9 + 0.205 * (len(rows) + 1) + 0.13 * len(legend_lines))
            fig, ax = plt.subplots(figsize=(11.7, fig_h))
            ax.axis("off")
            if caption:
                ax.set_title(caption, fontsize=7.2, wrap=True, loc="left", pad=14)

            cells = [[clean(c)[0] for c in r] for r in rows]
            tbl = ax.table(cellText=cells, colLabels=[clean(h)[0] for h in header],
                           cellLoc="left", colLoc="left", loc="upper center")
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(6.4)
            tbl.scale(1, 1.18)

            for (ri, ci), cell in tbl.get_celld().items():
                cell.set_linewidth(0.3)
                cell.set_edgecolor("#cccccc")
                if ri == 0:
                    cell.set_text_props(weight="bold")
                    cell.set_facecolor("#ececec")
                else:
                    raw = rows[ri - 1][ci] if ci < len(rows[ri - 1]) else ""
                    if clean(raw)[1]:
                        cell.set_text_props(weight="bold")
                    # a peptide's first sub-row starts a group: rule above it to keep the
                    # three-pose blocks legible, which is the whole point of the restructure
                    if ci == 0 and raw.strip():
                        cell.set_edgecolor("#999999")
            tbl.auto_set_column_width(col=list(range(len(header))))

            if legend_lines:
                # Place the legend just under the table by MEASURING where the table actually
                # ends. Estimating it from the row count left an inch of whitespace on the long
                # tables, because the drawn row height is not the nominal one once the table has
                # been scaled and its columns auto-sized.
                fig.canvas.draw()
                bbox = tbl.get_window_extent(fig.canvas.get_renderer())
                y0 = bbox.transformed(fig.transFigure.inverted()).y0
                fig.text(0.012, max(0.03, y0 - 0.22 / fig_h), "\n".join(legend_lines),
                         fontsize=6.0, va="top", ha="left", linespacing=1.45, color="#222222")

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"  page: {len(rows)} rows x {len(header)} cols, "
                  f"{len(legend_lines)} legend lines  {(caption or '')[:50]}")
    print(f"wrote {os.path.relpath(a.out, REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
