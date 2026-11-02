#!/usr/bin/env python
"""Colour-coded PDF of the accounting workbook's Pearson correlation matrices.

Reads paper_info/PB_runs_accounting_corr_<mol>.csv (+ the _n.csv with the
pairwise-complete sample counts per cell) and writes
paper_info/correlation_matrices.pdf — one page per molecule, each page carrying
the correlation heatmap (values annotated, diverging blue-gray-red on a fixed
-1..+1 scale) and a small companion heatmap of the cell n's (sequential blue).

    .venv/bin/python code/plot_corr_matrices.py
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, Normalize

REPO = Path(__file__).resolve().parent.parent
INFO = REPO / "paper_info"
OUT = INFO / "correlation_matrices.pdf"

# Diverging pair: blue <-> red with neutral gray midpoint (light surface).
BLUE, MID, RED = "#2a78d6", "#f0efec", "#e34948"
INK, INK_LIGHT = "#0b0b0b", "#ffffff"

DIV_CMAP = LinearSegmentedColormap.from_list("bluered", [BLUE, MID, RED])
SEQ_CMAP = LinearSegmentedColormap.from_list("seqblue", ["#e7e9f2", BLUE])


def wrap(label: str) -> str:
    words = label.split()
    lines, cur = [], ""
    for w in words:
        if cur and len(cur) + 1 + len(w) > 12:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return "\n".join(lines)


def read_matrix(path: Path):
    with path.open() as fh:
        rows = list(csv.reader(fh))
    labels = [c.strip() for c in rows[0][1:]]
    data = np.array(
        [[float(c) if c.strip() else np.nan for c in r[1:]] for r in rows[1:]]
    )
    return labels, data


def lum(hexcolor: str) -> float:
    r, g, b = (int(hexcolor[i : i + 2], 16) / 255 for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def draw(ax, data, labels, cmap, norm, fmt="{:.2f}", xrot=0):
    im = ax.imshow(data, cmap=cmap, norm=norm)
    n_rows, n_cols = data.shape
    # thin white gap between cells
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    for i in range(n_rows):
        for j in range(n_cols):
            if np.isnan(data[i, j]):
                continue
            fill = matplotlib.colors.to_hex(cmap(norm(data[i, j])))
            ink = INK_LIGHT if lum(fill) < 0.45 else INK
            if i == j:  # trivial self-correlation / same-column n
                ax.text(j, i, "—", ha="center", va="center", fontsize=6, color=ink)
            else:
                ax.text(j, i, fmt.format(data[i, j]), ha="center", va="center",
                        fontsize=7, color=ink)
    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_xticklabels([wrap(l) for l in labels], fontsize=7,
                       rotation=xrot, ha="left" if xrot else "center")
    ax.set_yticklabels([wrap(l) for l in labels], fontsize=7)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def main() -> None:
    with PdfPages(OUT) as pdf:
        for mol in ("octinoxate", "oxybenzone"):
            labels, r = read_matrix(INFO / f"PB_runs_accounting_corr_{mol}.csv")
            _, n = read_matrix(INFO / f"PB_runs_accounting_corr_{mol}_n.csv")
            fig, axes = plt.subplots(
                1, 2, figsize=(10.5, 6.2),
                gridspec_kw={"width_ratios": [1.5, 1]}, layout="constrained",
            )
            fig.suptitle(f"{mol.capitalize()} — Pearson correlation matrix "
                         f"(PB_runs_accounting.xlsx, {mol.capitalize()} sheet)",
                         fontsize=12)
            norm_r = Normalize(vmin=-1, vmax=1)
            im = draw(axes[0], r, labels, DIV_CMAP, norm_r)
            cb = fig.colorbar(im, ax=axes[0], shrink=0.8, pad=0.02)
            cb.set_label("Pearson r", fontsize=8)
            cb.ax.tick_params(labelsize=7)
            axes[0].set_title("r", fontsize=10, pad=28)

            norm_n = Normalize(vmin=0, vmax=max(120, np.nanmax(n)))
            im_n = draw(axes[1], n, labels, SEQ_CMAP, norm_n, fmt="{:.0f}",
                        xrot=45)
            cb2 = fig.colorbar(im_n, ax=axes[1], shrink=0.8, pad=0.02)
            cb2.set_label("n per cell", fontsize=8)
            cb2.ax.tick_params(labelsize=7)
            axes[1].set_title("pairwise-complete n", fontsize=10, pad=28)

            pdf.savefig(fig)
            plt.close(fig)
    print(f"{OUT}: {OUT.stat().st_size / 1e3:.0f} kB")


if __name__ == "__main__":
    main()