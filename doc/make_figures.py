#!/usr/bin/env python3
"""Draw the two figures of the user guide from the runs of doc/run_examples.sh.

    python3 doc/make_figures.py

reads   doc/out/eijgenraam/E1.txt         (examples/fig1Eijgenraam.mk)
        doc/out/eg1_eps0.1/E0.txt, E1.txt  (Eg1 at eps = 0.1, H = 4)
writes  doc/fig/eijgenraam.pdf, doc/fig/eg1_cover.pdf

Requires matplotlib.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from boxio import read_boxes  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

FILL = "#3987e5"     # boxes
EDGE = "#1c5cab"
INK = "#0b0b0b"      # the exact set
INK2 = "#52514e"     # the hull B1
MUTED = "#898781"    # axes

plt.rcParams.update({
    "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def hull(boxes):
    n = len(boxes[0])
    return [(min(b[i][0] for b in boxes), max(b[i][1] for b in boxes)) for i in range(n)]


def draw_boxes(ax, boxes):
    for b in boxes:
        (x0, x1), (y0, y1) = b[0], b[1]
        ax.add_patch(Rectangle((x0, y0), max(x1 - x0, 1e-9), max(y1 - y0, 1e-9),
                               facecolor=FILL, alpha=0.28, edgecolor=EDGE,
                               linewidth=0.35))


def draw_hull(ax, B):
    (x0, x1), (y0, y1) = B[0], B[1]
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                           edgecolor=INK2, linewidth=0.8, linestyle=(0, (4, 2))))


def frame(ax, B, frac=0.06):
    (x0, x1), (y0, y1) = B[0], B[1]
    dx, dy = (x1 - x0) * frac or 0.05, (y1 - y0) * frac or 0.05
    ax.set_xlim(x0 - dx, x1 + dx)
    ax.set_ylim(y0 - dy, y1 + dy)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def eijgenraam(outdir):
    C = read_boxes(os.path.join(HERE, "out", "eijgenraam", "E1.txt"))
    B1 = hull(C)
    fig, ax = plt.subplots(figsize=(4.0, 2.6))
    draw_boxes(ax, C)
    draw_hull(ax, B1)
    xs = [-1 + 2 * k / 400 for k in range(401)]
    ax.plot(xs, [x * x for x in xs], color=INK, linewidth=0.9)
    ax.plot([0], [1], marker="o", markersize=4, markerfacecolor="white",
            markeredgecolor=INK, markeredgewidth=0.8)
    ax.annotate("(0, 1)", (0, 1), xytext=(6, -10), textcoords="offset points",
                color=INK, fontsize=7.5)
    frame(ax, B1)
    fig.savefig(os.path.join(outdir, "eijgenraam.pdf"))
    print("wrote doc/fig/eijgenraam.pdf  (%d boxes)" % len(C))


def eg1(outdir):
    E0 = read_boxes(os.path.join(HERE, "out", "eg1_eps0.1", "E0.txt"))
    C = read_boxes(os.path.join(HERE, "out", "eg1_eps0.1", "E1.txt"))
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(6.2, 2.7))
    draw_boxes(a0, E0)
    frame(a0, hull(E0))
    draw_boxes(a1, C)
    draw_hull(a1, hull(C))
    frame(a1, hull(C))
    fig.tight_layout(w_pad=2.0)
    fig.savefig(os.path.join(outdir, "eg1_cover.pdf"))
    print("wrote doc/fig/eg1_cover.pdf  (%d sub-boxes, %d boxes in C)" % (len(E0), len(C)))


def main():
    outdir = os.path.join(HERE, "fig")
    os.makedirs(outdir, exist_ok=True)
    eijgenraam(outdir)
    eg1(outdir)


if __name__ == "__main__":
    main()
