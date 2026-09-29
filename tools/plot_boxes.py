#!/usr/bin/env python3
"""Draw an end-cover the way Figure 4 of the paper does.

    python tools/plot_boxes.py <dir> [--title "Eg1, H=4"] [--save fig.png]

<dir> must contain E0.txt (the initial sub-boxes) and E1.txt (the cover C
at time H), i.e. the output of `EnCorelib` at output level 3 or above --
see repro/figure4.sh.

n = 2   two panels: the subdivision of B0 on top, the end-cover below.
n = 3   three panels: the xy-, xz- and yz-projections of the cover, with
        the corresponding projections of B0 inset.

Boxes belonging to the same initial sub-box get the same colour, so the
two panels can be read together.

Requires matplotlib.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boxio import read_boxes, bounds, pad  # noqa: E402

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
except ImportError:  # pragma: no cover
    sys.exit("this script needs matplotlib:  pip install matplotlib")


def draw(ax, boxes, i, j, colors, lw=0.6, alpha=0.55):
    for k, b in enumerate(boxes):
        (x0, x1), (y0, y1) = b[i], b[j]
        w, h = x1 - x0, y1 - y0
        ax.add_patch(Rectangle((x0, y0), w if w > 0 else 1e-12,
                               h if h > 0 else 1e-12,
                               facecolor=colors[k % len(colors)], alpha=alpha,
                               edgecolor="0.25", linewidth=lw))
    lo, hi = pad(*bounds(boxes, i))
    ax.set_xlim(lo, hi)
    lo, hi = pad(*bounds(boxes, j))
    ax.set_ylim(lo, hi)
    ax.set_aspect("auto")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", help="directory holding E0.txt and E1.txt")
    ap.add_argument("--title", default=None)
    ap.add_argument("--save", default=None, help="write a PNG instead of showing")
    args = ap.parse_args()

    e0p = os.path.join(args.dir, "E0.txt")
    e1p = os.path.join(args.dir, "E1.txt")
    for p in (e0p, e1p):
        if not os.path.exists(p):
            sys.exit("missing %s -- run repro/figure4.sh first" % p)

    E0 = read_boxes(e0p)
    E1 = read_boxes(e1p)
    if not E1:
        sys.exit("E1.txt is empty")
    n = len(E1[0])

    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % 20) for i in range(max(len(E0), len(E1)))]

    title = args.title or os.path.basename(os.path.normpath(args.dir))
    info = os.path.join(args.dir, "run.info")
    if os.path.exists(info):
        with open(info) as fh:
            title += "   (" + fh.readline().strip() + ")"

    if n == 2:
        fig, axes = plt.subplots(2, 1, figsize=(5.2, 8.0))
        draw(axes[0], E0, 0, 1, colors)
        axes[0].set_title("initial sub-boxes of $B_0$   (#%d)" % len(E0))
        axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
        draw(axes[1], E1, 0, 1, colors)
        axes[1].set_title("end-cover $C$ at time $H$   (#%d)" % len(E1))
        axes[1].set_xlabel("x"); axes[1].set_ylabel("y")
    else:
        pairs = [(0, 1, "xy"), (0, 2, "xz"), (1, 2, "yz")]
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
        for ax, (i, j, lab) in zip(axes, pairs):
            draw(ax, E1, i, j, colors)
            ax.set_title("%s-projection of $C$   (#%d)" % (lab, len(E1)))
            ax.set_xlabel(lab[0]); ax.set_ylabel(lab[1])

    fig.suptitle(title)
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=160)
        print("wrote", args.save)
    else:
        plt.show()


if __name__ == "__main__":
    main()
