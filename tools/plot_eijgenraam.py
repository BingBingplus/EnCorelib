#!/usr/bin/env python3
"""Reproduce the bottom row of Figure 1 (Eijgenraam's example).

    python tools/plot_eijgenraam.py repro/out/fig1/eijgenraam [--save fig1b.png]

The directory must contain sub-directories eps<value>/ each holding an
E1.txt, as produced by repro/figure1.sh.

Panel (a) shows why no single convex set can enclose the end set tightly:
the exact end set at T = 1 is the parabola  y = x^2,  |x| <= 1,  and every
convex enclosure of it contains the midpoint (0,1) of its two endpoints,
although (0,1) is not on the parabola.  The remaining panels show the
eps-covers computed by EnCorelib.

Requires matplotlib.
"""

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boxio import read_boxes  # noqa: E402

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Polygon
except ImportError:  # pragma: no cover
    sys.exit("this script needs matplotlib:  pip install matplotlib")

T = 1.0


def parabola(m=400):
    xs = [-1.0 + 2.0 * k / (m - 1) for k in range(m)]
    return xs, [T * x * x for x in xs]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", help="repro/out/fig1/eijgenraam")
    ap.add_argument("--save", default=None)
    args = ap.parse_args()

    runs = []
    for d in sorted(glob.glob(os.path.join(args.dir, "eps*"))):
        m = re.search(r"eps([0-9.]+)", os.path.basename(d))
        f = os.path.join(d, "E1.txt")
        if m and os.path.exists(f):
            runs.append((float(m.group(1)), read_boxes(f)))
    runs.sort(key=lambda r: -r[0])
    if not runs:
        sys.exit("no eps*/E1.txt found under %s -- run repro/figure1.sh" % args.dir)

    fig, axes = plt.subplots(1, 1 + len(runs), figsize=(4.0 * (1 + len(runs)), 4.0))
    if 1 + len(runs) == 1:
        axes = [axes]

    xs, ys = parabola()

    # (a) the convex-enclosure obstruction
    ax = axes[0]
    ax.add_patch(Polygon(list(zip(xs, ys)) + [(1.0, T), (-1.0, T)],
                         closed=True, facecolor="0.85", edgecolor="0.5",
                         label="any convex enclosure"))
    ax.plot(xs, ys, "k-", lw=2, label="exact end set $E$")
    ax.plot([0.0], [T], "rx", ms=10, mew=2, label="(0,1) forced in, not in $E$")
    ax.plot([-1.0, 1.0], [T, T], "r--", lw=1)
    ax.set_title("(a) no convex set is tight")
    ax.legend(fontsize=8, loc="lower center")

    # (b), (c), ... the eps-covers
    for ax, (eps, boxes) in zip(axes[1:], runs):
        for b in boxes:
            (x0, x1), (y0, y1) = b[0], b[1]
            ax.add_patch(Rectangle((x0, y0), max(x1 - x0, 1e-12), max(y1 - y0, 1e-12),
                                   facecolor="tab:blue", alpha=0.35,
                                   edgecolor="tab:blue", linewidth=0.5))
        ax.plot(xs, ys, "k-", lw=1.5)
        ax.set_title(r"$\varepsilon$ = %g,   #(C) = %d" % (eps, len(boxes)))

    for ax in axes:
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-0.15, 1.35)
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    fig.suptitle(r"Eijgenraam:  $x'=0,\ y'=x^2$,  $x(0)\in[-1,1]$, $y(0)=0$,  $T=1$")
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=160)
        print("wrote", args.save)
    else:
        plt.show()


if __name__ == "__main__":
    main()
