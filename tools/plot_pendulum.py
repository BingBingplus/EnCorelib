#!/usr/bin/env python3
"""Reproduce the top row of Figure 1 (undamped pendulum).

    python tools/plot_pendulum.py repro/out/fig1/pendulum [--save fig1a.png]

Left panel   a non-validated RK45 integration, the same family of method as
             MATLAB's ode45 with default tolerances.  The orbit started at
             (0, 2.0) makes one loop and then veers onto a different orbit,
             which is qualitatively wrong -- the exact orbit through (0,2)
             is below the separatrix (energy 1 - cos x + y^2/2 = 2 < 4) and
             must stay bounded in x.
Right panel  the validated enclosures produced by ./Encoretraj (see
             repro/figure1.sh), which stay on the correct orbit.

The left panel needs scipy; if MATLAB is available run tools/pendulum_ode45.m
instead, which is what the paper used.
"""

import argparse
import glob
import math
import os
import re
import sys

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    sys.exit("this script needs matplotlib:  pip install matplotlib")

Y0S = [-2.5, 1.5, 2.0, 2.5]
TEND = 20.0


def read_traj(path):
    ts, xs, ys = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            v = [float(t) for t in line.split()]
            ts.append(v[0])
            xs.append(0.5 * (v[1] + v[2]))
            ys.append(0.5 * (v[3] + v[4]))
    return ts, xs, ys


def rk45(y0):
    try:
        from scipy.integrate import solve_ivp
    except ImportError:
        return None
    f = lambda t, u: [u[1], -math.sin(u[0])]
    sol = solve_ivp(f, (0.0, TEND), [0.0, y0], method="RK45", dense_output=True,
                    rtol=1e-3, atol=1e-6, max_step=0.5)
    return sol.y[0], sol.y[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", nargs="?", default="repro/out/fig1/pendulum")
    ap.add_argument("--save", default=None)
    args = ap.parse_args()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    colors = {-2.5: "tab:blue", 1.5: "tab:green", 2.0: "tab:red", 2.5: "tab:purple"}

    ok = False
    for y0 in Y0S:
        r = rk45(y0)
        if r is None:
            continue
        ok = True
        axes[0].plot(r[0], r[1], lw=1.0, color=colors[y0], label="y(0)=%g" % y0)
    axes[0].set_title("RK45 (same family as ode45), no error control on the orbit")
    if not ok:
        axes[0].text(0.5, 0.5, "install scipy, or run tools/pendulum_ode45.m",
                     ha="center", transform=axes[0].transAxes)

    found = False
    for path in sorted(glob.glob(os.path.join(args.dir, "traj_y*.txt"))):
        m = re.search(r"traj_y(-?[0-9.]+)\.txt", os.path.basename(path))
        if not m:
            continue
        y0 = float(m.group(1))
        _, xs, ys = read_traj(path)
        axes[1].plot(xs, ys, lw=1.0, color=colors.get(y0, "k"), label="y(0)=%g" % y0)
        found = True
    axes[1].set_title("validated enclosures from ./Encoretraj")
    if not found:
        axes[1].text(0.5, 0.5, "run repro/figure1.sh first",
                     ha="center", transform=axes[1].transAxes)

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle(r"undamped pendulum  $x'=y,\ y'=-\sin x$,  $x(0)=0$")
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=160)
        print("wrote", args.save)
    else:
        plt.show()


if __name__ == "__main__":
    main()
