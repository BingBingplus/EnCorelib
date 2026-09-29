#!/usr/bin/env bash
# ===========================================================================
#  repro/figure1.sh  --  data for Figure 1
#
#  BOTTOM ROW -- Eijgenraam's example [9, Ex. 6.4.2, p.127]
#
#      x' = 0,  y' = x^2,   x(0) in [-1,1],  y(0) = 0
#
#  At T = 1 the exact end set is the parabola {(x, x^2) : |x| <= 1}.  Any
#  single convex enclosure of it must contain the midpoint (0,1) of the two
#  endpoints (-1,1) and (1,1), although (0,1) is not in the set -- panel (a).
#  Panels (b),(c),(d) are the eps-covers for eps = 0.3, 0.1, 0.02.
#
#  TOP ROW -- undamped pendulum  x' = y, y' = -sin x
#
#      Validated trajectories from (0, y0), y0 in {-2.5, 1.5, 2.0, 2.5}.
#      MATLAB's ode45 gets the y0 = 2.0 orbit qualitatively wrong; run
#      tools/pendulum_ode45.m (MATLAB) or tools/pendulum_compare.py (Python)
#      for the left-hand panel.
#
#  Output:  repro/out/fig1/eijgenraam/eps<value>/{E0,E1}.txt
#           repro/out/fig1/pendulum/traj_y<value>.txt
# ===========================================================================
set -u

cd "$(dirname "$0")/.." || exit 1
BIN=./EnCorelib
[ -x "$BIN" ] || BIN=./EnCorelib.exe
[ -x "$BIN" ] || BIN=./bin/EnCorelib.exe
if [ ! -x "$BIN" ]; then echo "build first (or use bin/):  make" >&2; exit 1; fi

# --- bottom row: Eijgenraam ------------------------------------------------
for e in 0.3 0.1 0.02; do
  dir=repro/out/fig1/eijgenraam/eps$e
  mkdir -p "$dir"
  echo "== Eijgenraam, eps=$e =="
  ( cd "$dir" && "../../../../../$BIN" 3 0 2 0 0 2 x y 0 "x^2" "$e" 20 1 0 -1 1 0 0 )
done

# --- top row: pendulum -----------------------------------------------------
TRAJ=./Encoretraj
[ -x "$TRAJ" ] || TRAJ=./Encoretraj.exe
if [ -x "$TRAJ" ]; then
  dir=repro/out/fig1/pendulum
  mkdir -p "$dir"
  for y0 in -2.5 1.5 2.0 2.5; do
    echo "== pendulum, y0=$y0 =="
    "$TRAJ" 2 x y y "-sin(x)" 20 20 0.02  0 0  "$y0" "$y0" > "$dir/traj_y$y0.txt"
  done
  echo "trajectories in $dir"
else
  echo
  echo "note: ./Encoretraj was not built, skipping the pendulum panel."
  echo "      build it with:  make Encoretraj"
fi

echo
echo "data written under repro/out/fig1/"
echo "plot with:  python tools/plot_eijgenraam.py repro/out/fig1/eijgenraam"
echo "            python tools/plot_pendulum.py   repro/out/fig1/pendulum"
