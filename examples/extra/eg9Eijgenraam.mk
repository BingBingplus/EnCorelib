# ══════════════════════════════════════════════════════════════════════════
# eg9Eijgenraam.mk – Eijgenraam's wrapping-effect example (Thesis, Ex. 6.4.2,
# pp. 127-128; also reproduced in Nedialkov's PhD thesis)
#
#   x' = 0
#   y' = x^2
#
# Initial set:  x in [-lambda, lambda],  y = 0.
# Exact solution set at time T:  the parabola  y = T*x^2,  x in [-lambda,lambda].
#
# Eijgenraam's point: any single convex enclosure (box, parallelepiped, ...)
# of this curved end set has excess area of quadratic order in lambda --
# it is *impossible* to tightly enclose the whole set with ONE convex region.
# Our method instead subdivides the initial box and covers the end set with
# many small (eps-small) boxes, so it does not hit this obstruction.
#
# Run (from version_2026/):
#   make run EG=examples/eg9Eijgenraam.mk           # EndCover mode (mode=0)
#   make run EG=examples/eg9Eijgenraam.mk output_mode=3
#
# CLI format (for reference):
#   ./EnCorelib output_mode mode method stepB tubedegree n vars... funs... eps order T debug lo1 hi1 ...
# ══════════════════════════════════════════════════════════════════════════

#── Problem ────────────────────────────────────────────────────────────────
n   = 2
var = x,y
ff  = 0,x^2

# Initial box: center ± width per coordinate.  x in [-lambda,lambda], y=0.
cen = 0.0,0.0
wid = 1.0,0.0

#── Enclosure parameters ───────────────────────────────────────────────────
eps   = 0.02       # max width of end-enclosure
order = 20         # Taylor order for stepB

#── Time ───────────────────────────────────────────────────────────────────
T = 1.0

#── Algorithm options ──────────────────────────────────────────────────────
# mode: 0 = EndCover (whole box), 1 = Boundary (edges/faces only)
mode       = 0

# method: 0 = tube (Bisect+EulerTube, the default), 1 = bisect, 2 = splitonly
method     = 2

# stepB:  0 = Lohner,  1 = +lognorm,  2 = Lohner-HO+lognorm,  3 = Lohner-HO
stepB      = 0

# tubedegree:  0 = default (order-1 = 19),  1 = Euler tube,  p>=2 = Taylor-p
tubedegree = 0

#── Output ─────────────────────────────────────────────────────────────────
# output_mode: 0=time 1=+hull 2=+counts 3=+E0/E1 files 4=+plot files
output_mode = 3

# debug: 0 = silent, 1 = print containment warnings
debug = 0
