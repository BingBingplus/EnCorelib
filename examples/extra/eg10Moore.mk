# ══════════════════════════════════════════════════════════════════════════
# eg10Moore.mk – Moore's rotation example (classical wrapping-effect test)
#
#   x' = y
#   y' = -x
#
# Initial box:  x in [0.95,1.05], y in [-0.05,0.05]  (a small square near (1,0))
#
# The flow is an exact rigid rotation:  (x(t),y(t)) = R(-t)*(x0,y0), a linear
# map with a skew-symmetric Jacobian [[0,1],[-1,0]] (purely imaginary
# eigenvalues, log-norm = 0). So the TRUE end set at any time T is just the
# initial square, rotated by angle -T about the origin -- same size, same
# shape, forever.  This is R. Moore's classical example showing that a naive
# interval method (re-enclosing the rotated box by an axis-aligned box every
# step) suffers *unbounded* wrapping-effect growth over many periods, even
# though the true reachable set never grows.  It is the standard benchmark
# for testing whether a validated integrator (QR-preconditioning / CR-Lohner
# doubleton / eps-cover with subdivision) avoids this blow-up.
#
# Run (from version_2026/):
#   make run EG=examples/eg10Moore.mk           # EndCover mode (mode=0)
#   make run EG=examples/eg10Moore.mk T=6.283185307179586   # 1 period
#
# CLI format (for reference):
#   ./EnCorelib output_mode mode method stepB tubedegree n vars... funs... eps order T debug lo1 hi1 ...
# ══════════════════════════════════════════════════════════════════════════

#── Problem ────────────────────────────────────────────────────────────────
n   = 2
var = x,y
ff  = y,-x

# Initial box: center ± width per coordinate.
cen = 1.0,0.0
wid = 0.05,0.05

#── Enclosure parameters ───────────────────────────────────────────────────
eps   = 0.02       # max width of end-enclosure
order = 20         # Taylor order for stepB

#── Time ───────────────────────────────────────────────────────────────────
T = 6.283185307179586   # one full period (2*pi)

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
output_mode = 2

# debug: 0 = silent, 1 = print containment warnings
debug = 0
