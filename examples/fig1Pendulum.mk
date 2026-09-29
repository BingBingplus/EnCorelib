# ==========================================================================
#  Figure 1 (top) -- undamped pendulum
#
#      x' = y
#      y' = -sin(x)
#
#  MATLAB's ode45 run from (x,y) = (0, 2.0) makes one loop and then veers
#  onto a different orbit, which is qualitatively wrong; the validated
#  computation converges to the saddle point.  See repro/figure1.sh and
#  tools/pendulum_ode45.m.
#
#  The .mk below covers the initial *point* (0,2) -- a degenerate box.
# ==========================================================================
n   = 2
var = x,y
ff  = y,-sin(x)

cen = 0.0,2.0
wid = 0.0,0.0

eps   = 1.0
order = 20
T     = 6.0

output_mode = 3
mode        = 0
method      = 0
stepB       = 0
tubedegree  = 0
debug       = 0
