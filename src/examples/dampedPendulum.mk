# ==========================================================================
#  A problem of your own -- the damped pendulum       (not in the paper)
#
#      x' = y
#      y' = -sin(x) - 0.2 y
#
#      B0 = Box(1,0)(0.1)         H = 10
#
#  The user guide (doc/userguide.pdf) builds this file step by step.
#      make run EG=src/examples/dampedPendulum.mk
#      make run EG=src/examples/dampedPendulum.mk eps=0.05 debug=1
# ==========================================================================
n   = 2
var = x,y
ff  = y,-sin(x)-0.2*y

cen = 1.0,0.0
wid = 0.1,0.1

eps   = 0.1
order = 20
T     = 10

output_mode = 2
mode        = 0
method      = 0
stepB       = 0
tubedegree  = 0
debug       = 0
