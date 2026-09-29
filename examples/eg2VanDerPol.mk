# ==========================================================================
#  Eg2 -- Van der Pol oscillator                      (paper, Table 1)
#
#      x' = y
#      y' = c (1 - x^2) y - x     c = 1
#
#      B0 = Box(-3,3)(0.1)        H in {1, 2}         Ref. [3, p.2]
# ==========================================================================
n   = 2
var = x,y
ff  = y,(1-x^2)*y-x

cen = -3.0,3.0
wid = 0.1,0.1

eps   = 1.0
order = 20
T     = 1.0

output_mode = 2
mode        = 0
method      = 0
stepB       = 0
tubedegree  = 0
debug       = 0
