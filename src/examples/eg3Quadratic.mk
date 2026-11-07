# ==========================================================================
#  Eg3 -- Quadratic system                            (paper, Table 1)
#
#      x' = y
#      y' = x^2
#
#      B0 = Box(1,-1)(0.05)       H in {1, 4}         Ref. [3, p.11]
# ==========================================================================
n   = 2
var = x,y
ff  = y,x^2

cen = 1.0,-1.0
wid = 0.05,0.05

eps   = 1.0
order = 20
T     = 1.0

output_mode = 2
mode        = 0
method      = 0
stepB       = 0
tubedegree  = 0
debug       = 0
