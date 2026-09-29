# ==========================================================================
#  Eg1 -- Lotka-Volterra predator-prey                (paper, Table 1)
#
#      x' = a x (1 - y)
#      y' = -b y (1 - x)          (a,b) = (2,1)
#
#      B0 = Box(1,3)(0.1)         H in {4, 5.5}       Ref. [19], [3, p.13]
#
#  Table 2 rows:   make table2-eg1
# ==========================================================================
n   = 2
var = x,y
ff  = 2*x-2*x*y,-y+x*y

cen = 1.0,3.0
wid = 0.1,0.1

eps   = 1.0
order = 20
T     = 4.0

output_mode = 2
mode        = 0
method      = 0
stepB       = 0
tubedegree  = 0
debug       = 0
