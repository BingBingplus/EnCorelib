# ==========================================================================
#  Eg4 -- FitzHugh-Nagumo neuron model                (paper, Table 1)
#
#      x' = x - x^3/3 - y + I
#      y' = e (x + a - b y)       (a,b,e,I) = (0.7, 0.8, 0.08, 0.5)
#
#      B0 = Box(1,0)(0.1)         H in {1, 4}         Ref. [29]
# ==========================================================================
n   = 2
var = x,y
ff  = x-x^3/3-y+0.5,0.08*(x+0.7-0.8*y)

cen = 1.0,0.0
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
