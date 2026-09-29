# ==========================================================================
#  Figure 1 (bottom) -- Eijgenraam's example    [9, Example 6.4.2, p.127]
#                                               (also [22, p.77-78])
#      x' = 0
#      y' = x^2        x(0) in [-1,1],  y(0) = 0
#
#  At time T the exact end set is the parabola  E = {(x, T x^2) : |x| <= 1},
#  so NO single convex set encloses E tightly: any convex enclosure must
#  contain the midpoint (0,1) of (-1,1) and (1,1), which is not in E.
#  An eps-cover has no such obstruction -- panels (b),(d),(c) of Figure 1
#  are this example with eps = 0.3, 0.1, 0.02.
#
#  Reproduce with:   make figure1
# ==========================================================================
n   = 2
var = x,y
ff  = 0,x^2

# x in [-1,1], y = 0
cen = 0.0,0.0
wid = 1.0,0.0

eps   = 0.1
order = 20
T     = 1.0

output_mode = 3
mode        = 0
method      = 0
stepB       = 0
tubedegree  = 0
debug       = 0
