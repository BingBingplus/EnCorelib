# The three competing solvers of Table 2

Table 2 of the paper compares `Ours(eps)` with three external validated
solvers.  Only `Ours(eps)` is produced by this repository
(`bash repro/table2.sh`); the other three need software that we cannot
redistribute.  This directory records exactly how they were driven.

| Row in Table 2 | Software | How it was run |
|---|---|---|
| `CAPD` | [CAPD](http://capd.ii.uj.edu.pl/) `CnRect2Set`, Cr-Lohner with r = 3 | `capd_compare.cpp` in this directory — build with `make capdref` |
| `VNODE-LP` | [VNODE-LP](https://www.cas.mcmaster.ca/~nedialk/vnodelp/) | `vnodelp_compare.cc` template below |
| `Taylor-model` | [INTLAB](https://www.tuhh.de/ti3/rump/intlab/) V14, `verifyode` with preconditioned Taylor models | `intlab_compare.m` template below |

All three were used **with their default configuration**, as stated in
Section 7.1 of the paper.  In particular the default minimum step sizes are
2^-20 (CAPD), 2^-52 (VNODE-LP) and 1e-4 (INTLAB); an entry `Invalid` in
Table 2 means the solver hit that threshold, and `No Output` means it
stopped without returning an enclosure.

## CAPD

```bash
make capdref
# ./capdref n  var_1..var_n  f_1..f_n  order H  lo_1 hi_1 .. lo_n hi_n
./capdref 2 x y 2*x-2*x*y -y+x*y 20 4    0.9 1.1  2.9 3.1      # Eg1, H=4
./capdref 2 x y 2*x-2*x*y -y+x*y 20 5.5  0.9 1.1  2.9 3.1      # Eg1, H=5.5 -> No Output
./capdref 3 x y z 10*\(y-x\) x*\(28-z\)-y x*y-8*z/3 20 4 \
          14.999 15.001 14.999 15.001 35.999 36.001            # Eg5, H=4
```

The `rho` column is then `wmax(B1)` of this run divided by `wmax(B1)` of
`Ours(1.0)` for the same example and horizon.

## VNODE-LP

`vnodelp_compare.cc` is a template; fill in the right-hand side and the
initial box and link against your VNODE-LP installation.

```cpp
#include "vnode.h"
using namespace vnodelp;

template<typename var_type>
void f(int n, var_type *yp, const var_type *y, var_type t, void *param) {
    yp[0] = 2*y[0] - 2*y[0]*y[1];
    yp[1] = -y[1] + y[0]*y[1];
}

int main() {
    const int n = 2;
    AD *ad = new FADBAD_AD(n, f, f);
    VNODE *solver = new VNODE(ad);

    interval t = 0.0, tend = 4.0;
    iVector y(n);
    y[0] = interval(0.9, 1.1);
    y[1] = interval(2.9, 3.1);

    solver->integrate(t, y, tend);
    if (!solver->successful()) { printf("Invalid\n"); return 1; }
    printVector(y);
    return 0;
}
```

## INTLAB Taylor-model

```matlab
startintlab
f = @(t,y) [2*y(1) - 2*y(1)*y(2); -y(2) + y(1)*y(2)];
y0 = [infsup(0.9,1.1); infsup(2.9,3.1)];
opt = verifyodeset('order', 20);          % defaults otherwise
tic
[T, Y] = verifyode(f, [0 4], y0, opt);
toc
B1 = Y(:, end);                            % end enclosure
disp([mid(B1) rad(B1)])
```
