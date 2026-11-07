# EndCover for the Initial Value Problem

Reference implementation of the algorithm in

> B. Zhang and C. Yap.
> **End Cover for Initial Value Problem: Complete Validated Algorithm with
> Complexity Analysis.**
> Condensed from arXiv [2502.00503](https://arxiv.org/abs/2502.00503) and
> [2602.00162](http://arxiv.org/abs/2602.00162).

Given an autonomous ODE `x' = f(x)`, a box `B0 ⊆ R^n`, a horizon `H > 0` and
an error bound `ε > 0`, the program computes a finite set `C` of boxes with

```
End(B0,H)  ⊆  ⋃ C  ⊆  End(B0,H) + [-ε,ε]^n
```

where `End(B0,H) = { x(H) : x(0) ∈ B0 }`.  `C` is called an **ε-cover** of the
reachable set.  Unlike the usual validated IVP solvers, the accuracy is a
*user-specified input*, not an outcome, and the algorithm does not need any tuning
parameters such as a minimum step size or an internal failure tolerance.

The tables and figures of the paper are regenerated from this directory:

```bash
make            # build
make check      # self-test
make table2     # Table 2   (the Ours(ε) rows)
make figure1    # Figure 1  (pendulum + Eijgenraam's example)
make figure4    # Figure 4  (geometry of end-covers)
make verify     # sample-point check of the computed covers
```

**On Windows you do not have to build anything.**  Ready-to-run executables
are shipped in [`bin/`](bin/); see section 2.1.

---

## 1. Requirements

| | |
|---|---|
| C++17 compiler | g++ ≥ 9, clang++ ≥ 10 |
| GNU make, bash | for the build and the scripts of `repro/` |
| [CAPD](http://capd.ii.uj.edu.pl/) | interval arithmetic; needed to build from source, see section 2.2 |
| Python 3 + matplotlib | *optional*, for the plots |
| MATLAB | *optional*, only for `tools/pendulum_ode45.m` |
| [SymEngine](https://github.com/symengine/symengine) | *optional*, for algebraic expansion of the right-hand side |

---

## 2. Installing and building

### 2.1 Ready-to-run Windows binaries

`bin/` holds prebuilt 64-bit Windows executables, so the program can be run
without installing a compiler or CAPD:

| file | what it is |
|---|---|
| `bin/EnCorelib.exe` | the solver — everything described in section 3 below |
| `bin/Encoretraj.exe` | this print the validated trajectory in Figure 1 (and `make verify`) |
| `bin/*.dll` | the runtime libraries the two executables need |

Keep the `.dll` files next to the `.exe` files and run it straight from the
repository root, from Command Prompt, PowerShell or Git Bash:
Input:
```bash
bin\EnCorelib.exe 2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 4 0  0.9 1.1  2.9 3.1
```
Output:
```
EnCorelib  mode=EndCover  method=tube  stepB=lohner  tubedegree=19  order=20  n=2
         f  = (2*x-2*x*y, -y+x*y)
         B0 = {[0.9, 1.1],[2.9, 3.1]}   H=4   eps=1
Hull(T=4.000000)={[1.093880, 1.853071],[0.146098, 0.243519]}
B1=(1.473475, 0.194809) +- (0.379596, 0.048710)
#(C)=5  #E0=5  max_{B in C} wmax(B)=0.551660
```

The following shell scripts are used to build the rows of Table 2 called `Ours(ε)`. 
```bash
bash repro/check.sh          # the self-test
bash repro/table2.sh         # Table 2
```

To build from source instead continue with 2.2.

### 2.2 Build CAPD

**Linux / macOS.**  Clone this to the folder of your choice; remember the path, you will
need it in step 2.3.

```bash
cd ~                                     # or any directory you prefer
git clone https://github.com/CAPDGroup/CAPD.git
cd CAPD && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCAPD_BUILD_EXAMPLES=OFF -DCAPD_BUILD_TESTS=OFF
make -j4
```

When this finishes you should have `~/CAPD/build/libcapd.a`.  `sudo make
install` is optional — the build tree alone is enough.

**Windows.**  Do the same inside [Cygwin](https://www.cygwin.com/) (install
the packages `gcc-g++`, `make`, `cmake`, `git`) and run everything below from
a Cygwin login shell.  Build and run from the same shell that built CAPD, so
that the include paths and the interval library match.  Inside Cygwin a
Windows path `C:\Users\you\CAPD` is written `/cygdrive/c/Users/you/CAPD`.

### 2.3 Tell this project where CAPD is

**Edit one file: `config.mk`.**  It does not exist yet — copy the template:

```bash
cd <this directory>          # the one containing Makefile and EnCorelib.cpp
cp config.mk.example config.mk
```

Open `config.mk` and change the single line

```make
CAPD_MASTER := /change/me/to/your/CAPD
```

to the **absolute path of the directory you that cloned CAPD into** — the one that
contains `capdDynSys/`, `capdAlg/`, `capdAux/`, `capdExt/` and `build/`.  For
the clone above, it is:

```make
CAPD_MASTER := /home/yourname/CAPD               # Linux
CAPD_MASTER := /Users/yourname/CAPD              # macOS
CAPD_MASTER := /cygdrive/c/Users/yourname/CAPD   # Windows + Cygwin
```

`config.mk` is git-ignored, so your local paths stay on your machine and you
never have to touch `Makefile` itself.

Alternatives, if you already have the CAPD library then link to it:

```bash
make CAPD_MASTER=/absolute/path/to/CAPD       # once per build
make CAPD_CONFIG=/usr/local/bin/capd-config   # an installed CAPD
```

and if `capd-config` is already on your `PATH`, plain `make` finds it with no
configuration at all.

### 2.4 Build

```bash
make
```

The first line of output tells you which CAPD it picked up:

```
CAPD source tree: /home/yourname/CAPD
Built EnCorelib
```

If instead you see `CAPD not found`, `CAPD_MASTER` is pointing at the wrong
directory — check that `$CAPD_MASTER/build/libcapd.a` exists.

This produces `./EnCorelib`.  Two companion programs are built on demand:

```bash
make Encoretraj   # validated trajectory printer  (Figure 1, top row)
make capdref      # reference integrator used by repro/competitors/
make tools        # both
```

If you modify the CAPD source, then you can use the following to recompile instead of above.

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build . --parallel
ctest                                   # runs repro/check.sh
```

Both build systems pass `-frounding-math`, which is required.

### 2.5 Check the build

```bash
make check
```

18 assertions; all of them should say `ok`.

---

## 3. Running

### 3.1 Command line

```
./EnCorelib  output_mode mode method stepB tubedegree n \
           var_1 ... var_n  f_1 ... f_n \
           eps order T debug \
           lo_1 hi_1 ... lo_n hi_n
```

| argument | values |
|---|---|
| `output_mode` | `0` time only · `1` +enclosure · `2` +counts · `3` +`E0/E1` files · `4` +plot files |
| `mode` | `0` EndCover (cover all of `B0`) · `1` Boundary (cover `∂B0` only) |
| `method` | which refinement operators Refine may use: `0` tube (default) · `1` bisect · `2` splitonly |
| `stepB` | `0` Lohner (default) · `1` Lohner+logNorm · `2` Lohner-HO+logNorm · `3` Lohner-HO |
| `tubedegree` | `0` = `order-1` · `1` Euler tube · `p ≥ 2` Taylor tube of degree `p` |
| `n` | dimension |
| `var_i`, `f_i` | variable names and right-hand sides, as strings |
| `eps` | the error bound `ε` |
| `order` | Taylor order `k` used by StepB (the paper uses `20`) |
| `T` | the horizon `H` |
| `debug` | `0` silent · `1` report containment checks |
| `lo_i hi_i` | the initial box `B0` |

Example — Lotka–Volterra, `B0 = Box(1,3)(0.1)`, `H = 4`, `ε = 1`:

```bash
./EnCorelib 2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 4 0  0.9 1.1  2.9 3.1
```

```
EnCorelib  mode=EndCover  method=tube  stepB=lohner  tubedegree=19  order=20  n=2
         f  = (2*x-2*x*y, -y+x*y)
         B0 = {[0.9, 1.1],[2.9, 3.1]}   H=4   eps=1
time(ms)=419
Hull(T=4.000000)={[1.093880, 1.853071],[0.146098, 0.243519]}
B1=(1.473475, 0.194809) +- (0.379596, 0.048710)
wmax(B1)=0.759191
#(C)=5  #E0=5  max_{B in C} wmax(B)=0.551660
```

`B1` is the `B1` column of Table 2 and `#(C)` its `#(C)` column.

### 3.2 The options

**`mode` — what gets covered.**
`mode=0` is EndCover proper: the initial box is subdivided until every piece
has an `ε`-small end-enclosure (Subsection 3.3 of the paper).
`mode=1` covers only the boundary of `B0` — the four edges in 2D, the `2n`
faces in `nD`.

**`method` — which refinement operators Refine may use.**

Section 5 gives Refine three operators: `Bisect`, `EulerTube` and `Split`.
`method` says which of them are switched on.

| value | name | operators |
|---|---|---|
| `0` | `tube` | `Bisect` + `EulerTube` + `Split`, driven by `StageFreeze()` (5.4), `NeedSplit()` (5.5) and `CanTerminate()` (5.6). This is `Refine` as described in Section 5, and the default. |
| `1` | `bisect` | the same loop with `EulerTube` switched off; stages are refined by `Bisect` alone. |
| `2` | `splitonly` | `Split` is the only operator; the accuracy comes from shrinking `E_0`. |

`method` is independent of `stepB`: `stepB` chooses how one step is enclosed,
`method` chooses how the scaffold is refined once those steps exist.

**`stepB` — which end-enclosure step.**

| value | name | representation |
|---|---|---|
| `0` | `lohner` | Lohner's method: the centred term is kept as a doubleton with a QR-preconditioned frame |
| `1` | `lohner+lognorm` | the same, intersected with the ball of Lemma 2.1 |
| `2` | `lohnerHO+lognorm` | Lohner's method with a higher-order a-priori enclosure, intersected with the ball of Lemma 2.1 |
| `3` | `lohnerHO` | as `2` without the ball |

**`tubedegree` — the tube threshold of Subsection 2.6.**
`1` is the Euler threshold `h_euler` of equation (2.6); `p ≥ 2` is the
degree-`p` Taylor threshold; `0` means `order-1`.

**`order`** is the Taylor order `k` of StepB, and **`eps`** is the error bound
`ε` of the cover.  There is nothing else to tune.

### 3.3 Output levels

Output is progressive — level `k` prints everything below `k` as well.

| level | adds |
|---|---|
| `0` | `time(ms)=` |
| `1` | `Hull(T)=`, `B1=(mid) +- (half-width)`, `wmax(B1)=`, and a header echoing the problem |
| `2` | `#(C)=`, `#E0=`, and the largest `wmax` over the individual boxes of `C` |
| `3` | writes `E0.txt`, `E1.txt` and (in 2D) `convex_hull.txt` |
| `4` | writes `E_0.txt`, `E_1.txt` for plotting |

| output file | contents |
|---|---|
| `out.txt` | a cumulative log, one tab-separated row per run, appended; |
| `E0.txt` | the initial sub-boxes that were integrated, one per box of `C` |
| `E1.txt` | the cover `C` at time `H` |
| `convex_hull.txt` | the 2D convex hull of `C`, as a closed polygon |
| `E_0.txt`, `E_1.txt` | as above plus the images of every `E0` box at times 0.1 / 0.4 / 0.7 and, in 2D, the trajectories of the four corners of `B0` |


### 3.4 Example files

`examples/*.mk` hold the problems of Table 1; `examples/_template.txt` lists
every field.  Run one with:

```bash
make run-eg1Volterra                       # Eg1 with the defaults in the file
make run EG=examples/eg5Lorenz.mk          # the same, spelled out
make run EG=examples/eg5Lorenz.mk eps=0.1 T=4 method=0 output_mode=4
make run-boundary-eg1Volterra              # force mode=1
make show-args EG=examples/eg5Lorenz.mk    # print the command without running it
```

Any field can be overridden on the command line.

| file | Table 1 entry |
|---|---|
| `examples/eg1Volterra.mk` | Eg1 Lotka–Volterra, `(a,b)=(2,1)`, `Box(1,3)(0.1)`, `H ∈ {4, 5.5}` |
| `examples/eg2VanDerPol.mk` | Eg2 Van der Pol, `c=1`, `Box(-3,3)(0.1)`, `H ∈ {1, 2}` |
| `examples/eg3Quadratic.mk` | Eg3 `(y, x²)`, `Box(1,-1)(0.05)`, `H ∈ {1, 4}` |
| `examples/eg4FitzHughNagumo.mk` | Eg4 FitzHugh–Nagumo, `(a,b,ε,I)=(0.7,0.8,0.08,0.5)`, `Box(1,0)(0.1)`, `H ∈ {1, 4}` |
| `examples/eg5Lorenz.mk` | Eg5 Lorenz, `(σ,ρ,β)=(10,28,8/3)`, `Box(15,15,36)(0.001)`, `H ∈ {1, 4}` |
| `examples/eg6Rossler.mk` | Eg6 Rössler, `(a,b,c)=(0.2,0.2,5.7)`, `Box(1,2,3)(0.1)`, `H ∈ {1, 4}` |
| `examples/fig1Eijgenraam.mk` | Figure 1 (bottom), `x'=0, y'=x²` |
| `examples/fig1Pendulum.mk` | Figure 1 (top), `x'=y, y'=-sin x` |

---

## 4. The tables and figures of the paper

### Table 2

```bash
make table2                                  # the Ours(ε) rows
make compare                                 # diff that run against the paper
bash repro/table2.sh eg2 eg4 eg6             # a subset of the examples
EPS="1.0 0.5 0.1 0.02" bash repro/table2.sh eg5
METHOD=2 bash repro/table2.sh                # a different refinement strategy
```

`make table2` writes `repro/out/table2.txt` and `repro/out/table2.tex`.
`make compare` (`repro/compare.py`) prints that run cell by cell next to
`repro/table2_paper.txt`, which holds Table 2 as printed.

Timings depend on the machine; the paper used a 13th Gen Intel Core
i7-13700HX (2.10 GHz) with 16 GB RAM.

The three competing solvers of Table 2 (INTLAB Taylor-model, VNODE-LP, CAPD)
are not part of this repository.  `repro/competitors/` records how each was
driven and contains a ready-to-build reference integrator:

```bash
make capdref
./capdref 2 x y 2*x-2*x*y -y+x*y 20 4  0.9 1.1  2.9 3.1
```

### Figure 1

```bash
make figure1
python tools/plot_eijgenraam.py repro/out/fig1/eijgenraam --save fig1_bottom.png
python tools/plot_pendulum.py   repro/out/fig1/pendulum   --save fig1_top.png
```

*Bottom row* — Eijgenraam's example `x'=0, y'=x²`, `x(0) ∈ [-1,1]`, `y(0)=0`.
At `T=1` the exact end set is the parabola `y = x²`.  Any convex enclosure of
it must contain the midpoint `(0,1)` of its two endpoints although `(0,1)` is
not in the set; an ε-cover does not.  `ε = 0.3, 0.1, 0.02` give covers of 23,
79 and 351 boxes.

*Top row* — the undamped pendulum, `x' = y, y' = -sin x`, started at
`(0, 2.0)`.  The MATLAB panel of the paper is produced by
`tools/pendulum_ode45.m`; `tools/plot_pendulum.py` draws the same panel with
SciPy.

### Figure 4

```bash
make figure4
python tools/plot_boxes.py repro/out/fig4/eg1 --save fig4_eg1.png
python tools/plot_boxes.py repro/out/fig4/eg5 --save fig4_eg5.png
```

For `n = 2` this draws the subdivision of `B0` and the corresponding
end-cover, one panel each; for `n = 3` the xy-, xz- and yz-projections of the
cover.  Boxes coming from the same initial sub-box share a colour.

### Checking the covers

```bash
make verify                                       # all examples, all methods
GRID=7 EPS=0.05 bash repro/verify_cover.sh eg2    # a finer test
```

`repro/verify_cover.sh` tests `End(B0,H) ⊆ ⋃ C` from the outside: for every
example it takes a grid of sample points of `B0`, integrates each of them
separately with `./Encoretraj`, and checks that each resulting enclosure lies
inside some box of `C`.

```
ok    eg1  method=0  #(C)=21     all 16 sample images covered
ok    eg1  method=1  #(C)=21     all 16 sample images covered
ok    eg1  method=2  #(C)=21     all 16 sample images covered
...
```

---

## 5. Code map

A file-by-file walkthrough is in **[CODE_SUMMARY.md](CODE_SUMMARY.md)**: it
covers every struct and function and the paper equation each one implements.

| file | contents | paper |
|---|---|---|
| `types.h` | boxes, `wmax`, midpoints, logNorm, `Scaffold`, `MiniScaffold` | Sec. 2.1, 2.4, 3.2, 5.1 |
| `stepAB.h` | `StepA` (both entry points), the four `StepB` variants, the tube thresholds, `SolutionChain` | Sec. 4.1, 4.2, eq. (2.2), (2.6) |
| `extend.h` | `Extend` — append one stage to the scaffold; `ExtendToHorizon`, `RebuildScaffold` | Sec. 4, 3.3 |
| `refine.h` | `Evaluate`, `Bisect`, `EulerTube`, `Split`, `StageFreeze`, `NeedSplit`, `CanTerminate`, `Refine` | Sec. 5 |
| `endcover.h` | `EndCover` — the queue of scaffolds over sub-boxes of `B0` | Sec. 3.3 |
| `boundary.h` | `runScaffold` (the Extend/Refine driver) and the Boundary mode | Sec. 3.3 |
| `symparse.h` | right-hand-side string pre-processing | — |
| `EnCorelib.cpp` | command line, reporting, output files | Sec. 7 |

### Notation

| paper | code |
|---|---|
| `B`, `B0` | `IVector` |
| `w(B)`, `wmax(B)`, `m(B)` | `width`, `wmax(B)`, `midpoint(B)` |
| `Box(S)` | `hullBox(A,B)` |
| `μ(A) = μ₂(A)` | `logNorm(J, n)` |
| `S = (t, E, F, G)` | `struct Scaffold { t, E, F, G }`, `S.m()`, `S.T()` |
| `G[i] = ((ℓᵢ, Eⁱ, Fⁱ), μⁱ, (δᵢ, h^i_euler))` | `struct MiniScaffold { ell, E, F, mu, delta, hTube }` |
| `StepA(E0,H,δ)`, `StepB(E0,h,F1)` | `StepA(...)`, `StepB(...)` |
| `S.Extend(H,ε₀)`, `S.Refine(ε₀)` | `Extend(...)`, `Refine(...)` |
| `S.Bisect(i)`, `S.EulerTube(i)`, `S.Split()` | `Bisect`, `EulerTube`, `Split` |
| `StageFreeze(i)` — eq. (5.4) | `StageFreeze(S, i, eps0)` |
| `NeedSplit()` — eq. (5.5) | `NeedSplit(S, eps0, n)` |
| `CanTerminate()` — eq. (5.6) | `CanTerminate(freezeFlag, needSplit)` |
| `N^tube(S)`, `mu_max` | `nTubeOf(S)`, `muMaxOf(S)` |
| the chain `E_{i-1} -> Quad[i,1] -> ... -> E_i` | one `SolutionChain`, walked by `Extend` and `Evaluate` |
| admissibility test — eq. (2.2) | `StepAatStep(...)` |
| `h_euler(H,M,μ,δ)` — eq. (2.6) | `hEuler(H,M,mu,delta)` |
| `EndCover_f(B0,H,ε₀) → C` | `endcover::EndCover(...)` |
| `ε₀` | `eps0` (CLI `eps`) |
| `C`, `B1 = Box(C)` | `CoverResult::C`, `CoverResult::hull` |

Everything runs in IEEE double precision.

---

## 6. Citing

```bibtex
@misc{zhang-yap-endcover,
  author = {Bingwei Zhang and Chee Yap},
  title  = {End Cover for Initial Value Problem: Complete Validated
            Algorithm with Complexity Analysis},
  note   = {arXiv:2502.00503, arXiv:2602.00162}
}
```

---

## 7. Acknowledgements

Interval arithmetic with directed rounding and the normalised Taylor
coefficients `f^[i]` are provided by **[CAPD](http://capd.ii.uj.edu.pl/)**
(the CAPD Group, Jagiellonian University).  We are grateful for it.

This work was funded by NSF Grant #CCF-2212462.
