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
*user-specified input*, not an outcome, and the algorithm has no tuning
parameters such as a minimum step size or an internal failure tolerance.

The tables and figures of the paper are regenerated from this directory:

```bash
make            # build
make check      # self-test
make table2     # Table 2   (the Ours(ε) rows)
make figure1    # Figure 1  (pendulum + Eijgenraam's example)
make figure4    # Figure 4  (geometry of end-covers)
make verify     # sample-point check of the computed covers
make audit      # box-by-box check of the computed covers
```

The **user guide**, [doc/userguide.pdf](doc/userguide.pdf), walks through
installation, a first run, writing your own problem, following a long run,
choosing the options and reading the error messages, with the real output of
every command it shows.  [CHANGES.md](CHANGES.md) lists what changed since the
first public version.

---

## 1. Requirements

| | |
|---|---|
| C++17 compiler | g++ ≥ 9, clang++ ≥ 10 |
| GNU make, bash | for the build and the scripts of `repro/` |
| [CAPD](http://capd.ii.uj.edu.pl/) | interval arithmetic; see section 2.1 |
| Python 3 + matplotlib | *optional*, for the plots |
| MATLAB | *optional*, only for `tools/pendulum_ode45.m` |
| [SymEngine](https://github.com/symengine/symengine) | *optional*, for algebraic expansion of the right-hand side |

---

## 2. Installing and building

### 2.1 Build CAPD

**Linux / macOS.**  Clone it wherever you like; remember the path, you will
need it in step 2.2.

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

### 2.2 Tell this project where CAPD is

**Edit one file: `config.mk`.**  It does not exist yet — copy the template:

```bash
cd <this directory>          # the one containing Makefile and EnCorelib.cpp
cp config.mk.example config.mk
```

Open `config.mk` and change the single line

```make
CAPD_MASTER := /change/me/to/your/CAPD
```

to the **absolute path of the directory you cloned CAPD into** — the one that
contains `capdDynSys/`, `capdAlg/`, `capdAux/`, `capdExt/` and `build/`.  For
the clone above that is:

```make
CAPD_MASTER := /home/yourname/CAPD               # Linux
CAPD_MASTER := /Users/yourname/CAPD              # macOS
CAPD_MASTER := /cygdrive/c/Users/yourname/CAPD   # Windows + Cygwin
```

`config.mk` is git-ignored, so your local paths stay on your machine and you
never have to touch `Makefile` itself.

Two alternatives, if you prefer not to create a file:

```bash
make CAPD_MASTER=/absolute/path/to/CAPD       # once per build
make CAPD_CONFIG=/usr/local/bin/capd-config   # an installed CAPD
```

and if `capd-config` is already on your `PATH`, plain `make` finds it with no
configuration at all.

### 2.3 Build

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

This produces `./EnCorelib`.  Three companion programs are built on demand:

```bash
make Encoretraj   # validated trajectory printer  (Figure 1, top row)
make capdref      # reference integrator used by repro/competitors/
make Encoreaudit  # box-by-box check of a cover   (make audit)
make tools        # all three
```

A CMake build is also provided; it can fetch CAPD itself, in which case
nothing has to be configured:

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build . --parallel
ctest                                   # runs repro/check.sh
```

Both build systems pass `-frounding-math`, which is required.

### 2.4 Check the build

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
| `tubedegree` | `0` or `1`: the Euler tube of Section 5.2 (no other degree is accepted) |
| `n` | dimension |
| `var_i`, `f_i` | variable names and right-hand sides, as strings |
| `eps` | the error bound `ε` |
| `order` | Taylor order `k` used by StepB (the paper uses `20`) |
| `T` | the horizon `H` |
| `debug` | progress on standard error: `0` silent · `1` one line per sub-box of `B0` · `2` + Refine passes · `3` + every stage and the containment checks (section 3.6) |
| `lo_i hi_i` | the initial box `B0` |

Example — Lotka–Volterra, `B0 = Box(1,3)(0.1)`, `H = 4`, `ε = 1`:

```bash
./EnCorelib 2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 4 0  0.9 1.1  2.9 3.1
```

```
EnCorelib  mode=EndCover  method=tube  stepB=lohner  tubedegree=1  order=20  n=2
         f  = (2*x-2*x*y, -y+x*y)
         B0 = {[0.9, 1.1],[2.9, 3.1]}   H=4   eps=1
time(ms)=671
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

**`tubedegree` — the tube of Subsection 5.2.**
`EulerTube` follows the Euler polygon (2.5), and its step-size threshold is
`h_euler` of equation (2.6).  `0` and `1` both select it; other values are
rejected, because a degree-`p` threshold is valid only for a degree-`p`
Taylor polygon, which is not implemented.

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

| file | contents |
|---|---|
| `E0.txt` | the initial sub-boxes that were integrated, one per box of `C` |
| `E1.txt` | the cover `C` at time `H` |
| `convex_hull.txt` | the 2D convex hull of `C`, as a closed polygon |
| `E_0.txt`, `E_1.txt` | as above plus the images of every `E0` box at times 0.1 / 0.4 / 0.7 and, in 2D, the trajectories of the four corners of `B0` |
| `out.txt` | one tab-separated row per run, appended; a cumulative log |

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
| `examples/dampedPendulum.mk` | not in the paper: `x'=y, y'=-sin x - 0.2y`, used in the user guide |

### 3.5 Writing the right-hand side

An `f_i` may use `+ - * / ^`, parentheses and the functions `sin cos exp log
sqrt sqr atan asin acos`.  `2x` means `2*x` and `(x+1)(y-1)` means
`(x+1)*(y-1)`; digits inside a name (`x1`) belong to the name.  In a `.mk`
file commas separate the `f_i`, so an `f_i` cannot contain one.

A decimal such as `0.2` is read as the nearest double, so the cover is
rigorous for the system with that coefficient.  A quotient such as `1/5` or
`8/3` is evaluated in interval arithmetic and encloses the exact value, and
scientific notation is turned into such a quotient (`2e-1` becomes
`(2/10^1)`).

Interval evaluation depends on how an expression is written.  When SymEngine
is found at build time every `f_i` is first multiplied out into monomials;
otherwise it is evaluated as written.  The `f =` line of the output shows the
form used.  On Lorenz the difference is large: to `H = 5` with `ε = 1` and
`stepB=0`, written as `10*(y-x), x*(28-z)-y, x*y-8*z/3`, it needed 585
sub-boxes of `B0` when evaluated as written and 73 multiplied out.  The files in
`examples/` and the scripts in `repro/` therefore use the multiplied-out form,
which is also the form the paper's numbers come from, and give the same
results with and without SymEngine.

### 3.6 Following a long run

The `debug` argument switches on progress reports, written to standard error
so that the results on standard output stay clean.  `debug=1` prints one line
per sub-box of `B0`, under a header:

```
$ ./EnCorelib 2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 0.5 20 4 1  0.9 1.1  2.9 3.1
     time    run depth  wmax(B) wmax(E_m) queue    done  outcome
     0.4s      1     0  2.0e-01   4.1e-01     4    0.0%  split, E_0 -> 1.0e-01
     0.5s      2     1  1.0e-01   4.6e-01     3   25.0%  resolved
     0.6s      3     1  1.0e-01   4.8e-01     2   50.0%  resolved
     0.8s      4     1  1.0e-01   1.6e-01     5   50.0%  split, E_0 -> 5.0e-02
   ...
     1.6s     12     2  5.0e-02   1.8e-01     1   93.8%  resolved
     1.7s     13     2  5.0e-02   1.9e-01     0  100.0%  resolved
     1.7s  EndCover finished: 13 sub-boxes processed, 10 resolved, deepest level 2, #(C)=13
```

`depth` says how often `B0` was subdivided to get the sub-box, `queue` how
many sub-boxes are still waiting and `done` how much of `B0` (by volume) is
finished.  `debug=2` adds one line per scaffold and per pass of Refine, and
`debug=3` one line per stage (`t`, the stage length, the enclosure) together
with the containment self-checks.  `2> progress.log` keeps the reports in a
file you can follow with `tail -f`.

### 3.7 When a run fails

A sub-box whose scaffold cannot be carried to `H` is subdivided and its
children are tried; with `debug=1` it shows up as `failed:` with the time
reached.  The run itself fails only when a sub-box that is a single point
fails.  Then subdividing cannot help, and the message says so:

```
$ make run EG=examples/eg5Lorenz.mk T=100 stepB=3 debug=1
     0.4s      1     0  2.0e-03         -     8    0.0%  failed: Extend: StepB cannot advance at t=3.83076 (H=100), split into 8
     0.8s      2     1  1.0e-03         -    15    0.0%  failed: Extend: StepB cannot advance at t=4.60846 (H=100), split into 8
   ...
    57.1s     39    38  7.1e-15         -   270    0.0%  failed: Extend: StepB cannot advance at t=26.5727 (H=100), split into 4
    59.3s     40    39  7.1e-15         -   273    0.0%  failed: Extend: StepB cannot advance at t=26.5367 (H=100), split into 4
    61.5s     41    40  7.1e-15         -   272    0.0%  failed: Extend: StepB cannot advance at t=26.583 (H=100)
[EndCover failed] Extend: StepB cannot advance at t=26.583 (H=100)
The sub-box of B0 on which this happened is a single point (depth 40 of the subdivision),
so subdividing B0 further cannot help: the solution through that point
cannot be enclosed up to the horizon.
Try a smaller T.
```

Here each level of subdivision gained about 0.76 time units, and the
solution through a single point could not be enclosed past t ≈ 26.6.
Arguments that are out of range and right-hand sides that cannot be read
end with exit status 2; failed computations with exit status 3.

### 3.8 Which options to use

Measured on the six problems of the paper (user guide, section 6):

- `eps` and `T` decide the result and most of the cost.  From `ε = 0.5` to
  `ε = 0.1` the number of boxes grew 4 to 12 times on five of the problems.
  For a chaotic system every level of subdivision gains a fixed amount of
  time and costs a factor `2^n` (user guide, section 7).
- `debug=1` is worth switching on for anything that runs longer than a minute.
- `method` 0, 1 and 2 gave covers with the same number of boxes and the same
  `B1` everywhere; `method=2` was the fastest, with 20% to 83% of the time of
  `method=0` at `ε = 0.1`.
- `stepB=1` gave the results of `stepB=0`, and `stepB=2` those of `stepB=3`,
  each at a higher cost.  `stepB=3` was the fastest variant in 12 of 13 runs.
- `order=20` was never the slowest; 15 or 30 is sometimes faster.
- `mode=1` is for when only the bounding box `B1` is needed.
- `method=2` with `stepB=3`, measured one run at a time, was 3 to 25 times
  faster than the settings of the example files and gave the same or a
  narrower `B1`, except on FitzHugh-Nagumo (6% wider).
- The example files keep the settings of the paper's experiments
  (`method=0 stepB=0 order=20`), which are also the ones `make table2` uses.

The user guide (section 6) has the measurements behind this advice.

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
inside some box of `C`.  It tests the left inclusion only: a box of `C` that
lies far from the end set would pass.

```
ok    eg1  method=0  #(C)=21     all 16 sample images covered
ok    eg1  method=1  #(C)=21     all 16 sample images covered
ok    eg1  method=2  #(C)=21     all 16 sample images covered
...
```

`make audit` (`repro/audit_cover.sh`, `tools/audit_cover.cpp`) checks the
cover box by box: every box of `C` must meet the validated images of the
centre and the corners of its own sub-box of `B0`, which `E0.txt` and
`E1.txt` pair line by line.

```
$ make audit
eps=0.1 stepB=0 order=20
eij  H=1   method=0  ok    79 boxes, 395 images: 0 boxes miss an image of their sub-box
eij  H=1   method=1  ok    79 boxes, 395 images: 0 boxes miss an image of their sub-box
eij  H=1   method=2  ok    79 boxes, 395 images: 0 boxes miss an image of their sub-box
eg1  H=4   method=0  ok    85 boxes, 425 images: 0 boxes miss an image of their sub-box
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
| `progress.h` | the progress reports selected by `debug` | — |
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
(the CAPD Group).  We are grateful for it.


