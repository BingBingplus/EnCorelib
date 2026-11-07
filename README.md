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
make audit      # box-by-box check of the computed covers
```

The **user guide**, [doc/userguide.pdf](doc/userguide.pdf), walks through
installation, a first run, writing your own problem, following a long run,
choosing the options and reading the error messages, with the real output of
every command it shows.  [CHANGES.md](CHANGES.md) lists what changed since the
first public version.

**On Windows you do not have to build anything.**  Ready-to-run executables
are shipped in [`bin/`](bin/); see section 2.1.

---

## 1. Requirements

| | |
|---|---|
| C++17 compiler | g++ ≥ 9, clang++ ≥ 10 |
| GNU make, bash | for the build and the scripts in `src/repro/` and `src/doc/` |
| [CAPD](http://capd.ii.uj.edu.pl/) | interval arithmetic; needed to build from source, see section 2.2 |
| Python 3 + matplotlib | *optional*, for the plots |
| MATLAB | *optional*, only for `src/tools/pendulum_ode45.m` |
| [SymEngine](https://github.com/symengine/symengine) | *optional*, for algebraic expansion of the right-hand side |

---

## 2. Installing and building

### 2.1 Ready-to-run Windows binaries

`bin/` holds prebuilt 64-bit Windows executables, so the program can be run
without installing a compiler or CAPD:

| file | what it is |
|---|---|
| `bin/EndCorelib.exe` | the solver — everything described in section 3 below |
| `bin/EndCoretraj.exe` | prints validated trajectories (Figure 1, top row, and `make verify`) |
| `bin/EndCoreaudit.exe` | the box-by-box check of a cover (`make audit`) |
| `bin/capdref.exe` | the CAPD reference integrator of `repro/competitors/` |
| `bin/*.dll` | the runtime libraries the executables need |

Keep the `.dll` files next to the `.exe` files and run them straight from the
repository root, from Command Prompt, PowerShell or Git Bash (in Git Bash
write `bin/EndCorelib.exe`):
Input:
```bash
bin\EndCorelib.exe 2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 4 0  0.9 1.1  2.9 3.1
```
Output:
```
EndCorelib  mode=EndCover  method=tube  stepB=lohner  tubedegree=1  order=20  n=2
         f  = (2*x-2*x*y, -y+x*y)
         B0 = {[0.9, 1.1],[2.9, 3.1]}   H=4   eps=1
time(ms)=338
Hull(T=4.000000)={[1.013875, 1.945561],[0.131115, 0.263816]}
B1=(1.479718, 0.197466) +- (0.465843, 0.066350)
wmax(B1)=0.931686
#(C)=2  #E0=2  max_{B in C} wmax(B)=0.798371
```

The following shell scripts are used to build the rows of Table 2 called
`Ours(ε)`.  Run from Git Bash, they use `bin/EndCorelib.exe` when no program
has been built:
```bash
bash src/repro/check.sh          # the self-test
bash src/repro/table2.sh         # Table 2
```

To build from source instead continue with 2.2.

### 2.2 Build CAPD

**Linux / macOS.**  Clone CAPD to the folder of your choice; remember the
path, you will need it in step 2.3.

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
cd <this directory>          # the one containing Makefile and src/
cp config.mk.example config.mk
```

Open `config.mk` and change the single line

```make
CAPD_MASTER := /change/me/to/your/CAPD
```

to the **absolute path of the directory that you cloned CAPD into** — the one
that contains `capdDynSys/`, `capdAlg/`, `capdAux/`, `capdExt/` and `build/`.
For the clone above, it is:

```make
CAPD_MASTER := /home/yourname/CAPD               # Linux
CAPD_MASTER := /Users/yourname/CAPD              # macOS
CAPD_MASTER := /cygdrive/c/Users/yourname/CAPD   # Windows + Cygwin
```

`config.mk` is git-ignored, so your local paths stay on your machine and you
never have to touch `Makefile` itself.

Alternatives, if you prefer not to create the file or already have the CAPD
library installed:

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
Built EndCorelib
```

If instead you see `CAPD not found`, `CAPD_MASTER` is pointing at the wrong
directory — check that `$CAPD_MASTER/build/libcapd.a` exists.

This produces `./EndCorelib`.  Three companion programs are built on demand:

```bash
make EndCoretraj   # validated trajectory printer  (Figure 1, top row)
make capdref       # reference integrator used by repro/competitors/
make EndCoreaudit  # box-by-box check of a cover   (make audit)
make tools         # all three
```

A CMake build is also provided.  It uses an installed CAPD, or downloads
CAPD and builds it from source first, so nothing has to be configured:

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build . --parallel
ctest                                   # runs src/repro/check.sh
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
./EndCorelib  output_mode mode method stepB tubedegree n \
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
| `debug` | progress on standard error: `0` silent · `1` one line per scaffold taken from the queue · `2` + Refine passes · `3` + every stage and the containment checks (section 3.6) |
| `lo_i hi_i` | the initial box `B0` |

Example — Lotka–Volterra, `B0 = Box(1,3)(0.1)`, `H = 4`, `ε = 1`:

```bash
./EndCorelib 2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 4 0  0.9 1.1  2.9 3.1
```

```
EndCorelib  mode=EndCover  method=tube  stepB=lohner  tubedegree=1  order=20  n=2
         f  = (2*x-2*x*y, -y+x*y)
         B0 = {[0.9, 1.1],[2.9, 3.1]}   H=4   eps=1
time(ms)=249
Hull(T=4.000000)={[1.013875, 1.945561],[0.131115, 0.263816]}
B1=(1.479718, 0.197466) +- (0.465843, 0.066350)
wmax(B1)=0.931686
#(C)=2  #E0=2  max_{B in C} wmax(B)=0.798371
```

`B1` is the `B1` column of Table 2 and `#(C)` its `#(C)` column.

### 3.2 The options

**`mode` — what gets covered.**
`mode=0` is EndCover proper: `Split` halves the initial box until every piece
has an `ε`-small end-enclosure (Subsection 3.3 of the paper).
`mode=1` runs EndCover on the boundary of `B0` only — the four edges in 2D,
the `2n` faces in `nD`.

**`method` — which refinement operators Refine may use.**

Section 5 gives Refine three operators: `Bisect`, `EulerTube` and `Split`.
`method` says which of them are switched on.

| value | name | operators |
|---|---|---|
| `0` | `tube` | `Bisect` + `EulerTube` + `Split`, driven by `StageFreeze()` (5.4), `NeedSplit()` (5.5) and `CanTerminate()` (5.6). This is `Refine` as described in Section 5, and the default. |
| `1` | `bisect` | the same loop with `EulerTube` switched off; stages are refined by `Bisect` alone. |
| `2` | `splitonly` | `Split` is the only operator; the accuracy comes from splitting `E_0`. |

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

| output file | contents |
|---|---|
| `out.txt` | a cumulative log, one tab-separated row per run, appended |
| `E0.txt` | the initial sub-boxes that were integrated, one per box of `C` |
| `E1.txt` | the cover `C` at time `H` |
| `convex_hull.txt` | the 2D convex hull of `C`, as a closed polygon |
| `E_0.txt`, `E_1.txt` | as above plus the images of every `E0` box at times 0.1 / 0.4 / 0.7 and, in 2D, the trajectories of the four corners of `B0` |

### 3.4 Example files

`src/examples/*.mk` hold the problems of Table 1; `src/examples/_template.txt` lists
every field.  Run one with:

```bash
make run-eg1Volterra                       # Eg1 with the defaults in the file
make run EG=src/examples/eg5Lorenz.mk          # the same, spelled out
make run EG=src/examples/eg5Lorenz.mk eps=0.1 T=4 method=0 output_mode=4
make run-boundary-eg1Volterra              # force mode=1
make show-args EG=src/examples/eg5Lorenz.mk    # print the command without running it
```

Any field can be overridden on the command line.

| file | Table 1 entry |
|---|---|
| `src/examples/eg1Volterra.mk` | Eg1 Lotka–Volterra, `(a,b)=(2,1)`, `Box(1,3)(0.1)`, `H ∈ {4, 5.5}` |
| `src/examples/eg2VanDerPol.mk` | Eg2 Van der Pol, `c=1`, `Box(-3,3)(0.1)`, `H ∈ {1, 2}` |
| `src/examples/eg3Quadratic.mk` | Eg3 `(y, x²)`, `Box(1,-1)(0.05)`, `H ∈ {1, 4}` |
| `src/examples/eg4FitzHughNagumo.mk` | Eg4 FitzHugh–Nagumo, `(a,b,ε,I)=(0.7,0.8,0.08,0.5)`, `Box(1,0)(0.1)`, `H ∈ {1, 4}` |
| `src/examples/eg5Lorenz.mk` | Eg5 Lorenz, `(σ,ρ,β)=(10,28,8/3)`, `Box(15,15,36)(0.001)`, `H ∈ {1, 4}` |
| `src/examples/eg6Rossler.mk` | Eg6 Rössler, `(a,b,c)=(0.2,0.2,5.7)`, `Box(1,2,3)(0.1)`, `H ∈ {1, 4}` |
| `src/examples/fig1Eijgenraam.mk` | Figure 1 (bottom), `x'=0, y'=x²` |
| `src/examples/fig1Pendulum.mk` | Figure 1 (top), `x'=y, y'=-sin x` |
| `src/examples/dampedPendulum.mk` | not in the paper: `x'=y, y'=-sin x - 0.2y`, used in the user guide |

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
form used.  On Lorenz the difference shows in the work: to `H = 5` with
`ε = 1` and `stepB=0`, written as `10*(y-x), x*(28-z)-y, x*y-8*z/3`, it took
116 scaffolds and 78 boxes when evaluated as written, and 95 scaffolds and
64 boxes multiplied out.  The files in
`src/examples/` and the scripts in `src/repro/` therefore use the multiplied-out form,
which is also the form the paper's numbers come from, and give the same
results with and without SymEngine.

### 3.6 Following a long run

The `debug` argument switches on progress reports, written to standard error
so that the results on standard output stay clean.  `debug=1` prints one line
per scaffold taken from the queue, under a header:

```
$ ./EndCorelib 2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 0.5 20 4 1  0.9 1.1  2.9 3.1
     time    run depth  wmax(B) wmax(E_m) queue    done  outcome
     0.4s      1     0  2.0e-01   4.1e-01     2   25.0%  resolved, E_0 -> 1.0e-01, 2 split off
     0.6s      2     2  1.0e-01   4.0e-01     1   50.0%  resolved
     0.9s      3     1  2.0e-01   3.7e-01     1   75.0%  resolved, E_0 -> 1.0e-01, 1 split off
     1.0s      4     2  1.0e-01   3.6e-01     0  100.0%  resolved
     1.0s  EndCover finished: 4 scaffolds processed, deepest level 2, #(C)=4
```

`depth` says how often `B0` was split to get the initial box of the
scaffold, `queue` how many scaffolds are still waiting and `done` how much of
`B0` (by volume) is covered.  `split off` counts the halves that `Split`
gave to new scaffolds.  `debug=2` adds one line per scaffold and per pass of Refine, and
`debug=3` one line per stage (`t`, the stage length, the enclosure) together
with the containment self-checks.  `2> progress.log` keeps the reports in a
file you can follow with `tail -f`.

### 3.7 When a run stops

A scaffold that cannot be carried to `H` has its initial box split in two,
and the halves are tried; with `debug=1` it shows up as `split:` with the
time reached.  The run stops only when double precision is exhausted, that
is, when such a box can no longer be halved, and the message says so:

```
$ make run EG=src/examples/eg5Lorenz.mk T=100 stepB=3 debug=1
     0.4s      1     0  2.0e-03         -     2    0.0%  split: Extend: StepB cannot advance at t=3.83076 (H=100)
     0.8s      2     1  2.0e-03         -     3    0.0%  split: Extend: StepB cannot advance at t=4.01253 (H=100)
   ...
   375.4s    121   120  7.1e-15         -   122    0.0%  split: Extend: StepB cannot advance at t=26.4072 (H=100)
   380.7s    122   121  7.1e-15         -   121    0.0%  stopped: Extend: StepB cannot advance at t=26.4064 (H=100)
[EndCover stopped] double precision exhausted
Extend: StepB cannot advance at t=26.4064 (H=100)
The sub-box of B0 on which this happened (depth 121 of the subdivision) cannot be halved in double precision.
Try a smaller T.
```

Here every three levels, which halve each side of the box once, gained
about 0.72 time units, and from depth 96 on the scaffolds stopped between
t = 26.40 and 26.78.
Arguments that are out of range and right-hand sides that cannot be read
end with exit status 2; computations that stop end with exit status 3.

### 3.8 Which options to use

Measured on the six problems of the paper (user guide, section 6):

- `eps` and `T` decide the result and most of the cost.  From `ε = 0.5` to
  `ε = 0.1` the number of boxes grew 2 to 12 times where the runs finished;
  with the settings of the example files, Lotka-Volterra, Lorenz and Rössler
  did not finish within ten minutes at `ε = 0.1`.  For a chaotic system every
  level of splitting gains a fixed amount of time and doubles the work (user
  guide, section 7).
- `debug=1` is worth switching on for anything that runs longer than a minute.
- `method=2` was the fastest strategy, up to 120 times faster than
  `method=0` at `ε = 0.1`, and finished every run; its `B1` was up to 21%
  wider at `ε = 0.5` and at most 0.5% wider at `ε = 0.1`.
- `stepB=1` gave the results of `stepB=0`, and `stepB=2` those of `stepB=3`,
  each at a higher cost.  `stepB=3` was the fastest variant wherever both it
  and `stepB=0` finished.
- `order=20` was the fastest order; 10, 15 or 30 sometimes gave a hull a few
  per cent narrower.
- `mode=1` is for when only the bounding box `B1` is needed.
- `method=2` with `stepB=3` finished every problem within a minute, 2.4 to
  86 times faster than the settings of the example files where those
  finished within ten minutes; `B1` was within 0.6% at `ε = 0.1` and 43% wider
  on Lorenz to `H = 5` at `ε = 1`.
- The example files keep the settings of the paper's experiments
  (`method=0 stepB=0 order=20`), which are also the ones `make table2` uses.

The user guide (section 6) has the measurements behind this advice.

---

## 4. The tables and figures of the paper

### Table 2

```bash
make table2                                  # the Ours(ε) rows
make compare                                 # diff that run against the paper
bash src/repro/table2.sh eg2 eg4 eg6             # a subset of the examples
EPS="1.0 0.5 0.1 0.02" bash src/repro/table2.sh eg5
METHOD=2 bash src/repro/table2.sh                # a different refinement strategy
```

`make table2` writes `repro/out/table2.txt` and `repro/out/table2.tex`.
`make compare` (`src/repro/compare.py`) prints that run cell by cell next to
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
python src/tools/plot_eijgenraam.py repro/out/fig1/eijgenraam --save fig1_bottom.png
python src/tools/plot_pendulum.py   repro/out/fig1/pendulum   --save fig1_top.png
```

*Bottom row* — Eijgenraam's example `x'=0, y'=x²`, `x(0) ∈ [-1,1]`, `y(0)=0`.
At `T=1` the exact end set is the parabola `y = x²`.  Any convex enclosure of
it must contain the midpoint `(0,1)` of its two endpoints although `(0,1)` is
not in the set; an ε-cover does not.  `ε = 0.3, 0.1, 0.02` give covers of 12,
40 and 174 boxes.

*Top row* — the undamped pendulum, `x' = y, y' = -sin x`, started at
`(0, 2.0)`.  The MATLAB panel of the paper is produced by
`src/tools/pendulum_ode45.m`; `src/tools/plot_pendulum.py` draws the same panel with
SciPy.

### Figure 4

```bash
make figure4
python src/tools/plot_boxes.py repro/out/fig4/eg1 --save fig4_eg1.png
python src/tools/plot_boxes.py repro/out/fig4/eg5 --save fig4_eg5.png
```

For `n = 2` this draws the subdivision of `B0` and the corresponding
end-cover, one panel each; for `n = 3` the xy-, xz- and yz-projections of the
cover.  Boxes coming from the same initial sub-box share a colour.

### Checking the covers

```bash
make verify                                       # all examples, all methods
GRID=7 EPS=0.05 bash src/repro/verify_cover.sh eg2    # a finer test
```

`src/repro/verify_cover.sh` tests `End(B0,H) ⊆ ⋃ C` from the outside: for every
example it takes a grid of sample points of `B0`, integrates each of them
separately with `./EndCoretraj`, and checks that each resulting enclosure lies
inside some box of `C`.  It tests the left inclusion only: a box of `C` that
lies far from the end set would pass.

```
ok    eg1  method=0  #(C)=10     all 16 sample images covered
ok    eg1  method=1  #(C)=10     all 16 sample images covered
ok    eg1  method=2  #(C)=12     all 16 sample images covered
...
```

`make audit` (`src/repro/audit_cover.sh`, `src/tools/audit_cover.cpp`) checks the
cover box by box: every box of `C` must meet the validated images of the
centre and the corners of its own sub-box of `B0`, which `E0.txt` and
`E1.txt` pair line by line.

```
$ make audit
eps=0.1 stepB=0 order=20
eij  H=1   method=0  ok    40 boxes, 200 images: 0 boxes miss an image of their sub-box
eij  H=1   method=1  ok    40 boxes, 200 images: 0 boxes miss an image of their sub-box
eij  H=1   method=2  ok    40 boxes, 200 images: 0 boxes miss an image of their sub-box
eg1  H=4   method=0  ok    48 boxes, 240 images: 0 boxes miss an image of their sub-box
...
```

---

## 5. Code map

A file-by-file walkthrough is in **[CODE_SUMMARY.md](CODE_SUMMARY.md)**: it
covers every struct and function and the paper equation each one implements.

| file | contents | paper |
|---|---|---|
| `src/types.h` | boxes, `wmax`, midpoints, logNorm, `Scaffold`, `MiniScaffold`, `Init` | Sec. 2.1, 2.4, 3.2, 5.1 |
| `src/stepAB.h` | `StepA` (both entry points), the four `StepB` variants, the tube thresholds, `SolutionChain` | Sec. 4.1, 4.2, eq. (2.2), (2.6) |
| `src/extend.h` | `Extend` — append one stage to the scaffold; `ExtendToHorizon`, `RebuildScaffold` | Sec. 4, 3.3 |
| `src/refine.h` | `Evaluate`, `Bisect`, `EulerTube`, `Split`, `StageFreeze`, `NeedSplit`, `CanTerminate`, `Refine` | Sec. 5 |
| `src/endcover.h` | `EndCover` — the queue of scaffolds | Sec. 3.3 |
| `src/boundary.h` | `runScaffold` (the Extend/Refine driver) and the Boundary mode | Sec. 3.3 |
| `src/progress.h` | the progress reports selected by `debug` | — |
| `src/symparse.h` | right-hand-side string pre-processing | — |
| `src/EndCorelib.cpp` | command line, reporting, output files | Sec. 7 |

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
| `Init(B)` | `Init(B, G0)` |
| `S.Bisect(i)`, `S.EulerTube(i)`, `S.Split()` | `Bisect`, `EulerTube`, `Split(S, Qout)` |
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
