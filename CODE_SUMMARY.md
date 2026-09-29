# EnCorelib — structure of the code

> Paper: B. Zhang and C. Yap, *End Cover for Initial Value Problem: Complete
> Validated Algorithm with Complexity Analysis* (arXiv:2502.00503 /
> arXiv:2602.00162).  This document walks through the headers in dependency
> order and gives, for every struct and function, what it does and which part
> of the paper it implements.

**The problem.**  Given an autonomous ODE `x' = f(x)`, an initial box `B0`, a
horizon `H` and an error bound `ε > 0`, compute a finite set of boxes `C` with

```
End(B0,H)  ⊆  ⋃ C  ⊆  End(B0,H) + [-ε,ε]^n
```

where `End(B0,H) = { x(H) : x(0) ∈ B0 }`.  `C` is an **ε-cover** of the
reachable set.

---

## 📁 Layout

```
version_2026/
├── types.h                # box utilities + the two core structs (Scaffold / MiniScaffold)
├── symparse.h             # right-hand-side string pre-processing (optional SymEngine)
│
├── stepAB.h               # ★ StepA / the four StepB variants / tube thresholds / SolutionChain
│
├── extend.h               # S.Extend(H, ε₀) — append one stage to the scaffold
├── refine.h               # S.Refine(ε₀) — Evaluate / Bisect / EulerTube /
│                          #   Split / StageFreeze / NeedSplit / CanTerminate
│
├── endcover.h             # the EndCover loop (the subdivision queue over B0)
├── boundary.h             # runScaffold (the Extend/Refine driver) + Boundary mode
│
├── EnCorelib.cpp          # command line, output levels, output files
│
├── examples/              # the six problems of Table 1 + the two cases of Figure 1 (.mk)
├── repro/                 # scripts: table2 / figure1 / figure4 / check / verify
│   └── competitors/       # reference integrator + VNODE-LP / INTLAB drivers
├── tools/                 # plot scripts (Python/MATLAB), Encoretraj, run.sh
└── extras/                # additional diagnostic programs
```

**Dependency order** (each layer uses only the ones above it):

```
types.h  →  stepAB.h  →  extend.h  →  refine.h
                                        ↓
                         endcover.h ←  boundary.h  →  EnCorelib.cpp
```

---

## 📚 File by file

### ⭐ Layer 1: infrastructure

---

## 📄 **types.h** — box utilities and the core structs


### Box utilities

| function | paper | description | in | out |
|------|---------|------|------|------|
| `hullBox(B1, B2)` | `Box(B1 ∪ B2)` | smallest axis-aligned box containing both | two IVector | IVector |
| `IntersectB(B1, B2)` | `B1 ∩ B2` | coordinatewise intersection; keeps `B1[i]` when a coordinate comes out empty | two IVector | IVector |
| `wmax(B)` | `wmax(B)` | largest edge width `max_i w(B)_i` | IVector | double |
| `wmin(B)` | `wmin(B)` | smallest edge width | IVector | double |
| `midpoint(B)` | `m(B)` | midpoint, as a plain vector of doubles | IVector | vector\<double\> |
| `euclideanRadius(B)` | `½·wmax(Ball(B))` | radius of the smallest **Euclidean ball** containing B, `‖w(B)‖₂ / 2` | IVector | double |

**Empty intersections in `IntersectB`.**  Both arguments enclose the same
solution set, so an empty coordinate can only come from outward rounding;
keeping `B1[i]`, itself a valid enclosure, is sound.

**`euclideanRadius`.**  Lemma 2.1 is stated in the Euclidean norm, so every
ball radius is the radius `|w(B)|₂/2` of the circumscribed ball; all ball
radii in the code go through this function.

### Derivatives and the logarithmic norm

| function | paper | description | in | out |
|------|---------|------|------|------|
| `computeJacobian(f, B)` | `J_f(B)` | Jacobian over a box, from a degree-2 jet | IMap, IVector | IMatrix |
| `logNorm(J, n)` | `μ(A) = μ₂(A)` | upper bound on the logarithmic norm (Subsection 2.4) | IMatrix, dimension | double |

**Two paths in `logNorm`:**

- `n == 2`: closed form.  The largest eigenvalue of the symmetric part
  `(J+Jᵀ)/2` is `(a+d)/2 + sqrt(((a-d)/2)² + b²)` with
  `a=J₀₀, d=J₁₁, b=(J₀₁+J₁₀)/2`.
- `n > 2`: the Euclidean logarithmic norm.

**Why μ may be negative.**  Unlike a p-norm, `μ(A)` can be negative, which is
why the paper uses it in place of a Lipschitz constant: a contracting flow
gets a **contracting** bound instead of an exponentially growing one.

### String utilities

| function | description | example |
|------|------|------|
| `Convert_to_IMap(SVar, SFun)` | assembles the map string the parser expects | `{x,y}, {y, -x}` → `"var:x,y;fun:y,-x;"` |

### Struct 1: `MiniScaffold` — `G[i]` (Subsection 5.1, [G1]–[G5])

```cpp
struct MiniScaffold {
    std::vector<double>  mu;      // μ^i : logNorm bound of each mini-step     [G4]
    double               delta;   // δ_i : current tube tolerance              [G5]
    double               hTube;   // h^i_euler : tube step-size threshold      [G5]
    int                  ell;     // ℓ_i : this stage holds 2^ℓ_i mini-steps   [G2]
    std::vector<IVector> E;       // **E**^i : 2^ℓ_i + 1 mini end-enclosures   [G3]
    std::vector<IVector> F;       // **F**^i : 2^ℓ_i mini full enclosures      [G3]
    int                  nTube;   // tube segments made by the last EulerTube(i)
};
```



### Struct 2: `Scaffold` — `S = (t, E, F, G)` (Subsection 3.2)

```cpp
struct Scaffold {
    std::vector<double>       t;  // S.t = (t₀ < t₁ < ... < t_m)
    std::vector<IVector>      E;  // S.E = (E₀, ..., E_m)
    std::vector<IVector>      F;  // S.F = (F₁, ..., F_m)
    std::vector<MiniScaffold> G;  // S.G = (G₁, ..., G_m)

    int    m() const;             // S.m — the number of stages
    double T() const;             // S.T = t.back() — the current final time
};
```

**Stage i is the admissible quad** `(E_{i-1}, Δt_i, F_i, E_i)` with
`Δt_i = t_i - t_{i-1}`, and every `E_i` contains `End(E₀, t_i)`.

### Struct 3: `CoverResult`

```cpp
struct CoverResult {
    std::vector<IVector> E0;    // the initial sub-boxes that were integrated
    std::vector<IVector> C;     // the ε-cover (the paper's C)
    IVector              hull;  // Box(C), the B1 column of Table 2
};
```

---

## 📄 **symparse.h** — right-hand-side string pre-processing

**What it does.**  Puts a user-written expression into the form the
right-hand-side parser accepts.

| function | description | in | out |
|------|------|------|------|
| `basic_sanitize(s)` | strips whitespace, inserts implicit products | `"2x + 3(y)"` → `"2*x+3*(y)"` | string |
| `expand_expression(s)` | entry point: optional SymEngine expansion, otherwise `basic_sanitize` | string | string |

**Compile switch.**  With `-DHAVE_SYMENGINE` the SymEngine path is used
(`^` ↔ `**` conversion plus `expand`); otherwise only the light clean-up runs.
Both paths give the same result on the built-in examples; SymEngine is needed
only where an expression such as `(x+y)^20` has to be expanded for real.

**Note.**  SymEngine rewrites expressions, e.g.
`x-x^3/3-y+0.5` → `0.5+x-y+(-1/3)*x^3`.  The rewritten form is what the header
line prints.

---

### ⭐ Layer 2: the one-step algorithms

---

## 📄 **stepAB.h** — StepA / StepB / tube thresholds / SolutionChain

**What it does.**  The whole one-step machinery of the paper: Subsection 4.1
(StepA), equation (2.2) (admissibility), Subsection 4.2 (StepB), Subsection
2.4 (the ball of Lemma 2.1), equation (2.6) (`h_euler`), and **how a chain of
mini-quads is evaluated** (`SolutionChain`).

### Taylor polynomial of a step (Subsection 2.2)

| function | description | in | out |
|------|------|------|------|
| `taylorPolyOverStep(f, E, h, k)` | `P(E,h) = Σ_{i<k} [0,h]^i f^[i](E)`, evaluated by Horner as Subsection 4.1 suggests | IMap, start box, step, order | IVector |

### StepA — two entry points, one algorithm (Subsection 4.1, Lemma 4.1)

`StepA(E₀, H, δ) → (h, F₁)`, where `(E₀, h, F₁)` is a δ-admissible triple with
`h ≤ H`.  Admissibility is the inclusion of equation (2.2):

```
Σ_{i<k} [0,h]^i f^[i](E₀) + [0,h]^k f^[k](F₁)  ⊆  F₁
```

| function | description | in | out |
|------|------|------|------|
| `StepA(f, E0, H, delta, k=STEPA_ORDER)` | **chooses** h as large as Lemma 4.1 allows, then certifies F₁ | IMap, E₀, upper bound H, δ, order | `pair<double h, IVector F1>` |
| `StepAatStep(f, E, h, delta, k=STEPA_ORDER)` | **h is already fixed** (by a scaffold stage, or by StepA itself); only the full enclosure is wanted | as above with h fixed | IVector F₁ |

**Pseudo-code** (line by line as in the paper):

```
h ← 0
While (H > h)
    F₁ ← Σ_{i<k} [0,H]^i f^[i](E₀) + Box(δ)          ← Horner
    M  ← ( ‖f^[k](F₁)‖₁, ..., ‖f^[k](F₁)‖ₙ )
    h  ← min{ H, minᵢ (δᵢ / Mᵢ)^{1/k} }               ← Lemma 4.1
    H  ← H / 2
F₁ ← Σ_{i<k} [0,h]^i f^[i](E₀) + [0,h]^k f^[k](F₁)   ← last line: the contraction
Return (h, F₁)
```

**★ The last line is the point.**  It replaces the inflated candidate `F₁` by
the (narrower) set the inclusion actually certifies.  Both entry points end
with it, so wherever the scaffold needs an admissible full enclosure it goes
through the same StepA; there is no second mechanism.

### The four StepB variants (Subsection 4.2)

StepB upgrades an admissible triple `(E₀, h, F₁)` to a quad
`(E₀, h, F₁, E₁)`.  **StepB is Lohner's method**: the centred term

```
( Σ_{i<k} h^i J_{f^[i]}(E₀) ) · (E₀ − m(E₀))
```

is a **matrix acting on E₀**, and it is kept in that form — a doubleton

```
E₁ = q₁ + C₁·r₀ + B₁·r₁
```

with a QR-preconditioned frame — instead of being re-enclosed in an
axis-aligned box.  Re-enclosing once per step is exactly the wrapping effect.

| value | function | description |
|----|------|------|
| 0 | `stepBlohner(f, degree, h, B0)` | Lohner's method |
| 1 | `stepBlohnerWithMu(f, degree, h, B0, mu)` | Lohner's method ∩ the ball of Lemma 2.1 |
| 2 | `stepBlohnerHOWithMu(f, degree, h, B0, mu)` | Lohner's method (higher-order a-priori enclosure) ∩ the ball |
| 3 | `stepBlohnerHO(f, degree, h, B0)` | as 2, without the ball |

| helper | description |
|---------|------|
| `StepB(f, degree, h, B0, mu, stepBtype)` | dispatcher over the four |

The four values of `stepBtype` are listed in the options summary at the end;
the default is 0.

### The ball of Lemma 2.1 (Subsection 2.4)

| function | description | in | out |
|------|------|------|------|
| `logNormBallTighten(f, order, h, E0, mu, dim, E)` | tightens E with `Ball_{x₁(h)}(r₀e^{μh})` | `r₀ = euclideanRadius(E₀)`, `x₁` the trajectory through `m(E₀)` | IVector |

When the ball cannot be computed, E is returned unchanged, i.e. not tightened.

### Tube step-size thresholds (equation (2.6), Lemma 2.2)

| function | description | in | out |
|------|------|------|------|
| `hEuler(H, M, mu, delta)` | the Euler tube threshold | stage length, `‖f^[2]‖`, μ, δ | double |
| `taylorCoeffNorm(f, B, p, degreeAvail)` | `‖f^[p+1]‖` | — | double |
| `hTaylorStep(H, Cbar, mu, delta, p)` | the degree-p threshold (`p=1` is `hEuler`) | — | double |
| `computeHTube(f, F, H, delta, mu, tubedegree, degreeAvail)` | dispatcher: `tubedegree==1` → Euler, `≥2` → Taylor-p | — | double |

```
                  ⎧ min{H, 2μδ / (M(e^{μH}−1) − μ²δ)}   μ > 0
h_euler(H,M,μ,δ) = ⎨ min{H, δ / (MH)}                    μ = 0
                  ⎩ min{H, 2μδ / (M(e^{μH}−1))}         μ < 0
```

While `0 < h ≤ h_euler` the Euler polygon stays inside the δ-tube of the exact
trajectory.  A return value of `0` means "no usable tube at this stage", which
is what an overflowing `e^{μH}` amounts to.

### ★ `SolutionChain` — how a chain of mini-quads is evaluated

A stage of the scaffold is a chain of admissible mini-quads (equation (5.1)):

```
E_{i−1} → Quad[i,1] → ... → Quad[i,2^ℓ_i] → E_i
```

The paper prescribes what each mini-quad **contains**, not how the chain is
**evaluated**.  Evaluating it one box at a time — take StepB, store the hull,
start the next mini-step from that hull — re-encloses the parallelepiped of
the centred term at every mini-step, and that is the wrapping effect.

`SolutionChain` instead carries a single doubleton

```
E_j = q_j + C_j·r₀ + B_j·r_j
```

anchored at `E₀` throughout, so `E_j` is the **accumulated affine map `A_j`
applied to the original `E₀`**.  A hull is taken only when some `E_j` is
**reported**, and a reported hull is never fed back into the propagation.

| member | description | in | out |
|---------|------|------|------|
| `SolutionChain(f, order, E0, stepBtype)` | construct, anchored at E₀; `stepBtype` picks the representation (as in StepB) | — | — |
| `restart(E0)` | start again at a new (smaller) E₀ from `t = 0`; used after `Split` | IVector | — |
| `proposedStep()` | the step StepB **would like** to take next | — | double |
| `step(hmax)` | advance by at most hmax: the step it wants when that fits, exactly hmax when it does not | double | bool |
| `marchTo(t)` | reach t exactly, halving any step that fails | double | bool |
| `hull()` / `time()` | the current hull / the current time | — | — |

**"Take one step" is the only primitive.**  There is no time map on the chain,
no "integrate straight to some instant".  The scaffold is laid out one step at
a time — each step is a stage, or a mini-step of one — so `step()` is the only
action needed.

**A step that fails is rolled back.**  `step()` does not throw; it puts the
doubleton back exactly as it was and returns `false`, and the caller
(`Extend` / `marchTo`) halves the step and tries again.  The doubleton itself
is never re-enclosed.

**Landing exactly.**  The reported box must be the solution at the target
instant, so the final piece shorter than one step is covered by
`step(remaining)` **exactly**.

**Not copyable.**  The chain owns `m_f` (declared before the propagator) and
`= delete`s the copy operations.

---

### ⭐ Layer 3: the two scaffold methods

---

## 📄 **extend.h** — `S.Extend(H, ε₀)` (Section 4)

**What it does.**  Appends one stage to the scaffold.

```cpp
inline void Extend(IMap f, Scaffold& S, SolutionChain& chain,
                   double eps0, double delta, int degree, double H,
                   int stepBtype, int tubedegree, int dim, int debuglevel);

inline void ExtendToHorizon(...);   // While (S.T < H) Extend(...)
inline void RebuildScaffold(...);   // drop the time grid and lay a new one from E₀
```

**StepA bounds the stage length from above.**


**Seven steps** (line by line against the paper's pseudo-code):

| step | action | paper |
|------|------|------|
| 1 | `(h, F_A) ← StepA(E_m, H - t_m, ε₀)`, the **upper bound** on the step | Subsection 4.1 |
| 2 | `chain.step(h)`: StepB takes one step, no longer than h; halved if it fails | Subsection 4.2 |
| 3 | `F₁ ← StepAatStep(E_m, Δt, ε₀)` (`F_A` reused when `Δt = h`) | eq. (2.2) |
| 4 | `μ ← logNorm(J_f(F₁))` | eq. (2.3) |
| 5 | `E₁ ← chain.hull() ∩ F₁` (stepB=1/2 also intersect the ball of Lemma 2.1) | Subsection 4.2 |
| 6 | `h_tube ← computeHTube(F₁, Δt, δ, μ, tubedegree)` | eq. (2.6) |
| 7 | append `t_{m+1}, F_{m+1}, E_{m+1}` and initialise `G_{m+1}` | Section 4 |

**Initial values of the new mini-scaffold:** `ℓ = 0` (one mini-step),
`E = {E₀, E₁}`, `F = {F₁, F₁}`, `mu = {μ, μ}`, `nTube = 1`.

**Checks under `debuglevel = 1`** (these should never fire): `E₀ ⊆ F₁`,
`E₁ ⊆ F₁`.

---

## 📄 **refine.h** — `S.Refine(ε₀)` (Section 5, **the main algorithm**)

**What it does.**  Refines `S.E` and `S.F`, keeping `S.t` fixed, until
`wmax(S.E[m]) ≤ ε₀` (the scaffold is then **ε₀-small**).

### Constants and enum

```cpp
enum RefineMethod {
    REFINE_TUBE      = 0,   // Bisect + EulerTube + Split (the Refine of Section 5, **default**)
    REFINE_BISECT    = 1,   // Bisect + Split, tube switched off
    REFINE_SPLITONLY = 2    // Split only; the stages are the steps StepB takes
};


```

`ℓ_i` controls how fine a grid `EulerTube` works on, and the bound of Lemma 6.1
is reached well inside this limit.

### `Evaluate` — walk the scaffold again, on the grid it already has

```cpp
inline void Evaluate(Scaffold& S, IMap f, SolutionChain& chain, int degree,
                     int stepBtype, int dim, double stepAdelta);
```

Refine changes two things about a scaffold: `Split` shrinks `E₀`, and `Bisect`
asks for stage i to be crossed in `2^{ℓ_i}` mini-steps instead of one.
`Evaluate` is what makes either of them take effect.

**It fills** `G[i].E[j]`, `G[i].F[j]`, `G[i].mu[j]` for every stage i and
mini-step j, plus `S.E[i]` and `S.F[i]`.

**One chain, all the way.**  After `chain.restart(S.E[0])` it uses
`marchTo(t_j)` to land exactly on each mini-step time.  No time map, no second
chain, and the doubleton is never re-enclosed.

**Each mini-step is certified on its own:**
`F[j] = StepAatStep(f, E[j-1], h₁, δ)`.  Halving the step narrows the full
enclosure `F` and with it `μ`, so the criteria of Subsection 5.3 are met
sooner and the tube of Theorem 5.1 has a finer grid to lie on.  The stage's own
`S.F[i]` was certified over the whole stage from a box no smaller than this
one, so it is equally valid for a single mini-step.

**Refining only tightens.**  The time grid does not move.  What the scaffold
already holds encloses the **same solution sets** under an `E₀` that was no
smaller, so the new values are intersected into the old ones.

### `Bisect(i)` (Subsection 5.1)

```cpp
inline bool Bisect(Scaffold& S, int idx, int levels = 1);   // false = ELL_MAX reached
```

**★ For the Lohner propagation, Bisect is a step-size halving with the
intermediate point inserted.**  The span `[t_{i-1}, t_i]` of stage i does not
change, but it is crossed in **twice as many, half as long** steps, and each of
them is certified by StepA on its own, so the full enclosures `F^i` shrink with
the step — that, rather than the end-enclosure of the stage, is what Bisect
buys.

Bisect only refines the **grid** (`ℓ_i += levels`); the finer chain of
mini-quads is produced by the next `Evaluate`, with the same doubleton, rather
than by re-integrating box by box — the latter would re-introduce the wrapping
effect.

**The `levels` argument.**  Refine already knows how fine the grid has to
become for Lemma 2.2 to hold (`ceil(log2(h / h_tube))`), so it gets there in
one move instead of one level per phase.

### `EulerTube(i)` (Subsection 5.2, Theorem 5.1)

```cpp
inline int EulerTube(Scaffold& S, int idx, int dim, IMap f, int debuglevel);
                                                    // returns the number of segments
```

With `E^i[0] ⊆ Ball_{q₀}(r₀)`, μ a logNorm bound on the stage, and
`h₁ = Δt_i / 2^ℓ_i`:

| Theorem 5.1 | formula | used as |
|-------------|------|------|
| (a) | `Ball_{q_j}(r₀ e^{jμh₁} + δ)` | end-enclosure → intersected into `G.E[j]` |
| (b) | `Box(Ball_{q_{j-1}}(r), Ball_{q_j}(r))`, `r = δ + max(r_{j-1}e^{μh₁}, r_{j-1})` | full enclosure → intersected into `G.F[j]` |

where `q₀, q₁, ...` is the Euler polygon of step size `h₁`.

**Condition (2.7) and tube segments.**  Theorem 5.1 applies only while the
Euler point `q_{j-1}` lies in the current end-enclosure `E^i[j-1]`
(Subsection 2.6).  When that test fails the tube is **restarted** at the
midpoint of `E^i[j-1]`; each restart begins a new **tube segment**, the count
goes into `G.nTube`, and they add up to `N^tube(S)` in the stopping criteria.

**Which μ.**  Theorem 5.1 needs **one** μ for the whole tube, so the largest of
the stage's mini-step bounds is taken (skipping the sentinel at index 0).

`r₀` is `euclideanRadius(box)`: Lemma 2.1 and Theorem 5.1 are both in the
Euclidean norm, so every ball radius is that of the circumscribed ball.

### `Split()` (Subsection 5.3)

```cpp
inline void Split(Scaffold& S, int dim, const std::vector<double>& target);
```

Halves `S.E[0]` towards `target` (towards the midpoint when no target is
given).  **The complement is not stored here**: the caller — EndCover, or the
boundary sweep — sees that `S.E[0]` has shrunk and re-queues the remaining
sub-boxes.  Those are the "split-off scaffolds returned by Refine" of the
pseudo-code.

### The three stopping criteria (Subsection 5.3)

| function | paper | formula | meaning |
|------|------|------|------|
| `muMaxOf(S)` | `μ_max` | | the global logNorm bound of the scaffold, item [G4] |
| `nTubeOf(S)` | `N^tube(S)` | `Σ_i max(1, G[i].nTube)` | total tube segments over all stages |
| `StageFreeze(S, i, eps0)` | eq. (5.4) | `e^{μ_max·T}·δ_i·N^tube < ε₀/8` | once true, stage i is not refined further |
| `NeedSplit(S, eps0, dim)` | eq. (5.5) | `e^{μ_max·T}·wmax(E₀) > ε₀/(2√n)` | true means `E₀` has to be split |
| `CanTerminate(freezeFlag, needSplit)` | eq. (5.6) | `(∧_i FreezeFlag[i]) ∧ ¬NeedSplit()` | Refine may wrap up |

### `Refine`

```cpp
inline void Refine(Scaffold& S, int dim, IMap f, SolutionChain& chain,
                   int method, int stepBtype, int degree,
                   double eps0, int tubedegree, int debuglevel,
                   const std::vector<double>& target = {});
```

**Early return.**  If `wmax(S.E.back()) ≤ ε₀` already holds, Refine returns at
once, without a single `Evaluate`.

**Handing the chain back.**  When Refine returns, `chain` is anchored at the
current `E₀` and standing at `S.T`, so the next `Extend` can carry straight on.

**The main loop** (parts A / B / C of Section 5):

```
tubes()                                    ← EulerTube on every stage fine enough
While (S is not ε₀-small)
    ── Part A: refine every stage that is not frozen        (method 0 / 1)
       If (StageFreeze(i))              FreezeFlag[i] ← 1, skip        (5.4)
       If (h_i ≤ h^i_euler and tube on) δ_i ← δ_i/2; update h^i_euler  ([G5])
       Bisect(i, levels)                                              (5.1)
         — freeze the stage when ELL_MAX is out of reach
       Bisected → Evaluate()   ← refining takes effect only on a new walk
       tubes()
    ── Part B: If (NeedSplit()) { Split(); reevaluate() }   (unconditional for method 2)
    ── Part C: If (CanTerminate() and still too wide) { Split(); reevaluate() }
```

**Why `Evaluate()` comes before Part B.**  Part A only changes `ℓ_i`; nothing
moves until the scaffold is walked again.  And `NeedSplit()` in Part B reads
`μ_max`, which has to be the refined value — otherwise Bisect would have no
influence on the decision to split.

**`reevaluate()` — re-evaluating after a Split:**

| case | what happens | why |
|------|------|------|
| some stage has been bisected (`ℓ_i > 0`) | `Evaluate()`, **keeping the time grid** | the grid carries Bisect's work |
| nothing has been refined | `RebuildScaffold()`, **a new grid** | nothing to keep, and a smaller box admits longer steps |

**Progress guard.**  A pass that has neither split nor narrowed the scaffold
forces a `Split` and re-evaluates; when `E₀` has degenerated to a point the
loop breaks.

---

### ⭐ Layer 4: the covering algorithm

---

## 📄 **endcover.h** — `EndCover_f(B0, H, ε₀) → C` (Subsection 3.3)

**What it does.**  The global algorithm: a work queue that subdivides `B0`.

```
C  ← {}                            the ε₀-cover
Q0 ← {Init(B0)}                    the scaffold queue
While (Q0 ≠ {})
    S ← Q0.pop()
    Q0.push(S.Refine(ε₀))
    While (S.t.back() < H)
        S.Extend(H, ε₀)
        Q0.push(S.Refine(ε₀))
    C.push(S.E.back())
Return C
```

### Types and enum

```cpp
enum CoverSet { COVER_ALL = 0, COVER_LEAVES = 1 };

struct EndCoverResult {
    std::vector<IVector> cover;      // C
    std::vector<IVector> initBoxes;  // the sub-boxes of B0 that were integrated
};

using EndEncFn = std::function<std::pair<IVector, IVector>(
                     const IVector&, double, const std::vector<double>&, double)>;
```

| `CoverSet` | meaning |
|-----------|------|
| `COVER_ALL` (default) | every box that is processed puts its end-enclosure into C |
| `COVER_LEAVES` | only boxes that needed no further splitting — the ε₀-cover in the strict sense of eq. (1.3), every box satisfying `wmax ≤ ε₀` |

### Functions

| function | description | in | out |
|------|------|------|------|
| `midpoint_vec(box)` | midpoint, as a vector of doubles | IVector | vector\<double\> |
| `boxesEqual(a, b, tol=0)` | endpoint-by-endpoint comparison | two IVector | bool |
| `splitAllDims(box)` | bisects along **every coordinate of positive width** | IVector | `vector<IVector>` (2^d of them, d = number of such coordinates) |
| `EndCover(B0, eps0, H, endEnc, maxSplits=2e6, coverSet=COVER_ALL)` | the main loop | — | EndCoverResult |

**★ `splitAllDims` skips degenerate coordinates.**  Splitting a coordinate of
zero width would only produce duplicate children.  When every coordinate is
degenerate (the box is a point) an empty vector is returned.

### The `endEnc` protocol

The callback returns `(ulB, olB)`:

| value | meaning |
|--------|------|
| `ulB` | `S.E[0]` after `Refine` — the part that was actually integrated |
| `olB` | `S.E.back()` — the end-enclosure at time H |

**The decision rule.**  `ulB == B` means all of B is covered (resolved);
`ulB ⊊ B` means `Refine` called `Split`, so the complement of B is still
uncovered → `splitAllDims(B)` and re-queue.  Those children are the
"split-off scaffolds returned by Refine" of the pseudo-code.

**Exceptions.**  When `endEnc` throws, B is **not** dropped: it is subdivided
and re-queued.  Only a box that is already a point, and so cannot be
subdivided, propagates the exception upwards.

---

## 📄 **boundary.h** — the scaffold driver and Boundary mode

### `runScaffold` — the engine shared by both modes

```cpp
inline Scaffold runScaffold(const IVector& B, IMap F,
                            const MiniScaffold& templateG,
                            double eps0, double delta, int degree, double H,
                            int method, int stepBtype, int tubedegree,
                            int debuglevel, int dim,
                            const std::vector<double>& target);
```

**The inner loop of EndCover:**

```
While (S.T < H)  S.Extend(H, ε₀)
S.Refine(ε₀)
```

`Refine` reads the whole scaffold and refines **all** of its stages, so it is
called once the scaffold has reached the horizon.

**The carried doubleton.**  `runScaffold` owns one `SolutionChain`; `Extend`
advances it one step per stage, and after `Split` / `Bisect` `Refine` walks it
again from `E₀`.  From `E₀` to `E_m` there is only this one doubleton, and it is
never re-enclosed on the way.

**Finiteness check.**  After each stage the endpoints of `S.E.back()` are
checked for finiteness; if they are not finite an exception sends the box back
to EndCover for subdivision.

### Boundary mode (`mode = 1`)

Covers only the **boundary** of `B0`: the four edges in 2D, the 2n faces in nD.

| function | dimension | how |
|------|------|------|
| `TwoDimEncAlgo(...)` | n = 2 | sweeps edge by edge, bisecting a segment where needed, down to a single point |
| `ThreeDimEncAlgo(...)` | n ≥ 3 | one `endcover::EndCover` run per face, 2n of them |

For a flow the image of the boundary bounds the image of the interior, so when
`End(B0,H)` is simply connected the hull of the boundary cover is an enclosure
of `End(B0,H)`.

### 2D convex hull helpers (for plotting)

| function | description |
|------|------|
| `struct Point2D { double x, y; }` | a point in the plane |
| `cross2D(o, a, b)` | cross product |
| `convexHull2D(pts)` | Andrew monotone chain |
| `saveConvexHull2D(boxes, fname, verbose)` | takes the four corners of every box, computes the hull, writes it as a closed polygon |

---

### ⭐ Layer 5: the driver

---

## 📄 **EnCorelib.cpp** — command line and output

### Command line

```
./EnCorelib  output_mode mode method stepB tubedegree n
           var_1 ... var_n  f_1 ... f_n
           eps order T debug
           lo_1 hi_1 ... lo_n hi_n
```

| argument | values |
|------|------|
| `output_mode` | 0 time / 1 +enclosure / 2 +counts / 3 +E0,E1 files / 4 +plot files |
| `mode` | 0 = EndCover, 1 = Boundary |
| `method` | 0 = tube, 1 = bisect, 2 = splitonly |
| `stepB` | 0 = lohner, 1 = +logNorm, 2 = lohnerHO+logNorm, 3 = lohnerHO |
| `tubedegree` | 0 = order−1, 1 = Euler, p ≥ 2 = Taylor-p |
| `eps` | the ε of eq. (1.3) |
| `order` | the Taylor order k of StepB (the paper uses 20) |
| `T` | the horizon H |
| `debug` | 0 silent, 1 print the containment checks |

### Internal functions

| function | description | in | out |
|------|------|------|------|
| `modeName / methodName / stepBName` | enum value → readable name | int | `const char*` |
| `saveBoxes(boxes, fname, verbose)` | writes a box file (scientific, 17 digits) | — | — |
| `midRadString(B)` | `"(mid,...) +- (half,...)"`, the B1 column of Table 2 | IVector | string |
| `appendLog(...)` | appends one tab-separated row to `out.txt` | — | — |
| `Init(B, delta, H)` | `Init(B0)` — a scaffold with no stages (Subsection 3.3) | — | Scaffold |
| `usage()` | prints the usage message | — | — |

**Sentinels set by `Init`:** `G.mu = {1000.0}`, `G.hTube = H`, `ℓ = 0`,
`nTube = 1`.

**`std::cout << std::unitbuf`** on the first line of main, so that output
already produced is not left sitting in the buffer.

### Output levels (progressive: level k also prints every lower level)

| level | adds |
|------|-----------|
| 0 | `time(ms)=` |
| 1 | `Hull(T)=`, `B1=(mid) +- (half-width)`, `wmax(B1)=`, and a header echoing the problem |
| 2 | `#(C)=`, `#E0=`, `max_{B in C} wmax(B)` |
| 3 | writes `E0.txt`, `E1.txt`, `convex_hull.txt` (2D) |
| 4 | writes `E_0.txt`, `E_1.txt` (with snapshots at t = 0.1/0.4/0.7 and, in 2D, the corner trajectories) |

| file | contents |
|------|------|
| `E0.txt` | the initial sub-boxes that were integrated, one per box of C |
| `E1.txt` | the cover C at time H — **and nothing else** |
| `convex_hull.txt` | the 2D convex hull of C, as a closed polygon |
| `E_0.txt` / `E_1.txt` | the two above, plus the intermediate images and the 2D corner trajectories |
| `out.txt` | a cumulative log, one row appended per run |

---

### ⭐ Supporting directories

---

## 📁 **repro/** — the tables and figures of the paper

| script | description |
|------|------|
| `table2.sh` | regenerates the `Ours(ε)` rows of Table 2 (`make table2`) |
| `table2_paper.txt` | Table 2 as printed, transcribed verbatim |
| `figure1.sh` | the data behind Figure 1: the Eijgenraam ε-cover and the pendulum orbit |
| `figure4.sh` | the data behind Figure 4: E0/E1 for the six examples |
| `check.sh` | the self-test, 18 assertions (`make check`); `BASELINE=1` re-takes the baseline |
| `verify_cover.sh` | integrates a grid of sample points of B0 separately with `Encoretraj` and confirms that every image lies inside some box of C; runs all three `method` values |
| `competitors/` | the reference integrator (`make capdref`) and the VNODE-LP / INTLAB driver templates |

## 📁 **tools/**

| file | description |
|------|------|
| `trajectory.cpp` → `Encoretraj` | validated trajectory printer (Figure 1, top row; also the independent reference for `verify_cover.sh`) |
| `run.sh` | invokes EnCorelib from the Makefile; splits `var`/`ff` on commas here, so the shell never sees the parentheses in `x*(28-z)-y` |
| `boxio.py` | box-file reader |
| `plot_boxes.py` | Figure 4 style plots (two panels in 2D, three projections in 3D) |
| `plot_eijgenraam.py` | Figure 1, bottom row |
| `plot_pendulum.py` / `pendulum_ode45.m` | Figure 1, top row (Python RK45 / MATLAB ode45) |

---

## 📊 Data flow

```
command line / examples/*.mk
    ↓
IMap F  ←  symparse.h (optional SymEngine expansion)
    ↓
┌──────────────────────────────────────────────────────┐
│ mode = 0: EndCover (endcover.h)                      │
│   work queue Q0 of sub-boxes of B0                   │
│   for each B:                                        │
│     ┌────────────────────────────────────────────┐   │
│     │ runScaffold (boundary.h)                    │   │
│     │   carries one doubleton (stepAB.h)          │   │
│     │   While (S.T < H)                           │   │
│     │     Extend (extend.h)                       │   │
│     │       ├ StepA   (stepAB.h) ← step bound     │   │
│     │       ├ StepB   (stepAB.h → Lohner form)    │   │
│     │       └ h_euler (stepAB.h)                  │   │
│     │   Refine (refine.h)                         │   │
│     │       ├ Evaluate ← the same SolutionChain   │   │
│     │       ├ Part A: Bisect(ℓ++) / EulerTube     │   │
│     │       ├ Part B: NeedSplit → Split           │   │
│     │       └ Part C: CanTerminate → wrap up      │   │
│     └────────────────────────────────────────────┘   │
│   returns (ulB, olB)                                 │
│   ulB ⊊ B ?  → splitAllDims(B) and re-queue          │
└──────────────────────────────────────────────────────┘
    ↓
CoverResult { E0, C, hull }
    ↓
output: time / Hull(T) / B1 / #(C) / E0.txt / E1.txt / ...
```



## 📋 Call graph

### One EndCover run

```
main (EnCorelib.cpp)
  └─ endcover::EndCover(B0, ε, H, endEnc)
       └─ endEnc(B, ε, m(B), H)                       ← lambda in main
            └─ runScaffold(B, ...)                     (boundary.h)
                 ├─ SolutionChain (carried, anchored at E₀)   (stepAB.h)
                 ├─ ExtendToHorizon(...)               (extend.h)
                 │    └─ Extend(...)  once per stage
                 │         ├─ StepA(f, E₀, H-t, ε)      (stepAB.h) ← step bound
                 │         ├─ chain.step(h)             (stepAB.h) ← one StepB step
                 │         ├─ StepAatStep(...)          (stepAB.h) ← certify at the actual Δt
                 │         ├─ logNorm(computeJacobian(...)) (types.h)
                 │         └─ computeHTube(...)         (stepAB.h)
                 └─ Refine(...)                        (refine.h)
                      ├─ Evaluate(...)            ← the same SolutionChain
                      │    ├─ StepAatStep(...)      (stepAB.h, per mini-step)
                      │    └─ logNormBallTighten()  (stepAB.h)
                      ├─ EulerTube(i)             (Theorem 5.1, method 0)
                      ├─ Bisect(i, levels)        (ℓ_i += levels, method 0/1)
                      ├─ NeedSplit / StageFreeze / CanTerminate
                      └─ Split(S, dim, target)
```

### Where each enclosure comes from

```
E_i (stage endpoint)    ← the chain's hull at t_i ∩ F_i ∩ the previous E_i
E^i[j] (mini-step)      ← the same chain's hull at t_j ∩ F^i[j]
F_i (stage enclosure)   ← StepAatStep (eq. 2.2), certified at the Δt actually taken
F^i[j]                  ← StepAatStep (per mini-step), then tightened by Thm 5.1(b)
μ^i[j]                  ← logNorm(J_f(F^i[j]))
```

---

## 📌 Conventions

1. **Time is absolute.**  The argument of `SolutionChain::marchTo(t)` is an
   absolute time and must not decrease; after `restart(E₀)` time starts at 0
   again.
2. **`mu[0]` is the sentinel 1000**, not a real logNorm bound.  Every loop
   starts at index 1.
3. **Every ball radius is `euclideanRadius`**, not `wmax/2` (Lemma 2.1 is in
   the Euclidean norm).
4. **`IntersectB` returns `B1[i]` on an empty coordinate** — both arguments
   enclose the same solution set, so an empty intersection can only come from
   rounding.
5. **One set of names.**  Every quantity in the code carries the paper's
   notation; there are no compatibility aliases.
6. **`-frounding-math` is required**; the directed roundings every interval
   bound depends on rely on it, and both build systems pass it.
7. **How `#(C)` is counted.**  `COVER_ALL` by default, i.e. every box that was
   **processed** counts, including those later cut by `Split`; `COVER_LEAVES`
   keeps only the boxes that needed no further splitting.  `output_mode ≥ 2`
   also prints `max_{B in C} wmax(B)`.

---

## 🎯 Suggested reading order

**To understand the algorithm:**
1. `types.h` — box utilities, the `Scaffold` and `MiniScaffold` structs
2. `stepAB.h` — StepA / StepB / `h_euler` (Section 4 and eq. 2.6)
3. `stepAB.h::SolutionChain` — how a chain of mini-quads is evaluated (Subsection 5.1)
4. `extend.h` — how a stage comes into being
5. `refine.h` — the main algorithm (all of Section 5)
6. `endcover.h` + `boundary.h` — the subdivision and the driver

**To just build and run:**
`README.md` §2 to build → `make check` → `make table2` → compare with
`repro/table2_paper.txt`

**For the numerics:**
1. `types.h::euclideanRadius` (the norm of Lemma 2.1)
2. `stepAB.h::StepAatStep` (the admissibility test of eq. 2.2)
3. `refine.h::EulerTube` (Theorem 5.1 and condition 2.7)
4. `repro/verify_cover.sh` (the sample-point check)

---

## 🔤 Notation (paper ↔ code)

| paper | code |
|------|------|
| `B`, `B₀` | `IVector` |
| `w(B)`, `wmax(B)`, `m(B)` | `width`, `wmax(B)`, `midpoint(B)` |
| `Box(S)` | `hullBox(A,B)` |
| `μ(A) = μ₂(A)` | `logNorm(J, n)` |
| `S = (t, E, F, G)` | `struct Scaffold { t, E, F, G }`, `S.m()`, `S.T()` |
| `G[i] = ((ℓ_i, E^i, F^i), μ^i, (δ_i, h^i_euler))` | `struct MiniScaffold { ell, E, F, mu, delta, hTube, nTube }` |
| `StepA(E₀,H,δ)`, `StepB(E₀,h,F₁)` | `StepA(...)`, `StepB(...)` |
| `S.Extend(H,ε₀)`, `S.Refine(ε₀)` | `Extend(...)`, `Refine(...)` |
| `S.Bisect(i)`, `S.EulerTube(i)`, `S.Split()` | `Bisect`, `EulerTube`, `Split` |
| `StageFreeze(i)` — eq. (5.4) | `StageFreeze(S, i, eps0)` |
| `NeedSplit()` — eq. (5.5) | `NeedSplit(S, eps0, n)` |
| `CanTerminate()` — eq. (5.6) | `CanTerminate(freezeFlag, needSplit)` |
| `N^tube(S)`, `μ_max` | `nTubeOf(S)`, `muMaxOf(S)` |
| admissibility test — eq. (2.2) | `StepAatStep(...)` |
| `h_euler(H,M,μ,δ)` — eq. (2.6) | `hEuler(H,M,mu,delta)` |
| the mini-quad chain `E_{i-1} → ... → E_i` | one `SolutionChain`, walked by `Extend` / `Evaluate` |
| `EndCover_f(B₀,H,ε₀) → C` | `endcover::EndCover(...)` |
| `ε₀` | `eps0` (`eps` on the command line) |
| `C`, `B₁ = Box(C)` | `CoverResult::C`, `CoverResult::hull` |

---

## 📝 The options

### The three Refine strategies

`method` and `stepB` are **orthogonal**: `stepB` chooses how **one step** is
enclosed (all four variants are Lohner's method), `method` chooses **which
operators Refine uses** to refine the scaffold once those steps exist.

| method | name | contents |
|--------|------|------|
| 0 | `tube` | Bisect + EulerTube + Split + (5.4)–(5.6); the Refine written out in Section 5, **the default** |
| 1 | `bisect` | the same with EulerTube switched off |
| 2 | `splitonly` | Split alone: the stages stay the steps Lohner's method took and are not subdivided, so the accuracy comes from shrinking E₀ |

### The four StepB representations

| stepB | representation |
|-------|------|
| 0 | Lohner's method (doubleton + QR preconditioning), **the default** |
| 1 | the same, plus the ball of Lemma 2.1 |
| 2 | Lohner's method with a higher-order a-priori enclosure, plus the ball |
| 3 | as 2, without the ball |

### The tube threshold `tubedegree`

| tubedegree | meaning |
|------------|------|
| 1 | the Euler threshold `h_euler`, equation (2.6) |
| p ≥ 2 | the degree-p Taylor threshold |
| 0 | `order - 1` |

### The rest

| argument | meaning |
|------|------|
| `eps` | the error bound ε of the cover |
| `order` | the Taylor order k of StepB (the paper uses 20) |
| `T` | the horizon H |
| `mode` | 0 = EndCover, 1 = Boundary |
| `output_mode` | 0..4, see the output levels above |

Everything runs in IEEE double precision.

---

**Covers:** every header of EnCorelib, the driver, and the `repro/` and
`tools/` directories.
**Size:** 7 headers, the driver, and the scripts and tools.
**Self-test:** `make check` and `make verify` both pass.
