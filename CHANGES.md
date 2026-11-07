# Changes

## November 2026

Changes since the first public version (September 2026).

### Enclosures

- A step of `SolutionChain` (stepAB.h) could end later than the caller
  allowed. When the error control of StepB chose a longer step than the
  previous one, the chain took it, so `marchTo(t)` could stop past `t`, and
  `Extend` could make a stage longer than the step StepA had certified.
  `Evaluate` then reported the enclosure of a later time at a mini-step
  node. Every step is now capped at the length the caller allows. The error
  showed most in sub-boxes that Refine had shrunk: in Eijgenraam's example
  with eps = 0.1, 34 of the 79 boxes of the cover missed the image of their
  own sub-box, and `B1` reached y = 1.21 although the end set ends at y = 1.
- The step-size threshold of the Euler tube is now equation (2.6) of the
  paper. `hEuler` had both cases mu > 0 and mu < 0 in a form about twice as
  permissive, `computeHTube` returned H for mu = 0, and the default
  `tubedegree = 0` meant `order - 1`, a degree-19 Taylor threshold applied
  to the Euler polygon that `EulerTube` follows. With that default,
  `method = 0` could return boxes that missed the solution set (Eg3 at
  H = 1 with eps = 0.1: 4 of 21 boxes). `tubedegree` now accepts 0 or 1,
  both the Euler tube.

A per-box check, `make audit` (new), in which every box of C must meet the
validated images of the centre and the corners of its own sub-box of B0,
finds no violation in 35 settings (the six examples of Table 1 and
Eijgenraam's example, with method 0, 1 and 2 at stepB 0 and with method 2 at
stepB 1 and 3). Before these two changes it found violations in 5 of the 7
problems it was run on, with method 0 and with method 1.

With the corrected threshold `method = 0` bisects more and takes longer;
section 6 of the user guide has the timings.

### Building

- `stepAB.h` called `sup()`, which only the filib interval type of CAPD
  provides. With CAPD built on its own interval type, as on macOS with
  clang, the build stopped at "use of undeclared identifier 'sup'". It now
  calls `rightBound()` and compiles with either type (checked with g++ 13
  against filib and with clang 14 and libc++ against CAPD's own intervals).

### Input

- Constants may be written in scientific notation, `2.5e-3`; they are
  evaluated as quotients such as `(2.5/10^3)`, which the parser encloses in
  interval arithmetic. Before, `1e-3` was mangled into `1*e-3`.
- A number directly followed by `(` is multiplied by it, `2(x+1)`, and digits
  inside a name (`x1`, `k2a`) stay part of the name.
- A right-hand side that cannot be read is reported with the parser's
  message and exit status 2; before, the program aborted. A failure in
  Boundary mode is reported with exit status 3 instead of aborting.
- `src/examples/*.mk` and the scripts in `src/repro/` write the right-hand sides
  multiplied out, the form in which they are evaluated when SymEngine is
  installed and the form the numbers of the paper come from. Results no
  longer depend on whether SymEngine was found: written as
  `10*(y-x), x*(28-z)-y, x*y-8*z/3`, Lorenz to H = 5 with stepB = 0 took
  116 scaffolds and 78 boxes in a build without SymEngine, and 95 scaffolds
  and 64 boxes multiplied out.

### Progress reports and messages

- `debug` now selects progress reports on standard error: 1, one line per
  scaffold taken from the queue (depth, widths, queue, share of B0 covered,
  outcome); 2, also one line per scaffold laid out and per pass of Refine; 3,
  also one line per stage and the containment self-checks that `debug = 1`
  used to switch on.
- Messages name the time reached, for example
  `Extend: StepB cannot advance at t=26.4064 (H=100)`. When a sub-box of B0
  can no longer be halved in double precision, the program stops with a
  message saying so and suggests a smaller T.

### Documentation

- `doc/userguide.pdf`, with the scripts that regenerate its examples
  (`src/doc/run_examples.sh`), its tables (`src/doc/compare_options.sh`,
  `src/doc/make_tables.py`) and its figures (`src/doc/make_figures.py`).
- `src/examples/dampedPendulum.mk`, a problem that is not in the paper.
- The program is called `EndCorelib`, its companion programs `EndCoretraj`
  and `EndCoreaudit`; the sources are in `src/`.
