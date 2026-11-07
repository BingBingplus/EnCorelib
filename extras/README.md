# extras/ — material that is **not** part of the paper

Nothing in this directory is referenced by the paper.  It is kept because the
two demonstrations below are useful when reasoning about the wrapping effect,
and because the recorded runs document how the figures in an appendix draft
were produced.

Build the two programs with

```bash
make extras-all
```

The sources are in `src/extras/`; the programs are written into this
directory.  They link against the same headers as `EndCorelib`, so they follow
any change to `src/refine.h` / `src/stepAB.h`.

| | |
|---|---|
| `src/extras/moore_direct_test.cpp` | Moore's rotation `x' = y, y' = -x` on `Box(1,0)(0.05)`. Compares `stepBlohnerHO` (stepB = 3) against `stepBlohner` (stepB = 0) step by step, with the enclosure put back into a box after every step, so that the raw wrapping-effect growth of each stepper is visible. The exact flow is a rigid rotation, so the true reachable set never grows; a naive interval stepper grows by `sqrt(2)` per step. |
| `src/extras/moore_endcover_direct.cpp` | the same problem driven through `runScaffold` / `EndCover` with stepB = 3, writing the cover under a label. |
| `moore_runs/`, `moore_direct_runs/`, `moore_endcover_direct_runs/` | recorded outputs and figures; the MATLAB scripts that draw them are in `src/extras/` under the same directory names and are run from the data directory, e.g. `cd extras/moore_runs` and then `run ../../src/extras/moore_runs/moore_figure.m` in MATLAB |
| `eijgenraam_runs/` | earlier runs of Eijgenraam's example; the version used for Figure 1 of the paper is produced by `src/repro/figure1.sh` instead |
| `reviewer_runs/` | scratch output from a reviewer response |

If you are packaging a minimal release of the code, this directory can be
deleted; `make`, `make check`, `make table2`, `make figure1` and
`make figure4` do not use it.
