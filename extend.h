#pragma once
// ===========================================================================
//  extend.h  --  S.Extend(H, eps0)                 (paper, Section 4)
//
//  Extend appends one stage to the scaffold S and initialises its
//  mini-scaffold:
//
//      (h, F_{m+1})  <- StepA(S.E.back(), H - S.T, delta)      Subsection 4.1
//      E_{m+1}       <- StepB(S.E.back(), h, F_{m+1})          Subsection 4.2
//      G_{m+1}       <- ((0, (S.E.back(), E_{m+1}), (F_{m+1})),
//                        mu, (delta, h_euler))
//
//  with mu >= mu_2(J_f(F_{m+1})) and h_euler = h_euler(h, ||f^[2](F)||, mu, delta).
//
//  StepA fixes the upper bound on the stage length.  StepB -- Lohner's method
//  -- has a step of its own, the one its error control asks for, so the stage
//  runs for
//
//      Dt = min( h , the step Lohner's method wants )
//
//  and is certified by StepA at that length: F is recomputed for Dt, and a
//  full enclosure certified for the longer h encloses the shorter stage as
//  well.  A step that cannot be taken is halved and tried again, so the stage
//  is shorter where the flow demands it.
//
//  ExtendToHorizon is the Extend loop of the pseudo-code of Subsection 3.3,
//  run until t_m reaches the horizon.
// ===========================================================================

#include "stepAB.h"

// Reaching this many stages means the step sizes have collapsed; the box is
// subdivided instead.
static constexpr int MAX_STAGES = 20000;

// ---------------------------------------------------------------------------
//  Extend
// ---------------------------------------------------------------------------
//  f           the vector field
//  S           scaffold, modified in place
//  chain       the doubleton carried along the scaffold (StepB, stepAB.h)
//  eps0        target error bound of the enclosing EndCover call; it is also
//              StepA's delta, the componentwise bound on the k-th remainder
//  delta       tube tolerance (delta_i of the new mini-scaffold)
//  degree      Taylor order k used by StepB
//  H           the horizon; the stage stops there at the latest
//  stepBtype   StepB variant, see stepAB.h
//  tubedegree  1 = Euler tube, p >= 2 = Taylor tube of degree p
//  dim         n
//  debuglevel  1 = report violated containments (these should never fire)
inline void Extend(IMap f, Scaffold& S, SolutionChain& chain,
                   double eps0, double delta, int degree, double H,
                   int stepBtype, int tubedegree, int dim, int debuglevel)
{
    const IVector E0    = S.E.back();
    const double  tPrev = S.t.back();
    const double  left  = H - tPrev;
    if (!(left > 0.0)) return;

    // -- 1. StepA: the largest delta-admissible step, and its full enclosure
    double  hA     = left;
    IVector FA;
    bool    haveFA = false;
    try {
        std::pair<double, IVector> pr = StepA(f, E0, left, eps0);
        if (pr.first > 0.0) {
            hA     = std::min(pr.first, left);
            FA     = pr.second;
            haveFA = true;
        }
    } catch (const std::exception&) {
        // No admissible step at this box size; what is left of the horizon
        // bounds StepB's step instead.
    }

    // -- 2. StepB: one step of Lohner's method, no longer than hA ----------
    double bound = hA;
    bool   moved = false;
    for (int retry = 0; retry < 40 && !moved; ++retry) {
        moved = chain.step(bound);
        if (!moved) {
            bound *= 0.5;
            if (!(bound > 1e-13 * std::max(1.0, H))) break;
        }
    }
    if (!moved) throw std::runtime_error("Extend: StepB cannot advance");

    const double tNow = std::min(chain.time(), H);
    const double h    = tNow - tPrev;
    if (!(h > 0.0)) throw std::runtime_error("Extend: stage of zero length");

    // -- 3. Full enclosure of the stage, at the length actually taken ------
    IVector F1;
    if (haveFA && !(h < hA)) {
        F1 = FA;
    } else {
        try {
            F1 = StepAatStep(f, E0, h, eps0);
        } catch (const std::exception&) {
            if (!haveFA) throw;
            F1 = FA;               // certified for hA >= h, so valid here too
        }
    }

    if (debuglevel == 1) {
        for (int i = 0; i < dim; ++i) {
            if (F1[i].leftBound() > E0[i].leftBound() ||
                F1[i].rightBound() < E0[i].rightBound()) {
                std::cout << "[Extend warn] E0 not contained in F1 at coord "
                          << i << ", t=" << tPrev << "\n";
            }
        }
    }

    // -- 4. logNorm bound on F1 -------------------------------------------
    const double mu = logNorm(computeJacobian(f, F1), dim);

    // -- 5. End-enclosure of the stage ------------------------------------
    IVector E1 = chain.hull();
    if (stepBtype == 1 || stepBtype == 2)
        E1 = logNormBallTighten(f, degree, h, E0, mu, dim, E1);
    E1 = IntersectB(E1, F1);

    if (debuglevel == 1) {
        for (int i = 0; i < dim; ++i) {
            if (E1[i].leftBound() < F1[i].leftBound() ||
                E1[i].rightBound() > F1[i].rightBound()) {
                std::cout << "[Extend warn] E1 not contained in F1 at coord "
                          << i << ", t=" << tNow << "\n";
            }
        }
    }

    // -- 6. Tube step-size threshold for the new stage ---------------------
    const double hTube = computeHTube(f, F1, h, delta, mu, tubedegree, degree);

    // -- 7. Append the stage ----------------------------------------------
    S.t.push_back(tNow);
    S.F.push_back(F1);
    S.E.push_back(E1);

    S.G.push_back(S.G[0]);                   // template, then override
    MiniScaffold& G = S.G.back();
    G.mu    = {mu, mu};   // one entry per mini-step endpoint (ell = 0 -> 1 step)
    G.delta = delta;
    G.hTube = hTube;
    G.ell   = 0;
    G.E     = {E0, E1};
    G.F     = {F1, F1};
    G.nTube = 1;
}

// ---------------------------------------------------------------------------
//  ExtendToHorizon : While (S.T < H) S.Extend(H, eps0)
// ---------------------------------------------------------------------------
inline void ExtendToHorizon(IMap f, Scaffold& S, SolutionChain& chain,
                            double eps0, double delta, int degree, double H,
                            int stepBtype, int tubedegree, int dim, int debuglevel)
{
    int guard = 0;
    while (S.t.back() < H && guard++ < MAX_STAGES)
        Extend(f, S, chain, eps0, delta, degree, H, stepBtype, tubedegree,
               dim, debuglevel);
    if (S.t.back() < H)
        throw std::runtime_error("Extend: stage limit reached before the horizon");
}

// ---------------------------------------------------------------------------
//  RebuildScaffold : throw the stage grid away and lay it out again
// ---------------------------------------------------------------------------
//  Split shrinks E_0, and a smaller box admits longer steps.  The grid is
//  chosen afresh from the new E_0; the scaffold keeps the final time it had.
inline void RebuildScaffold(IMap f, Scaffold& S, SolutionChain& chain,
                            double eps0, double delta, int degree,
                            int stepBtype, int tubedegree, int dim, int debuglevel)
{
    const double Hcur = S.t.back();
    MiniScaffold G0   = S.G[0];
    G0.E.assign(1, S.E[0]);
    G0.F.assign(1, S.E[0]);
    G0.mu.assign(1, 1000.0);              // sentinel: no bound known yet
    G0.ell = 0;

    S.t = std::vector<double>(1, 0.0);
    S.E = std::vector<IVector>(1, S.E[0]);
    S.F = std::vector<IVector>(1, S.E[0]);
    S.G = std::vector<MiniScaffold>(1, G0);

    chain.restart(S.E[0]);
    ExtendToHorizon(f, S, chain, eps0, delta, degree, Hcur, stepBtype,
                    tubedegree, dim, debuglevel);
}
