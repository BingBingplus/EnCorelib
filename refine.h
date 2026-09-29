#pragma once
// ===========================================================================
//  refine.h  --  S.Refine(eps0)                    (paper, Section 5)
//
//  Refine turns an m-stage scaffold into an eps0-small one: it refines S.E
//  and S.F, keeping S.t fixed, until wmax(S.E[m]) <= eps0.  The subroutines
//  and the stopping criteria are those of Section 5:
//
//      Bisect(i)        Subsection 5.1  -- halve every mini-step of stage i
//      EulerTube(i)     Subsection 5.2  -- Theorem 5.1 tube enclosures
//      Split()          Subsection 5.3  -- shrink S.E[0]
//      StageFreeze(i)   equation (5.4)
//      NeedSplit()      equation (5.5)
//      CanTerminate()   equation (5.6)
//
//  and the main loop has the three parts A (refine every stage), B (split
//  check) and C (terminate check) of the pseudo-code.
//
//  The chain of admissible mini-quads
//
//      E_{i-1} -> Quad[i,1] -> ... -> Quad[i,2^l_i] -> E_i
//
//  is evaluated with a single Lohner doubleton carried from S.E[0] through
//  every mini-step of every stage, so the matrix term of StepB acts on the
//  original E_0 throughout.  SolutionChain in stepAB.h has the details.
//
//  `method` says which of the three operators Refine is allowed to use.  It
//  is independent of `stepB`: `stepB` chooses the one-step enclosure (always
//  Lohner's method, in one of four flavours), `method` chooses how the
//  scaffold is refined once those steps exist.
//
//    method = 0 "tube"       Bisect + EulerTube + Split.  Stage i is cut into
//                            2^l_i mini-steps and the tube of Theorem 5.1 is
//                            laid along them.  This is Refine as it is
//                            written out in Section 5, and the default.
//    method = 1 "bisect"     Bisect + Split; the tube is switched off.
//    method = 2 "splitonly"  Split only.  The stages stay the steps Lohner's
//                            method took and are not subdivided any further,
//                            so the accuracy comes from shrinking E_0.
// ===========================================================================

#include <stdexcept>
#include "extend.h"   // Extend, ExtendToHorizon, RebuildScaffold, stepAB.h

// Refinement strategies (`method` on the command line).
//  Which of the three refinement operators of Section 5 Refine may use.
//  Orthogonal to `stepB`, which chooses the one-step enclosure.
enum RefineMethod {
    REFINE_TUBE      = 0,   // Bisect + EulerTube + Split   (Refine of Section 5)
    REFINE_BISECT    = 1,   // Bisect + Split
    REFINE_SPLITONLY = 2    // Split only; the stages are the steps StepB takes
};

// A stage holds at most 2^ELL_MAX mini-steps; l_i controls how fine a grid
// EulerTube works on, and the bound of Lemma 6.1 is reached well inside this.
static constexpr int ELL_MAX = 6;

// ═══════════════════════════════════════════════════════════════════════════
//  Evaluate : walk the scaffold again, on the grid it already has
// ═══════════════════════════════════════════════════════════════════════════
//  Refine changes two things about a scaffold: Split shrinks E_0, and Bisect
//  asks for stage i to be crossed in 2^{l_i} mini-steps instead of one.
//  Evaluate is what makes either of them take effect.  It carries one
//  doubleton from E_0 along the whole grid and fills, for every stage i and
//  every mini-step j,
//
//      S.G[i].E[j]   end-enclosure of the j-th mini-quad
//      S.G[i].F[j]   full enclosure of the j-th mini-quad   (equation (2.2))
//      S.G[i].mu[j]  a logNorm bound on that full enclosure       (item [G4])
//      S.E[i], S.F[i]
//
//  Every mini-step carries its own StepA certificate: over half the step the
//  full enclosure is narrower, so mu is smaller, the criteria of Subsection
//  5.3 are met sooner, and the tube of Theorem 5.1 has a finer grid to lie on.
//
//  The time grid itself is left alone.  Whatever the scaffold already holds
//  encloses the same solution sets for a box that was no smaller, so the new
//  values are intersected into the old ones.
inline void Evaluate(Scaffold& S, IMap f, SolutionChain& chain, int degree,
                     int stepBtype, int dim, double stepAdelta)
{
    chain.restart(S.E[0]);
    const int m = S.m();

    S.G[0].E[0] = S.E[0];
    S.G[0].F[0] = S.E[0];
    S.F[0]      = S.E[0];

    for (int i = 1; i <= m; ++i) {
        MiniScaffold& G  = S.G[static_cast<size_t>(i)];
        const int     N  = 1 << G.ell;
        const double  t0 = S.t[static_cast<size_t>(i - 1)];
        const double  t1 = S.t[static_cast<size_t>(i)];
        const double  h1 = (t1 - t0) / static_cast<double>(N);

        std::vector<IVector> E(static_cast<size_t>(N) + 1, IVector(dim));
        std::vector<IVector> F(static_cast<size_t>(N) + 1, IVector(dim));
        std::vector<double>  mu(static_cast<size_t>(N) + 1, 0.0);

        E[0]  = (i == 1) ? S.E[0] : S.G[static_cast<size_t>(i - 1)].E.back();
        F[0]  = E[0];
        mu[0] = 1000.0;                      // sentinel: no bound known yet

        IVector stageF;
        for (int j = 1; j <= N; ++j) {
            const double tj = (j == N) ? t1 : t0 + static_cast<double>(j) * h1;

            // equation (2.2) for this mini-step.  The stage's own full
            // enclosure bounds it as well, having been certified over the
            // whole stage from a box no smaller than this one.
            IVector Fj;
            try {
                Fj = StepAatStep(f, E[static_cast<size_t>(j - 1)], h1, stepAdelta);
            } catch (const std::exception&) {
                Fj = S.F[static_cast<size_t>(i)];
            }

            if (!chain.marchTo(tj))
                throw std::runtime_error("Refine: StepB cannot advance");

            const double muj = logNorm(computeJacobian(f, Fj), dim);
            IVector      Ej  = chain.hull();
            if (stepBtype == 1 || stepBtype == 2)
                Ej = logNormBallTighten(f, degree, h1,
                                        E[static_cast<size_t>(j - 1)], muj, dim, Ej);

            E[static_cast<size_t>(j)]  = IntersectB(Ej, Fj);
            F[static_cast<size_t>(j)]  = Fj;
            mu[static_cast<size_t>(j)] = muj;
            stageF = (j == 1) ? Fj : hullBox(stageF, Fj);
        }

        // Both the old and the new values enclose End(E_0, t_i), so keep the
        // better of the two.
        E[static_cast<size_t>(N)] = IntersectB(E[static_cast<size_t>(N)],
                                               S.E[static_cast<size_t>(i)]);
        stageF = IntersectB(stageF, S.F[static_cast<size_t>(i)]);

        G.E  = E;
        G.F  = F;
        G.mu = mu;
        S.E[static_cast<size_t>(i)] = E.back();
        S.F[static_cast<size_t>(i)] = stageF;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Bisect(i)                                       (paper, Subsection 5.1)
// ═══════════════════════════════════════════════════════════════════════════
//  Subdivides every mini-step of stage i into two, so that stage i goes from
//  2^l_i to 2^{l_i+1} uniform admissible mini-quads.
//
//  For the Lohner propagation this is a step-size halving with the
//  intermediate point inserted: the stage keeps the span [t_{i-1}, t_i] it
//  was given, but it is crossed in twice as many, half as long, steps, and
//  each of them gets its own StepA certificate, so the full enclosures F^i
//  shrink with the step.  The finer chain of mini-quads is produced by the
//  next Evaluate(), which walks the whole scaffold with one doubleton.
//
//  `levels` bisections are applied at once, so that a grid Lemma 2.2 needs
//  several levels finer is reached in a single pass.
//
//  Returns false when the grid has reached ELL_MAX and cannot be refined.
inline bool Bisect(Scaffold& S, int idx, int levels = 1) {
    MiniScaffold& G = S.G[static_cast<size_t>(idx)];
    if (G.ell >= ELL_MAX) return false;
    G.ell = std::min(ELL_MAX, G.ell + std::max(1, levels));
    return true;
}

// ═══════════════════════════════════════════════════════════════════════════
//  EulerTube(i)                                    (paper, Subsection 5.2)
// ═══════════════════════════════════════════════════════════════════════════
//  Applies Theorem 5.1 to every mini-step of stage idx.  With E^i[0] inside
//  Ball_{q0}(r0) and mu a logNorm bound on the stage:
//
//    (a) Ball_{q_j}(r0 e^{j mu h1} + delta)               is an end-enclosure
//    (b) Box(Ball_{q_{j-1}}(r), Ball_{q_j}(r)),
//        r = delta + max(r_{j-1} e^{mu h1}, r_{j-1})      is a full enclosure
//
//  where q_0, q_1, ... is the Euler polygon of step size h1 = Dt_i / 2^l_i.
//
//  Theorem 5.1 applies only while the Euler point q_{j-1} lies in the current
//  end-enclosure E^i[j-1] -- condition (2.7), via Subsection 2.6.  When that
//  test fails the tube is restarted at the midpoint of E^i[j-1]; each restart
//  begins a new tube segment.  The number of segments is stored in G.nTube
//  and feeds N^tube(S) in the stopping criteria.
inline int EulerTube(Scaffold& S, int idx, int dim, IMap f, int debuglevel) {
    MiniScaffold& G     = S.G[static_cast<size_t>(idx)];
    const double  delta = G.delta;
    const int     N     = 1 << G.ell;
    const double  h1    = (S.t[static_cast<size_t>(idx)] - S.t[static_cast<size_t>(idx - 1)])
                          / static_cast<double>(N);
    if (static_cast<int>(G.E.size()) < N + 1) return G.nTube;

    // Theorem 5.1 needs one mu for the whole tube, so take the largest of the
    // mini-step bounds of this stage (entry 0 is the sentinel).
    double mu = 0.0;
    if (G.mu.size() > 1) {
        mu = G.mu[1];
        for (size_t j = 2; j < G.mu.size(); ++j) mu = std::max(mu, G.mu[j]);
    }

    auto ballBox = [&](const IVector& q, double r) {
        IVector b(dim);
        for (int k = 0; k < dim; ++k) b[k] = q[k] + interval(-r, r);
        return b;
    };
    auto contains = [&](const IVector& outer, const IVector& inner) {
        for (int k = 0; k < dim; ++k)
            if (inner[k].leftBound()  < outer[k].leftBound() ||
                inner[k].rightBound() > outer[k].rightBound()) return false;
        return true;
    };
    auto startSegment = [&](const IVector& box, IVector& q, double& r0) {
        const std::vector<double> c = midpoint(box);
        q = IVector(dim);
        for (int k = 0; k < dim; ++k)
            q[k] = interval(c[static_cast<size_t>(k)], c[static_cast<size_t>(k)]);
        r0 = euclideanRadius(box);
    };

    IVector q(dim);
    double  r0 = 0.0;
    int     step = 0, segments = 1;
    startSegment(G.E[0], q, r0);

    IVector stageF = G.F[1];

    for (int j = 1; j <= N; ++j) {
        // condition (2.7): the Euler point must lie in the current enclosure
        if (!contains(G.E[static_cast<size_t>(j - 1)], q)) {
            startSegment(G.E[static_cast<size_t>(j - 1)], q, r0);
            step = 0;
            ++segments;
            if (debuglevel == 1)
                std::cout << "[EulerTube] restart at stage=" << idx << " j=" << j << "\n";
        }

        const double rPrev = r0 * std::exp(static_cast<double>(step)     * mu * h1);
        const double rNext = r0 * std::exp(static_cast<double>(step + 1) * mu * h1);

        IVector qPrev = q;
        q = q + f(q) * interval(h1);        // Euler step of the centre trajectory

        // (b) full enclosure of mini-step j
        const double rFull = delta + std::max(rNext, rPrev);
        IVector tubeF = hullBox(ballBox(qPrev, rFull), ballBox(q, rFull));
        G.F[static_cast<size_t>(j)] = IntersectB(G.F[static_cast<size_t>(j)], tubeF);
        stageF = hullBox(stageF, G.F[static_cast<size_t>(j)]);
        G.mu[static_cast<size_t>(j)] = logNorm(computeJacobian(f, G.F[static_cast<size_t>(j)]), dim);

        // (a) end enclosure of mini-step j
        G.E[static_cast<size_t>(j)] =
            IntersectB(G.E[static_cast<size_t>(j)], ballBox(q, rNext + delta));

        ++step;
    }

    G.nTube = segments;
    S.E[static_cast<size_t>(idx)] = G.E.back();
    S.F[static_cast<size_t>(idx)] = stageF;
    return segments;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Split()                                         (paper, Subsection 5.3)
// ═══════════════════════════════════════════════════════════════════════════
//  Halves S.E[0] towards `target` (its midpoint when no target is given).
//  The complementary part is not stored here: the caller -- EndCover, or the
//  boundary sweep -- sees that S.E[0] has shrunk and re-queues the remaining
//  sub-boxes.  Those are the "split-off scaffolds returned by Refine".
inline void Split(Scaffold& S, int dim, const std::vector<double>& target) {
    const bool hasTarget = !target.empty() && static_cast<int>(target.size()) == dim;
    IVector E0(dim);
    for (int i = 0; i < dim; ++i) {
        const double lo = S.E[0][i].leftBound();
        const double hi = S.E[0][i].rightBound();
        if (hasTarget) {
            const double t = target[static_cast<size_t>(i)];
            E0[i] = interval(lo + (t - lo) * 0.5, hi - (hi - t) * 0.5);
        } else {
            E0[i] = interval(lo + (hi - lo) / 4.0, hi - (hi - lo) / 4.0);
        }
    }
    S.E[0]      = E0;
    S.F[0]      = E0;
    S.G[0].E[0] = E0;
    S.G[0].F[0] = E0;
}

// ═══════════════════════════════════════════════════════════════════════════
//  The three stopping criteria                    (paper, Subsection 5.3)
// ═══════════════════════════════════════════════════════════════════════════

// mu_max : the global logNorm bound of the scaffold              (item [G4]).
// Index 0 of each mu^i is the sentinel planted by Extend, and is skipped.
inline double muMaxOf(const Scaffold& S) {
    double muMax = -1e300;
    for (const auto& G : S.G)
        for (size_t j = 1; j < G.mu.size(); ++j)
            if (G.mu[j] < 100.0) muMax = std::max(muMax, G.mu[j]);
    return (muMax < -1e299) ? 0.0 : muMax;
}

// N^tube(S) : the total number of tube segments over all stages.
inline int nTubeOf(const Scaffold& S) {
    int n = 0;
    for (size_t i = 1; i < S.G.size(); ++i) n += std::max(1, S.G[i].nTube);
    return std::max(1, n);
}

// StageFreeze(i)                                   (paper, equation (5.4))
//
//      e^{mu_max T} delta_i N^tube  <  eps0 / 8
//
// Once this holds, refining stage i further cannot change the outcome, so the
// stage is skipped for the rest of Refine.
inline bool StageFreeze(const Scaffold& S, int i, double eps0) {
    const double T = S.t.back() - S.t.front();
    const double lhs = std::exp(muMaxOf(S) * T)
                     * S.G[static_cast<size_t>(i)].delta
                     * static_cast<double>(nTubeOf(S));
    return lhs < eps0 / 8.0;
}

// NeedSplit()                                      (paper, equation (5.5))
//
//      e^{mu_max T} wmax(S.E[0])  >  eps0 / (2 sqrt(n))
//
// If this fails, the initial box is already small enough that refining the
// stages alone can reach eps0, and S.E[0] must not be split further.
inline bool NeedSplit(const Scaffold& S, double eps0, int dim) {
    const double T = S.t.back() - S.t.front();
    const double lhs = std::exp(muMaxOf(S) * T) * wmax(S.E[0]);
    return lhs > eps0 / (2.0 * std::sqrt(static_cast<double>(dim)));
}

// CanTerminate()                                   (paper, equation (5.6))
//
//      ( AND_i FreezeFlag[i] )  AND  ( NeedSplit() = False )
inline bool CanTerminate(const std::vector<char>& freezeFlag, bool needSplit) {
    if (needSplit) return false;
    for (size_t i = 1; i < freezeFlag.size(); ++i)
        if (!freezeFlag[i]) return false;
    return true;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Refine
// ═══════════════════════════════════════════════════════════════════════════
//  S           scaffold, refined in place
//  dim         n
//  chain       the doubleton carried along the scaffold.  Refine leaves it
//              anchored at the current E_0 and standing at S.T, so that the
//              next Extend can carry on from there.
//  method      RefineMethod, see the enum above
//  eps0        Refine returns once wmax(S.E.back()) <= eps0
//  target      contraction point for Split()
inline void Refine(Scaffold& S,
                   int dim,
                   IMap f,
                   SolutionChain& chain,
                   int method,
                   int stepBtype,
                   int degree,
                   double eps0,
                   int tubedegree,
                   int debuglevel,
                   const std::vector<double>& target = {})
{
    const bool   useBisect = (method == REFINE_TUBE || method == REFINE_BISECT);
    const bool   useTube   = (method == REFINE_TUBE);
    const bool   splitOnly = !useBisect && !useTube;
    const double delta     = S.G[0].delta;

    if (wmax(S.E.back()) <= eps0) return;

    int m = S.m();
    std::vector<char> freezeFlag(static_cast<size_t>(m) + 1, 0);

    // Lay the tube of Theorem 5.1 along every stage whose grid is already
    // fine enough for Lemma 2.2 to apply.
    auto tubes = [&]() {
        if (!useTube) return;
        for (int i = 1; i <= m; ++i) {
            const MiniScaffold& G = S.G[static_cast<size_t>(i)];
            const double Dt = S.t[static_cast<size_t>(i)] - S.t[static_cast<size_t>(i - 1)];
            if (Dt / static_cast<double>(1 << G.ell) <= G.hTube)
                EulerTube(S, i, dim, f, debuglevel);
        }
    };

    // Split has shrunk E_0, so the scaffold is evaluated once more.  Bisect's
    // work lives in the time grid, so the grid is kept as soon as there is any
    // refinement in it; otherwise the smaller box gets a grid of its own.
    auto reevaluate = [&]() {
        bool refined = false;
        for (int i = 1; i <= m; ++i)
            if (S.G[static_cast<size_t>(i)].ell > 0) { refined = true; break; }
        if (refined) {
            Evaluate(S, f, chain, degree, stepBtype, dim, eps0);
        } else {
            RebuildScaffold(f, S, chain, eps0, delta, degree, stepBtype,
                            tubedegree, dim, debuglevel);
            m = S.m();
            freezeFlag.assign(static_cast<size_t>(m) + 1, 0);
        }
        tubes();
    };

    tubes();

    while (eps0 < wmax(S.E.back())) {
        const double before = wmax(S.E.back());

        // ---- Part A : refine every stage that is not frozen --------------
        bool bisected = false;
        if (useBisect) {
            for (int i = 1; i <= m; ++i) {
                if (freezeFlag[static_cast<size_t>(i)]) continue;
                MiniScaffold& G = S.G[static_cast<size_t>(i)];
                const double Dt = S.t[static_cast<size_t>(i)] - S.t[static_cast<size_t>(i - 1)];
                const double h  = Dt / static_cast<double>(1 << G.ell);

                // equation (5.4): once the stage's own tolerance is small
                // enough, refining it further cannot change the outcome.
                if (StageFreeze(S, i, eps0)) {
                    freezeFlag[static_cast<size_t>(i)] = 1;
                    continue;
                }

                // The grid is fine enough for Lemma 2.2, so the tube of
                // Theorem 5.1 is already lying along this stage.  Ask for a
                // narrower one next time.                        (item [G5])
                if (useTube && G.hTube > 0.0 && h <= G.hTube) {
                    G.delta /= 2.0;
                    double mu = 0.0;
                    for (size_t j = 1; j < G.mu.size(); ++j) mu = std::max(mu, G.mu[j]);
                    G.hTube = std::min(G.hTube,
                                       computeHTube(f, S.F[static_cast<size_t>(i)], Dt,
                                                    G.delta, mu, tubedegree, degree));
                }

                // Bisect(i): halve the step of this stage, in one move when
                // the tube asks for a grid several levels finer.      (5.1)
                int levels = 1;
                if (useTube && G.hTube > 0.0 && h > G.hTube)
                    levels = static_cast<int>(std::ceil(std::log2(h / G.hTube)));
                if (G.ell + levels > ELL_MAX || !Bisect(S, i, levels))
                    freezeFlag[static_cast<size_t>(i)] = 1;
                else
                    bisected = true;
            }
        }

        // The refinement takes effect when the scaffold is walked again, and
        // the split check below reads mu_max off the refined scaffold.
        if (bisected) Evaluate(S, f, chain, degree, stepBtype, dim, eps0);
        tubes();

        // ---- Part B : split check, equation (5.5) ------------------------
        // With Split as the only operator the criterion does not gate it.
        const bool needSplit = splitOnly ? true : NeedSplit(S, eps0, dim);
        if (needSplit) {
            if (wmax(S.E[0]) == 0.0) break;     // E_0 is a point
            Split(S, dim, target);
            reevaluate();
        }

        // ---- Part C : terminate check, equation (5.6) --------------------
        if (CanTerminate(freezeFlag, needSplit) && eps0 < wmax(S.E.back())) {
            if (wmax(S.E[0]) == 0.0) break;
            Split(S, dim, target);
            reevaluate();
        }

        // (5.4)-(5.6) are sufficient conditions; Split also runs when a pass
        // has neither split nor narrowed the scaffold.
        if (!needSplit && !(wmax(S.E.back()) < before)) {
            if (wmax(S.E[0]) == 0.0) break;
            Split(S, dim, target);
            reevaluate();
        }
    }
}
