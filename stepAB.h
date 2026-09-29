#pragma once
// ===========================================================================
//  stepAB.h  --  StepA, StepB, the tube thresholds, and the solution chain
//
//  This file implements the one-step machinery of the paper:
//
//      StepA                    Subsection 4.1, Lemma 4.1
//      admissibility test       equation (2.2)
//      StepB                    Subsection 4.2
//      Lemma 2.1 ball           Subsection 2.4
//      h_euler(H,M,mu,delta)    equation (2.6), Lemma 2.2
//      Taylor tube              the degree-p generalisation of the Euler tube
//      SolutionChain            how a stage chain is evaluated (see below)
//
//  Interval arithmetic and the normalised Taylor coefficients f^[i] are taken
//  from the underlying validated-arithmetic library; everything above is
//  built here.
// ===========================================================================

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include "types.h"

// Taylor order k used inside StepA, i.e. the order of the remainder that
// equation (2.2) certifies.  StepB runs at the user's order (the paper uses 20).
static constexpr int STEPA_ORDER = 5;

// ═══════════════════════════════════════════════════════════════════════════
//  Taylor polynomial of a step                     (paper, Subsection 2.2)
// ═══════════════════════════════════════════════════════════════════════════
//  P(E, h) = sum_{i<k} [0,h]^i f^[i](E),  evaluated by Horner's rule as
//  recommended in Subsection 4.1.
inline IVector taylorPolyOverStep(IMap f, const IVector& E, double h, int k) {
    const int n = E.dimension();
    ICnOdeSolver s(f, k);
    s.computeCoefficients(E, k);

    IVector P(n);
    for (int j = 0; j < n; ++j) {
        P[j] = s.coefficient(j, k - 1);
        for (int d = k - 2; d >= 1; --d)
            P[j] = P[j] * interval(0.0, h) + s.coefficient(j, d);
        P[j] = P[j] * interval(0.0, h) + E[j];
    }
    return P;
}

// ═══════════════════════════════════════════════════════════════════════════
//  StepA                                       (paper, Subsection 4.1)
// ═══════════════════════════════════════════════════════════════════════════
//
//  StepA(E0, H, delta) -> (h, F1) with (E0, h, F1) a delta-admissible triple
//  and h <= H.  Admissibility is the inclusion of equation (2.2):
//
//      sum_{i<k} [0,h]^i f^[i](E0) + [0,h]^k f^[k](F1)  subset  F1.
//
//  Two entry points, one algorithm:
//
//    StepA(f, E0, H, delta)          chooses h as large as Lemma 4.1 allows
//    StepAatStep(f, E0, h, delta)    h is already fixed -- by a scaffold
//                                    stage, or by StepA itself -- and only
//                                    the full enclosure is wanted
//
//  Both finish with the contraction of the last line of the StepA pseudo-code,
//
//      F1  <-  sum_{i<k} [0,h]^i f^[i](E0) + [0,h]^k f^[k](F1),
//
//  which replaces the inflated candidate by the (narrower) set the inclusion
//  actually certifies.

// StepA with the step size prescribed.  Returns the contracted F1.
inline IVector StepAatStep(IMap f, const IVector& E, double h,
                           double delta, int k = STEPA_ORDER) {
    const int n = E.dimension();
    const IVector P = taylorPolyOverStep(f, E, h, k);
    const interval hk(0.0, std::pow(h, k));

    // Candidate F1 = P + Box(inflation), inflated until (2.2) holds.
    double inflation = std::max(delta, 1e-14);
    const double inflationMax = std::max(1.0, 1e3 * (wmax(E) + 1.0));

    for (int attempt = 0; attempt < 40 && inflation <= inflationMax;
         ++attempt, inflation *= 2.0) {
        IVector F(n);
        for (int j = 0; j < n; ++j) F[j] = P[j] + interval(-inflation, inflation);

        IVector Fc(n);
        try {
            ICnOdeSolver s(f, k);
            s.computeCoefficients(F, k);
            for (int j = 0; j < n; ++j)
                Fc[j] = P[j] + hk * s.coefficient(j, k);   // the contraction
        } catch (const std::exception&) {
            continue;                                      // F too wide for f
        }

        bool admissible = true;
        for (int j = 0; j < n; ++j)
            if (!(Fc[j].leftBound()  > F[j].leftBound() &&
                  Fc[j].rightBound() < F[j].rightBound())) { admissible = false; break; }
        if (admissible) return Fc;
    }
    throw std::runtime_error("StepA: no admissible full enclosure at this step size");
}

// StepA proper: pick h, then certify and contract F1.
//
//   h  <-  0
//   While (H > h)
//       F1 <- sum_{i<k} [0,H]^i f^[i](E0) + Box(delta)
//       M  <- ( ||f^[k](F1)||_1, ..., ||f^[k](F1)||_n )
//       h  <- min{ H, min_i (delta_i / M_i)^{1/k} }        (Lemma 4.1)
//       H  <- H/2
//   F1 <- sum_{i<k} [0,h]^i f^[i](E0) + [0,h]^k f^[k](F1)
//
//  `delta` bounds the k-th order remainder componentwise; it is the delta of
//  Subsection 4.1, not the tube tolerance delta_i of item [G5].
inline std::pair<double, IVector> StepA(IMap f, IVector E0, double H,
                                        double delta, int k = STEPA_ORDER) {
    const int n = E0.dimension();
    ICnOdeSolver solver(f, k);
    solver.computeCoefficients(E0, k);

    double trialH = H;
    double h      = 0.0;
    IVector F1(n);

    while (h < trialH / 2.0) {
        IVector P(n);
        for (int j = 0; j < n; ++j) {
            P[j] = solver.coefficient(j, k - 1);
            for (int d = k - 2; d >= 1; --d)
                P[j] = P[j] * interval(0.0, trialH) + solver.coefficient(j, d);
            P[j] = P[j] * interval(0.0, trialH) + E0[j];
        }
        for (int j = 0; j < n; ++j) F1[j] = P[j] + interval(-delta, delta);

        ICnOdeSolver probe(f, k);
        probe.computeCoefficients(F1, k);
        IVector dk(n);
        for (int j = 0; j < n; ++j) dk[j] = probe.coefficient(j, k);

        capd::vectalg::EuclNorm<IVector, IMatrix> euclNorm;
        const double M = sup(euclNorm(dk));
        h = std::pow(delta / M, 1.0 / k);
        if (h >= trialH) { h = trialH; break; }
        trialH /= 2.0;
    }

    // Last line of the pseudo-code: replace the inflated F1 by what (2.2)
    // certifies at the chosen h.  The loop's F1 is delta-admissible already,
    // by the loop invariant.
    try {
        F1 = StepAatStep(f, E0, h, delta, k);
    } catch (const std::exception&) {
    }
    return {h, F1};
}

// ═══════════════════════════════════════════════════════════════════════════
//  StepB                                        (paper, Subsection 4.2)
// ═══════════════════════════════════════════════════════════════════════════
//
//  StepB(E0, h, F1) -> E1 upgrades the admissible triple (E0, h, F1) to an
//  admissible quad (E0, h, F1, E1) by computing an end-enclosure of
//  End(E0, h).  StepB is Lohner's method: the centred term
//
//      ( sum_{i<k} h^i J_{f^[i]}(E0) ) . (E0 - m(E0))
//
//  is a *matrix acting on E0*, and it is kept as such -- a doubleton
//
//      E1 = q1 + C1 . r0 + B1 . r1
//
//  with a QR-preconditioned frame -- instead of being re-enclosed in an
//  axis-aligned box.  Re-enclosing it once per step is precisely what
//  produces the wrapping effect; see SolutionChain at the end of this file.
//
//  Four variants:
//
//    stepB = 0   Lohner's method
//    stepB = 1   Lohner's method intersected with the ball of Lemma 2.1
//    stepB = 2   Lohner's method with a higher-order a-priori enclosure,
//                intersected with the ball of Lemma 2.1
//    stepB = 3   as stepB = 2 without the ball
//
//  Variants 2 and 3 carry a tighter a-priori enclosure and are sharper on
//  mild problems; variants 0 and 1 hold up better when the flow expands
//  strongly.

// Lohner's method: one step of size h from the box B0.
inline IVector stepBlohner(IMap f, int degree, double h, IVector B0) {
    ICnOdeSolver solver(f, degree);
    ICnTimeMap timeMap(solver);
    CnRect2Set s(B0, 3.0);
    return timeMap(interval(h), s);
}

// Lohner's method with a higher-order a-priori enclosure.
inline IVector stepBlohnerHO(IMap f, int degree, double h, IVector B0) {
    IOdeSolver solver(f, degree);
    ITimeMap timeMap(solver);
    C0HORect2Set s(B0);
    return timeMap(interval(h), s);
}

// ── Lemma 2.1 ──────────────────────────────────────────────────────────────
//  If x1 is the trajectory through m(E0) and mu is a logNorm bound on a full
//  enclosure of (E0, h), then
//
//      End(E0, h)  subset  Ball_{x1(h)}( r0 e^{mu h} ),
//      r0 = (1/2) wmax(Ball(E0)),
//
//  so the box around that ball can only tighten an end-enclosure.
//
//  r0 is the radius of the smallest *Euclidean* ball containing E0: Lemma 2.1
//  is stated in the Euclidean norm, and the max half-width would understate
//  the spread.  The bound is skipped whenever e^{mu h} carries no information.
inline IVector logNormBallTighten(IMap f, int order, double h, const IVector& E0,
                                  double mu, int dim, const IVector& E) {
    const double r = euclideanRadius(E0) * std::exp(mu * h);
    if (!std::isfinite(r) || r <= 0.0 || r >= 1e12) return E;

    IVector centre(dim);
    const std::vector<double> c = midpoint(E0);
    for (int j = 0; j < dim; ++j)
        centre[j] = interval(c[static_cast<size_t>(j)], c[static_cast<size_t>(j)]);

    IVector q;
    try {
        q = stepBlohner(f, order, h, centre);            // thin: x1(h)
    } catch (const std::exception&) {
        return E;
    }

    IVector ball(dim);
    for (int j = 0; j < dim; ++j) {
        ball[j] = q[j] + interval(-r, r);
        if (!std::isfinite(ball[j].leftBound()) || !std::isfinite(ball[j].rightBound()))
            return E;
    }
    return IntersectB(E, ball);
}

inline IVector stepBlohnerWithMu(IMap f, int degree, double h, IVector B0, double mu) {
    IVector E = stepBlohner(f, degree, h, B0);
    return logNormBallTighten(f, degree, h, B0, mu, B0.dimension(), E);
}

inline IVector stepBlohnerHOWithMu(IMap f, int degree, double h, IVector B0, double mu) {
    IVector E = stepBlohnerHO(f, degree, h, B0);
    return logNormBallTighten(f, degree, h, B0, mu, B0.dimension(), E);
}

// StepB dispatcher.
inline IVector StepB(IMap f, int degree, double h, const IVector& B0,
                     double mu, int stepBtype) {
    switch (stepBtype) {
    case 0: return stepBlohner(f, degree, h, B0);
    case 1: return stepBlohnerWithMu(f, degree, h, B0, mu);
    case 2: return stepBlohnerHOWithMu(f, degree, h, B0, mu);
    case 3: return stepBlohnerHO(f, degree, h, B0);
    default: return stepBlohner(f, degree, h, B0);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tube thresholds                            (paper, equation (2.6))
// ═══════════════════════════════════════════════════════════════════════════
//
//                       / min{H, 2 mu delta / (M (e^{mu H} - 1) - mu^2 delta)}  mu > 0
//  h_euler(H,M,mu,delta) = | min{H, delta / (M H)}                              mu = 0
//                       \ min{H, 2 mu delta / (M (e^{mu H} - 1))}               mu < 0
//
//  H     : stage length
//  M     : bound on || f^[2] || over the full enclosure
//  mu    : logNorm bound on the full enclosure
//  delta : target tube radius at the end of the stage
//
//  If 0 < h <= h_euler then the Euler polygon stays inside the delta-tube of
//  the exact trajectory (Lemma 2.2).  A return value of 0 means "no usable
//  tube at this stage", which is what an overflowing e^{mu H} amounts to.
inline double hEuler(double H, double M, double mu, double delta) {
    if (M == 0.0 || delta == 0.0) return H;
    const double e = std::exp(mu * H);
    if (!std::isfinite(e)) return 0.0;
    const double h = (mu > 0.0)
        ? (2.0 * mu * delta) / (M * (e - 1.0) - mu * mu * delta)
        : (2.0 * mu * delta) / (M * (e - 1.0));
    if (!std::isfinite(h) || h <= 0.0) return 0.0;
    return std::min(H, h);
}

// || f^[p+1] || over the box B.
inline double taylorCoeffNorm(IMap f, const IVector& B, int p, int degreeAvail) {
    const int n = B.dimension();
    const int solverOrder = std::max(p + 1, degreeAvail);
    ICnOdeSolver solver(f, solverOrder);
    solver.computeCoefficients(B, solverOrder);

    IVector coeffs(n);
    for (int i = 0; i < n; ++i)
        coeffs[i] = solver.coefficient(i, p + 1);

    capd::vectalg::EuclNorm<IVector, IMatrix> euclNorm;
    return std::max(0.0, euclNorm(coeffs).rightBound());
}

// Degree-p tube threshold.  With Cbar = || f^[p+1] ||:
//
//   h = (mu delta / (Cbar (e^{mu H} - 1)))^{1/p}                  (mu > 0)
//   h = (delta / (Cbar H))^{1/p}                                  (mu = 0)
//   h = (2 mu delta / (Cbar (e^{mu H} - 1) - mu^2 delta))^{1/p}   (mu < 0)
//
// For p = 1 this is h_euler above.
inline double hTaylorStep(double H, double Cbar, double mu, double delta, int p) {
    if (H <= 0.0 || delta <= 0.0 || p <= 0) return 0.0;
    if (Cbar <= 0.0) return H;

    double candidate = H;
    const double tiny = 1e-12;

    if (std::abs(mu) < tiny) {
        const double ratio = delta / (Cbar * H);
        if (ratio > 0.0) candidate = std::pow(ratio, 1.0 / static_cast<double>(p));
    } else if (mu >= 0.0) {
        const double den = Cbar * (std::exp(mu * H) - 1.0);
        const double ratio = (den > 0.0) ? (mu * delta / den) : -1.0;
        if (ratio > 0.0) candidate = std::pow(ratio, 1.0 / static_cast<double>(p));
    } else {
        const double den = Cbar * (std::exp(mu * H) - 1.0) - (mu * mu) * delta;
        const double ratio = (std::abs(den) > tiny) ? (2.0 * mu * delta / den) : -1.0;
        if (ratio > 0.0) candidate = std::pow(ratio, 1.0 / static_cast<double>(p));
    }

    if (!std::isfinite(candidate) || candidate <= 0.0) return H;
    return std::min(H, candidate);
}

//  tubedegree == 1 : Euler tube      (equation (2.6))
//  tubedegree >= 2 : Taylor tube of that degree
inline double computeHTube(IMap f, const IVector& F_enclosure,
                            double H, double delta, double mu,
                            int tubedegree, int degreeAvail) {
    if (tubedegree == 1) {
        // M = || J_f(F) . f(F) || / 2  bounds || f^[2] || on F.
        IMatrix J      = computeJacobian(f, F_enclosure);
        IVector fval   = f(F_enclosure);
        IVector coeffs = J * fval;
        capd::vectalg::EuclNorm<IVector, IMatrix> euclNorm;
        const double M = sup(euclNorm(coeffs)) / 2.0;
        if (M == 0.0 || std::abs(mu) < 1e-14) return H;
        return hEuler(H, M, mu, delta);
    }
    const double Cbar = taylorCoeffNorm(f, F_enclosure, tubedegree, degreeAvail);
    return hTaylorStep(H, Cbar, mu, delta, tubedegree);
}

// ═══════════════════════════════════════════════════════════════════════════
//  SolutionChain : how a chain of mini-quads is evaluated
// ═══════════════════════════════════════════════════════════════════════════
//
//  A stage of the scaffold is a chain of admissible mini-quads
//
//      E_{i-1} -> Quad[i,1] -> ... -> Quad[i,2^l_i] -> E_i        (eq. (5.1))
//
//  The centred term of StepB is a parallelepiped, and storing it as an
//  axis-aligned box at every mini-step is what produces the wrapping effect.
//  SolutionChain therefore carries a single doubleton
//
//      E_j = q_j + C_j . r0 + B_j . r_j
//
//  anchored at E_0 through the whole chain, so that E_j is the accumulated
//  affine map A_j applied to the original E_0.  A hull is taken only when a
//  value of E_j is reported, and a reported hull is never fed back into the
//  propagation.
//
//  `stepBtype` selects the same representation as the StepB variant above.
class SolutionChain {
public:
    SolutionChain(const IMap& f, int order, const IVector& E0, int stepBtype)
        : m_f(f), m_order(order), m_plain(stepBtype <= 1)
    {
        // The propagator refers to the map, so the chain owns m_f (declared
        // first) and is not copyable.
        restart(E0);
    }

    SolutionChain(const SolutionChain&)            = delete;
    SolutionChain& operator=(const SolutionChain&) = delete;

    // Start the chain again at E_0, at time 0.  Refine's Split shrinks E_0,
    // and the scaffold is then evaluated once more from the smaller box.
    void restart(const IVector& E0) {
        m_hull       = E0;
        m_t          = 0.0;
        m_widthLimit = 1e4 * (wmax(E0) + 1.0);
        if (m_plain) {
            m_plainSolver = std::make_unique<ICnOdeSolver>(m_f, m_order);
            m_plainSet    = std::make_unique<CnRect2Set>(E0, 3.0);
        } else {
            m_hoSolver = std::make_unique<IOdeSolver>(m_f, m_order);
            m_hoSet    = std::make_unique<C0HORect2Set>(E0);
        }
    }

    // ── one step at a time ──────────────────────────────────────────────
    //
    //  The scaffold is laid out one step at a time, every step becoming a
    //  stage or a mini-step of one.  Three primitives are enough:
    //
    //      proposedStep()  the step StepB would like to take next
    //      step(hmax)      take it when it is no longer than hmax, and a step
    //                      of exactly hmax when it is not
    //      marchTo(t)      reach t exactly, halving any step that fails
    //
    //  step() returns false and leaves the doubleton exactly as it was when
    //  the step cannot be taken, so the caller halves it and tries again: a
    //  halved step is a shorter stage, or one more mini-step.  The doubleton
    //  itself is never re-enclosed in a box.

    double proposedStep() const {
        return m_plain ? m_plainSolver->getStep().rightBound()
                       : m_hoSolver->getStep().rightBound();
    }

    bool step(double hmax) {
        if (!(hmax > 0.0)) return false;
        const double hp = proposedStep();
        return doStep((hp > 0.0 && hp <= hmax) ? 0.0 : hmax);
    }

    bool marchTo(double t) {
        const double hMin = 1e-13 * std::max(1.0, std::fabs(t));
        int guard = 0;
        while (m_t < t && guard++ < 200000) {
            double h = t - m_t;
            while (!step(h)) {
                h *= 0.5;
                if (!(h > hMin)) return false;
            }
        }
        return m_t >= t;
    }

    const IVector& hull() const { return m_hull; }
    double time() const { return m_t; }

private:
    //  One step: h = 0 means "the step StepB wants", h > 0 means exactly h.
    //  The doubleton is put back as it was when the step does not go through.
    bool doStep(double h) {
        if (m_plain) {
            CnRect2Set snap(*m_plainSet);
            bool ok = true;
            try {
                if (h > 0.0) { m_plainSolver->turnOffStepControl(); m_plainSolver->setStep(h); }
                m_plainSet->move(*m_plainSolver);
            } catch (const std::exception&) { ok = false; }
            if (h > 0.0) m_plainSolver->turnOnStepControl();
            if (ok) {
                const IVector next(*m_plainSet);
                if (usable(next)) {
                    m_hull = next;
                    m_t    = m_plainSet->getCurrentTime().rightBound();
                    return true;
                }
            }
            *m_plainSet = snap;
            return false;
        }
        C0HORect2Set snap(*m_hoSet);
        bool ok = true;
        try {
            if (h > 0.0) { m_hoSolver->turnOffStepControl(); m_hoSolver->setStep(h); }
            m_hoSet->move(*m_hoSolver);
        } catch (const std::exception&) { ok = false; }
        if (h > 0.0) m_hoSolver->turnOnStepControl();
        if (ok) {
            const IVector next(*m_hoSet);
            if (usable(next)) {
                m_hull = next;
                m_t    = m_hoSet->getCurrentTime().rightBound();
                return true;
            }
        }
        *m_hoSet = snap;
        return false;
    }

    // The step is accepted only when it returns a finite box; otherwise the
    // caller halves it.
    bool usable(const IVector& B) const {
        for (int k = 0; k < B.dimension(); ++k)
            if (!std::isfinite(B[k].leftBound()) || !std::isfinite(B[k].rightBound()))
                return false;
        return wmax(B) <= m_widthLimit;
    }

    IMap m_f;                       // owned; the propagator refers to it
    int  m_order;
    bool m_plain;
    std::unique_ptr<ICnOdeSolver> m_plainSolver;
    std::unique_ptr<CnRect2Set>   m_plainSet;
    std::unique_ptr<IOdeSolver>   m_hoSolver;
    std::unique_ptr<C0HORect2Set> m_hoSet;
    IVector m_hull;
    double  m_t;                    // time of m_hull
    double  m_widthLimit;
};
