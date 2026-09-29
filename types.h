#pragma once
// ===========================================================================
//  types.h  --  EndCover for IVP (reference implementation)
//
//  Basic interval/box utilities and the two data structures of the paper:
//
//      B. Zhang and C. Yap,
//      "End Cover for Initial Value Problem: Complete Validated Algorithm
//       with Complexity Analysis".
//
//  NOTATION (paper  <->  code)
//  ---------------------------------------------------------------------
//      B, B0            IVector                  axis-aligned box
//      w(B), wmax(B)    width(), wmax()          width / max-width  (Sec. 2.1)
//      m(B)             midpoint()               midpoint           (Sec. 2.1)
//      Box(S)           hullBox()                smallest enclosing box
//      mu(A) = mu_2(A)  logNorm()                logarithmic norm   (Sec. 2.4)
//      f^[i]            normalised Taylor coefficients (Sec. 2.2)
//
//      S = (t, E, F, G) struct Scaffold          scaffold           (Sec. 3.2)
//        S.t                Scaffold::t          t0 < t1 < ... < tm
//        S.E                Scaffold::E          end-enclosures  E0..Em
//        S.F                Scaffold::F          full-enclosures F1..Fm
//        S.G                Scaffold::G          refinement data G1..Gm
//        S.m                Scaffold::m()        number of stages
//        S.T                Scaffold::T()        current final time = t.back()
//
//      G[i] = ((l_i, E^i, F^i), mu^i, (delta_i, h^i_euler))
//                         struct MiniScaffold    mini-scaffold      (Sec. 5.1)
//        l_i                MiniScaffold::ell    refinement level
//        E^i                MiniScaffold::E      2^l_i + 1 mini end-enclosures
//        F^i                MiniScaffold::F      2^l_i mini full-enclosures
//        mu^i               MiniScaffold::mu     per-mini-step logNorm bounds
//        delta_i            MiniScaffold::delta  current tube tolerance
//        h^i_euler          MiniScaffold::hTube  tube step-size threshold
//
//      C                  struct CoverResult     the eps-cover      (Sec. 1)
// ===========================================================================

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

// ═══════════════════════════════════════════════════════════════════════════
//  Box utilities   (paper, Subsection 2.1)
// ═══════════════════════════════════════════════════════════════════════════

// Box(B1 u B2): smallest axis-aligned box containing B1 and B2.
inline IVector hullBox(const IVector& B1, const IVector& B2) {
    int n = B1.dimension();
    IVector B(n);
    for (int i = 0; i < n; ++i)
        B[i] = interval(std::min(B1[i].leftBound(), B2[i].leftBound()),
                        std::max(B1[i].rightBound(), B2[i].rightBound()));
    return B;
}

// B1 n B2, coordinatewise.  If the intersection is empty in coordinate i we
// keep B1[i]; this only happens through outward rounding of two enclosures of
// the same solution set, so keeping the (valid) enclosure B1[i] is sound.
inline IVector IntersectB(const IVector& B1, const IVector& B2) {
    int n = B1.dimension();
    IVector B(n);
    for (int i = 0; i < n; ++i) {
        double lo = std::max(B1[i].leftBound(), B2[i].leftBound());
        double hi = std::min(B1[i].rightBound(), B2[i].rightBound());
        B[i] = (lo <= hi) ? interval(lo, hi) : B1[i];
    }
    return B;
}

// wmax(B) = max_i w(B)_i          (paper, Subsection 2.1)
inline double wmax(const IVector& B) {
    double m = 0.0;
    for (int i = 0; i < B.dimension(); ++i)
        m = std::max(m, B[i].rightBound() - B[i].leftBound());
    return m;
}

// wmin(B) = min_i w(B)_i
inline double wmin(const IVector& B) {
    if (B.dimension() == 0) return 0.0;
    double m = B[0].rightBound() - B[0].leftBound();
    for (int i = 1; i < B.dimension(); ++i)
        m = std::min(m, B[i].rightBound() - B[i].leftBound());
    return m;
}

// m(B): midpoint of a box, as a plain vector of doubles.
inline std::vector<double> midpoint(const IVector& B) {
    std::vector<double> c(static_cast<size_t>(B.dimension()));
    for (int i = 0; i < B.dimension(); ++i)
        c[static_cast<size_t>(i)] = (B[i].leftBound() + B[i].rightBound()) / 2.0;
    return c;
}

// Radius of the smallest ball Ball(B) containing the box B:
//   |w(B)|_2 / 2.
inline double euclideanRadius(const IVector& B) {
    double s = 0.0;
    for (int i = 0; i < B.dimension(); ++i) {
        double w = B[i].rightBound() - B[i].leftBound();
        s += w * w;
    }
    return std::sqrt(s) / 2.0;
}

// J_f(B): the Jacobian of f evaluated on a box.
inline IMatrix computeJacobian(IMap f, IVector B) {
    const int n = B.dimension();
    IJet jet(n, n, 2);
    f(B, jet);
    IMatrix J(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            J[i][j] = jet(i, j);
    return J;
}

// Upper bound on the logarithmic norm mu_2(J)        (paper, Subsection 2.4)
//
//   mu_2(J) = largest eigenvalue of the symmetric part (J + J^T)/2.
//
// For n == 2 the eigenvalue is available in closed form; for n > 2 we fall
// back on the Euclidean logarithmic norm.
inline double logNorm(const IMatrix& J, int n) {
    if (n == 2) {
        interval a  = J[0][0];
        interval d  = J[1][1];
        interval b  = (J[0][1] + J[1][0]) / 2;
        interval ad = a - d;
        interval D  = ad * ad / 4 + b * b;
        double D_sup = D.rightBound();
        if (D_sup > 0.0)
            return (a + d).rightBound() / 2.0 + std::sqrt(D_sup);
        else
            return (a + d).rightBound() / 2.0;
    }
    capd::vectalg::EuclLNorm<IVector, IMatrix> lnorm;
    return lnorm(J).rightBound();
}

// Build the map string "var:x,y,...;fun:f1,f2,...;" the parser expects.
inline std::string Convert_to_IMap(const std::vector<std::string>& SVar,
                                    const std::vector<std::string>& SFun) {
    std::string s = "var:";
    for (size_t i = 0; i < SVar.size(); ++i) {
        s += SVar[i];
        if (i + 1 < SVar.size()) s += ",";
    }
    s += ";fun:";
    for (size_t i = 0; i < SFun.size(); ++i) {
        s += SFun[i];
        if (i + 1 < SFun.size()) s += ",";
    }
    s += ";";
    return s;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Data structures
// ═══════════════════════════════════════════════════════════════════════════

// Mini-scaffold  G[i] = ((l_i, E^i, F^i), mu^i, (delta_i, h^i_euler))
// (paper, Subsection 5.1, items [G1]-[G5]).
struct MiniScaffold {
    std::vector<double>  mu;      // mu^i : logNorm bound of each mini-step  [G4]
    double               delta;   // delta_i : current tube tolerance        [G5]
    double               hTube;   // h^i_euler : tube step-size threshold    [G5]
    int                  ell;     // l_i : stage i holds 2^l_i mini-steps    [G2]
    std::vector<IVector> E;       // E^i : 2^l_i + 1 mini end-enclosures     [G3]
    std::vector<IVector> F;       // F^i : 2^l_i mini full-enclosures        [G3]
    int                  nTube;   // tube segments made by the last EulerTube(i)
};

// Scaffold  S = (t, E, F, G)                      (paper, Subsection 3.2)
//
//   t = (t_0 < t_1 < ... < t_m),  E = (E_0,...,E_m),  F = (F_1,...,F_m)
//
// Stage i is the admissible quad (E_{i-1}, Dt_i, F_i, E_i) with
// Dt_i = t_i - t_{i-1}.  Each E_i contains End(E_0, t_i).
struct Scaffold {
    std::vector<double>       t;  // S.t
    std::vector<IVector>      E;  // S.E
    std::vector<IVector>      F;  // S.F
    std::vector<MiniScaffold> G;  // S.G

    int    m() const { return static_cast<int>(t.size()) - 1; }  // S.m
    double T() const { return t.back(); }                        // S.T
};

// Result of a cover computation.
//
//   E0   : the initial sub-boxes that were actually integrated
//   C    : the eps-cover of End(B0,H) -- one box per element of C
//   hull : Box(C), the "B1" column of Table 2 in the paper
struct CoverResult {
    std::vector<IVector> E0;
    std::vector<IVector> C;
    IVector              hull;
};

// Pretty-print helper for vectors of doubles.
inline std::ostream& operator<<(std::ostream& os, const std::vector<double>& v) {
    os << "{";
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i + 1 < v.size()) os << ",";
    }
    return os << "}";
}
