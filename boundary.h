#pragma once
// ===========================================================================
//  boundary.h  --  scaffold driver + the Boundary covering mode
//
//  runScaffold()  runs one scaffold (Extend / Refine) for a single initial
//                 box up to time H.  It is the common engine of both modes.
//
//  mode = 0 (EndCover)  see EnCorelib.cpp: EndCover subdivides all of B0.
//  mode = 1 (Boundary)  only the boundary of B0 is covered -- the four edges
//                       in 2D, the 2n faces in nD.  For a flow, the image of
//                       the boundary bounds the image of the interior, so the
//                       hull of the boundary cover is an enclosure of
//                       End(B0,H) whenever End(B0,H) is simply connected.
// ===========================================================================

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <string>
#include <algorithm>
#include <memory>

#include "refine.h"    // Bisect, EulerTube, Split, Refine, Extend
#include "endcover.h"  // EndCover (used for the nD faces)

// ═══════════════════════════════════════════════════════════════════════════
//  runScaffold : Init(B) -> Extend/Refine until S.T == H
// ═══════════════════════════════════════════════════════════════════════════
//  This is the inner loop of EndCover (Subsection 3.3):
//
//      While (S.t.back() < H)  S.Extend(H, eps0)
//      S.Refine(eps0)
//
//  Refine reads the whole scaffold and refines all of its stages, so it is
//  applied once the scaffold has reached the horizon.
//
//  `target` is the contraction point used by Split() (the midpoint of the
//  sub-box in EndCover mode, a boundary point in Boundary mode).
inline Scaffold runScaffold(
    const IVector& B,
    IMap F,
    const MiniScaffold& templateG,
    double eps0, double delta, int degree,
    double H,
    int method, int stepBtype, int tubedegree, int debuglevel, int dim,
    const std::vector<double>& target)
{
    MiniScaffold G = templateG;
    if (!G.E.empty()) G.E[0] = B;
    if (!G.F.empty()) G.F[0] = B;

    Scaffold S = {{0.0}, {B}, {B}, {G}};

    // The doubleton carried along the scaffold.  Extend advances it one step
    // per stage and Refine walks it again whenever Split or Bisect changes
    // what the scaffold has to enclose, so from E_0 to E_m the propagation
    // goes through a single doubleton and is never re-enclosed in a box.
    SolutionChain chain(F, degree, S.E[0], stepBtype);

    //  While (S.T < H)  S.Extend(H, eps0)
    ExtendToHorizon(F, S, chain, eps0, delta, degree, H, stepBtype,
                    tubedegree, dim, debuglevel);

    for (int d = 0; d < S.E.back().dimension(); ++d) {
        if (!std::isfinite(S.E.back()[d].leftBound()) ||
            !std::isfinite(S.E.back()[d].rightBound()))
            throw std::runtime_error("runScaffold: non-finite enclosure");
    }

    //  S.Refine(eps0)
    Refine(S, dim, F, chain, method, stepBtype, degree, eps0,
           tubedegree, debuglevel, target);
    return S;
}

// ═══════════════════════════════════════════════════════════════════════════
//  TwoDimEncAlgo  --  Boundary mode, n = 2 (the four edges of B0)
// ═══════════════════════════════════════════════════════════════════════════
inline void TwoDimEncAlgo(
    CoverResult& Sdata,
    IMap F, Scaffold& S,
    double eps0, double delta, int degree,
    double H,
    int method, int stepBtype, int tubedegree, int debuglevel, int dim,
    bool verbose = true)
{
    IVector Sbound(dim);
    bool hasBound = false;
    const int    maxDepth = 60;
    const double minWidth = 0.0;

    using MkTgt = std::function<std::vector<double>(const IVector&)>;
    std::function<void(int, interval, int, double, double, MkTgt, int)> processEdge;

    processEdge = [&](int fixedDim, interval fixedVal, int varDim,
                      double a, double b, MkTgt makeTarget, int depth) {
        if (!(a < b)) return;
        double cur = a;
        while (cur < b) {
            IVector edgeBox = S.E[0];
            edgeBox[fixedDim] = fixedVal;
            edgeBox[varDim]   = interval(cur, b);
            std::vector<double> tgt = makeTarget(edgeBox);
            try {
                Scaffold st = runScaffold(edgeBox, F, S.G[0],
                                          eps0, delta, degree, H,
                                          method, stepBtype, tubedegree, debuglevel, dim, tgt);
                Sdata.E0.push_back(st.E[0]);
                Sdata.C.push_back(st.E.back());
                hasBound = hasBound ? (Sbound = hullBox(Sbound, st.E.back()), true)
                                    : (Sbound = st.E.back(), true);

                double nxt = st.E[0][varDim].rightBound();
                cur = (nxt > cur) ? nxt : (cur + b) / 2.0;
            } catch (const std::exception& e) {
                if ((b - cur) <= minWidth || depth >= maxDepth) {
                    IVector pt = edgeBox;
                    pt[varDim] = interval(cur, cur);
                    try {
                        Scaffold st = runScaffold(pt, F, S.G[0],
                                                  eps0, delta, degree, H,
                                                  method, stepBtype, tubedegree, debuglevel, dim,
                                                  makeTarget(pt));
                        Sdata.E0.push_back(st.E[0]);
                        Sdata.C.push_back(st.E.back());
                        hasBound = hasBound ? (Sbound = hullBox(Sbound, st.E.back()), true)
                                            : (Sbound = st.E.back(), true);
                        cur = b;
                        continue;
                    } catch (...) {
                        throw std::runtime_error(
                            std::string("TwoDimEncAlgo: edge segment [") +
                            std::to_string(cur) + "," + std::to_string(b) + "] failed: " + e.what());
                    }
                }
                double mid = (cur + b) / 2.0;
                processEdge(fixedDim, fixedVal, varDim, cur, mid, makeTarget, depth + 1);
                processEdge(fixedDim, fixedVal, varDim, mid, b,   makeTarget, depth + 1);
                return;
            }
        }
    };

    const double xL = S.E[0][0].leftBound(),  xR = S.E[0][0].rightBound();
    const double yL = S.E[0][1].leftBound(),  yR = S.E[0][1].rightBound();

    processEdge(0, interval(xL, xL), 1, yL, yR,
        [&](const IVector& box){ return std::vector<double>{box[0].leftBound(),  box[1].leftBound()}; }, 0);
    if (verbose) std::cout << "  [Boundary 2D] left edge done\n";

    processEdge(0, interval(xR, xR), 1, yL, yR,
        [&](const IVector& box){ return std::vector<double>{box[0].rightBound(), box[1].leftBound()}; }, 0);
    if (verbose) std::cout << "  [Boundary 2D] right edge done\n";

    processEdge(1, interval(yR, yR), 0, xL, xR,
        [&](const IVector& box){ return std::vector<double>{box[0].leftBound(),  box[1].rightBound()}; }, 0);
    if (verbose) std::cout << "  [Boundary 2D] top edge done\n";

    processEdge(1, interval(yL, yL), 0, xL, xR,
        [&](const IVector& box){ return std::vector<double>{box[0].leftBound(),  box[1].leftBound()}; }, 0);
    if (verbose) std::cout << "  [Boundary 2D] bottom edge done\n";

    Sdata.hull = hasBound ? Sbound : S.E[0];
}

// ═══════════════════════════════════════════════════════════════════════════
//  ThreeDimEncAlgo  --  Boundary mode, n >= 3 (EndCover applied to each face)
// ═══════════════════════════════════════════════════════════════════════════
inline void ThreeDimEncAlgo(
    CoverResult& Sdata,
    IMap F, const Scaffold& S,
    double eps0, double delta, int degree,
    double H,
    int method, int stepBtype, int tubedegree, int debuglevel, int dim,
    bool verbose = true)
{
    IVector Sbound(dim);
    bool hasBound = false;

    auto endEnc = [&](const IVector& B, double epsLocal,
                      const std::vector<double>& p, double Hlocal)
        -> std::pair<IVector, IVector>
    {
        Scaffold st = runScaffold(B, F, S.G[0],
                                  epsLocal, delta, degree, Hlocal,
                                  method, stepBtype, tubedegree, debuglevel, dim, p);
        Sdata.E0.push_back(B);
        const IVector& olB = st.E.back();
        Sdata.C.push_back(olB);
        hasBound = hasBound ? (Sbound = hullBox(Sbound, olB), true) : (Sbound = olB, true);
        return {st.E.front(), olB};
    };

    const IVector& B0 = S.E[0];
    for (int fixDim = 0; fixDim < dim; ++fixDim) {
        for (int side = 0; side < 2; ++side) {
            IVector face = B0;
            double fv = side == 0 ? B0[fixDim].leftBound() : B0[fixDim].rightBound();
            face[fixDim] = interval(fv, fv);
            if (verbose)
                std::cout << "  [Boundary nD] face dim=" << fixDim << " side=" << side << "\n";
            endcover::EndCover(face, eps0, H, endEnc, 2000000);
        }
    }
    Sdata.hull = hasBound ? Sbound : B0;
}

// ═══════════════════════════════════════════════════════════════════════════
//  2-D convex hull of a box list (for plotting)
// ═══════════════════════════════════════════════════════════════════════════
struct Point2D { double x, y; };

inline double cross2D(const Point2D& o, const Point2D& a, const Point2D& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

inline std::vector<Point2D> convexHull2D(std::vector<Point2D> pts) {
    const int n = static_cast<int>(pts.size());
    if (n <= 1) return pts;
    std::sort(pts.begin(), pts.end(), [](const Point2D& a, const Point2D& b){
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    pts.erase(std::unique(pts.begin(), pts.end(), [](const Point2D& a, const Point2D& b){
        return std::abs(a.x - b.x) < 1e-14 && std::abs(a.y - b.y) < 1e-14;
    }), pts.end());
    const int m = static_cast<int>(pts.size());
    std::vector<Point2D> h(static_cast<size_t>(2 * m));
    int k = 0;
    for (int i = 0; i < m; ++i) {
        while (k >= 2 && cross2D(h[static_cast<size_t>(k-2)], h[static_cast<size_t>(k-1)],
                                 pts[static_cast<size_t>(i)]) <= 0) --k;
        h[static_cast<size_t>(k++)] = pts[static_cast<size_t>(i)];
    }
    for (int i = m-2, t = k+1; i >= 0; --i) {
        while (k >= t && cross2D(h[static_cast<size_t>(k-2)], h[static_cast<size_t>(k-1)],
                                 pts[static_cast<size_t>(i)]) <= 0) --k;
        h[static_cast<size_t>(k++)] = pts[static_cast<size_t>(i)];
    }
    h.resize(static_cast<size_t>(k - 1));
    return h;
}

inline void saveConvexHull2D(const std::vector<IVector>& boxes, const std::string& fname,
                             bool verbose = true) {
    std::vector<Point2D> pts;
    for (const auto& box : boxes) {
        if (box.dimension() < 2) continue;
        double x0 = box[0].leftBound(), x1 = box[0].rightBound();
        double y0 = box[1].leftBound(), y1 = box[1].rightBound();
        pts.push_back({x0, y0}); pts.push_back({x0, y1});
        pts.push_back({x1, y0}); pts.push_back({x1, y1});
    }
    auto hull = convexHull2D(pts);
    std::ofstream fout(fname);
    fout << std::fixed << std::setprecision(15);
    for (const auto& p : hull) fout << p.x << " " << p.y << "\n";
    if (!hull.empty()) fout << hull[0].x << " " << hull[0].y << "\n";
    if (verbose)
        std::cout << "convex hull -> " << fname << " (" << hull.size() << " pts)\n";
}
