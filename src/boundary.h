#pragma once
// ===========================================================================
//  boundary.h  --  scaffold driver + the Boundary covering mode
//
//  runScaffold()  brings one scaffold to the horizon H and refines it
//                 (Extend / Refine).  It is the common engine of both modes.
//
//  mode = 0 (EndCover)  see EndCorelib.cpp: EndCover covers all of B0.
//  mode = 1 (Boundary)  EndCover is applied to each of the 2n faces of B0 --
//                       the four edges in 2D.  For a flow, the image of the
//                       boundary bounds the image of the interior, so the
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
#include "endcover.h"  // EndCover

// ═══════════════════════════════════════════════════════════════════════════
//  runScaffold : bring S to the horizon, then S.Refine(eps0)
// ═══════════════════════════════════════════════════════════════════════════
//  This is the inner loop of EndCover (Subsection 3.3):
//
//      While (S.t.back() < H)  S.Extend(H, eps0)
//      Q0.push(S.Refine(eps0))
//
//  A new scaffold Init(B) is extended to H.  A scaffold split off by Split()
//  already reaches H; it is a clone of the scaffold it came from, and is
//  evaluated from its own E_0 first.  Refine reads the whole scaffold and
//  refines all of its stages, so it is applied once the scaffold has reached
//  the horizon.  The scaffolds Refine splits off are returned.
inline std::vector<Scaffold> runScaffold(
    Scaffold& S,
    IMap F,
    double eps0, double delta, int degree,
    double H,
    int method, int stepBtype, int tubedegree, int debuglevel, int dim)
{
    // The doubleton carried along the scaffold.  Extend advances it one step
    // per stage and Refine walks it again whenever Split or Bisect changes
    // what the scaffold has to enclose, so from E_0 to E_m the propagation
    // goes through a single doubleton and is never re-enclosed in a box.
    SolutionChain chain(F, degree, S.E[0], stepBtype);

    if (S.m() == 0) {
        //  While (S.T < H)  S.Extend(H, eps0)
        ExtendToHorizon(F, S, chain, eps0, delta, degree, H, stepBtype,
                        tubedegree, dim, debuglevel);
    } else if (!Reevaluate(S, F, chain, eps0, degree, stepBtype, tubedegree,
                           dim, debuglevel) && debuglevel >= progress::REFINE) {
        progress::report("    split-off scaffold evaluated: " + std::to_string(S.m())
                         + " stages, wmax(E_0)=" + progress::sci(wmax(S.E[0]))
                         + " wmax(E_m)=" + progress::sci(wmax(S.E.back())));
    }

    for (int d = 0; d < S.E.back().dimension(); ++d) {
        if (!std::isfinite(S.E.back()[d].leftBound()) ||
            !std::isfinite(S.E.back()[d].rightBound()))
            throw std::runtime_error("runScaffold: non-finite enclosure");
    }

    //  Q0.push(S.Refine(eps0))
    return Refine(S, dim, F, chain, method, stepBtype, degree, eps0,
                  tubedegree, debuglevel);
}

// ═══════════════════════════════════════════════════════════════════════════
//  BoundaryCover  --  Boundary mode: EndCover on each face of B0
// ═══════════════════════════════════════════════════════════════════════════
//  The faces are taken in the order x_1 = lo, x_1 = hi, x_2 = lo, ...  Each
//  face is a box of B0 with one coordinate fixed, and EndCover splits it
//  along its other coordinates.
inline void BoundaryCover(CoverResult& Sdata, const IVector& B0,
                          const MiniScaffold& G0, const endcover::ScaffoldRun& run,
                          int debuglevel, bool verbose = true)
{
    const int dim = B0.dimension();
    for (int fixDim = 0; fixDim < dim; ++fixDim) {
        for (int side = 0; side < 2; ++side) {
            IVector face = B0;
            const double fv = side == 0 ? B0[fixDim].leftBound() : B0[fixDim].rightBound();
            face[fixDim] = interval(fv, fv);
            const std::string name = "face x" + std::to_string(fixDim + 1) + "=" + progress::num(fv);
            if (debuglevel >= progress::SUBBOX) progress::report(name);
            endcover::EndCoverResult res = endcover::EndCover(face, G0, run, 2000000, debuglevel);
            Sdata.C.insert(Sdata.C.end(), res.cover.begin(), res.cover.end());
            Sdata.E0.insert(Sdata.E0.end(), res.initBoxes.begin(), res.initBoxes.end());
            if (verbose)
                std::cout << "  [Boundary] " << name << ": " << res.cover.size() << " boxes\n";
        }
    }
    Sdata.hull = B0;
    for (size_t i = 0; i < Sdata.C.size(); ++i)
        Sdata.hull = (i == 0) ? Sdata.C[0] : hullBox(Sdata.hull, Sdata.C[i]);
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
