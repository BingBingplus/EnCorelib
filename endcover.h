#pragma once
// ===========================================================================
//  endcover.h  --  EndCover_f(B0, H, eps0) -> C    (paper, Subsection 3.3)
//
//      C <- {}                          initialise the eps0-cover
//      Q0 <- {Init(B0)}                 initialise the queue of scaffolds
//      While (Q0 != {})
//          S <- Q0.pop()
//          Q0.push(S.Refine(eps0))
//          While (S.t.back() < H)
//              S.Extend(H, eps0)
//              Q0.push(S.Refine(eps0))
//          C.push(S.E.back())
//      Return C
//
//  In this implementation the scaffold work for one sub-box is done by the
//  caller-supplied callable `endEnc`, which returns the pair
//
//      (ulB, olB) = (S.E[0] after Refine,  S.E.back())
//
//  Refine's Split() shrinks S.E[0]; whenever ulB is a proper sub-box of the
//  box B we handed in, the complement of ulB inside B has not been covered,
//  so B is subdivided into 2^n children and they are pushed back on the
//  queue.  Those children are exactly the "split-off scaffolds returned by
//  S.Refine(eps0)" of the pseudo-code above.
//
//  With debug >= 1 every processed sub-box is reported on one line; see
//  progress.h.
// ===========================================================================

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
#include "types.h"
#include "progress.h"

namespace endcover {

// Raised when the scaffold run fails on a sub-box of B0 that is a single
// point.  Subdividing cannot help then: the solution through that point
// cannot be enclosed up to the horizon.
struct PointFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// ---------------------------------------------------------------------------
//  Which boxes end up in C
// ---------------------------------------------------------------------------
//   COVER_ALL    every processed box contributes its end-enclosure.  This is
//                the setting used for the experiments in the paper.
//   COVER_LEAVES only boxes that were resolved without splitting contribute;
//                the result is the strict eps0-cover of equation (1.3), in
//                which every box satisfies wmax <= eps0.
enum CoverSet { COVER_ALL = 0, COVER_LEAVES = 1 };

inline std::vector<double> midpoint_vec(const IVector& box) {
    const int n = box.dimension();
    std::vector<double> m(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        m[static_cast<size_t>(i)] = (box[i].leftBound() + box[i].rightBound()) / 2.0;
    return m;
}

inline bool boxesEqual(const IVector& a, const IVector& b, double tol = 0.0) {
    if (a.dimension() != b.dimension()) return false;
    for (int i = 0; i < a.dimension(); ++i)
        if (std::abs(a[i].leftBound()  - b[i].leftBound())  > tol ||
            std::abs(a[i].rightBound() - b[i].rightBound()) > tol)
            return false;
    return true;
}

// Split a box into 2^d children by halving every coordinate that has
// positive width (d of them).  Degenerate coordinates are left alone:
// halving them would only produce duplicate children, which for a box with
// a flat face -- Eijgenraam's example has y(0) = 0 -- multiplies the work
// by 2 at every level of the subdivision without covering anything new.
// If no coordinate can be split the box is a single point and the returned
// vector is empty.
inline std::vector<IVector> splitAllDims(const IVector& box) {
    const int n = box.dimension();

    std::vector<int> splitDims;
    std::vector<double> a(static_cast<size_t>(n)), m(static_cast<size_t>(n)), b(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        a[static_cast<size_t>(i)] = box[i].leftBound();
        b[static_cast<size_t>(i)] = box[i].rightBound();
        m[static_cast<size_t>(i)] = (a[static_cast<size_t>(i)] + b[static_cast<size_t>(i)]) / 2.0;
        if (a[static_cast<size_t>(i)] < m[static_cast<size_t>(i)] &&
            m[static_cast<size_t>(i)] < b[static_cast<size_t>(i)])
            splitDims.push_back(i);
    }
    if (splitDims.empty()) return {};

    const size_t count = static_cast<size_t>(1) << splitDims.size();
    std::vector<IVector> out;
    out.reserve(count);
    for (size_t mask = 0; mask < count; ++mask) {
        IVector child = box;
        for (size_t k = 0; k < splitDims.size(); ++k) {
            const int i = splitDims[k];
            const bool upper = ((mask >> k) & 1u) != 0u;
            child[i] = upper ? interval(m[static_cast<size_t>(i)], b[static_cast<size_t>(i)])
                             : interval(a[static_cast<size_t>(i)], m[static_cast<size_t>(i)]);
        }
        out.push_back(child);
    }
    return out;
}

// The share of B0 that its sub-box B makes up, measured in the coordinates
// in which B0 has positive width (all of B0 when B0 is a point).
inline double shareOf(const IVector& B, const IVector& B0) {
    double s = 1.0;
    for (int i = 0; i < B0.dimension(); ++i) {
        const double w0 = B0[i].rightBound() - B0[i].leftBound();
        if (w0 > 0.0) s *= (B[i].rightBound() - B[i].leftBound()) / w0;
    }
    return s;
}

struct EndCoverResult {
    std::vector<IVector> cover;      // C
    std::vector<IVector> initBoxes;  // the sub-boxes of B0 that were integrated
};

// Signature of the per-box scaffold run:
//   endEnc(B, eps0, m(B), H) -> (ulB, olB)
using EndEncFn = std::function<std::pair<IVector, IVector>(
    const IVector&, double, const std::vector<double>&, double)>;

inline EndCoverResult EndCover(const IVector& B0, double eps0, double H,
                               const EndEncFn& endEnc,
                               int maxSplits  = 2000000,
                               int coverSet   = COVER_ALL,
                               int debuglevel = 0) {
    if (maxSplits <= 0)
        throw std::invalid_argument("EndCover: maxSplits must be positive");

    // The working queue of scaffolds.  `depth` is the number of
    // subdivisions between B0 and the sub-box; it is only reported.
    struct Pending { IVector box; int depth; };
    std::vector<Pending> Q0 = {{B0, 0}};
    EndCoverResult result;
    result.cover.reserve(256);
    result.initBoxes.reserve(256);
    int splits = 0;

    long   runs = 0, nResolved = 0;
    int    deepest = 0;
    double done = 0.0;                  // share of B0 resolved so far

    // One line per processed sub-box, in columns under a header.
    const auto col = progress::column;
    auto subBoxReport = [&](int depth, const IVector& B, const std::string& wEnd,
                            const std::string& outcome) {
        progress::report(col(std::to_string(runs), 5) + col(std::to_string(depth), 6)
                         + col(progress::sci(wmax(B)), 9) + col(wEnd, 10)
                         + col(std::to_string(Q0.size()), 6)
                         + col(progress::percent(done), 8) + "  " + outcome);
    };
    if (debuglevel >= progress::SUBBOX)
        progress::header(col("run", 5) + col("depth", 6) + col("wmax(B)", 9)
                         + col("wmax(E_m)", 10) + col("queue", 6) + col("done", 8)
                         + "  outcome");

    while (!Q0.empty()) {
        const IVector B     = Q0.back().box;
        const int     depth = Q0.back().depth;
        Q0.pop_back();
        ++runs;
        deepest = std::max(deepest, depth);

        const std::vector<double> mid = midpoint_vec(B);
        std::vector<IVector> children;
        std::string outcome, wEnd = "-";   // for the progress report
        try {
            auto [ulB, olB] = endEnc(B, eps0, mid, H);
            const bool resolved = boxesEqual(ulB, B);
            // A box that cannot be subdivided any further (B is a point) is
            // treated as resolved: dropping it would break the cover.
            if (!resolved) children = splitAllDims(B);

            if (resolved || children.empty() || coverSet == COVER_ALL) {
                result.cover.push_back(olB);
                result.initBoxes.push_back(ulB);
            }
            if (children.empty()) {
                ++nResolved;
                done += shareOf(B, B0);
                outcome = "resolved";
            } else {
                // Refine had to Split: the rest of B is still uncovered.
                outcome = "split, E_0 -> " + progress::sci(wmax(ulB));
            }
            wEnd = progress::sci(wmax(olB));
        } catch (const std::exception& ex) {
            // The scaffold run failed (typically an enclosure blow-up).
            // Do not drop B: subdivide and retry on the children.
            children = splitAllDims(B);
            outcome  = std::string("failed: ") + ex.what();
            if (children.empty()) {
                if (debuglevel >= progress::SUBBOX) subBoxReport(depth, B, wEnd, outcome);
                throw PointFailure(std::string(ex.what())
                    + "\nThe sub-box of B0 on which this happened is a single point"
                      " (depth " + std::to_string(depth) + " of the subdivision),"
                      "\nso subdividing B0 further cannot help: the solution through"
                      " that point\ncannot be enclosed up to the horizon.");
            }
            outcome += ", split into " + std::to_string(children.size());
        }

        if (!children.empty()) {
            if (++splits > maxSplits)
                throw std::runtime_error("EndCover: exceeded maxSplits");
            for (const IVector& c : children) Q0.push_back({c, depth + 1});
        }
        if (debuglevel >= progress::SUBBOX) subBoxReport(depth, B, wEnd, outcome);
    }

    if (debuglevel >= progress::SUBBOX)
        progress::report("EndCover finished: " + std::to_string(runs)
                         + (runs == 1 ? " sub-box" : " sub-boxes") + " processed, "
                         + std::to_string(nResolved) + " resolved, deepest level "
                         + std::to_string(deepest) + ", #(C)=" + std::to_string(result.cover.size()));
    return result;
}

} // namespace endcover
