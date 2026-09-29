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
// ===========================================================================

#include <functional>
#include <stdexcept>
#include <vector>
#include "types.h"

namespace endcover {

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
                               int maxSplits = 2000000,
                               int coverSet  = COVER_ALL) {
    if (maxSplits <= 0)
        throw std::invalid_argument("EndCover: maxSplits must be positive");

    std::vector<IVector> Q0 = {B0};     // the working queue of scaffolds
    EndCoverResult result;
    result.cover.reserve(256);
    result.initBoxes.reserve(256);
    int splits = 0;

    while (!Q0.empty()) {
        IVector B = Q0.back();
        Q0.pop_back();

        const std::vector<double> mid = midpoint_vec(B);
        try {
            auto [ulB, olB] = endEnc(B, eps0, mid, H);
            const bool resolved = boxesEqual(ulB, B);
            // A box that cannot be subdivided any further (B is a point) is
            // treated as resolved: dropping it would break the cover.
            const std::vector<IVector> children = resolved ? std::vector<IVector>()
                                                           : splitAllDims(B);

            if (resolved || children.empty() || coverSet == COVER_ALL) {
                result.cover.push_back(olB);
                result.initBoxes.push_back(ulB);
            }

            if (!children.empty()) {
                // Refine had to Split: the rest of B is still uncovered.
                if (++splits > maxSplits)
                    throw std::runtime_error("EndCover: exceeded maxSplits");
                Q0.insert(Q0.end(), children.begin(), children.end());
            }
        } catch (const std::exception&) {
            // The scaffold run failed (typically an enclosure blow-up).
            // Do not drop B: subdivide and retry on the children.
            if (++splits > maxSplits) throw;
            auto children = splitAllDims(B);
            if (children.empty()) throw;       // a point we cannot enclose
            Q0.insert(Q0.end(), children.begin(), children.end());
        }
    }
    return result;
}

} // namespace endcover
