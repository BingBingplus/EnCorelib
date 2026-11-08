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
//  The scaffold work for one popped scaffold -- Extend to H, then Refine --
//  is done by the caller-supplied callable `run` (runScaffold in boundary.h),
//  which returns the initial boxes of the scaffolds that Refine's Split()
//  split off.  Their scaffolds Init(E0-) are pushed on the queue, and
//  S.E.back() goes into C.
//
//  When the scaffold of a sub-box cannot be carried to the horizon, its
//  initial box is split in two by halving a longest side, and both halves go
//  back on the queue as new scaffolds.
//
//  With debug >= 1 every processed scaffold is reported on one line; see
//  progress.h.
// ===========================================================================

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
#include "types.h"
#include "progress.h"

namespace endcover {

// Raised when the scaffold of a sub-box of B0 cannot be carried to the
// horizon and the sub-box cannot be halved in double precision.
struct PrecisionLimit : std::runtime_error {
    using std::runtime_error::runtime_error;
};

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
    std::vector<IVector> initBoxes;  // S.E[0] of the scaffold behind each box of C
};

// Signature of the scaffold run:  run(S) brings S to the horizon, refines it
// and returns the initial boxes of the scaffolds Refine split off.
using ScaffoldRun = std::function<std::vector<IVector>(Scaffold&)>;

//  B0          the initial box
//  G0          the mini-scaffold every new scaffold starts from (Init)
//  run         Extend and Refine, see above
//  maxSplits   bound on the number of splits
//  debuglevel  1 = one line per processed scaffold
inline EndCoverResult EndCover(const IVector& B0, const MiniScaffold& G0,
                               const ScaffoldRun& run,
                               int maxSplits  = 2000000,
                               int debuglevel = 0) {
    if (maxSplits <= 0)
        throw std::invalid_argument("EndCover: maxSplits must be positive");

    // The working queue of scaffolds, each held as Init(B) of its initial box
    // B.  `depth` is the number of splits between B0 and B; it is only
    // reported.
    struct Pending { IVector B; int depth; };
    std::vector<Pending> Q0;
    Q0.push_back({B0, 0});
    EndCoverResult result;
    result.cover.reserve(256);
    result.initBoxes.reserve(256);
    long splits = 0;

    long   runs = 0;
    int    deepest = 0;
    double done = 0.0;                  // share of B0 covered so far

    // One line per processed scaffold, in columns under a header.
    const auto col = progress::column;
    auto report = [&](int depth, const IVector& B, const std::string& wEnd,
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
        const Pending P = Q0.back();
        Q0.pop_back();
        const IVector& B = P.B;
        Scaffold S = Init(B, G0);
        ++runs;
        deepest = std::max(deepest, P.depth);

        std::string outcome, wEnd = "-";   // for the progress report
        try {
            const std::vector<IVector> Qout = run(S);
            const int k = static_cast<int>(Qout.size());
            // The j-th split-off is j splits below B; S.E[0] is k splits below.
            for (int j = 0; j < k; ++j)
                Q0.push_back({Qout[static_cast<size_t>(j)], P.depth + j + 1});
            splits += k;

            result.cover.push_back(S.E.back());
            result.initBoxes.push_back(S.E[0]);
            done += shareOf(S.E[0], B0);
            deepest = std::max(deepest, P.depth + k);
            wEnd = progress::sci(wmax(S.E.back()));
            outcome = "resolved";
            if (k > 0)
                outcome += ", E_0 -> " + progress::sci(wmax(S.E[0])) + ", "
                         + std::to_string(k) + " split off";
        } catch (const std::exception& ex) {
            // The scaffold could not be carried to the horizon: split B.
            IVector lower, upper;
            if (!halveLongestSide(B, lower, upper)) {
                if (debuglevel >= progress::SUBBOX)
                    report(P.depth, B, wEnd, std::string("stopped: ") + ex.what());
                throw PrecisionLimit(std::string(ex.what())
                    + "\nThe sub-box of B0 on which this happened (depth "
                    + std::to_string(P.depth) + " of the subdivision) cannot be"
                      " halved in double precision.");
            }
            Q0.push_back({upper, P.depth + 1});
            Q0.push_back({lower, P.depth + 1});
            ++splits;
            outcome = std::string("split: ") + ex.what();
        }

        if (splits > maxSplits)
            throw std::runtime_error("EndCover: exceeded maxSplits");
        if (debuglevel >= progress::SUBBOX) report(P.depth, B, wEnd, outcome);
    }

    if (debuglevel >= progress::SUBBOX)
        progress::report("EndCover finished: " + std::to_string(runs)
                         + (runs == 1 ? " scaffold" : " scaffolds") + " processed, deepest level "
                         + std::to_string(deepest) + ", #(C)=" + std::to_string(result.cover.size()));
    return result;
}

} // namespace endcover
