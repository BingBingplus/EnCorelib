#pragma once
// ===========================================================================
//  progress.h  --  progress reports while EndCover runs
//
//  The `debug` argument of EndCorelib chooses how much is reported.  Each
//  level adds to the one below it:
//
//      0   nothing
//      1   one line for every scaffold that EndCover has taken from its
//          queue, in columns: its number, the depth of its initial box B in
//          the subdivision of B0, wmax(B), wmax of its end-enclosure, the
//          number of scaffolds still queued, the share of B0 covered so far,
//          and the outcome (resolved, or split and at which time)
//      2   one line for every scaffold laid out to the horizon and for
//          every pass of Refine
//      3   one line for every stage appended by Extend (t_i, the stage
//          length, E_i), and the containment self-checks, which should
//          never fire
//
//  Every line starts with the wall-clock time since the program started.
//  The reports go to standard error; standard output carries the results
//  only, so it can be redirected or parsed exactly as without them.
// ===========================================================================

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "types.h"

namespace progress {

// The `debug` levels described above.
enum Level { SUBBOX = 1, REFINE = 2, STAGE = 3 };

// Wall-clock seconds since the first call.  main() makes that call when the
// program starts.
inline double elapsed() {
    using Clock = std::chrono::steady_clock;
    static const Clock::time_point start = Clock::now();
    return std::chrono::duration<double>(Clock::now() - start).count();
}

// x in short scientific notation, e.g. 1.2e-04.
inline std::string sci(double x, int digits = 1) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(digits) << x;
    return os.str();
}

// x with `digits` significant digits, e.g. 4.96028.
inline std::string num(double x, int digits = 6) {
    std::ostringstream os;
    os << std::setprecision(digits) << x;
    return os.str();
}

// A share in [0,1] as a percentage with one decimal, e.g. 56.3%.
inline std::string percent(double share) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << 100.0 * share << '%';
    return os.str();
}

// A box as {[a, b],[c, d],...} with `digits` significant digits.
inline std::string box(const IVector& B, int digits = 6) {
    std::ostringstream os;
    os << std::setprecision(digits) << '{';
    for (int i = 0; i < B.dimension(); ++i) {
        if (i > 0) os << ',';
        os << '[' << B[i].leftBound() << ", " << B[i].rightBound() << ']';
    }
    os << '}';
    return os.str();
}

// s right-aligned in a field of `width` characters.
inline std::string column(const std::string& s, int width) {
    const int pad = width - static_cast<int>(s.size());
    return (pad > 0 ? std::string(static_cast<size_t>(pad), ' ') : std::string()) + s;
}

// Write one report line, prefixed with the elapsed time.
inline void report(const std::string& text) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << std::setw(8) << elapsed()
       << "s  " << text << '\n';
    std::cerr << os.str() << std::flush;
}

// A header line for columns of report lines, aligned with their text.
inline void header(const std::string& text) {
    std::ostringstream os;
    os << std::setw(9) << "time" << "  " << text << '\n';
    std::cerr << os.str() << std::flush;
}

} // namespace progress
