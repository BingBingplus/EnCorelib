// ===========================================================================
//  tools/audit_cover.cpp  --  Encoreaudit
//
//  A box-by-box check of a cover computed by EnCorelib.
//
//  EnCorelib (output_mode >= 3) writes the sub-boxes of B0 to E0.txt and
//  their end-enclosures, the boxes of the cover C, to E1.txt; line k of one
//  belongs to line k of the other.  For every k this tool takes the centre
//  and the corners of E0[k], encloses the solution through each of them up
//  to time H with a plain validated integrator, and checks that each such
//  enclosure meets E1[k].  An enclosure disjoint from E1[k] proves that
//  E1[k] is not an end-enclosure of E0[k].  repro/verify_cover.sh checks
//  that every sample image lies in some box of C; this check is stronger,
//  because it ties every box of C to its own sub-box.
//
//  Usage
//      ./Encoreaudit n  var_1..var_n  f_1..f_n  order H  E0.txt E1.txt
//
//  Prints one summary line; the exit status is 1 when some box fails.
// ===========================================================================

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

#include "../types.h"
#include "../symparse.h"

static std::vector<IVector> readBoxes(const std::string& path, int n) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    const std::regex pair(R"(\[\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\])");
    std::vector<IVector> out;
    std::string line;
    while (std::getline(in, line)) {
        IVector B(n);
        int i = 0;
        for (std::sregex_iterator it(line.begin(), line.end(), pair), end;
             it != end && i < n; ++it, ++i)
            B[i] = interval(std::stod((*it)[1]), std::stod((*it)[2]));
        if (i == n) out.push_back(B);
    }
    return out;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) throw std::runtime_error("too few arguments");
        const int n = std::stoi(argv[1]);
        if (n <= 0 || argc != 6 + 2 * n) throw std::runtime_error("wrong number of arguments");
        std::vector<std::string> V, F;
        for (int i = 0; i < n; ++i) V.push_back(argv[2 + i]);
        for (int i = 0; i < n; ++i) F.push_back(expand_expression(argv[2 + n + i]));
        const int    order = std::stoi(argv[2 + 2 * n]);
        const double H     = std::stod(argv[3 + 2 * n]);
        const auto E0 = readBoxes(argv[4 + 2 * n], n);
        const auto E1 = readBoxes(argv[5 + 2 * n], n);
        if (E0.size() != E1.size()) throw std::runtime_error("E0 and E1 differ in length");
        IMap f(Convert_to_IMap(V, F), 3.0);

        int bad = 0, samples = 0, skipped = 0;
        double worst = 0.0;
        for (size_t k = 0; k < E0.size(); ++k) {
            std::vector<IVector> pts;
            IVector c(n);
            for (int i = 0; i < n; ++i)
                c[i] = interval((E0[k][i].leftBound() + E0[k][i].rightBound()) / 2.0);
            pts.push_back(c);
            for (int mask = 0; mask < (1 << n); ++mask) {
                IVector p(n);
                for (int i = 0; i < n; ++i)
                    p[i] = interval((mask >> i) & 1 ? E0[k][i].rightBound()
                                                    : E0[k][i].leftBound());
                pts.push_back(p);
            }
            bool boxBad = false;
            for (const IVector& p : pts) {
                IVector Y;
                try {
                    IOdeSolver solver(f, order);
                    ITimeMap timeMap(solver);
                    C0HORect2Set s(p);
                    Y = timeMap(interval(H), s);
                } catch (const std::exception&) { ++skipped; continue; }
                ++samples;
                double gap = 0.0;                    // distance of Y from E1[k]
                for (int i = 0; i < n; ++i) {
                    gap = std::max(gap, E1[k][i].leftBound() - Y[i].rightBound());
                    gap = std::max(gap, Y[i].leftBound() - E1[k][i].rightBound());
                }
                if (gap > 0.0) { boxBad = true; worst = std::max(worst, gap); }
            }
            if (boxBad) ++bad;
        }
        std::printf("%s  %zu boxes, %d images: %d boxes miss an image of their sub-box",
                    bad ? "FAIL" : "ok  ", E1.size(), samples, bad);
        if (bad) std::printf(" (largest gap %.3g)", worst);
        if (skipped) std::printf(", %d samples could not be integrated", skipped);
        std::printf("\n");
        return bad ? 1 : 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Encoreaudit: %s\n\nUsage: Encoreaudit n var_1..var_n f_1..f_n"
                             " order H E0.txt E1.txt\n", e.what());
        return 2;
    }
}
