// ===========================================================================
//  tools/trajectory.cpp  --  Encoretraj
//
//  A validated trajectory printer, used for the top row of Figure 1.
//
//  It repeatedly applies StepB (the same routine EnCorelib uses)
//  on a uniform time grid and prints one line per grid point:
//
//      t   lo_1 hi_1   lo_2 hi_2   ...   lo_n hi_n
//
//  Every printed box is a rigorous enclosure of x(t).  Unlike EnCorelib this
//  tool makes no attempt to control the width: it simply shows where the
//  validated solution goes, which is exactly what Figure 1 (top right)
//  displays.  If an enclosure blows up the program stops and says so.
//
//  Usage
//      ./Encoretraj n  var_1..var_n  f_1..f_n  order T dt  lo_1 hi_1 .. lo_n hi_n
//
//  Example (undamped pendulum from the point (0,2), to t = 20, step 0.02):
//      ./Encoretraj 2 x y y "-sin(x)" 20 20 0.02  0 0  2 2
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "../types.h"
#include "../stepAB.h"
#include "../symparse.h"

int main(int argc, char* argv[]) {
    try {
        if (argc < 6) throw std::runtime_error("too few arguments");
        const int n = std::stoi(argv[1]);
        if (n <= 0) throw std::runtime_error("n must be > 0");
        if (argc < 5 + 4 * n) throw std::runtime_error("wrong number of arguments");

        std::vector<std::string> SVar, SFun;
        for (int i = 0; i < n; ++i) SVar.push_back(argv[2 + i]);
        for (int i = 0; i < n; ++i) SFun.push_back(expand_expression(argv[2 + n + i]));

        const int    order = std::stoi(argv[2 + 2 * n]);
        const double T     = std::stod(argv[3 + 2 * n]);
        const double dt    = std::stod(argv[4 + 2 * n]);
        if (order < 3 || T <= 0.0 || dt <= 0.0)
            throw std::runtime_error("need order >= 3, T > 0, dt > 0");

        IVector B(n);
        for (int i = 0; i < n; ++i) {
            double lo = std::stod(argv[5 + 2 * n + 2 * i]);
            double hi = std::stod(argv[5 + 2 * n + 2 * i + 1]);
            B[i] = capd::interval(lo, hi);
        }

        IMap F(Convert_to_IMap(SVar, SFun), 3.0);

        std::cout << "# t";
        for (int i = 0; i < n; ++i)
            std::cout << "  " << SVar[i] << "_lo  " << SVar[i] << "_hi";
        std::cout << "\n" << std::setprecision(15);

        const int steps = static_cast<int>(std::ceil(T / dt));
        IVector cur = B;
        double  t   = 0.0;
        for (int k = 0; k <= steps; ++k) {
            std::cout << t;
            for (int i = 0; i < n; ++i)
                std::cout << ' ' << cur[i].leftBound() << ' ' << cur[i].rightBound();
            std::cout << "\n";
            if (k == steps) break;

            const double h = std::min(dt, T - t);
            cur = stepBlohner(F, order, h, cur);
            t  += h;

            for (int i = 0; i < n; ++i) {
                if (!std::isfinite(cur[i].leftBound()) || !std::isfinite(cur[i].rightBound())) {
                    std::cerr << "enclosure blew up at t = " << t << "\n";
                    return 1;
                }
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Encoretraj: " << e.what() << "\n\n"
                  << "Usage: Encoretraj n var_1..var_n f_1..f_n order T dt "
                     "lo_1 hi_1 .. lo_n hi_n\n"
                     "e.g.:  Encoretraj 2 x y y \"-sin(x)\" 20 20 0.02  0 0  2 2\n";
        return 2;
    }
}
