// ===========================================================================
//  repro/competitors/capd_compare.cpp  --  the "CAPD" rows of Table 2
//
//  Runs CAPD's Cr-Lohner method (CnRect2Set with r = 3, the setting used in
//  the paper) once from B0 to time H and prints the resulting end-enclosure
//  in the same format as EnCorelib, so the rho column of Table 2 can be
//  recomputed:
//
//      rho(CAPD) = wmax(CAPD) / wmax(Ours(1.0)).
//
//  When CAPD aborts -- which happens on several of the long-horizon rows --
//  the program prints "No Output" and exits with status 1, matching the
//  entry in Table 2.
//
//  Build:   make capdref
//  Usage:   ./capdref n var_1..var_n f_1..f_n order T lo_1 hi_1 .. lo_n hi_n
//  Example: ./capdref 2 x y 2*x-2*x*y -y+x*y 20 4  0.9 1.1  2.9 3.1
// ===========================================================================

#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>
#include <stdexcept>

#include "../../types.h"
#include "../../symparse.h"

int main(int argc, char* argv[]) {
    try {
        if (argc < 5) throw std::runtime_error("too few arguments");
        const int n = std::stoi(argv[1]);
        if (n <= 0) throw std::runtime_error("n must be > 0");
        if (argc < 4 + 4 * n) throw std::runtime_error("wrong number of arguments");

        std::vector<std::string> SVar, SFun;
        for (int i = 0; i < n; ++i) SVar.push_back(argv[2 + i]);
        for (int i = 0; i < n; ++i) SFun.push_back(expand_expression(argv[2 + n + i]));

        const int    order = std::stoi(argv[2 + 2 * n]);
        const double H     = std::stod(argv[3 + 2 * n]);

        IVector B(n);
        for (int i = 0; i < n; ++i)
            B[i] = capd::interval(std::stod(argv[4 + 2 * n + 2 * i]),
                                  std::stod(argv[4 + 2 * n + 2 * i + 1]));

        IMap F(Convert_to_IMap(SVar, SFun), 3.0);

        auto t0 = std::chrono::high_resolution_clock::now();
        ICnOdeSolver solver(F, order);
        ICnTimeMap   timeMap(solver);
        CnRect2Set   s(B, 3.0);                 // Cr-Lohner with r = 3
        IVector      E = timeMap(capd::interval(H), s);
        auto t1 = std::chrono::high_resolution_clock::now();

        const double secs =
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "time(ms)=" << static_cast<long long>(secs * 1000.0 + 0.5) << "\n";
        std::cout << "Hull(T=" << H << ")=" << E << "\n";
        std::cout << "B1=(";
        for (int i = 0; i < n; ++i) {
            std::cout << (E[i].leftBound() + E[i].rightBound()) / 2.0;
            if (i + 1 < n) std::cout << ", ";
        }
        std::cout << ") +- (";
        for (int i = 0; i < n; ++i) {
            std::cout << (E[i].rightBound() - E[i].leftBound()) / 2.0;
            if (i + 1 < n) std::cout << ", ";
        }
        std::cout << ")\nwmax(B1)=" << wmax(E) << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cout << "No Output\n";
        std::cerr << "capdref: " << e.what() << "\n";
        return 1;
    }
}
