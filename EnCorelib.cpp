// ===========================================================================
//  EnCorelib.cpp  --  EndCover for the Initial Value Problem
//
//  Reference implementation for
//
//      B. Zhang and C. Yap,
//      "End Cover for Initial Value Problem: Complete Validated Algorithm
//       with Complexity Analysis".
//
//  Given an autonomous ODE  x' = f(x),  a box B0, a horizon H and an error
//  bound eps > 0, the program computes a finite set C of boxes with
//
//      End(B0,H)  subset  Union(C)  subset  End(B0,H) + [-eps,eps]^n .
//
//  ------------------------------------------------------------------------
//  USAGE
//
//    ./EnCorelib output_mode mode method stepB tubedegree n
//              var_1 ... var_n  f_1 ... f_n
//              eps order T debug
//              lo_1 hi_1 ... lo_n hi_n
//
//    output_mode  progressive output level, 0..4 (see below)
//    mode         0 = EndCover (cover all of B0)   1 = Boundary (cover dB0)
//    method       which refinement operators Refine may use (see refine.h)
//                   0 = Bisect + EulerTube + Split   (Refine of Section 5)
//                   1 = Bisect + Split
//                   2 = Split only
//    stepB        0 = Lohner              1 = Lohner + logNorm
//                 2 = Lohner-HO + logNorm 3 = Lohner-HO
//    tubedegree   0 = order-1             1 = Euler tube   p>=2 = Taylor-p
//    n            dimension
//    eps          the error bound eps of equation (1.3)
//    order        Taylor order k used by StepB (the paper uses 20)
//    T            the horizon H
//    debug        0 = silent, 1 = report containment checks
//    lo_i hi_i    the initial box B0
//
//  OUTPUT LEVELS (progressive: level k also prints everything below k)
//
//    0   time(ms)
//    1   + Hull(T), B1 = (mid) +- (half-width), wmax
//    2   + #(C), #E0, wmax over the individual boxes of C
//    3   + writes E0.txt, E1.txt (E1 = the cover C only), convex_hull.txt (2D)
//    4   + writes E_0.txt, E_1.txt for plotting: E_1 additionally contains
//          the images of every E0 box at times 0.1/0.4/0.7 and, in 2D, the
//          trajectories of the four corners of B0
//
//  Every run appends one row to the cumulative log "out.txt".
// ===========================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "boundary.h"   // runScaffold, Boundary modes, convex hull
#include "endcover.h"   // EndCover
#include "symparse.h"

using namespace std;

// ---------------------------------------------------------------------------
//  Small helpers
// ---------------------------------------------------------------------------

static const char* modeName(int m)   { return m == 0 ? "EndCover" : "Boundary"; }
static const char* methodName(int m) {
    switch (m) {
    case REFINE_TUBE:   return "tube";
    case REFINE_BISECT: return "bisect";
    case REFINE_SPLITONLY: return "splitonly";
    default:            return "?";
    }
}
static const char* stepBName(int m) {
    switch (m) {
    case 0: return "lohner";
    case 1: return "lohner+lognorm";
    case 2: return "lohnerHO+lognorm";
    case 3: return "lohnerHO";
    default: return "?";
    }
}

static void saveBoxes(const std::vector<IVector>& boxes, const std::string& fname,
                      bool verbose) {
    std::ofstream fout(fname);
    if (!fout) { std::cerr << "Cannot open " << fname << "\n"; return; }
    fout << std::scientific << std::setprecision(17);
    for (size_t i = 0; i < boxes.size(); ++i) {
        fout << "Box " << i << ": ";
        for (int d = 0; d < boxes[i].dimension(); ++d) {
            fout << "[" << boxes[i][d].leftBound() << ", " << boxes[i][d].rightBound() << "]";
            if (d + 1 < boxes[i].dimension()) fout << " x ";
        }
        fout << "\n";
    }
    if (verbose)
        std::cout << "written -> " << fname << " (" << boxes.size() << " boxes)\n";
}

// "(m1, m2, ...) +- (r1, r2, ...)" -- the B1 column of Table 2.
static std::string midRadString(const IVector& B) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6);
    os << "(";
    for (int i = 0; i < B.dimension(); ++i) {
        os << (B[i].leftBound() + B[i].rightBound()) / 2.0;
        if (i + 1 < B.dimension()) os << ", ";
    }
    os << ") +- (";
    for (int i = 0; i < B.dimension(); ++i) {
        os << (B[i].rightBound() - B[i].leftBound()) / 2.0;
        if (i + 1 < B.dimension()) os << ", ";
    }
    os << ")";
    return os.str();
}

// One tab-separated row per run, appended to out.txt.
static void appendLog(const std::string& tag, int mode, int method, int stepBtype,
                      double eps, int order, double H, int n,
                      const IVector& B0, const IVector& B1,
                      size_t nCover, double timeSec) {
    const bool fresh = !std::ifstream("out.txt").good();
    std::ofstream f("out.txt", std::ios::app);
    if (!f) return;
    if (fresh) {
        f << "# EndCover run log.  One row per run; columns are tab separated.\n"
             "# tag\tmode\tmethod\tstepB\teps\torder\tH\tn\tB1_mid\tB1_halfwidth"
             "\twmax(B1)\ttime_s\t#C\tB0\n";
    }
    f << std::setprecision(12);
    f << tag << '\t' << modeName(mode) << '\t' << methodName(method) << '\t'
      << stepBtype << '\t' << eps << '\t' << order << '\t' << H << '\t' << n << '\t';
    f << '(';
    for (int i = 0; i < n; ++i) { f << (B1[i].leftBound()+B1[i].rightBound())/2.0; if (i+1<n) f << ','; }
    f << ")\t(";
    for (int i = 0; i < n; ++i) { f << (B1[i].rightBound()-B1[i].leftBound())/2.0; if (i+1<n) f << ','; }
    f << ")\t" << wmax(B1) << '\t' << timeSec << '\t' << nCover << '\t' << B0 << '\n';
}

// Init(B): the 0-stage scaffold of Subsection 3.3.
static Scaffold Init(const IVector& B, double delta, double H) {
    MiniScaffold G;
    G.mu    = {1000.0};   // sentinel: no logNorm bound known yet
    G.delta = delta;
    G.hTube = H;
    G.ell   = 0;
    G.E     = {B};
    G.F     = {B};
    G.nTube = 1;
    return Scaffold{{0.0}, {B}, {B}, {G}};
}

static void usage() {
    std::cerr <<
      "Usage: EnCorelib output_mode mode method stepB tubedegree n \\\n"
      "               var1 ... varn  f1 ... fn \\\n"
      "               eps order T debug  lo1 hi1 ... lon hin\n\n"
      "  output_mode 0=time 1=+hull 2=+counts 3=+E0/E1 files 4=+plot files\n"
      "  mode        0=EndCover 1=Boundary\n"
      "  method      0=tube (default) 1=bisect 2=splitonly\n"
      "  stepB       0=lohner 1=lohner+lognorm 2=lohnerHO+lognorm 3=lohnerHO\n"
      "  tubedegree  0=order-1 1=Euler p>=2=Taylor-p\n";
}

// ═══════════════════════════════════════════════════════════════════════════
//  main
// ═══════════════════════════════════════════════════════════════════════════
int main(int argc, char* argv[]) {
    std::cout << std::unitbuf;

    int    output_mode = 2;
    int    mode        = 0;   // 0 = EndCover, 1 = Boundary
    int    method      = REFINE_TUBE;   // see refine.h / README S3.2
    int    stepB_type  = 0;             // Lohner's method
    int    tubedegree  = 0;
    int    n           = 0;
    double eps         = 1.0;
    int    order       = 20;
    double T           = 1.0;
    int    debuglevel  = 0;
    std::vector<std::string> SVar, SFun;
    IVector B;

    // ── Argument parsing ────────────────────────────────────────────────
    try {
        if (argc < 11) throw std::runtime_error("too few arguments");

        output_mode = std::stoi(argv[1]);
        mode        = std::stoi(argv[2]);
        method      = std::stoi(argv[3]);
        stepB_type  = std::stoi(argv[4]);
        tubedegree  = std::stoi(argv[5]);
        n           = std::stoi(argv[6]);
        if (n <= 0) throw std::runtime_error("n must be > 0");
        if (mode   < 0 || mode   > 1) throw std::runtime_error("mode must be 0 or 1");
        if (method < 0 || method > 2) throw std::runtime_error("method must be 0, 1 or 2");
        if (stepB_type < 0 || stepB_type > 3) throw std::runtime_error("stepB must be 0..3");

        for (int i = 0; i < n; ++i) SVar.push_back(argv[7 + i]);
        for (int i = 0; i < n; ++i) SFun.push_back(argv[7 + n + i]);

        eps        = std::stod(argv[7 + 2 * n]);
        order      = std::stoi(argv[8 + 2 * n]);
        T          = std::stod(argv[9 + 2 * n]);
        debuglevel = std::stoi(argv[10 + 2 * n]);
        if (eps   <= 0.0) throw std::runtime_error("eps must be > 0");
        if (T     <= 0.0) throw std::runtime_error("T must be > 0");
        if (order <  3  ) throw std::runtime_error("order must be >= 3");

        if (argc < 11 + 3 * n) throw std::runtime_error("missing interval bounds");
        B.resize(n);
        for (int i = 0; i < n; ++i) {
            double lo = std::stod(argv[11 + 2 * n + 2 * i]);
            double hi = std::stod(argv[11 + 2 * n + 2 * i + 1]);
            if (lo > hi) throw std::runtime_error("lo > hi in the initial box");
            B[i] = capd::interval(lo, hi);
        }
    } catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << "\n\n";
        usage();
        return 2;
    }

    if (tubedegree <= 0) tubedegree = std::max(1, order - 1);

    // ── Build the vector field ──────────────────────────────────────────
    for (int i = 0; i < n; ++i) SFun[i] = expand_expression(SFun[i]);
    const std::string mapStr = Convert_to_IMap(SVar, SFun);
    const int degree = order;
    IMap F(mapStr, 3.0);       // degree-3 jet: enough for Jacobians

    const double eps0  = eps;
    const double delta = eps0 / 10.0;
    const bool   chat  = (output_mode >= 3);   // progress chatter

    if (output_mode >= 1) {
        std::cout << "EnCorelib  mode=" << modeName(mode)
                  << "  method=" << methodName(method)
                  << "  stepB=" << stepBName(stepB_type)
                  << "  tubedegree=" << tubedegree
                  << "  order=" << order << "  n=" << n << "\n"
                  << "         f  = (";
        for (int i = 0; i < n; ++i) { std::cout << SFun[i]; if (i+1<n) std::cout << ", "; }
        std::cout << ")\n         B0 = " << B << "   H=" << T << "   eps=" << eps0 << "\n";
    }

    CoverResult Dat;
    Dat.hull = B;

    auto t0 = std::chrono::high_resolution_clock::now();

    // ─────────────────────────────────────────────────────────────────────
    //  mode 0 : EndCover                          (paper, Subsection 3.3)
    // ─────────────────────────────────────────────────────────────────────
    if (mode == 0) {
        MiniScaffold templateG = Init(B, delta, T).G[0];

        auto endEnc = [&](const IVector& box, double epsLocal,
                          const std::vector<double>& p, double Hlocal)
            -> std::pair<IVector, IVector>
        {
            Scaffold S = runScaffold(box, F, templateG,
                                     epsLocal, delta, degree, Hlocal,
                                     method, stepB_type, tubedegree, debuglevel, n, p);
            return {S.E.front(), S.E.back()};
        };

        endcover::EndCoverResult res;
        try {
            res = endcover::EndCover(B, eps0, T, endEnc, 2000000, endcover::COVER_ALL);
        } catch (const std::exception& ex) {
            std::cerr << "[EndCover failed] " << ex.what()
                      << "\nTry a larger eps or a smaller T.\n";
            return 3;
        }
        Dat.C  = res.cover;
        Dat.E0 = res.initBoxes;
    }
    // ─────────────────────────────────────────────────────────────────────
    //  mode 1 : Boundary
    // ─────────────────────────────────────────────────────────────────────
    else {
        Scaffold S = Init(B, delta, T);
        if (n == 2)
            TwoDimEncAlgo(Dat, F, S, eps0, delta, degree, T,
                          method, stepB_type, tubedegree, debuglevel, n, chat);
        else
            ThreeDimEncAlgo(Dat, F, S, eps0, delta, degree, T,
                            method, stepB_type, tubedegree, debuglevel, n, chat);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsedSec =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;

    // ── B1 = Box(C) : the bounding box of the cover ─────────────────────
    if (!Dat.C.empty()) {
        IVector hull = Dat.C.front();
        for (size_t i = 1; i < Dat.C.size(); ++i) hull = hullBox(hull, Dat.C[i]);
        Dat.hull = hull;
    }
    double wmaxC = 0.0;
    for (const auto& b : Dat.C) wmaxC = std::max(wmaxC, wmax(b));

    // ── Reporting (progressive) ─────────────────────────────────────────
    std::cout << "time(ms)=" << static_cast<long long>(elapsedSec * 1000.0 + 0.5) << "\n";

    if (output_mode >= 1) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Hull(T=" << T << ")=" << Dat.hull << "\n"
                  << "B1=" << midRadString(Dat.hull) << "\n"
                  << "wmax(B1)=" << wmax(Dat.hull) << "\n";
    }
    if (output_mode >= 2) {
        std::cout << "#(C)=" << Dat.C.size() << "  #E0=" << Dat.E0.size()
                  << "  max_{B in C} wmax(B)=" << wmaxC << "\n";
        if (mode == 0 && wmaxC > eps0)
            std::cout << "note: C contains boxes wider than eps; these come from"
                         " sub-boxes that were split further (see endcover.h).\n";
    }

    // ── Files ────────────────────────────────────────────────────────────
    if (output_mode >= 3) {
        saveBoxes(Dat.E0, "E0.txt", chat);
        saveBoxes(Dat.C,  "E1.txt", chat);
        if (n == 2) saveConvexHull2D(Dat.C, "convex_hull.txt", chat);
    }

    if (output_mode >= 4) {
        // Plotting extras: images of every E0 box at intermediate times and,
        // in 2D, the trajectories of the four corners of B0.
        std::vector<double> snap = {0.1, 0.4, 0.7};
        snap.erase(std::remove_if(snap.begin(), snap.end(),
                    [&](double s){ return s > T + 1e-14; }), snap.end());

        std::vector<IVector> plotE0 = Dat.E0;
        std::vector<IVector> plotE1 = Dat.C;
        for (const auto& e0 : Dat.E0)
            for (double s : snap)
                plotE1.push_back(stepBlohner(F, degree, s, e0));

        if (n == 2) {
            std::vector<double> ct = snap;
            ct.push_back(T);
            for (int i = 0; i < 4; ++i) {
                IVector corner(2);
                corner[0] = (i & 1) ? interval(B[0].rightBound()) : interval(B[0].leftBound());
                corner[1] = (i & 2) ? interval(B[1].rightBound()) : interval(B[1].leftBound());
                plotE0.push_back(corner);
                for (double s : ct) plotE1.push_back(stepBlohner(F, degree, s, corner));
            }
            saveConvexHull2D(Dat.E0, "convex_hull_E0.txt", chat);
        }
        saveBoxes(plotE0, "E_0.txt", chat);
        saveBoxes(plotE1, "E_1.txt", chat);
    }

    appendLog(mode == 0 ? "EndCover" : "Boundary", mode, method, stepB_type,
              eps0, order, T, n, B, Dat.hull, Dat.C.size(), elapsedSec);
    return 0;
}
