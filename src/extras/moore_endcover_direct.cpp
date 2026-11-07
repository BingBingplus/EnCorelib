// moore_endcover_direct.cpp  --  standalone diagnostic (not part of EnCorelib)
//
// EndCover on Moore's rotation example x' = y, y' = -x with stepB = 3, the
// StepB variant with the higher-order a-priori enclosure and without the
// ball of Lemma 2.1.  It runs the same scaffold driver as EnCorelib
// (runScaffold: Extend to the horizon, then Refine) and writes the cover
// under a label, so that runs for several end times T can be compared.
//
// Usage: moore_endcover_direct <T> <eps> <degree> <label>
// Output: E0_<label>.txt, E1_<label>.txt (the "Box n: [..] x [..]" format
// of EnCorelib; line k of E0 is the sub-box whose end-enclosure is line k
// of E1).

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include "../boundary.h"
#include "../endcover.h"

static void saveBoxes(const std::vector<IVector>& boxes, const std::string& fname) {
    std::ofstream fout(fname);
    fout << std::fixed << std::setprecision(17);
    for (size_t i = 0; i < boxes.size(); ++i) {
        fout << "Box " << i << ": ";
        for (int d = 0; d < boxes[i].dimension(); ++d) {
            fout << "[" << boxes[i][d].leftBound() << ", " << boxes[i][d].rightBound() << "]";
            if (d + 1 < boxes[i].dimension()) fout << " x ";
        }
        fout << "\n";
    }
}

static Scaffold makeInitialStage(const IVector& B, double delta, double T_hint) {
    MiniScaffold G;
    G.mu   = {1000.0};
    G.delta = delta;
    G.hTube = T_hint;
    G.ell   = 0;
    G.E   = {B};
    G.F   = {B};
    G.nTube = 1;
    Scaffold S = {{0.0}, {B}, {B}, {G}};
    return S;
}

int main(int argc, char** argv) {
    double T      = (argc > 1) ? std::atof(argv[1]) : 0.0;
    double veps   = (argc > 2) ? std::atof(argv[2]) : 0.02;
    int    degree = (argc > 3) ? std::atoi(argv[3]) : 20;
    std::string label = (argc > 4) ? argv[4] : "out";

    const int stepB_type  = 3;
    const int debuglevel  = 0;
    const int n = 2;

    IMap F("var:x,y;fun:y,-x;", 3.0);
    IVector B(2);
    B[0] = interval(0.95, 1.05);
    B[1] = interval(-0.05, 0.05);

    const double delta = veps / 10.0;
    MiniScaffold G0 = makeInitialStage(B, delta, T).G[0];

    auto endEnc = [&](const IVector& box, double epsLocal,
                      const std::vector<double>& target, double Hlocal)
        -> std::pair<IVector, IVector>
    {
        Scaffold S = runScaffold(box, F, G0, epsLocal, delta, degree, Hlocal,
                                 REFINE_TUBE, stepB_type, /*tubedegree*/ 1,
                                 debuglevel, n, target);
        return {S.E.front(), S.E.back()};
    };

    auto t0 = std::chrono::high_resolution_clock::now();
    endcover::EndCoverResult coverAtT;
    try {
        coverAtT = endcover::EndCover(B, veps, T, endEnc, 4000000);
    } catch (const std::exception& ex) {
        std::cerr << "[EndCover failed] T=" << T << " : " << ex.what() << "\n";
        return 3;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsedSec =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0;

    IVector unionBox = coverAtT.cover.empty() ? B : coverAtT.cover[0];
    for (size_t i = 1; i < coverAtT.cover.size(); ++i)
        unionBox = hullBox(unionBox, coverAtT.cover[i]);

    std::cout << "T=" << T << "  eps=" << veps
              << "  boxes=" << coverAtT.cover.size()
              << "  wmax(union)=" << wmax(unionBox)
              << "  time(s)=" << elapsedSec << "\n";

    saveBoxes(coverAtT.cover,     "E1_" + label + ".txt");
    saveBoxes(coverAtT.initBoxes, "E0_" + label + ".txt");
    return 0;
}
