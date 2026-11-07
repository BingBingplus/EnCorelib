// moore_direct_test.cpp  --  standalone diagnostic (not part of EnCorelib)
//
// Moore's rotation example x' = y, y' = -x, propagated with one StepB
// variant at a time from a single box, with the enclosure put back into an
// axis-aligned box after every step.  This shows the wrapping effect that
// EnCorelib avoids by carrying one doubleton along a scaffold
// (SolutionChain in stepAB.h).
//
// This program repeatedly applies ONE step method to a SINGLE box (no
// EndCover subdivision) over many fixed-size steps, and records the
// enclosure at each step for both:
//   - stepBlohnerHO   (stepB=3 : plain Taylor + interval remainder,
//                          re-enclosed in an axis-aligned box every step
//                          -- exactly the classical "naive interval
//                          method" that suffers wrapping effect)
//   - stepBlohner       (stepB=0 : CR-Lohner / doubleton, for reference)
//
// Output: moore_direct_runs/direct_wmax.txt, moore_direct_runs/lohner_wmax.txt
//   columns: t  wmax  xlo xhi ylo yhi

#include "../types.h"
#include "../stepAB.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

static void dumpRow(std::ofstream& f, double t, const IVector& B) {
    f << std::setprecision(17)
      << t << " " << wmax(B) << " "
      << B[0].leftBound() << " " << B[0].rightBound() << " "
      << B[1].leftBound() << " " << B[1].rightBound() << "\n";
}

int main(int argc, char** argv) {
    int nPeriods       = (argc > 1) ? std::atoi(argv[1]) : 3;
    int stepsPerPeriod = (argc > 2) ? std::atoi(argv[2]) : 64;
    int degree         = (argc > 3) ? std::atoi(argv[3]) : 20;

    const double pi = std::acos(-1.0);
    const double period = 2.0 * pi;
    const double h       = period / stepsPerPeriod;
    const int totalSteps = nPeriods * stepsPerPeriod;

    IMap F("var:x,y;fun:y,-x;", 3.0);

    IVector B0(2);
    B0[0] = interval(0.95, 1.05);
    B0[1] = interval(-0.05, 0.05);

    std::ofstream outDirect("moore_direct_runs/direct_wmax.txt");
    std::ofstream outLohner("moore_direct_runs/lohner_wmax.txt");

    IVector Bd = B0, Bl = B0;
    dumpRow(outDirect, 0.0, Bd);
    dumpRow(outLohner, 0.0, Bl);

    std::cout << "h=" << h << "  totalSteps=" << totalSteps
              << "  degree=" << degree << "\n";

    bool directAlive = true;
    for (int s = 1; s <= totalSteps; ++s) {
        double t = s * h;

        if (directAlive) {
            try {
                Bd = stepBlohnerHO(F, degree, h, Bd);
                if (!std::isfinite(Bd[0].leftBound()) || !std::isfinite(Bd[0].rightBound()) ||
                    !std::isfinite(Bd[1].leftBound()) || !std::isfinite(Bd[1].rightBound()) ||
                    wmax(Bd) > 1e6) {
                    std::cerr << "[direct] blew up / non-finite at step " << s
                              << " (t=" << t << "), stopping direct-method updates.\n";
                    directAlive = false;
                }
            } catch (const std::exception& ex) {
                std::cerr << "[direct] exception at step " << s << ": " << ex.what()
                          << " -- stopping direct-method updates.\n";
                directAlive = false;
            }
        }
        dumpRow(outDirect, t, Bd);

        Bl = stepBlohner(F, degree, h, Bl);
        dumpRow(outLohner, t, Bl);

        if (s % stepsPerPeriod == 0) {
            std::cout << "t=" << t << "  wmax(direct)=" << wmax(Bd)
                      << "  wmax(lohner)=" << wmax(Bl) << "\n";
        }
    }
    return 0;
}
