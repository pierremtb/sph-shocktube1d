#include <vector>
#include <cmath>

#include "density.hpp"
#include "direct-find.hpp"
#include "eos.hpp"
#include "linear-momentum.hpp"
#include "energy-equation.hpp"

using namespace std;

void doSingleStep(
    int particlesCount,
    vector<double>& position,
    vector<double>& velocity,
    vector<double>& density,
    vector<double>& energy,
    vector<double>& pressure,
    vector<double>& L,
    vector<double>& H,
    double mass,
    double gamma,
    double h) {

        // Preparing
        vector<int> pairI, pairJ;
        vector<double> w, dw;
        int niac = 0;
        double w0 = 0;

        // 1
        doDirectFind(particlesCount, niac, pairI, pairJ, w0, position, w, dw, h);

        // 2
        computeDensity(particlesCount, density, mass, w, pairI, pairJ, niac, w0, h);

        // 3
        computePressure(particlesCount, density, pressure, energy, gamma);

        // 4
        computeL(particlesCount, L, density, pressure, dw, pairI, pairJ, niac, mass);

        // 5
        computeH(particlesCount, H, density, velocity, pressure, dw, pairI, pairJ, niac, mass);

    }

