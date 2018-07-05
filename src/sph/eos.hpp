#include <vector>
#include <cmath>
using namespace std;

void computePressure(int particlesCount, vector<double>& density, vector<double>& pressure, vector<double>& energy, double gamma) {
    for (int i = 0; i < particlesCount; i++) {

        // Equation of state p = rho.u.(gamma - 1)
        pressure[i] = density[i] * energy[i] * (gamma - 1.0);
    }
}
