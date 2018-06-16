#include <vector>
#include <cmath>
using namespace std;

void computeDensity(int particlesCount, vector<double>& density, double mass, vector<double>& w, vector<int>& pairI, vector<int>& pairJ, int niac, double w0, double h) {
    
    // Considering the self-effect
    for (int i = 0; i < particlesCount; i++) {
        density[i] = mass * w0;
    }

    // Computing with the neighbours
    for (int k = 0; k < niac; k++) {
        density[pairI[k]] += mass * w[k];
        density[pairJ[k]] += mass * w[k];
    }
}