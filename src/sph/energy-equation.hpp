#include <vector>
#include <cmath>
using namespace std;

void computeH(int particlesCount, vector<double>& H, vector<double>& density, vector<double> & velocity, vector<double>& pressure, vector<double>& dw, vector<int>& pairI, vector<int>& pairJ, int niac, double mass) {
    
    // Zero-initialization, since dwdx(0) = 0
    for (int i = 0; i < particlesCount; i++) {
        H[i] = 0;
    }

    // Computing with the neighbours
    for (int k = 0; k < niac; k++) {

        // Getting pairs
        int i = pairI[k];
        int j = pairJ[k];

        // Calculating constant terms
        double PVij = -0.5 * (pressure[i]/(density[i]*density[i]) + pressure[j]/(density[j]*density[j])) * (velocity[i] - velocity[j]);
        
        // Wrapping up for both particles, remembering dwdx_ij = -dwdx_ji
        H[i] += dw[k] * PVij * mass;
        H[j] += dw[k] * PVij * mass;
    }
}
