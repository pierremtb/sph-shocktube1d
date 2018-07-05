#include <vector>
#include <cmath>
using namespace std;

void computeL(int particlesCount, vector<double>& L, vector<double>& density, vector<double>& pressure, vector<double>& dw, vector<int>& pairI, vector<int>& pairJ, int niac, double mass) {
    
    // Zero-initialization, since dwdx(0) = 0
    for (int i = 0; i < particlesCount; i++) {
        L[i] = 0;
    }

    // Computing with the neighbours
    for (int k = 0; k < niac; k++) {

        // Getting pairs
        int i = pairI[k];
        int j = pairJ[k];

        // Calculating constant terms
        double Pij = (pressure[i]/(density[i]*density[i]) + pressure[j]/(density[j]*density[j]));
        
        // Wrapping up for both particles, remembering dwdx_ij = -dwdx_ji
        L[i] += dw[k] * Pij * mass;
        L[j] -= dw[k] * Pij * mass;
    }
}
