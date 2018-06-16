#include <vector>
#include <cmath>
#include "kernels.hpp"
using namespace std;

void doDirectFind(int particlesCount, int& niac, vector<int>& pairI, vector<int>& pairJ, double& w0, vector<double>& position, vector<double>& w, vector<double>& dw, double h) {

        for (int i = 0; i < particlesCount - 1; i++) {
            for (int j = i + 1; j < particlesCount; j++) {

                // Checking smoothing length
                double r = position[i] - position[j];

                if (fabs(r) < 2 * h) {
                    
                    // Increasing the number of interacting pairs
                    niac++;

                    // Storing the 1st and 2nd particles index of the pair
                    pairI.push_back(i);
                    pairJ.push_back(j);

                    // Computing smoothing kernel and derivative
                    kernels::Result1D kernel = kernels::splineCubic(r, h);
                    w.push_back(kernel.w);
                    dw.push_back(kernel.dw);
                }
            }
        }
                    
        // Getting the w0 term
        w0 = kernels::splineCubic(0, h).w;
}