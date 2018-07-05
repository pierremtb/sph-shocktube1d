// g++ -o run src/main.cpp -lm --std=c++11

#include <iostream>
#include <vector>
#include <cmath>

#include "sph/new-single-step.hpp"
#include "utils/output.hpp"

using namespace std;

const int LEFT_PARTICLES_COUNT = 320;
const int RIGHT_PARTICLES_COUNT = 80;
const int PARTICLES_COUNT = LEFT_PARTICLES_COUNT + RIGHT_PARTICLES_COUNT;

const double RIGHT_MIN = 0;
const double RIGHT_MAX = 0.6;
const double LEFT_MIN = -0.6;
const double LEFT_MAX = RIGHT_MIN - (RIGHT_MIN - LEFT_MIN) / (LEFT_PARTICLES_COUNT + 2);

const double LEFT_DENSITY = 1.0;
const double LEFT_VELOCITY= 0.0;
const double LEFT_ENERGY = 2.5;
const double LEFT_PRESSURE = 1.0;

const double RIGHT_DENSITY = 0.25;
const double RIGHT_VELOCITY= 0.0;
const double RIGHT_ENERGY = 1.795;
const double RIGHT_PRESSURE = 0.1795;

const double PARTICLE_MASS = 0.001875;
const double GAMMA = 1.4;

const double DX = 0.6/80.0;
const double HSML = 2.0 * DX;
//const double KAPPA = 2;

const unsigned int TIME_STEPS = 40;
const double TIME_DELTA = 0.005;

int main()
{
    // TESTING KERNEL
    
    // cout << "x,w,dwdx" <<endl;
    // for (int i = 0; i<PARTICLES_COUNT; i++) {
    //     double dx = 4*HSML/PARTICLES_COUNT;
    //     double r = -2*HSML+i*dx;
    //     kernels::Result1D res = kernels::splineCubic(r, HSML);
    //     cout << r << "," << res.w << "," << res.dw << endl;
    // }
    // return 0;


    // VARIABLES 

    int iteration = 0;
    double time = 0;

    vector<double> position(PARTICLES_COUNT);
    vector<double> velocity(PARTICLES_COUNT);
    vector<double> density(PARTICLES_COUNT);
    vector<double> energy(PARTICLES_COUNT);
    vector<double> pressure(PARTICLES_COUNT);

    vector<double> L(PARTICLES_COUNT);
    vector<double> H(PARTICLES_COUNT);


    // INITIALIZATION

    for (int i = 0; i < PARTICLES_COUNT; i++) {
        if (i < 320) {
            position[i] = i * (DX/4.0) - 0.6 + DX/4.0;
            density[i] = LEFT_DENSITY;
            velocity[i] = LEFT_VELOCITY;
            energy[i] = LEFT_ENERGY;
            pressure[i] = LEFT_PRESSURE;
        } else {
            position[i] = (i-320) * DX + position[319] + 0.5 * DX;
            density[i] = RIGHT_DENSITY;
            velocity[i] = RIGHT_VELOCITY;
            energy[i] = RIGHT_ENERGY;
            pressure[i] = RIGHT_PRESSURE;
        }

        L[i] = 0.0;
        H[i] = 0.0;
    }

    output(PARTICLES_COUNT, position, velocity, density, energy, pressure, iteration);
    
    iteration++;
    time += TIME_DELTA;


    // SIMULATION STEPS

    int N_block = 20; // fixed particles.

    while (iteration <= TIME_STEPS) {
        vector<double> velocityMin(PARTICLES_COUNT);
        vector<double> energyMin(PARTICLES_COUNT);

        // First LF integration
        if (iteration > 1) {
            //velocityMin = velocity;
            //energyMin = energy;

            for (int i = N_block; i < PARTICLES_COUNT - N_block; i++) {
                velocityMin[i] = velocity[i];
                energyMin[i]=energy[i];
                velocity[i] += TIME_DELTA * L[i]/2.0;
                energy[i] += TIME_DELTA * H[i]/2.0;
            }
        }

        // SPH Single Step
        doSingleStep(PARTICLES_COUNT, position, velocity, density, energy, pressure, L, H, PARTICLE_MASS, GAMMA, HSML);

        // Second LF integration
        for (int i = N_block; i < PARTICLES_COUNT - N_block; i++) {
            if (iteration == 1) {
                velocity[i] += TIME_DELTA * L[i]/2.0;
                energy[i] += TIME_DELTA * H[i]/2.0;

                position[i] += TIME_DELTA * velocity[i];
            } else {
                velocity[i] = velocityMin[i] + TIME_DELTA * L[i];
                energy[i] = energyMin[i] + TIME_DELTA * H[i];

                position[i] += TIME_DELTA * velocity[i];
            }
        }

        // Iteration data export
        output(PARTICLES_COUNT, position, velocity, density, energy, pressure, iteration);
        cout << time << endl;

        // Next iteration
        time += TIME_DELTA;
        iteration++;
    }

    return 0;
    
}
