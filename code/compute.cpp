#include "compute.hpp"
#include "simulation.hpp"
#include "parameter.hpp"
#include <algorithm>
#include <cmath>

Compute::Compute(Simulation *simu) : Parameter(simu) { }
Compute::~Compute() { }

double Compute::gyrationRadius(double** r) 
{
    double rg = 0;
    double rcm[DIM];
    std::fill(&rcm[0], &rcm[0] + DIM, 0);
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            rcm[j] += r[i][j];
        }
    }

    for (int i = 0; i < DIM; ++i) {
        rcm[i] = rcm[i] / nBead;
    }

    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            rg += (r[i][j] - rcm[j]) * 
                (r[i][j] - rcm[j]);
        }
    }

    return sqrt(rg / nBead);
}
