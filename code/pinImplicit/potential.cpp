#include "potential.hpp"
#include "parameter.hpp"
#include "simulation.hpp"
#include "bead.hpp"
#include <cmath>

Potential::Potential(Simulation *simu) : Parameter(simu) { }
Potential::~Potential() { }

double Potential::LennardJones() 
{	
    double **r = bead->r;
    double plj = 0.0;
    double r0 = 0.75;
    double eps = 1.0;
    for (int i = 0; i < nBead; ++i) {
        for (int j = i+1; j < nBead; ++j) {
            double rsd = 0;
            for (int k = 0; k < DIM; ++k) {
                rsd = rsd + (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
            }
            if (pow(rsd, 3) <= 2*pow(r0,6)) {
                double r6;
                r6 = pow(r0*r0/rsd, 3);
                plj += 4*eps*(r6*r6 - r6) + eps;
            }
        }
    }
    return plj;
}
