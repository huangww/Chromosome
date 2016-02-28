#include "potential.hpp"
#include "input.hpp"
#include "constant.hpp"
#include "bead.hpp"
#include <cmath>

Potential::Potential()  { }
Potential::~Potential() { }


double Potential::LennardJones(int N, double **r) 
{	
    double plj = 0.0;
    double r0 = 0.75;
    double eps = 1.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
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
