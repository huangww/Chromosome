#include "force.hpp"
#include "input.hpp"
#include "random.hpp"
#include "particle.hpp"
#include "bead.hpp"
#include "constant.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

Force::Force()  { }
Force::~Force() { }

void Force::setParameter(Input *input) 
{
    nSite = int(input->parameter["nSite"]);
    nPar = int(input->parameter["nPar"]);

    nBead = int(input->parameter["nBead"]);
    tempEff = input->parameter["tempEff"];
    dt = input->parameter["dt"];
}

void Force::print(double* f) 
{
    for (int i = 0; i < nPar; ++i) {
        std::cout << i << '\t' << f[i] << std::endl;
    }
}
void Force::print(double** f) 
{
    for (int i = 0; i < nBead; ++i) {
        std::cout << i << '\t';
        for (int j = 0; j < DIM; ++j) {
            std::cout << f[i][j] << '\t';
        }
        std::cout << std::endl;
    }
}


double* Force::repulsive(double *x, double* f)
{
    std::fill(&f[0], &f[0] + nPar, 0);

    double r0 = 1.0;
    double eps = 1.0;
    for (int i = 1; i < nPar; ++i) {
        int j = i - 1;
        double d0 = fabs(x[i] - x[j]);
        if (pow(d0, 6) <= 2*pow(r0, 6)) {
            double d6 = pow(r0/d0, 6);
            f[i] = f[i] + 48 * eps * (d6 - 0.5)
                * d6 * (x[i] - x[j]) / (d0 * d0);
            f[j] = f[j] - 48 * eps * (d6 - 0.5)
                * d6 * (x[i] - x[j]) / (d0 * d0);
        }

    }

    return f;
}

double** Force::repulsive(double** r, double** f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

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
                for (int k = 0; k < DIM; ++k) {
                    f[i][k] = f[i][k] + 48 * eps *
                        (r6 - 0.5) * r6 *
                        (r[i][k]-r[j][k])/rsd;
                    f[j][k] = f[j][k] - 48 * eps *
                        (r6 - 0.5) * r6 *
                        (r[i][k]-r[j][k])/rsd;
                }

            }
        }

    }

    return f;
}



double* Force::brownian(double* f)
{
    double temp = 1.0;
    for (int i = 0; i < nPar; i++) {
        f[i] = sqrt(2.0*temp/dt) * GaussRan(seed);
    }
    return f;
}

double** Force::brownian(double** f)
{
    double temp = 1.0;
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            f[i][j] = sqrt(2.0*temp/dt) * GaussRan(seed);
        }
    }
    return f;
}


double* Force::external(double* f)
{
    for (int i = 0; i < nPar; ++i) {
        f[i] = - 1.0 / tempEff;
    }

    return f;
}

double** Force::external(double** f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            f[i][0] = 1.0 / tempEff;
        }
    }
    return f;
}


double* Force::boundary(double*x, double* f)
{
    std::fill(&f[0], &f[0] + nPar, 0);

    double r0 = 0.5;
    double eps = 1.0;

    double d0 = fabs(x[0]);
    if (pow(d0, 6) <= 2*pow(r0, 6)) {
        double d6 = pow(r0/d0, 6);
        f[0] = 48 * eps * (d6 - 0.5)
            * d6 * fabs(x[0] - 0) / (d0 * d0);
    }
    d0 = fabs(x[nPar-1] - nSite);
    if (pow(d0, 6) <= 2*pow(r0, 6)) {
        double d6 = pow(r0/d0, 6);
        f[nPar-1] =  - 48 * eps * (d6 - 0.5)
            * d6 * fabs(x[nPar-1] - nSite) / (d0 * d0);
    }

    return f;
}
