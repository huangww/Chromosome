#include "compute.hpp"
#include "constant.hpp"
#include <algorithm>
#include <cmath>

Compute::Compute()  { }
Compute::~Compute() { }

double Compute::gyrationRadius(int N, double** r) 
{
    
    double rg = 0;
    double rcm[DIM];
    std::fill(&rcm[0], &rcm[0] + DIM, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < DIM; ++j) {
            rcm[j] += r[i][j];
        }
    }

    for (int i = 0; i < DIM; ++i) {
        rcm[i] = rcm[i] / N;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < DIM; ++j) {
            rg += (r[i][j] - rcm[j]) * 
                (r[i][j] - rcm[j]);
        }
    }

    return sqrt(rg / N);
}

double Compute::gyrationRadius(int N, double* r) 
{
    double rcm = 0;
    for (int i = 0; i < N; ++i) {
        rcm += r[i];
    }
    rcm = rcm / N;

    double rg = 0;
    for (int i = 0; i < N; ++i) {
        rg += (r[i] - rcm) * (r[i] - rcm);
    }

    return sqrt(rg / N);
}

double Compute::theta(double *r1, double *r2)
{
    double rs1=0, rs2=0, rdot12=0;
    for (int i = 0; i < DIM; ++i) {
       rs1 += r1[i]*r2[i]; 
       rs2 += r2[i]*r2[i]; 
       rdot12 += r1[i]*r2[i];
    }

    return acos(rdot12 / sqrt(rs1*rs2));
}

