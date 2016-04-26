#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "solver.h"
#include "random.h"
#include "utilities.h"
#include "output.h"


void CalculateFc(double r[numUnit][DIM],
        double rs[numUnit][DIM],
        int link[numLink][2],
        int g[numLink][numLink],
        double fc[numUnit][DIM])
{
    double u[numLink][DIM];
    CalculateVectorU(r, link, u);

    double b[numLink][DIM];
    CalculateVectorB(rs, link, b);

    double tension[numLink];
    memset(tension, 0, sizeof(tension[0])*numLink);
    Picard(b, u, g, tension);

    memset(fc, 0, sizeof(fc[0][0]) * numUnit * DIM);
    for (int i = 0; i < numUnit; i++) {
        for (int j = 0; j < numLink; j++) {
            if (link[j][0] == i) {
                for (int k = 0; k < DIM; k++) {
                    fc[i][k] = fc[i][k] + tension[j] * u[j][k];
                }
            }
            if (link[j][1] == i) {
                for (int k = 0; k < DIM; k++) {
                    fc[i][k] = fc[i][k] - tension[j] * u[j][k];
                }
            }
        }
    }

}

void GenerateFb(double fb[numUnit][DIM],
        unsigned long seed)
{
    double temp = 1;
    for (int i = 0; i < numUnit; i++) {
        for (int j = 0; j < DIM; j++) {
            fb[i][j] = sqrt(2.0*temp/dt) * GaussRan(seed);
        }
    }

}

double LennardJonesPotential(double r[numUnit][DIM])
{	
    double plj = 0.0;
    double r0 = 0.75;
    double eps = 1.0;
    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = i+1; j < numUnit; ++j)
        {
            double rsd = 0;
            for (int k = 0; k < DIM; ++k)
            {
                rsd = rsd + (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
            }
            if (pow(rsd, 3) <= 2*pow(r0,6))
            {
                double r6;
                r6 = pow(r0*r0/rsd, 3);
                plj += 4*eps*(r6*r6 - r6) + eps;
            }
        }

    }
    return plj;
}

void LennardJonesForce(double r[numUnit][DIM],
        double f[numUnit][DIM])
{
    memset(f, 0, sizeof(f[0][0])*numUnit*DIM);
    double r0 = 0.75;
    double eps = 1.0;
    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = i+1; j < numUnit; ++j)
        {
            double rsd = 0;
            for (int k = 0; k < DIM; ++k)
            {
                rsd = rsd + (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
            }
            if (pow(rsd, 3) <= 2*pow(r0,6))
            {
                double r6;
                r6 = pow(r0*r0/rsd, 3);
                for (int k = 0; k < DIM; ++k)
                {
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

}


extern int dgetrf_(int* m, int* n, double* A, 
        int* lda, int* iPIv, int* info);

extern int dgetri_(int* m, double* a, int* lda,
        int* iPIv, double* work, int* lwork, int* info);


void PseudoForce(double r[numUnit][DIM],
        int link[numLink][2],
        int g[numLink][numLink],
        double f[numUnit][DIM])
{
    memset(f, 0, sizeof(f[0][0])*numUnit*DIM);
    double u[numLink][DIM];
    CalculateVectorU(r, link, u);
    double metric[numLink*numLink];
    CalculateAij(u, u, g, metric); 
    for (int i = 0; i < numLink*numLink; ++i)
    {
        metric[i] = -metric[i];
    }

    /* calculate the inverse of the metric matrix */
    int n = numLink;
    int iPIv[numLink];
    int info;
    double work[numLink];
    dgetrf_(&n, &n, metric, &n, iPIv, &info);
    dgetri_(&n, metric, &n, iPIv, work, &n, &info);

    for (int k = 0; k < numUnit; ++k)
    {
        for (int i = 0; i < numLink; ++i)
        {
            for (int j = i+1; j < numLink; ++j)
            {
                if (abs(g[i][j]) == 1)
                {
                    int i0,i1,j0,j1;
                    i0 = link[i][0];
                    i1 = link[i][1];
                    j0 = link[j][0];
                    j1 = link[j][1];
                    if ((k-i0)*(k-i1)*(k-j0)*(k-j1)==0)
                    {
                        double uij;
                        uij = ddot(&u[i][0], &u[j][0]); 
                        double pgr[DIM];
                        for (int m = 0; m < DIM; ++m)
                        {
                            pgr[m] = (delta(k,i1) - delta(k,i0)) *
                                (u[j][m] - uij * u[i][m]) +
                                (delta(k,j1) - delta(k, j0)) *
                                (u[i][m] - uij * u[j][m]);			
                            f[k][m] = f[k][m] + g[j][i] * metric[i*numLink+j] * pgr[m];
                        }
                    }


                }
            }
        }
    }

}
