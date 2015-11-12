#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "utilities.h"
/* #include <omp.h> */

void CalculateAij(double b[numLink][DIM], 
        double u[numLink][DIM],
        int g[numLink][numLink], 
        double A[numLink*numLink])
{
    memset(A, 0, sizeof(A[0])*numLink*numLink);
    for (int i = 0; i < numLink; i++) 
    {
        for (int j = 0; j < numLink; j++) 
        {
            /* note the order of A indices, C or F, A[j*rodnumber+i] = A[i][j] */
            A[j*numLink+i] = g[i][j]*ddot(&b[i][0],&u[j][0]);
        }
    }
}

void CalculateB(double b[numLink][DIM],
        double u[numLink][DIM],
        int g[numLink][numLink], 
        double x[numLink],
        double B[numLink])
{
    for (int i = 0; i < numLink; i++)
    {
        double temp[DIM];
        memset(temp, 0, sizeof(temp[0])*DIM);
        for (int j = 0; j < numLink; j++)
        {
            if (g[i][j] != 0)
            {
                for (int k = 0; k < DIM; k++)
                {
                    temp[k] += g[i][j]*x[j]*u[j][k];
                }
            }
        }
        double dotBiBi = 0.0;
        double dotTiTi = 0.0;
        for (int k = 0; k < DIM; k++)
        {
            dotBiBi += b[i][k]*b[i][k];
            dotTiTi += temp[k]*temp[k];
        }
        B[i] = (1.0 - dotBiBi)/(2.0*dt) -
            dt * dotTiTi/2.0;
    }
}



extern int dgetrf_(int* m, int* n, double* A, 
        int* lda, int* iPIv, int* info);

extern int dgetrs_(char* s, int* n, int* nrhs, 
        double* A, int* lda, int* iPIv, 
        double* B, int* ldb, int* info);

void Picard(double b[numLink][DIM], 
        double u[numLink][DIM],
        int g[numLink][numLink],
        double x[numLink])
{
    double A[numLink*numLink];
    CalculateAij(b, u, g, A);

    int n = numLink;
    int lda = numLink;
    int iPIv[numLink];
    int info;
    dgetrf_(&n, &n, A, &lda, iPIv, &info);

    double xold[numLink];
    for (int step = 0; step < 200; step++)
    {
        memcpy(xold, x, sizeof(xold));
        double B[numLink];
        CalculateB(b, u, g, x, B);

        char s = 'N';
        int nrhs = 1;
        dgetrs_(&s, &n, &nrhs, A, &n, iPIv, B, &n, &info);

        memcpy(x, B, sizeof(B));
        double maxDiff = fabs(xold[0] - x[0]);
        for (int i = 1; i < numLink; i++)
        {
            if (fabs(xold[i] - x[i]) > maxDiff)
            {
                maxDiff = fabs(xold[i] - x[i]);
            }
        }
        if (maxDiff < 1e-6) return;
    }

    printf("MaxStep exceeded in Picard\n");

}
