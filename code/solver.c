#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void CalculateAij(double b[rodNumber][dimension], 
		double u[rodNumber][dimension],
		int g[rodNumber][rodNumber], 
		double A[rodNumber*rodNumber])
{
	for (int i = 0; i < rodNumber; i++) 
	{
		for (int j = 0; j < rodNumber; j++) 
		{
			double dotBiUj = 0.0;
			for (int k = 0; k < dimension; k++)
		       	{
				dotBiUj = dotBiUj + b[i][k]*u[j][k];
			}
			/*note the order of A indices, C or F*/
			A[j*rodNumber+i] = g[i][j]*dotBiUj;
		}
	}
}

void CalculateB(double b[rodNumber][dimension],
		double u[rodNumber][dimension],
		int g[rodNumber][rodNumber], 
		double x[rodNumber],
		double B[rodNumber])
{
	for (int i = 0; i < rodNumber; i++)
       	{
		double temp[dimension] = {0};
		for (int j = 0; j < rodNumber; j++)
	       	{
			for (int k = 0; k < dimension; k++)
		       	{
				temp[k] = temp[k] + g[i][j]*x[j]*u[j][k];
			}
		}
		double dotBiBi = 0.0;
		double dotTiTi = 0.0;
		for (int k = 0; k < dimension; k++)
	       	{
			dotBiBi = dotBiBi + b[i][k]*b[i][k];
			dotTiTi = dotTiTi + temp[k]*temp[k];
		}
		B[i] = (1.0 - dotBiBi)/(2.0*dt) -
			dt * dotTiTi/2.0;
	}
}



extern int dgetrf_(int* m, int* n, double* A, 
		int* lda, int* ipiv, int* info);

extern int dgetrs_(char* s, int* n, int* nrhs, 
		double* A, int* lda, int* ipiv, 
		double* B, int* ldb, int* info);

void Picard(double b[rodNumber][dimension], 
		double u[rodNumber][dimension],
		int g[rodNumber][rodNumber],
		double x[rodNumber])
{
	double A[rodNumber*rodNumber] = {0};
	CalculateAij(b, u, g, A);

	int n = rodNumber;
	int lda = rodNumber;
	int ipiv[rodNumber];
	int info;
	dgetrf_(&n, &n, A, &lda, ipiv, &info);
	
	double xold[rodNumber];
	for (int step = 0; step < maxStep; step++)
       	{
		memcpy(xold, x, sizeof(xold));
		double B[rodNumber] = {0};
		CalculateB(b, u, g, x, B);

		char s = 'N';
		int nrhs = 1;
		dgetrs_(&s, &n, &nrhs, A, &n, ipiv, B, &n, &info);

		memcpy(x, B, sizeof(B));
		double maxDiff = fabs(xold[0] - x[0]);
		for (int i = 1; i < rodNumber; i++)
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
