#include <math.h>
#include <string.h>
#include "solver.h"
#include "random.h"

void CalculateVectorU(double r[beadNumber][dimension], 
		int link[rodNumber][2], 
		double u[rodNumber][dimension])
{
	/*Calculate rod unit vector u*/
	for (int i = 0; i < rodNumber; i++) {
		double uLength = 0;
		for (int j = 0; j < dimension; j++) {
			u[i][j] = r[link[i][1]][j] - r[link[i][0]][j]; 
			uLength = uLength + u[i][j]*u[i][j];
		}
		uLength = sqrt(uLength);
		for (int j = 0; j < dimension; j++) {
			u[i][j] = u[i][j]/uLength;
		}
	}
}

void CalculateVectorB(double rs[beadNumber][dimension], 
		int link[rodNumber][2], 
		double b[rodNumber][dimension])
{
	/*Calculate rod unit vector B*/
	for (int i = 0; i < rodNumber; i++) {
		for (int j = 0; j < dimension; j++) {
			b[i][j] = rs[link[i][1]][j] - rs[link[i][0]][j]; 
		}
	}
}

void CalculateFc(double r[beadNumber][dimension],
		double rs[beadNumber][dimension],
		int link[rodNumber][2],
		int g[rodNumber][rodNumber],
		double fc[beadNumber][dimension])
{
	double u[rodNumber][dimension];
	CalculateVectorU(r, link, u);

	double b[rodNumber][dimension];
	CalculateVectorB(rs, link, b);

	double tension[rodNumber]= {0};
	Picard(b, u, g, tension);
	memset(fc, 0, sizeof(fc[0][0]) * beadNumber * dimension);
	for (int i = 0; i < beadNumber; i++) {
		for (int j = 0; j < rodNumber; j++) {
			if (link[j][0] == i) {
				for (int k = 0; k < dimension; k++) {
					fc[i][k] = fc[i][k] + tension[j] * u[j][k];
				}
			}
			if (link[j][1] == i) {
				for (int k = 0; k < dimension; k++) {
					fc[i][k] = fc[i][k] - tension[j] * u[j][k];
				}
			}
		}
	}

	
}

void GenerateFb(double fb[beadNumber][dimension],
		long* seed)
{
	for (int i = 0; i < beadNumber; i++) {
		for (int j = 0; j < dimension; j++) {
			fb[i][j] = sqrt(2.0/dt) * gasdev(seed);
		}
	}
	
}

void LennardJones(double r[beadNumber][dimension],
		double f[beadNumber][dimension])
{
	memset(f, 0, sizeof(f[0][0])*beadNumber*dimension);
	double r0 = 0.75;
	double eps = 6.0;
	for (int i = 0; i < beadNumber; ++i)
	{
		for (int j = 0; j < beadNumber; ++j)
		{
			if (i != j)
			{
				double rd = 0;
				for (int k = 0; k < dimension; ++k)
				{
					rd = rd + (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
				}
				if (pow(rd, 3) <= 2*pow(r0,6))
				{

					for (int k = 0; k < dimension; ++k)
					{
						double r6;
						r6 = pow(r0*r0/rd, 3);
						f[j][k] = f[j][k] + 4 * eps *
							(12*r6*r6 - 6*r6) *
							(r[i][k]-r[j][k])/rd;
					}
					
				}
			}
		}
	
	}
	
}

extern double ddot_(int* n, double* dx, int* incx, 
		double* dy, int* incy);

extern int daxpy_(int* n, double* da,
		double* dx, int* incx,
		double* dy, int* incy);

extern int dgetrf_(int* m, int* n, double* A, 
		int* lda, int* ipiv, int* info);

extern int dgetri_(int* m, double* a, int* lda,
		int* ipiv, double* work, int* lwork, int* info);

void PseudoForce(double u[beadNumber][dimension],
		int link[rodNumber][2],
		int g[rodNumber][rodNumber],
		double f[beadNumber][dimension])
{
	memset(f, 0, sizeof(f[0][0])*beadNumber*dimension);
	double metric[rodNumber*rodNumber];
	CalculateAij(u, u, g, metric);

	/* calculate the inverse of the metric matrix */
	int n = rodNumber;
	int ipiv[rodNumber];
	int info;
	double work[rodNumber];
	dgetrf_(&n, &n, metric, &n, ipiv, &info);
	dgetri_(&n, metric, &n, ipiv, work, &n, &info);

	for (int k = 0; k < beadNumber; ++k)
	{
		for (int i = 0; i < rodNumber; ++i)
		{
			for (int j = i+1; j < rodNumber; ++j)
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
						int nd = dimension;
						int inc = 1;
						uij = ddot_(&nd, *(u+i), &inc, *(u+j), &inc); 
					}


				}
			}
		}
	}



}
