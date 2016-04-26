#ifndef UTILITIES_H_QP9WGDEO
#define UTILITIES_H_QP9WGDEO

#include "main.h"

double ddot(double* dx, double *dy);

double Distance(double *point1, double *point2);

void PrintMatrix(double matrix[numUnit][DIM]);

void CalculateVectorU(double r[numUnit][DIM], 
		int link[numLink][2], 
		double u[numLink][DIM]);

void CalculateVectorB(double rs[numUnit][DIM], 
		int link[numLink][2], 
		double b[numLink][DIM]);

void MatrixMulVector(double matrix[DIM][DIM],
		double vector[DIM]);

static inline int delta(int i, int j)
{
	return i==j ? 1 : 0;
}

#endif /* end of include guard: UTILITIES_H_QP9WGDEO */
