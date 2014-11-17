#ifndef UTILITIES_H_QP9WGDEO
#define UTILITIES_H_QP9WGDEO

#include "main.h"

double ddot(double* dx, double *dy);

void PrintMatrix(double matrix[beadNumber][dimension]);

void CalculateVectorU(double r[beadNumber][dimension], 
		int link[rodNumber][2], 
		double u[rodNumber][dimension]);

void CalculateVectorB(double rs[beadNumber][dimension], 
		int link[rodNumber][2], 
		double b[rodNumber][dimension]);

static inline int delta(int i, int j)
{
	return i==j ? 1 : 0;
}

#endif /* end of include guard: UTILITIES_H_QP9WGDEO */
