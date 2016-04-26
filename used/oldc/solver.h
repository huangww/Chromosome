#ifndef SOLVER_H_FBEAFDEZ
#define SOLVER_H_FBEAFDEZ

#include "main.h"


void Picard(double b[numLink][DIM], 
		double u[numLink][DIM],
		int g[numLink][numLink],
		double x[numLink]);

void CalculateAij(double b[numLink][DIM], 
		double u[numLink][DIM],
		int g[numLink][numLink], 
		double A[numLink*numLink]);


#endif /* end of include guard: SOLVER_H_FBEAFDEZ */
