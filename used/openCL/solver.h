#ifndef SOLVER_H_FBEAFDEZ
#define SOLVER_H_FBEAFDEZ

#include "main.h"


void Picard(double b[rodNumber][dimension], 
		double u[rodNumber][dimension],
		int g[rodNumber][rodNumber],
		double x[rodNumber]);

void CalculateAij(double b[rodNumber][dimension], 
		double u[rodNumber][dimension],
		int g[rodNumber][rodNumber], 
		double A[rodNumber*rodNumber]);


#endif /* end of include guard: SOLVER_H_FBEAFDEZ */
