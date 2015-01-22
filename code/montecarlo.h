#ifndef MONTECARLO_H_PSRHUOLL
#define MONTECARLO_H_PSRHUOLL


#include "main.h"

int MonteCarloMove(int N,
		double r[N][dimension],
		int topolType,
		double Teff,
		unsigned long seed);

void Equilibration(double r[beadNumber][dimension], 
		int topolType,
		double Teff,
		unsigned long seed,
		int equilibrateSteps);

#endif /* end of include guard: MONTECARLO_H_PSRHUOLL */
