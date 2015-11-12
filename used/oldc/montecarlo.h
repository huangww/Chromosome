#ifndef MONTECARLO_H_PSRHUOLL
#define MONTECARLO_H_PSRHUOLL


#include "main.h"

int MonteCarloMove(int N,
		double r[N][DIM],
		int topologyType,
		double Teff,
		unsigned long seed);

void Equilibration(double r[numUnit][DIM], 
		int topologyType,
		double Teff,
		unsigned long seed,
		int equilibrateSteps);

void MoveRing(int N,
	double r[N][DIM],
	unsigned long seed);

#endif /* end of include guard: MONTECARLO_H_PSRHUOLL */
