#ifndef EQUILIBRATE_H_WDABYV54
#define EQUILIBRATE_H_WDABYV54

#include "main.h"

int MonteCarloMove(double r[beadNumber][dimension],
		double Teff,
		unsigned long seed);

void Equilibration(double r[beadNumber][dimension], 
		double Teff,
		unsigned long seed,
		int equilibrateSteps);

#endif /* end of include guard: EQUILIBRATE_H_WDABYV54 */
