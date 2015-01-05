#ifndef MDRUN_H_Q2M3FEG0
#define MDRUN_H_Q2M3FEG0

#include "main.h"

void MDRun(double r[beadNumber][dimension],
		int link[rodNumber][2],
		int g[rodNumber][rodNumber],
		double v0,
		unsigned long seed);

#endif /* end of include guard: MDRUN_H_Q2M3FEG0 */
