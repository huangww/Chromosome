#ifndef MDRUN_H_Q2M3FEG0
#define MDRUN_H_Q2M3FEG0

#include "main.h"

void MDRun(double r[numUnit][DIM],
		int link[numLink][2],
		int g[numLink][numLink],
		double v0,
		unsigned long seed);

#endif /* end of include guard: MDRUN_H_Q2M3FEG0 */
