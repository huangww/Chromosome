#ifndef FORCES_H

#define FORCES_H

#include "main.h"

void CalculateVectorU(double r[beadNumber][dimension], 
		int link[rodNumber][2], 
		double u[rodNumber][dimension]);

void CalculateVectorB(double rs[beadNumber][dimension], 
		int link[rodNumber][2], 
		double b[rodNumber][dimension]);

void CalculateFc(double r[beadNumber][dimension],
		double rs[beadNumber][dimension],
		int link[rodNumber][2],
		int g[rodNumber][rodNumber],
		double fc[beadNumber][dimension]);

void GenerateFb(double fb[beadNumber][dimension],
		long* seed);
#endif /* end of include guard: FORCES_H */
