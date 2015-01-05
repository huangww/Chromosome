#ifndef FORCES_H

#define FORCES_H

#include "main.h"


void CalculateFc(double r[beadNumber][dimension],
		double rs[beadNumber][dimension],
		int link[rodNumber][2],
		int g[rodNumber][rodNumber],
		double fc[beadNumber][dimension]);

void GenerateFb(double fb[beadNumber][dimension],
		unsigned long seed);

void PseudoForce(double r[beadNumber][dimension],
		int link[rodNumber][2],
		int g[rodNumber][rodNumber],
		double f[beadNumber][dimension]);

void LennardJones(double r[beadNumber][dimension],
		double f[beadNumber][dimension]);

#endif /* end of include guard: FORCES_H */
