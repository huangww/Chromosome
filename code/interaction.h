#ifndef INTERACTION_H_1DVCLSIK
#define INTERACTION_H_1DVCLSIK


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

void LennardJonesForce(double r[beadNumber][dimension],
		double f[beadNumber][dimension]);

double LennardJonesPotential(double r[beadNumber][dimension]);

#endif /* end of include guard: INTERACTION_H_1DVCLSIK */
