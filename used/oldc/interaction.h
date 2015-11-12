#ifndef INTERACTION_H_1DVCLSIK
#define INTERACTION_H_1DVCLSIK


#include "main.h"


void CalculateFc(double r[numUnit][DIM],
		double rs[numUnit][DIM],
		int link[numLink][2],
		int g[numLink][numLink],
		double fc[numUnit][DIM]);

void GenerateFb(double fb[numUnit][DIM],
		unsigned long seed);

void PseudoForce(double r[numUnit][DIM],
		int link[numLink][2],
		int g[numLink][numLink],
		double f[numUnit][DIM]);

void LennardJonesForce(double r[numUnit][DIM],
		double f[numUnit][DIM]);

double LennardJonesPotential(double r[numUnit][DIM]);

#endif /* end of include guard: INTERACTION_H_1DVCLSIK */
