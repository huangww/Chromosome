#ifndef INITIALIZATION_H

#define INITIALIZATION_H

#include "main.h"

int DetermineNumLink(int topologyType, int numUnit);

void ConfigFixedChain(int N, double r[N][DIM]);

void InitializeConfiguration(int topologyType,
		double r[numUnit][DIM]);

void RandomizeConfiguration(int topologyType,
        unsigned long seed,
        double r[numUnit][DIM]);

#endif /* end of include guard: INITIALIZATION_H */
