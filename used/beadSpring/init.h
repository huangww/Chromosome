#ifndef INIT_H_8YELSBKT
#define INIT_H_8YELSBKT

#include "main.h"

void ConfigFixedChain(int N, double r[N][DIM]);

void InitializeConfiguration(int topologyType,
		double r[][DIM]);

void InitializeVelocity(unsigned long seed,
        double v[numUnit][DIM]);

#endif /* end of include guard: INIT_H_8YELSBKT */
