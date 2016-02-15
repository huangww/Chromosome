#ifndef FORCES_H_NBVYOWYJ
#define FORCES_H_NBVYOWYJ


#include "main.h"
void ForceAll(double r[][DIM],
        double v[][DIM],
        int link[][2],
        int angleList[][3],
        unsigned long seed,
        int step,
        double f[][DIM]);

#endif /* end of include guard: FORCES_H_NBVYOWYJ */
