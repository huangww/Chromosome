#ifndef UTILITIES_H_RIVGLQSO
#define UTILITIES_H_RIVGLQSO

#include "main.h"

double Distance(double *r1, double *r2);
double Angle(double *r1, double *r2);
void DistanceTable(double r[numUnit][DIM],
                double table[numUnit][numUnit]);

#endif /* end of include guard: UTILITIES_H_RIVGLQSO */
