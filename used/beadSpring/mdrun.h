
#ifndef MDRUN_H_WLH9KP3X
#define MDRUN_H_WLH9KP3X

#include "main.h"

void MDRun(unsigned long seed,
        int step,
        int link[][2],
        int angleList[][3],
        double r[][DIM], 
        double v[][DIM],
        double f[][DIM]);

#endif /* end of include guard: MDRUN_H_WLH9KP3X */
