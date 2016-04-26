#include "main.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>


double Distance(double *r1, double *r2)
{
        double distance = 0;
        for (int i = 0; i < DIM; ++i)
        {
                distance += (r1[i] - r2[i])*(r1[i] - r2[i]);
        }
        return sqrt(distance);
}

