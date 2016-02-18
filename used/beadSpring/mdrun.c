#include "main.h"
#include "forces.h"
#include "output.h"

void MDRun(unsigned long seed,
        int step,
        int link[][2],
        int angleList[][3],
        double r[][DIM], 
        double v[][DIM],
        double f[][DIM])
{
    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = 0; j < DIM; ++j)
        {
            r[i][j] += v[i][j]*dt + 0.5*f[i][j]*dt*dt;
            v[i][j] += 0.5*f[i][j]*dt;
        }
    }

    ForceAll(r, v, link, angleList, seed, step, f);

    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = 0; j < DIM; ++j)
        {
            v[i][j] += 0.5*f[i][j]*dt;
        }
    }

}
