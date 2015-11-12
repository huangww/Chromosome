#include "main.h"
#include <string.h>

double GyrationRadiusSquare(double r[numUnit][DIM])
{
	double rCM[DIM];
        memset(rCM, 0, sizeof(rCM[0])*DIM);
	for (int i = 0; i < numUnit; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			rCM[j] += r[i][j];
		}
	}
	for (int i = 0; i < DIM; ++i)
	{
		rCM[i] = rCM[i] / numUnit;
	}

	double rg = 0;
	for (int i = 0; i < numUnit; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			rg += (r[i][j] - rCM[j]) *
				(r[i][j] - rCM[j]);  
		}
	}
		
	return rg / numUnit;
}

