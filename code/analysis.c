#include "main.h"

double GyrationRadiusSquare(double r[beadNumber][dimension])
{
	double rCM[dimension] = {0};
	for (int i = 0; i < beadNumber; ++i)
	{
		for (int j = 0; j < dimension; ++j)
		{
			rCM[j] += r[i][j];
		}
	}
	for (int i = 0; i < dimension; ++i)
	{
		rCM[i] = rCM[i] / beadNumber;
	}

	double rg = 0;
	for (int i = 0; i < beadNumber; ++i)
	{
		for (int j = 0; j < dimension; ++j)
		{
			rg += (r[i][j] - rCM[j]) *
				(r[i][j] - rCM[j]);  
		}
	}
		
	return rg / beadNumber;
}

