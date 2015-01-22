#include "main.h"
#include "random.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "output.h"
#include "interaction.h"
#include "utilities.h"
#include "initialization.h"


double CalculateEnergy(int N, double r[N][dimension])
{
	double energy = 0;
	for (int i = 0; i < N; ++i)
	{
		energy -= r[i][0];
	}
	/* energy += LennardJonesPotential(r); */
	return energy;
}

void RotateAxis(int index1,
		int index2,
		double r[][dimension],
		double axis[dimension])
{
	double dsum = 0;
	for (int i = 0; i < dimension; ++i)
	{
		axis[i] = r[index2][i] - r[index1][i];
		dsum += axis[i]*axis[i];
	}
	for (int i = 0; i < dimension; ++i)
	{
		axis[i] = axis[i]/sqrt(dsum);
	}

}


void RotateMatrix(double rotate[dimension][dimension],
		double axis[dimension],
		double theta)
{
	rotate[0][0] = cos(theta) + axis[0]*axis[0]*
		(1-cos(theta));
	rotate[0][1] = axis[0]*axis[1]*(1-cos(theta)) -
		axis[2]*sin(theta);
	rotate[0][2] = axis[0]*axis[2]*(1-cos(theta)) + 
		axis[1]*sin(theta);
	rotate[1][0] = axis[1]*axis[0]*(1-cos(theta)) +
		axis[2]*sin(theta);
	rotate[1][1] = cos(theta) + axis[1]*axis[1]*
		(1-cos(theta));
	rotate[1][2] = axis[1]*axis[2]*(1-cos(theta)) -
		axis[0]*sin(theta);
	rotate[2][0] = axis[2]*axis[0]*(1-cos(theta)) -
		axis[1]*sin(theta);
	rotate[2][1] = axis[2]*axis[1]*(1-cos(theta)) +
		axis[0]*sin(theta);
	rotate[2][2] = cos(theta) + axis[2]*axis[2]*
		(1-cos(theta));
}


void Move(int N,
	double r[N][dimension],
	int topolType,
	unsigned long seed)
{
	int index1, index2, temp;
	index1 = (int) (Ran(seed)*N);
	index2 = (int) (Ran(seed)*N);
	while(index2 == index1) {
		index2 = (int) (Ran(seed)*N);
	}
	if (index1 > index2)
	{
		temp = index2;
		index2 = index1;
		index1 = temp;
	}

	double axis[dimension];
	RotateAxis(index1, index2, r, axis);
	double theta;
	double phi = pi;
	theta = (2*Ran(seed)-1)*phi;
	double matrix[dimension][dimension];
	RotateMatrix(matrix, axis, theta);
	if (topolType == 0 && Ran(seed) > 0.5)
	{
		for (int i = index2+1; i < N+index1; ++i)
		{
			double point[dimension];
			memcpy(point, &r[i%N][0], sizeof(point));
			for (int j = 0; j < dimension; ++j)
			{
				point[j] = point[j] - r[index2][j];
			}
			MatrixMulVector(matrix, point);
			for (int j = 0; j < dimension; ++j)
			{
				point[j] = point[j] + r[index2][j];
			}
			memcpy(&r[i%N][0], point, sizeof(point));
		}
		for (int i = N-1; i >= 0; --i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				r[i][j] = r[i][j] - r[0][j];
			}
		}
	}
	else
	{
		for (int i = index1+1; i < index2; ++i)
		{
			double point[dimension];
			memcpy(point, &r[i][0], sizeof(point));
			for (int j = 0; j < dimension; ++j)
			{
				point[j] = point[j] - r[index1][j];
			}
			MatrixMulVector(matrix, point);
			for (int j = 0; j < dimension; ++j)
			{
				point[j] = point[j] + r[index1][j];
			}
			memcpy(&r[i][0], point, sizeof(point));
		}
	}
		
}


int MonteCarloMove(int N,
		double r[N][dimension],
		int topolType,
		double Teff,
		unsigned long seed)
{
	double rTry[N][dimension];
	double dE;
	memcpy(rTry, r, sizeof(r[0][0])*N*dimension);
	Move(N,rTry, topolType, seed);
	dE = CalculateEnergy(N,rTry) - CalculateEnergy(N,r);
	if (dE <= 0)
	{
		memcpy(r, rTry, sizeof(r[0][0])*N*dimension);
		return 1;
	}
	else if (Ran(seed) < exp(-dE/Teff))
	{
		memcpy(r, rTry, sizeof(r[0][0])*N*dimension);
		return 1;
	}
	return 0;
}

void EquilibrateCentromerePair(
		double r[beadNumber][dimension],
		double Teff,
		unsigned long seed,
		int equilibrateSteps)
{
	int ringSize = (beadNumber+2)/2;
	double ring[ringSize][dimension];
	memcpy(ring, r, sizeof(r[0][0])*dimension*ringSize);
	for (int step = 0; step < equilibrateSteps; ++step)
	{
		MonteCarloMove(ringSize, ring, 0, Teff, seed);
	}
	memcpy(r, ring, sizeof(r[0][0])*dimension*ringSize);

	double chain1[cm+1][dimension];
	memset(chain1, 0, sizeof(r[0][0])*dimension);
	memcpy(&chain1[cm][0], &r[cm][0], sizeof(r[0][0])*dimension);
	ConfigFixedChain(cm+1, chain1);

	double chain2[ringSize-cm+1][dimension];
	memset(chain2, 0, sizeof(r[0][0])*dimension);
	memcpy(&chain2[ringSize-cm][0], &r[cm][0], sizeof(r[0][0])*dimension);
	ConfigFixedChain(ringSize-cm+1, chain2);

	for (int step = 0; step < equilibrateSteps; ++step)
	{
		MonteCarloMove(cm+1, chain1, 1, Teff, seed);
		MonteCarloMove(ringSize-cm+1, chain2, 1, Teff, seed);
	}

	memcpy(&r[ringSize][0], &chain1[1][0], sizeof(r[0][0])*dimension*(cm-1));
	for (int i = 1; i < ringSize-cm; ++i)
	{
		memcpy(&r[beadNumber-i][0], &chain2[i][0], sizeof(r[0][0])*dimension);
	}

}

void Equilibration(double r[beadNumber][dimension], 
		int topolType,
		double Teff,
		unsigned long seed,
		int equilibrateSteps)
{
	if (topolType == 3)
	{
		EquilibrateCentromerePair(r, Teff, seed, equilibrateSteps);
	}
	else
	{
		for (int step = 0; step < equilibrateSteps; ++step)
		{
			MonteCarloMove(beadNumber, r, topolType, Teff, seed);
		}
		
	}
	
}
