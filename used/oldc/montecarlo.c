#include "main.h"
#include "random.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "output.h"
#include "interaction.h"
#include "utilities.h"
#include "initialization.h"


double CalculateEnergy(int N, double r[N][DIM])
{
	double energy = 0;
	for (int i = 0; i < N; ++i)
	{
		energy -= r[i][0];
	}
	energy += LennardJonesPotential(r);
	return energy;
}

void RotateAxis(int index1,
		int index2,
		double r[][DIM],
		double axis[DIM])
{
	double dsum = 0;
	for (int i = 0; i < DIM; ++i)
	{
		axis[i] = r[index2][i] - r[index1][i];
		dsum += axis[i]*axis[i];
	}
	for (int i = 0; i < DIM; ++i)
	{
		axis[i] = axis[i]/sqrt(dsum);
	}

}


void RotateMatrix(double rotate[DIM][DIM],
		double axis[DIM],
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


void MoveRing(int N,
	double r[N][DIM],
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

	double axis[DIM];
	RotateAxis(index1, index2, r, axis);
	double theta;
	double phi = PI;
	theta = (2*Ran(seed)-1)*phi;
	double matrix[DIM][DIM];
	RotateMatrix(matrix, axis, theta);
	if (Ran(seed) > 0.5)
	{
		for (int i = index2+1; i < N+index1; ++i)
		{
			double point[DIM];
			memcpy(point, &r[i%N][0], sizeof(point));
			for (int j = 0; j < DIM; ++j)
			{
				point[j] = point[j] - r[index2][j];
			}
			MatrixMulVector(matrix, point);
			for (int j = 0; j < DIM; ++j)
			{
				point[j] = point[j] + r[index2][j];
			}
			memcpy(&r[i%N][0], point, sizeof(point));
		}
		for (int i = N-1; i >= 0; --i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				r[i][j] = r[i][j] - r[0][j];
			}
		}
	}
	else
	{
		for (int i = index1+1; i < index2; ++i)
		{
			double point[DIM];
			memcpy(point, &r[i][0], sizeof(point));
			for (int j = 0; j < DIM; ++j)
			{
				point[j] = point[j] - r[index1][j];
			}
			MatrixMulVector(matrix, point);
			for (int j = 0; j < DIM; ++j)
			{
				point[j] = point[j] + r[index1][j];
			}
			memcpy(&r[i][0], point, sizeof(point));
		}
	}
		
}

void MoveChain(int N,
	double r[N][DIM],
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

	double axis[DIM];
	RotateAxis(index1, index2, r, axis);
	double theta;
	double phi = PI;
	theta = (2*Ran(seed)-1)*phi;
	double matrix[DIM][DIM];
	RotateMatrix(matrix, axis, theta);
	for (int i = index1+1; i < index2; ++i)
	{
		double point[DIM];
		memcpy(point, &r[i][0], sizeof(point));
		for (int j = 0; j < DIM; ++j)
		{
			point[j] = point[j] - r[index1][j];
		}
		MatrixMulVector(matrix, point);
		for (int j = 0; j < DIM; ++j)
		{
			point[j] = point[j] + r[index1][j];
		}
		memcpy(&r[i][0], point, sizeof(point));
	}
		
}

void MoveRingPair(int N,
		double r[N][DIM],
		unsigned long seed)
{
	int ringSize = (N+1)/2;
	double ring1[ringSize][DIM];
	double ring2[ringSize][DIM];
	memcpy(ring1, r, sizeof(r[0][0])*DIM*ringSize);
	memcpy(ring2, r, sizeof(r[0][0])*DIM);
	memcpy(&ring2[1][0], &r[ringSize][0], sizeof(r[0][0])*DIM*(ringSize-1));
	MoveRing(ringSize, ring1, seed);
	MoveRing(ringSize, ring2, seed);
	memcpy(r, ring1, sizeof(r[0][0])*DIM*ringSize);
	memcpy(&r[ringSize][0], &ring2[1][0], sizeof(r[0][0])*DIM*(ringSize-1));

}

void MoveCentromerePair(int N,
		double r[N][DIM],
		unsigned long seed)
{
	int ringSize = (N+2)/2;
	double ring[ringSize][DIM];
	memcpy(ring, r, sizeof(r[0][0])*DIM*ringSize);
	MoveRing(ringSize, ring, seed);
	memcpy(r, ring, sizeof(r[0][0])*DIM*ringSize);

	double chain1[cm+1][DIM];
	memset(chain1, 0, sizeof(r[0][0])*DIM);
	memcpy(&chain1[cm][0], &r[cm][0], sizeof(r[0][0])*DIM);
	ConfigFixedChain(cm+1, chain1);

	double chain2[ringSize-cm+1][DIM];
	memset(chain2, 0, sizeof(r[0][0])*DIM);
	memcpy(&chain2[ringSize-cm][0], &r[cm][0], sizeof(r[0][0])*DIM);
	ConfigFixedChain(ringSize-cm+1, chain2);

	MoveChain(cm+1, chain1, seed);
	MoveChain(ringSize-cm+1, chain2, seed);

	memcpy(&r[ringSize][0], &chain1[1][0], sizeof(r[0][0])*DIM*(cm-1));
	for (int i = 1; i < ringSize-cm; ++i)
	{
		memcpy(&r[numUnit-i][0], &chain2[i][0], sizeof(r[0][0])*DIM);
	}

}

void MoveThreePair(int N,
		double r[N][DIM],
		unsigned long seed)
{
	int pairSize1 = 2*monomer[0]-1;
	double ringPair1[pairSize1][DIM];
	int pairSize2 = 2*monomer[1]-1;
	double ringPair2[pairSize1][DIM];
	int pairSize3 = 2*monomer[2]-1;
	double ringPair3[pairSize1][DIM];

	memcpy(ringPair1, r, sizeof(r[0][0])*DIM*pairSize1);
	memcpy(ringPair2, r, sizeof(r[0][0])*DIM);
	memcpy(&ringPair2[1][0], &r[pairSize1][0], sizeof(r[0][0])*DIM*(pairSize2-1));
	memcpy(ringPair3, r, sizeof(r[0][0])*DIM);
	memcpy(&ringPair3[1][0], &r[pairSize1+pairSize2-1][0], sizeof(r[0][0])*DIM*(pairSize3-1));

	MoveRingPair(pairSize1, ringPair1, seed);
	MoveRingPair(pairSize2, ringPair2, seed);
	MoveRingPair(pairSize3, ringPair3, seed);

	memcpy(r, ringPair1, sizeof(r[0][0])*DIM*pairSize1);
	memcpy(&r[pairSize1][0], &ringPair2[1][0], sizeof(r[0][0])*DIM*(pairSize2-1));
	memcpy(&r[pairSize1+pairSize2-1][0], &ringPair3[1][0], sizeof(r[0][0])*DIM*(pairSize3-1));
}

void Move(int N, double r[N][DIM],
		int topologyType,
		unsigned long seed)
{
	switch (topologyType) {
		case 0:
			MoveRing(N, r, seed);
			break;
		case 1:
			MoveChain(N, r, seed);
			break;
		case 2:
			MoveRingPair(N, r, seed);
			break;
		case 3:
			MoveCentromerePair(N, r, seed);
			break;
		case 4:
			MoveThreePair(N, r, seed);
			break;
		default:
			MoveRing(N, r, seed);
			
	}
}

int MonteCarloMove(int N,
		double r[N][DIM],
		int topologyType,
		double Teff,
		unsigned long seed)
{
	double rTry[N][DIM];
	double dE;
	memcpy(rTry, r, sizeof(r[0][0])*N*DIM);
	Move(N, rTry, topologyType, seed);
	dE = CalculateEnergy(N,rTry) - CalculateEnergy(N,r);
	if (dE <= 0)
	{
		memcpy(r, rTry, sizeof(r[0][0])*N*DIM);
		return 1;
	}
	else if (Ran(seed) < exp(-dE/Teff))
	{
		memcpy(r, rTry, sizeof(r[0][0])*N*DIM);
		return 1;
	}
	return 0;
}


void Equilibration(double r[numUnit][DIM], 
		int topologyType,
		double Teff,
		unsigned long seed,
		int equilibrateSteps)
{
	for (int step = 0; step < equilibrateSteps; ++step)
	{
		MonteCarloMove(numUnit, r, topologyType, Teff, seed);
	}

}
