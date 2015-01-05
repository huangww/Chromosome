#include "main.h"
#include "random.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "output.h"


double CalculateEnergy(double r[beadNumber][dimension])
{
	double energy = 0;
	for (int i = 0; i < beadNumber; ++i)
	{
		energy -= r[i][0];
	}
	return energy;
}

void RotateAxis(int index1,
		int index2,
		double r[beadNumber][dimension],
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

void MatrixMulVector(double matrix[dimension][dimension],
		double vector[dimension])
{
	double result[dimension] = {0};
	for (int i = 0; i < dimension; ++i)
	{
		for (int j = 0; j < dimension; ++j)
		{
			result[i] += matrix[i][j]*vector[j];
		}
	}

	memcpy(vector, result,  sizeof(result));
}

void Move(double r[beadNumber][dimension],
		unsigned long seed)
{
	int index1, index2, temp;
	index1 = (int) (Ran(seed)*beadNumber);
	index2 = (int) (Ran(seed)*beadNumber);
	while(index2 == index1) {
		index2 = (int) (Ran(seed)*beadNumber);
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
	if (Ran(seed) > 0.5)
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
	else
	{
		for (int i = index2+1; i < beadNumber+index1; ++i)
		{
			double point[dimension];
			memcpy(point, &r[i%beadNumber][0], sizeof(point));
			for (int j = 0; j < dimension; ++j)
			{
				point[j] = point[j] - r[index2][j];
			}
			MatrixMulVector(matrix, point);
			for (int j = 0; j < dimension; ++j)
			{
				point[j] = point[j] + r[index2][j];
			}
			memcpy(&r[i%beadNumber][0], point, sizeof(point));
		}
		for (int i = beadNumber-1; i >= 0; --i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				r[i][j] = r[i][j] - r[0][j];
			}
		}
	}
		
}


int MonteCarloMove(double r[beadNumber][dimension],
		double Teff,
		unsigned long seed)
{
	double rTry[beadNumber][dimension];
	double dE;
	memcpy(rTry, r, sizeof(r[0][0])*beadNumber*dimension);
	Move(rTry, seed);
	dE = CalculateEnergy(rTry) - CalculateEnergy(r);
	if (dE <= 0)
	{
		memcpy(r, rTry, sizeof(r[0][0])*beadNumber*dimension);
		return 1;
	}
	else if (Ran(seed) < exp(-dE/Teff))
	{
		memcpy(r, rTry, sizeof(r[0][0])*beadNumber*dimension);
		return 1;
	}
	return 0;
}

void Equilibration(double r[beadNumber][dimension], 
		double Teff,
		unsigned long seed,
		int equilibrateSteps)
{
	for (int step = 0; step < equilibrateSteps; ++step)
	{
		MonteCarloMove(r, Teff, seed);

	}
}

