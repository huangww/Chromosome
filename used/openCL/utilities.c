#include <math.h>
#include <stdio.h>
#include "main.h"

double ddot(double* dx, double *dy)
{
	double result = 0;
	for (int i = 0; i < dimension; ++i)
	{
		result = result + dx[i] * dy[i];
	}
	return result;
}

void PrintMatrix(double matrix[beadNumber][dimension])
{
	printf("\n");
	for (int i = 0; i < beadNumber; i++) {
		for (int j = 0; j < dimension; j++) {
			printf("%lf  ", matrix[i][j]);
		}
		printf("\n");
	}
}


void CalculateVectorU(double r[beadNumber][dimension], 
		int link[rodNumber][2], 
		double u[rodNumber][dimension])
{
	/*Calculate rod unit vector u*/
	for (int i = 0; i < rodNumber; i++) {
		double uLength = 0;
		for (int j = 0; j < dimension; j++) {
			u[i][j] = r[link[i][1]][j] - r[link[i][0]][j]; 
			uLength = uLength + u[i][j]*u[i][j];
		}
		uLength = sqrt(uLength);
		for (int j = 0; j < dimension; j++) {
			u[i][j] = u[i][j]/uLength;
		}
	}
}

void CalculateVectorB(double rs[beadNumber][dimension], 
		int link[rodNumber][2], 
		double b[rodNumber][dimension])
{
	/*Calculate rod unit vector B*/
	for (int i = 0; i < rodNumber; i++) {
		for (int j = 0; j < dimension; j++) {
			b[i][j] = rs[link[i][1]][j] - rs[link[i][0]][j]; 
		}
	}
}

