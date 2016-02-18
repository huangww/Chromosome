#include <stdio.h>
#include <string.h>
#include <math.h>
#include "main.h"

double ddot(double* dx, double *dy)
{
	double result = 0;
	for (int i = 0; i < DIM; ++i)
	{
		result = result + dx[i] * dy[i];
	}
	return result;
}

double Distance(double *point1, double *point2)
{
	double result = 0;
	for (int i = 0; i < DIM; ++i)
	{
		result = result + 
			(point1[i] - point2[i]) *
			(point1[i] - point2[i]);
	}
	return sqrt(result);

}

void PrintMatrix(double matrix[numUnit][DIM])
{
	printf("\n");
	for (int i = 0; i < numUnit; i++) {
		for (int j = 0; j < DIM; j++) {
			printf("%lf  ", matrix[i][j]);
		}
		printf("\n");
	}
}

void CalculateVectorU(double r[numUnit][DIM], 
		int link[numLink][2], 
		double u[numLink][DIM])
{
	/*Calculate rod unit vector u*/
	for (int i = 0; i < numLink; i++) {
		double uLength = 0;
		for (int j = 0; j < DIM; j++) {
			u[i][j] = r[link[i][1]][j] - r[link[i][0]][j]; 
			uLength = uLength + u[i][j]*u[i][j];
		}
		uLength = sqrt(uLength);
		for (int j = 0; j < DIM; j++) {
			u[i][j] = u[i][j]/uLength;
		}
	}
}

void CalculateVectorB(double rs[numUnit][DIM], 
		int link[numLink][2], 
		double b[numLink][DIM])
{
	/*Calculate rod unit vector B*/
	for (int i = 0; i < numLink; i++) {
		for (int j = 0; j < DIM; j++) {
			b[i][j] = rs[link[i][1]][j] - rs[link[i][0]][j]; 
		}
	}
}

void MatrixMulVector(double matrix[DIM][DIM],
		double vector[DIM])
{
	double result[DIM];
        memset(result, 0, sizeof(result[0])*DIM);
	for (int i = 0; i < DIM; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			result[i] += matrix[i][j]*vector[j];
		}
	}

	memcpy(vector, result,  sizeof(result));
}


