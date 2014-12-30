#include <string.h>
#include <math.h>
#include <stdio.h>

#define dimension 3
#define pi 3.1415926

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

int main(void)
{

	double axis[dimension] = {1, 0, 0};
	double theta;
	theta = pi;
	double matrix[dimension][dimension];
	RotateMatrix(matrix, axis, theta);
	double point[dimension] = {3, 2, 0};
	MatrixMulVector(matrix, point);
	printf("%lf\t%lf\t%lf\n", point[0], point[1], point[2]);
}


