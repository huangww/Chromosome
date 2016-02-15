#include "ultilities.hpp"
#include <cmath>
#include <algorithm>

double Dot(double* dx, double *dy, int dim)
{
    double result = 0;
    for (int i = 0; i < dim; ++i) {
        result +=  dx[i] * dy[i];
    }
    return result;
}

double Distance(double *point1, double *point2, int dim)
{
    double result = 0;
    for (int i = 0; i < dim; ++i) {
        result += (point1[i] - point2[i]) 
            * (point1[i] - point2[i]);
    }
    return sqrt(result);
}

double Mean(int N, double* x)
{
    double sum = 0;
    for (int i = 0; i < N; ++i) {
        sum += x[i];
    }
    return sum/N;
}

double Variance(int N, double* x)
{
    double xm = Mean(N, x);
    double sum = 0;
    for (int i = 0; i < N; ++i) {
        sum += (x[i] - xm)*(x[i] - xm);
    }
    return sum/N;
}

void MatrixMulVector(double** matrix,
		double* vector, int dim)
{
    double result[dim];
    std::fill(&result[0], &result[0]+dim, 0);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            result[i] += matrix[i][j]*vector[j];
        }
    }
    std::copy(&result[0], &result[0] + dim, &vector[0]);
}

