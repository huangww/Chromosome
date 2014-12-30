#include "omp.h"
#include <stdio.h>

static long num_steps = 100000000;
double step;
int main(void)
{
	double pi = 0.0;
	double x, sum=0.0;
	step = 1.0/(double) num_steps;
#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < num_steps; ++i)
	{
		x = (i+0.5)*step;
		sum = 4.0/(1.0 +x*x);
	}
	pi = step * sum;
	printf("%lf\n", pi);
	return 0;
}
