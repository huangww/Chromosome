#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{
	int count = 1000;
	double r=0.001, f=0;
	FILE *outPutTest;
	outPutTest = fopen("outPutTest.dat","w");

	for (int i = 0; i < count; ++i)
	{
		double r0 = 0.75;
		double eps = 1.0;
		double rsd = r*r;
		double r6 ;
		r6 = pow(r0*r0/rsd, 3);
		f =  - 48 * eps *
		(r6 - 0.5) * r6  /r;
		fprintf(outPutTest, "%lf\t%lf\n", r, f);
		/* printf("%lf\t%lf\n", r, f); */
		r = r + 0.0001;
	}
	return 0;
	fclose(outPutTest);
}

