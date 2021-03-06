#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include "main.h"

void ConfigSingleRing(cl_double3* r)
{
	memset(r, 0, sizeof(cl_double3)*beadNumber);
	r[1].x = r[0][0] + cos(pi/6.0);
	r[1][1] = r[0][1] + sin(pi/6.0);
	r[beadNumber-1][0] = r[0][0] + cos(-pi/6.0);
	r[beadNumber-1][1] = r[0][1] + sin(-pi/6.0);
	for (int i = 2; i < beadNumber/2; i++) 
	{
		r[i][0] = r[1][0] + i - 1;	
		r[i][1] = r[1][1];	
		r[beadNumber-i][0] = r[beadNumber-1][0] + i - 1;	
		r[beadNumber-i][1] = r[beadNumber-1][1];	
	}

	if (beadNumber % 2 == 0)
       	{
		r[beadNumber/2][0] = beadNumber/2 + sqrt(3.0) - 2.0;
	}
}

void ConfigRingPair(double r[beadNumber][dimension])
{
	memset(r, 0, sizeof(r[0][0])*beadNumber*dimension);
	r[1][0] = r[0][0] + cos(pi/6.0);
	r[1][1] = r[0][1] + sin(pi/6.0);
	int ring = beadNumber/2;
	r[ring][0] = r[0][0] + cos(-pi/6.0);
	r[ring][1] = r[0][1] + sin(-pi/6.0);
	int halfRing = (beadNumber-1)/4;
	for (int i = 0; i < halfRing; ++i)
	{
		r[i+1][0] = r[1][0] + i; 
		r[i+1][1] = r[1][1];
		r[ring-i][0] = r[ring][0] + i;
		r[ring-i][1] = r[ring][1];
	}
	if (ring % 2 == 1)
	{
		r[halfRing+1][1] = r[0][0] + halfRing - 1 + sqrt(3.0);
	}
	for (int i = 1; i <= ring; ++i)
	{
		r[i+ring][0] = - r[i][0];
		r[i+ring][1] = r[i][1];
	}
	
}

void InitializeConfiguration(cl_double3* r)
{ 	
	struct stat sb;
	/*check if there is an input file*/
	if (stat("input.in", &sb) == 0) 
	{
		FILE *inputfile;
		inputfile = fopen("input.in", "r");
		for (int i = 0; i < beadNumber; i++) {
			fscanf(inputfile, "%lf\t%lf\t%lf\n", 
					&r[i].z, &r[i].z, &r[i].z);
		}
		fclose(inputfile);
	}
	else {
		ConfigSingleRing(r);
		/* ConfigRingPair(r); */
	}
}

