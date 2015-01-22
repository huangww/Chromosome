#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include "random.h"
#include "main.h"
#include "utilities.h"


void ConfigSingleChain(int N, double r[N][dimension])
{
	memset(r, 0, sizeof(r[0][0])*N*dimension);
	for (int i = 0; i < N; ++i)
	{
		r[i][0] = i;
	}
}

void ConfigFixedChain(int N, double r[N][dimension])
{
	double dcm;
	dcm= Distance(&r[0][0], &r[N-1][0]);
	double theta0, theta, alpha;
	theta0 = pi/2 - acos(r[N-1][2]/dcm);
	alpha= acos(r[N-1][0]/sqrt(r[N-1][0]*r[N-1][0]+
				r[N-1][1]*r[N-1][1]));

	if (N%2 == 0)
	{
		theta = acos((dcm-1)/(N-2));
		for (int i = 0; i < (N-2)/2; ++i)
		{
			r[i+1][0] = r[i][0] + cos(theta+theta0)*cos(alpha);
			r[i+1][1] = r[i][1] + cos(theta+theta0)*sin(alpha);
			r[i+1][2] = r[i][2] + sin(theta+theta0);
			r[N-2-i][0] = r[N-1-i][0] - cos(theta-theta0)*cos(alpha);
			r[N-2-i][1] = r[N-1-i][1] - cos(theta-theta0)*sin(alpha);
			r[N-2-i][2] = r[N-1-i][2] + sin(theta-theta0);
		}
	}
	else
	{
		theta = acos(dcm/(N-1));
		for (int i = 0; i < (N-1)/2; ++i)
		{
			r[i+1][0] = r[i][0] + cos(theta+theta0)*cos(alpha);
			r[i+1][1] = r[i][1] + cos(theta+theta0)*sin(alpha);
			r[i+1][2] = r[i][2] + sin(theta+theta0);
			r[N-2-i][0] = r[N-1-i][0] - cos(theta-theta0)*cos(alpha);
			r[N-2-i][1] = r[N-1-i][1] - cos(theta-theta0)*sin(alpha);
			r[N-2-i][2] = r[N-1-i][2] + sin(theta-theta0);
		}
	}
}

void ConfigSingleRing(int N, double r[N][dimension])
{
	memset(r, 0, sizeof(r[0][0])*N*dimension);

	r[1][0] = r[0][0] + cos(pi/6.0);
	r[1][1] = r[0][1] + sin(pi/6.0);
	r[N-1][0] = r[0][0] + cos(-pi/6.0);
	r[N-1][1] = r[0][1] + sin(-pi/6.0);

	if (N % 2 == 0)
       	{
		r[N/2][0] = N/2 + sqrt(3.0) - 2.0;
		for (int i = 2; i < N/2; i++) 
		{
			r[i][0] = r[1][0] + i - 1;	
			r[i][1] = r[1][1];	
			r[N-i][0] = r[N-1][0] + i - 1;	
			r[N-i][1] = r[N-1][1];	
		}
	}
	else
	{
		for (int i = 2; i < N/2+1; i++) 
		{
			r[i][0] = r[1][0] + i - 1;	
			r[i][1] = r[1][1];	
			r[N-i][0] = r[N-1][0] + i - 1;	
			r[N-i][1] = r[N-1][1];	
		}
	}
}

void ConfigRingPair(int N, double r[N][dimension])
{
	memset(r, 0, sizeof(r[0][0])*N*dimension);
	r[1][0] = r[0][0] + cos(pi/6.0);
	r[1][1] = r[0][1] + sin(pi/6.0);
	int ringSize = N/2;
	r[ringSize][0] = r[0][0] + cos(-pi/6.0);
	r[ringSize][1] = r[0][1] + sin(-pi/6.0);
	int halfRing = (N-1)/4;
	for (int i = 0; i < halfRing; ++i)
	{
		r[i+1][0] = r[1][0] + i; 
		r[i+1][1] = r[1][1];
		r[ringSize-i][0] = r[ringSize][0] + i;
		r[ringSize-i][1] = r[ringSize][1];
	}
	if (ringSize % 2 == 1)
	{
		r[halfRing+1][1] = r[0][0] + halfRing - 1 + sqrt(3.0);
	}
	for (int i = 1; i <= ringSize; ++i)
	{
		r[i+ringSize][0] = - r[i][0];
		r[i+ringSize][1] = r[i][1];
	}
	
}

void ConfigCentromerePair(int N, double r[N][dimension])
{
	memset(r, 0, sizeof(r[0][0])*N*dimension);
	int ringSize = (N+2)/2;
	if (cm > ringSize)
	{
		printf("%s\n", "Centromere position improper!");
	}
	ConfigSingleRing(ringSize, r);

	double chain1[cm+1][dimension];
	memset(chain1, 0, sizeof(r[0][0])*dimension);
	memcpy(&chain1[cm][0], &r[cm][0], sizeof(r[0][0])*dimension);
	ConfigFixedChain(cm+1, chain1);
	memcpy(&r[ringSize][0], &chain1[1][0], sizeof(r[0][0])*dimension*(cm-1));

	double chain2[ringSize-cm+1][dimension];
	memset(chain2, 0, sizeof(r[0][0])*dimension);
	memcpy(&chain2[ringSize-cm][0], &r[cm][0], sizeof(r[0][0])*dimension);
	ConfigFixedChain(ringSize-cm+1, chain2);
	for (int i = 1; i < ringSize-cm; ++i)
	{
		chain2[i][2] = -chain2[i][2];
		memcpy(&r[N-i][0], &chain2[i][0], sizeof(r[0][0])*dimension);
	}
		
}

void InitializeConfiguration(double r[beadNumber][dimension])
{ 	
	struct stat sb;
	/*check if there is an input file*/
	if (stat("input.in", &sb) == 0) 
	{
		FILE *inputfile;
		inputfile = fopen("input.in", "r");
		for (int i = 0; i < beadNumber; i++) {
			fscanf(inputfile, "%lf\t%lf\t%lf\n", 
					&r[i][0], &r[i][1], &r[i][2]);
		}
		fclose(inputfile);
	}
	else {
		/* ConfigSingleRing(beadNumber, r); */
		/* ConfigFixedChain(beadNumber, r); */
		ConfigCentromerePair(beadNumber, r);
	}
}

