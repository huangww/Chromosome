#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "main.h"
#include "initialization.h"
#include "topology.h"
#include "output.h"
#include "montecarlo.h"
#include "mdrun.h"
#include "analysis.h"
	
int main(int argc, char *argv[])
{
	/* Step I: initialization */
	clock_t startTime = clock();

	/* 1. Initialize parameters */
	/* unsigned long seed = 5489; */
	unsigned long seed = 2;
	if (argc > 1 && atoi(argv[1])>0)
	{
		seed = atoi(argv[1]);
	}
	double Teff = 10.0;
	double v0;
	v0 = 1.0/Teff;
	/* v0 = 0; */

	/* 2. Initializing the connecting topology */
	int link[rodNumber][2];
	int g[rodNumber][rodNumber];
	Topology(link, g);

	/* 3. Initialize the configuration and thermalise*/
	double r[beadNumber][dimension]= {{0}};
	InitializeConfiguration(r);
	Equilibration(r, 3, Teff, seed, 1e4);
	
	/* 4. Open files for data output */
	char *outputDir = "data/";
	char parameters[80];
	if (v0 ==0)
		sprintf(parameters, "N%d_Tinf_%ld", beadNumber, seed);
	else
		sprintf(parameters, "N%d_T%d_%ld", beadNumber, (int)Teff, seed);
	char fileName[80];

	/* a) Output topological information */
	FILE  *topolFile;
	memset(fileName, 0, sizeof(fileName));
	sprintf(fileName, "%stopol.dat", outputDir);
	topolFile = fopen(fileName,"w");
	OutputTopol(topolFile, link);
	fclose(topolFile);

	/* b) Main output: configuration information */
	FILE *outputFile;
	memset(fileName, 0, sizeof(fileName));
	sprintf(fileName, "%sr_%s.dat", outputDir, parameters);
	outputFile = fopen(fileName,"w");

	/* c) Output input script for LAMMPS */
	/* FILE *lammpsFile; */
	/* memset(fileName, 0, sizeof(fileName)); */
	/* sprintf(fileName, "%slammpsIn_N%d.dat", outputDir, beadNumber); */
	/* lammpsFile = fopen(fileName,"w"); */
	/* Output4lammps(lammpsFile, link, r); */
	/* fclose(lammpsFile); */

	/* d) Output gyration radius */
	/* FILE *gyration; */
	/* memset(fileName, 0, sizeof(fileName)); */
	/* sprintf(fileName, "%srg_N%d_T%d_%ld.dat", outputDir, beadNumber, (int)(Teff), seed); */
	/* gyration = fopen(fileName,"w"); */
		
	/* Step II: Run the simulation */
	/* test */
	OutputConfiguration(outputFile, r);

	/* Type A: Monte Carlo Simulation */
	for (int step = 0; step < 0; step++) 
	/* for (int step = 0; step < runSteps; step++)  */
	{
		MonteCarloMove(beadNumber, r, 3, Teff, seed);
		/* Output samples */
		/* if (step % (int)(1e4) == 0)  */
		{
			fprintf(outputFile, "# step = %d\n",step);
			OutputConfiguration(outputFile, r);
			/* fprintf(gyration, "%lf\n",  */
			/* 		GyrationRadiusSquare(r)); */
		}
		/* Output states to the screen */
		if (step % (int)(runSteps/100) == 0)
		{
			printf("%d %% done!\n", step/(int)(runSteps/100));
		}

	}

	/* Type B: Molecular Dynamics Simulation */
	/* for (int step = 0; step < 0; step++)  */
	for (int step = 0; step < runSteps; step++) 
	{
		MDRun(r, link, g, v0, seed);
		/* Output samples */
		if (step % (int)(1e2) == 0) 
		{
			fprintf(outputFile, "# t = %f\n",step*dt);
			OutputConfiguration(outputFile, r);
		}
		/* Output states to the screen */
		if (step % (int)(runSteps/100) == 0)
		{
			printf("%d %% done!\n", step/(int)(runSteps/100));
		}

	}
		
	fclose(outputFile);
	/* fclose(gyration); */
	clock_t endTime = clock();
	double elapsedTime = (double)(endTime - startTime)/CLOCKS_PER_SEC;
	printf("Timming: %f seconds\n", elapsedTime);
	return 0;
}
