#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "initialization.h"
#include "topology.h"
#include "forces.h"
#include "random.h"
#include "main.h"
#include "utilities.h"
#include "output.h"
#include "equilibrate.h"
	
int main(int argc, char *argv[])
{
	/* test function */
	/* Step I: initialization */
	clock_t startTime = clock();

	/* 1. Initialize parameters */
	unsigned long seed = 5489;
	if (argc > 1 && atoi(argv[1])>0)
	{
		seed = atoi(argv[1]);
	}
	double Teff = 5000.0;
	double v0;
	v0 = 1.0/Teff;

	/* 2. Initializing the connecting topology */
	int link[rodNumber][2];
	int g[rodNumber][rodNumber];
	Topology(link, g);

	/* 3. Open files for data output */
	char *outputDir = "data/";
	char fileName[80];
	/* FILE * temp; */
	/* temp = fopen("temp.txt","w+"); */

	/* Output topological information */
	FILE  *topolFile;
	memset(fileName, 0, sizeof(fileName));
	sprintf(fileName, "%stopol.dat", outputDir);
	topolFile = fopen(fileName,"w");
	OutputTopol(topolFile, link);
	fclose(topolFile);

	/* Main output: configuration information */
	FILE *outputFile;
	memset(fileName, 0, sizeof(fileName));
	sprintf(fileName, "%sr_N%d_T%d_%ld.dat", outputDir, beadNumber, (int)(Teff), seed);
	outputFile = fopen(fileName,"w");

	/* Output input script for LAMMPS */
	/* FILE *lammpsFile; */
	/* memset(fileName, 0, sizeof(fileName)); */
	/* sprintf(fileName, "%slammpsIn_N%d.dat", outputDir, beadNumber); */
	/* lammpsFile = fopen(fileName,"w"); */
	/* Output4lammps(lammpsFile, link, r); */
	/* fclose(lammpsFile); */

	/* 4. Initialize the configuration and thermalise*/
	double r[beadNumber][dimension];
	InitializeConfiguration(r);
	Equilibration(r, Teff, seed, 1e7);
	double rs[beadNumber][dimension];
	memcpy(rs, r, sizeof(r));

	/* Step II: Evolution of the dynamical system */
	for (int step = 0; step < runSteps; step++) 
	{
		/* Generate random force */
		double fb[beadNumber][dimension] = {{0}};
		GenerateFb(fb, seed);
	
		/* calculate pseudo-force */
		double fa[beadNumber][dimension] = {{0}};
		/* PseudoForce(r, link, g, fa); */

		/* calculate Lennard-Jones  */
		/* double flj[beadNumber][dimension] = {{0}}; */
		/* LennardJones(r, flj); */

		/* Add external force */
		/* 1) pinned SPB */
		for (int i = 0; i < dimension; ++i)
		{
			fa[0][i] = fa[0][i] - 2.0e3*r[0][i];
		}

		/* 2) external force field */
		for (int i = 1; i < beadNumber; ++i)
		{
			fa[i][0] = fa[i][0] + 1.0 * v0;
		}
		
		/* Predictive step */
		for (int i = 0; i < beadNumber; i++) 
		{
			for (int j = 0; j < dimension; j++) 
			{
				/* rs[i][j] = r[i][j] + (fa[i][j]+fb[i][j]+flj[i][j])*dt; */
				rs[i][j] = r[i][j] + (fa[i][j]+fb[i][j])*dt;
			}
		}

		/* Calculate the constraint force on rods */
		double fc[beadNumber][dimension];
		CalculateFc(r, rs, link, g, fc);

		/* Corrective step */
		for (int i = 0; i < beadNumber; i++) 
		{
			for (int j = 0; j < dimension; j++)
		       	{
				r[i][j] = rs[i][j] + fc[i][j]*dt;
			}
		}
		
		/* Output samples */
		if (step % (int)(1e4) == 0) 
		{
			fprintf(outputFile, "# t = %f\n",step*dt);
			OutputConfiguration(outputFile, r);
		}
		if (step % (int)(runSteps/100) == 0)
		{
			printf("%d %% done!\n", step/(int)(runSteps/100));
		}

	}
		
	fclose(outputFile);
	/* fclose(temp); */
	clock_t endTime = clock();
	double elapsedTime = (double)(endTime - startTime)/CLOCKS_PER_SEC;
	printf("Timming: %f seconds\n", elapsedTime);
	return 0;
}
