#include <stdio.h>
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

int main(void)
{
	/* test function */
	/* Step I: initialization */
	clock_t startTime = clock();
	/* 1. Initialize the configuration */
	double r[beadNumber][dimension];
	InitializeConfiguration(r);

	/* 2. Initialize parameters */
	double rs[beadNumber][dimension];
	memcpy(rs, r, sizeof(r));
	long seed = 3;

	/* 3. Initializing random generator and topology */
	int link[rodNumber][2];
	int g[rodNumber][rodNumber];
	Topology(link, g);

	/* 4. Open files for data output */
	Output4lammps(link, r);
	FILE *outputFile, *topolFile;
	char *outputDir = "data/";
	char fileName[80];
	sprintf(fileName, "%sr_N%d_test.dat", outputDir, beadNumber);
	outputFile = fopen(fileName,"w");

	/* Output topological information */
	sprintf(fileName, "%stopol.dat", outputDir);
	topolFile = fopen(fileName,"w");
	for (int i = 0; i < rodNumber; i++) 
	{
		fprintf(topolFile, "%d\t%d\n", link[i][0], link[i][1]);
	}
	fclose(topolFile);
	/* temp = fopen("temp.txt","w+"); */
	/* randomForce = fopen("randomForce.txt","w+"); */

	/* Step II: Evolution of the dynamical system */
	for (int step = 0; step < maxStep; step++) {
		double fa[beadNumber][dimension] = {{0}};
		/* Generate random force */
		GenerateFb(fa, &seed);

		/* calculate pseudoforce */
		double fb[beadNumber][dimension] = {{0}};
		PseudoForce(r, link, g, fb);

		/* Predictive step */
		for (int i = 0; i < beadNumber; i++) 
		{
			for (int j = 0; j < dimension; j++) {
				rs[i][j] = r[i][j] + (fa[i][j]+fb[i][j])*dt;
			}
		}
		/* Calculate the constraint force on rods */
		double fc[beadNumber][dimension];
		CalculateFc(r, rs, link, g, fc);

		/* Corrective step */
		for (int i = 0; i < beadNumber; i++) 
		{
			for (int j = 0; j < dimension; j++) {
				r[i][j] = rs[i][j] + fc[i][j]*dt;
			}
		}
		
		/* Output samples */
		if (step % (int)(1e3) == 0) 
		{
			for (int i = 0; i < beadNumber; i++) {
				for (int j = 0; j < dimension; j++) {
					fprintf(outputFile, "%lf\n", r[i][j]);
				}
			}
		}
		
	}
		

	fclose(outputFile);
	/* fclose(temp); */
	/* fclose(randomForce); */
	clock_t endTime = clock();
	double elapsedTime = (double)(endTime - startTime)/CLOCKS_PER_SEC;
	printf("Timming: %f seconds\n", elapsedTime);
	return 0;
}
