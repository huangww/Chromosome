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

	/* 3. Initializing the connecting topology */
	int link[rodNumber][2];
	int g[rodNumber][rodNumber];
	Topology(link, g);

	/* 4. Open files for data output */
	FILE *outputFile;
	char *outputDir = "data/";
	char fileName[80];
	sprintf(fileName, "%sr_N%d_test.dat", outputDir, beadNumber);
	outputFile = fopen(fileName,"w");

	/* Output topological information */
	OutputTopol(link, outputDir);
	Output4lammps(link, r, outputDir);

	/* temp = fopen("temp.txt","w+"); */
	/* FILE * randomForce; */
	/* randomForce = fopen("randomForce.txt","w+"); */

	/* Step II: Evolution of the dynamical system */
	for (int step = 0; step < maxStep; step++) {
		double fb[beadNumber][dimension] = {{0}};
		/* Generate random force */
		GenerateFb(fb, &seed);
	
		/* calculate pseudoforce */
		double fa[beadNumber][dimension] = {{0}};
		PseudoForce(r, link, g, fa);

		/* Add external force */
		/* 1) pinned SPB */
		for (int i = 0; i < dimension; ++i)
		{
			fa[0][i] = fa[0][i] - 1.0e3*r[0][i];
		}
		/* 2) external force feild */
		double v0 = 1.0;
		for (int i = 1; i < beadNumber; ++i)
		{
			fa[i][0] = fa[i][0] + 1.0 * v0;
		}
		
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
			for (int j = 0; j < dimension; j++)
		       	{
				r[i][j] = rs[i][j] + fc[i][j]*dt;
			}
		}
		
		/* Output samples */
		if (step % (int)(1e4) == 0 /* && step > maxStep/2 */) 
		{
			for (int i = 0; i < beadNumber; i++) 
			{
				for (int j = 0; j < dimension; j++) 
				{
					fprintf(outputFile, "%lf\n", r[i][j]);
				}
			}
			/* printf("%d steps done!\n", step); */
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
