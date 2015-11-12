#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "main.h"
#include "input.h"
#include "initialization.h"
#include "topology.h"
#include "output.h"
#include "montecarlo.h"
#include "mdrun.h"
#include "analysis.h"


/* Some running setups */
/* topologyType = 0, 1, 2, 3, 4 */
/* 0 : single ring */
/* 1 : single chain */
/* 2 : ring pair */
/* 3 : centremere pair */
/* 4 : three pair */
const int topologyType = 0; 
const int numUnit = 100;

int numLink;
const double dt = 1e-4;
double const Teff = 10.0;
int runSteps = 1.e6;

const int cm = 100;
int monomer[3] = {279,227,123};
/* numUnit = (279+227+123)*2-5 */
/* const int numUnit = 1253; */
/* int monomer[3] = {28,23,12}; */


int main(int argc, char *argv[])
{
    /* Step I: initialization */
    clock_t startTime = clock();

    /* 1. Initialize parameters */
    unsigned long seed = 5489;
    if (argc > 1 && atoi(argv[1])>0)
    {
        seed = atoi(argv[1]);
    }
    double v0 = 1.0/Teff;
    /* double v0 = 0.0; */
    numLink = DetermineNumLink(topologyType, 
            numUnit);

    /* 2. Initializing the connecting topology */
    int link[numUnit][2];
    int g[numLink][numLink];
    Topology(topologyType, link, g);

    /* 3. Initialize the configuration and thermalise*/
    double r[numUnit][DIM];
    InitializeConfiguration(topologyType, r);
    /* RandomizeConfiguration(topologyType, seed, r); */
    /* Equilibration(r, topologyType, Teff, seed, 1e6); */

    /* 4. Open files for data output */
    char *outputDir = "data/";
    char parameters[80];
    if (v0 ==0)
        sprintf(parameters, "N%d_Tinf_%ld", numUnit, seed);
    else
        sprintf(parameters, "N%d_T%d_%ld", numUnit, (int)Teff, seed);
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
    /* sprintf(fileName, "%slammpsIn_N%d.dat", outputDir, numUnit); */
    /* lammpsFile = fopen(fileName,"w"); */
    /* Output4lammps(lammpsFile, link, r); */
    /* fclose(lammpsFile); */

    /* d) Output gyration radius */
    /* FILE *gyration; */
    /* memset(fileName, 0, sizeof(fileName)); */
    /* sprintf(fileName, "%srg_N%d_T%d_%ld.dat", outputDir, numUnit, (int)(Teff), seed); */
    /* gyration = fopen(fileName,"w"); */

    /* Step II: Run the simulation */
    /* test */
    /* OutputConfiguration(outputFile, r); */

    /* Type A: Monte Carlo Simulation */
    for (int step = 0; step < 0; step++) 
        /* for (int step = 0; step < runSteps; step++)  */
    {
        MonteCarloMove(numUnit, r, topologyType, Teff, seed);
        /* Output samples */
        if (step % (int)(1e2) == 0 && step >= 1e5) 
        {
            fprintf(outputFile, "# step = %d\n",step);
            OutputConfiguration(outputFile, r);
            /* fprintf(gyration, "%lf\n",  */
            /* GyrationRadiusSquare(r)); */
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
        /* Output samples */
        if (step % (int)(1e4) == 0) 
        {
            fprintf(outputFile, "# t = %f\n",step*dt);
            OutputConfiguration(outputFile, r);
        }
        /* Output states to the screen */
        if (step % (int)(runSteps/100) == 0)
        {
            printf("%d %% done!\n", step/(int)(runSteps/100));
        }
        MDRun(r, link, g, v0, seed);
    }

    fclose(outputFile);
    /* fclose(gyration); */
    clock_t endTime = clock();
    double elapsedTime = (double)(endTime - startTime)/CLOCKS_PER_SEC;
    printf("Timming: %f seconds\n", elapsedTime);
    return 0;
}
