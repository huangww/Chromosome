#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "main.h"
#include "input.h"
#include "init.h"
#include "bond.h"
#include "output.h"
#include "mdrun.h"

/* 0: single chain */
/* 1: single chain */
/* 2: ring pair */
/* 3: pair with centromere */
/* 4: three pairs */
const int topologyType = 4;

/* int numUnit = 21; */
int numLink;
int numAngle;
double Teff = 10.0;
double dt = 0.001;
int runSteps = 1.e5;

/* int numUnit = 1253; */
/* int monomer[3] = {279,227,123}; */
int numUnit = 121;
int monomer[3] = {28,23,12};
int debug = 0;
const int cm = 100;

int main(int argc, char *argv[])
{
    clock_t startTime = clock();

    /* Step I: set parameters */
    unsigned long seed = 5489;
    if (argc > 1 && atoi(argv[1])>0)
    {
        seed = atoi(argv[1]);
    }
    numLink = DetermineNumLink(topologyType, numUnit);
    numAngle = DetermineNumAngle(topologyType, numUnit);

    /* allocate memory */
    double (*r)[DIM] = malloc(sizeof *r * numUnit);
    double (*v)[DIM] = malloc(sizeof *v * numUnit);
    double (*f)[DIM] = malloc(sizeof *f * numUnit);
    int (*link)[2] = malloc(sizeof *link * numLink);
    int (*angleList)[3] = malloc(sizeof *angleList * numLink);

    /* Step II: set up system */
    BondList(topologyType, link);
    if(AngleList(link, angleList) != numAngle) getchar();
    InitializeConfiguration(topologyType, r);
    InitializeVelocity(seed, v);

    /* Step III: open file for output */
    char *outputDir = "data/";
    char fileName[80];

    /* a) Output topological information */
    memset(fileName, 0, sizeof(fileName));
    sprintf(fileName, "%stopol.dat", outputDir);
    FILE *topolFile = fopen(fileName,"w");
    OutputBonds(topolFile, link);
    fclose(topolFile);

    memset(fileName, 0, sizeof(fileName));
    sprintf(fileName, "%sangles.dat", outputDir);
    FILE *angleFile = fopen(fileName,"w");
    OutputAngles(angleFile, angleList);
    fclose(angleFile);

    /* b) Main output: configuration information */
    memset(fileName, 0, sizeof(fileName));
    char suffix[80];
    sprintf(suffix, "N%d_T%d_%ld", numUnit, 
            (int)Teff, seed);
    sprintf(fileName, "%sr_%s.dat", outputDir, suffix);
    FILE *outputFile = fopen(fileName,"w");

    /* Step IV: Run and output */
    for (int step = 0; step < runSteps; ++step)
    {
        MDRun(seed, step, link, angleList, r, v, f); 
        if (step % (int)(1e3) == 0) 
        {
            fprintf(outputFile, "# t = %f\n",step*dt);
            OutputFrame(outputFile, r);
        }
        /* Output states to the screen */
        if (step % (int)(runSteps/100) == 0)
        {
            printf("%d %% done!\n", step/(int)(runSteps/100));
        }
    }

    /* Step V: release resources */
    free(r);
    free(v);
    free(f);
    free(link);
    free(angleList);
    fclose(outputFile);

    clock_t endTime = clock();
    float elapsedTime = (endTime - startTime)/CLOCKS_PER_SEC;
    printf("Timming: %f seconds\n", elapsedTime);
    return 0;
}



