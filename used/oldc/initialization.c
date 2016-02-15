#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include "random.h"
#include "main.h"
#include "utilities.h"
#include "montecarlo.h"

void ConfigSingleChain(int N, double r[N][DIM])
{
    memset(r, 0, sizeof(r[0][0])*N*DIM);
    for (int i = 0; i < N; ++i)
    {
        r[i][0] = i;
    }
}

void ConfigFixedChain(int N, double r[N][DIM])
{
    double dcm;
    dcm= Distance(&r[0][0], &r[N-1][0]);
    double theta0, theta, alpha;
    theta0 = PI/2 - acos(r[N-1][2]/dcm);
    alpha= atan2(r[N-1][1], r[N-1][0]);

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

void ConfigSingleRing(int N, double r[N][DIM])
{
    memset(r, 0, sizeof(r[0][0])*N*DIM);

    r[1][0] = r[0][0] + cos(PI/6.0);
    r[1][1] = r[0][1] + sin(PI/6.0);
    r[N-1][0] = r[0][0] + cos(-PI/6.0);
    r[N-1][1] = r[0][1] + sin(-PI/6.0);

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

void StretchedRingConfiguration(int N, double r[N][DIM])
{
    memset(r, 0, sizeof(r[0][0])*N*DIM);
    for (int i = 0; i < N; ++i) {
        if (i<=N/2) {
            r[i][0] = i;
        } else {
            r[i][0] = N - i;
        }
    }
}

void QuenchedRingConfiguration(int N, double r[N][DIM])
{
    memset(r, 0, sizeof(r[0][0])*N*DIM);
    for (int i = 0; i < N; ++i) {
        r[i][0] = i%2;
        if (i%4==3)  r[i][0] = -1;
    }

}

void ConfigRingPair(int N, double r[N][DIM])
{
    memset(r, 0, sizeof(r[0][0])*N*DIM);

    int ringSize = (N+1)/2;
    double ring1[ringSize][DIM];
    double ring2[ringSize][DIM];

    ConfigSingleRing(ringSize, ring1);
    ConfigSingleRing(ringSize, ring2);
    memcpy(r, ring1, sizeof(r[0][0])*DIM*ringSize);
    for (int i = 1; i < ringSize; ++i)
    {
        r[ringSize+i-1][0] = - ring2[i][0];	
        r[ringSize+i-1][1] = - ring2[i][1];	
    }
}

void ConfigCentromerePair(int N, double r[N][DIM])
{
    memset(r, 0, sizeof(r[0][0])*N*DIM);
    int ringSize = (N+2)/2;
    if (cm > ringSize)
    {
        printf("%s\n", "Improper centromere position!");
    }
    ConfigSingleRing(ringSize, r);

    double chain1[cm+1][DIM];
    memset(chain1, 0, sizeof(r[0][0])*DIM);
    memcpy(&chain1[cm][0], &r[cm][0], sizeof(r[0][0])*DIM);
    ConfigFixedChain(cm+1, chain1);
    memcpy(&r[ringSize][0], &chain1[1][0], sizeof(r[0][0])*DIM*(cm-1));

    double chain2[ringSize-cm+1][DIM];
    memset(chain2, 0, sizeof(r[0][0])*DIM);
    memcpy(&chain2[ringSize-cm][0], &r[cm][0], sizeof(r[0][0])*DIM);
    ConfigFixedChain(ringSize-cm+1, chain2);
    for (int i = 1; i < ringSize-cm; ++i)
    {
        chain2[i][2] = -chain2[i][2];
        memcpy(&r[N-i][0], &chain2[i][0], sizeof(r[0][0])*DIM);
    }

}

void ConfigThreePair(int N, double r[N][DIM])
{
    memset(r, 0, sizeof(r[0][0])*N*DIM);
    int *m = monomer;

    int pairSize1 = 2*m[0]-1;
    double pair1[pairSize1][DIM];
    ConfigRingPair(pairSize1, pair1);
    memcpy(r, pair1, sizeof(r[0][0])*DIM*pairSize1);

    int pairSize2 = 2*m[1]-1;
    double pair2[pairSize2][DIM];
    ConfigRingPair(pairSize2, pair2);
    for (int i = 1; i < pairSize2; ++i)
    {
        int i0 = 2*m[0]-2;
        r[i0+i][0] = pair2[i][1];
        r[i0+i][1] = pair2[i][0];
    }

    int pairSize3 = 2*m[2]-1;
    double pair3[pairSize3][DIM];
    ConfigRingPair(pairSize3, pair3);
    for (int i = 1; i < pairSize3; ++i)
    {
        int i0 = 2*(m[0]+m[1])-4;
        r[i0+i][0] = pair3[i][1];
        r[i0+i][2] = pair3[i][0];
    }

}

void InitializeConfiguration(int topologyType,
        double r[numUnit][DIM])
{ 	
    struct stat sb;
    /*check if there is an input file*/
    if (stat("input.in", &sb) == 0) 
    {
        FILE *inputfile;
        inputfile = fopen("input.in", "r");
        for (int i = 0; i < numUnit; i++) {
            fscanf(inputfile, "%lf\t%lf\t%lf\n", 
                    &r[i][0], &r[i][1], &r[i][2]);
        }
        fclose(inputfile);
    }
    else {
        switch (topologyType) {
            case 0:
                ConfigSingleRing(numUnit, r);
                /* StretchedRingConfiguration(numUnit, r); */
                /* QuenchedRingConfiguration(numUnit, r); */
                break;
            case 1:
                ConfigSingleChain(numUnit, r);
                break;
            case 2:
                ConfigRingPair(numUnit, r);
                break;
            case 3:
                ConfigCentromerePair(numUnit, r);
                break;
            case 4:
                ConfigThreePair(numUnit, r);
                break;
            default:
                ConfigSingleRing(numUnit, r);

        }
    }
}

void RandomizeConfiguration(int topologyType,
        unsigned long seed,
        double r[numUnit][DIM])
{
    int moveSteps = 1e6*Ran(seed) + 1e3;
    switch (topologyType) {
        case 0:
            for (int i = 0; i < moveSteps; ++i)
            {
                MoveRing(numUnit, r, seed);
            }
            break;
        default:
            MoveRing(numUnit, r, seed);
    }

}
