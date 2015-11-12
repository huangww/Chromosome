#include "config.hpp"
#include "simulation.hpp"
#include "parameter.hpp"
#include "random.hpp"
#include "bead.hpp"
#include "ultilities.hpp"
#include <cmath>
#include <iostream>

Config::Config(Simulation *simu) : Parameter(simu) { } 
Config::~Config() { }

double** Config::init() 
{
    double **r = bead->r;
    switch (topoType) {
        case 0:
            ring(nBead, r);
            break;
        case 1:
            chain(nBead, r);
            break;
        default:
            ring(nBead, r);
    }

    return r;

}


void Config::ring(int N, double **pos) 
{
    std::fill(&pos[0][0], &pos[0][0] + N * DIM, 0);
    
    pos[1][0] = pos[0][0] + cos(PI/6.0);
    pos[1][1] = pos[0][1] + sin(PI/6.0);
    pos[N-1][0] = pos[0][0] + cos(-PI/6.0);
    pos[N-1][1] = pos[0][1] + sin(-PI/6.0);

    if (N % 2 == 0) {
        pos[N/2][0] = N/2 + sqrt(3.0) - 2.0;
        for (int i = 2; i < N/2; i++) {
            pos[i][0] = pos[1][0] + i - 1;	
            pos[i][1] = pos[1][1];	
            pos[N-i][0] = pos[N-1][0] + i - 1;	
            pos[N-i][1] = pos[N-1][1];	
        }
    }
    else {
        for (int i = 2; i < N/2+1; i++) {
            pos[i][0] = pos[1][0] + i - 1;	
            pos[i][1] = pos[1][1];	
            pos[N-i][0] = pos[N-1][0] + i - 1;	
            pos[N-i][1] = pos[N-1][1];	
        }
    }
        
}

void Config::chain(int N, double **pos) 
{
    std::fill(&pos[0][0], &pos[0][0] + N * DIM, 0);

    for (int i = 0; i < N; ++i) {
        pos[i][0] = i;
    }
}

void Config::fixedChain(int N, double **pos) 
{
    double dcm;
    dcm= Distance(&pos[0][0], &pos[N-1][0], DIM);
    double theta0, theta, alpha;
    theta0 = PI/2 - acos(pos[N-1][2]/dcm);
    alpha= atan2(pos[N-1][1], pos[N-1][0]);

    if (N%2 == 0) {
        theta = acos((dcm-1)/(N-2));
        for (int i = 0; i < (N-2)/2; ++i) {
            pos[i+1][0] = pos[i][0] + cos(theta+theta0)*cos(alpha);
            pos[i+1][1] = pos[i][1] + cos(theta+theta0)*sin(alpha);
            pos[i+1][2] = pos[i][2] + sin(theta+theta0);
            pos[N-2-i][0] = pos[N-1-i][0] - cos(theta-theta0)*cos(alpha);
            pos[N-2-i][1] = pos[N-1-i][1] - cos(theta-theta0)*sin(alpha);
            pos[N-2-i][2] = pos[N-1-i][2] + sin(theta-theta0);
        }
    } else {
        theta = acos(dcm/(N-1));
        for (int i = 0; i < (N-1)/2; ++i) {
            pos[i+1][0] = pos[i][0] + cos(theta+theta0)*cos(alpha);
            pos[i+1][1] = pos[i][1] + cos(theta+theta0)*sin(alpha);
            pos[i+1][2] = pos[i][2] + sin(theta+theta0);
            pos[N-2-i][0] = pos[N-1-i][0] - cos(theta-theta0)*cos(alpha);
            pos[N-2-i][1] = pos[N-1-i][1] - cos(theta-theta0)*sin(alpha);
            pos[N-2-i][2] = pos[N-1-i][2] + sin(theta-theta0);
        }
    }
}

void Config::straightRing(int N, double **pos) 
{
    std::fill(&pos[0][0], &pos[0][0] + N * DIM, 0);

    for (int i = 0; i < N; ++i) {
        if (i<=N/2) {
            pos[i][0] = i;
        } else {
            pos[i][0] = N - i;
        }
    }
}

void Config::quenchedRing(int N, double **pos) 
{
    std::fill(&pos[0][0], &pos[0][0] + N * DIM, 0);

    for (int i = 0; i < N; ++i) {
        pos[i][0] = i%2;
    }
}

void Config::ringPair(int N, double **pos) 
{
    std::fill(&pos[0][0], &pos[0][0] + N * DIM, 0);

    int ringSize = (N+1)/2;
    double **ring1 = create2DArray<double>(ringSize, DIM);
    double **ring2 = create2DArray<double>(ringSize, DIM);

    ring(ringSize, ring1);
    ring(ringSize, ring2);
    std::copy(&ring1[0][0], &ring1[0][0] + ringSize * DIM, &pos[0][0]);
    for (int i = 1; i < ringSize; ++i) {
        pos[ringSize+i-1][0] = - ring2[i][0];	
        pos[ringSize+i-1][1] = - ring2[i][1];	
    }

    delete2DArray(ring1);
    delete2DArray(ring2);
}

void Config::ringPairWithCentromere(int N, double **pos) 
{
    std::fill(&pos[0][0], &pos[0][0] + N * DIM, 0);

    int ringSize = (N+2)/2;
    int cm = ringSize;
    if (cm > ringSize) { 
       std::cout << "Improper centromere position!" << std::endl;
    }
    ring(ringSize, pos);

    double **chain1 = create2DArray<double>(cm+1,DIM);
    std::fill(&chain1[0][0], &chain1[0][0] + DIM, 0);
    std::copy(&pos[cm][0], &pos[cm][0] + DIM, &chain1[cm][0]);
    fixedChain(cm+1, chain1);
    std::copy(&chain1[1][0], &chain1[1][0] + DIM * (cm-1), &pos[ringSize][0]);

    double **chain2 = create2DArray<double>(ringSize-cm+1,DIM);
    std::fill(&chain2[0][0], &chain2[0][0] + DIM, 0);
    std::copy(&pos[cm][0], &pos[cm][0] + DIM, &chain2[ringSize-cm][0]);
    fixedChain(ringSize-cm+1, chain2);
    for (int i = 1; i < ringSize-cm; ++i) {
        chain2[i][2] = -chain2[i][2];
        std::copy(&chain2[i][0], &chain2[i][0] + DIM, &pos[N-i][0]);
    }

    delete2DArray(chain1);
    delete2DArray(chain2);
}

void Config::threeRingPair(int N, double **pos) 
{
    std::fill(&pos[0][0], &pos[0][0] + N * DIM, 0);

    int m[3] = {10, 10, 10};

    int pairSize1 = 2*m[0]-1;
    double **pair1 = create2DArray<double>(pairSize1,DIM);
    ringPair(pairSize1, pair1);
    std::copy(&pair1[0][0], &pair1[0][0] + pairSize1 * DIM, &pos[0][0]);

    int pairSize2 = 2*m[1]-1;
    double **pair2 = create2DArray<double>(pairSize2,DIM);
    ringPair(pairSize2, pair2);
    for (int i = 1; i < pairSize2; ++i) {
        int i0 = 2*m[0]-2;
        pos[i0+i][0] = pair2[i][1];
        pos[i0+i][1] = pair2[i][0];
    }

    int pairSize3 = 2*m[2]-1;
    double **pair3 = create2DArray<double>(pairSize3,DIM);
    ringPair(pairSize3, pair3);
    for (int i = 1; i < pairSize3; ++i) {
        int i0 = 2*(m[0]+m[1])-4;
        pos[i0+i][0] = pair3[i][1];
        pos[i0+i][2] = pair3[i][0];
    }

    delete2DArray(pair1);
    delete2DArray(pair2);
    delete2DArray(pair3);

}
