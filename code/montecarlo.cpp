#include "montecarlo.hpp"
#include "input.hpp"
#include "potential.hpp"
#include "ultilities.hpp"
#include "random.hpp"
#include "config.hpp"
#include <cmath>
#include <algorithm>
#include <random>

Montecarlo::Montecarlo() { }
Montecarlo::~Montecarlo() { }

void Montecarlo::setParameter(Input *input) 
{
    if (input->parameter.count("seed") == 0) {
        std::random_device rd;
        seed = rd();
        std::cout << "seed = " << seed << std::endl;
    } else {
        seed = long(input->parameter["seed"]);
    }
    if (input->parameter.count("tempEff") == 0) {
        throw "Parameter \"tempEff\" is not specified!";
    }
    tempEff = input->parameter["tempEff"];
    nBead = int(input->parameter["nBead"]);
    topoType = int(input->parameter["topoType"]);
}

void Montecarlo::pivot(int N, int *ipivot) 
{
    ipivot[0] = int(Ran(seed)*N);
    ipivot[1] = int(Ran(seed)*N);
    while(ipivot[0] == ipivot[1]) {
        ipivot[1] = (int) (Ran(seed)*N);
    }
    if (ipivot[0] > ipivot[1]) {
        int temp = ipivot[0];
        ipivot[0] = ipivot[1];
        ipivot[1] = temp;
    }
}

void Montecarlo::rotateAxis(int* ipivot, double** r, 
        double* axis) 
{
    double dsum = 0;
    for (int i = 0; i < DIM; ++i) {
        axis[i] = r[ipivot[1]][i] - r[ipivot[0]][i];
        dsum += axis[i]*axis[i];
    }
    for (int i = 0; i < DIM; ++i) {
        axis[i] = axis[i]/sqrt(dsum);
    }
}

void Montecarlo::rotateMatrix(
        double theta,
        double* axis,
        double** matrix)
{
    matrix[0][0] = cos(theta) + axis[0]*axis[0]*
        (1-cos(theta));
    matrix[0][1] = axis[0]*axis[1]*(1-cos(theta)) -
        axis[2]*sin(theta);
    matrix[0][2] = axis[0]*axis[2]*(1-cos(theta)) + 
        axis[1]*sin(theta);
    matrix[1][0] = axis[1]*axis[0]*(1-cos(theta)) +
        axis[2]*sin(theta);
    matrix[1][1] = cos(theta) + axis[1]*axis[1]*
        (1-cos(theta));
    matrix[1][2] = axis[1]*axis[2]*(1-cos(theta)) -
        axis[0]*sin(theta);
    matrix[2][0] = axis[2]*axis[0]*(1-cos(theta)) -
        axis[1]*sin(theta);
    matrix[2][1] = axis[2]*axis[1]*(1-cos(theta)) +
        axis[0]*sin(theta);
    matrix[2][2] = cos(theta) + axis[2]*axis[2]*
		(1-cos(theta));
}
    
void Montecarlo::moveRing(int N, double** r) 
{
    // pre-move, randomly select pivot, calculate rotate axis and rotate matrix
    int ipivot[2];
    pivot(N, ipivot);
    
    double axis[DIM];
    rotateAxis(ipivot, r, axis);

    double theta;
    double phi = PI;
    theta = (2*Ran(seed)-1)*phi;
    double** matrix = create2DArray<double>(DIM, DIM);
    rotateMatrix(theta, axis, matrix);

    // rotate move, with equal probability to rotate left or right part
    if (Ran(seed) > 0.5) {
        for (int i = ipivot[1]+1; i < N+ipivot[0]; ++i) {
            double point[DIM];
            std::copy(&r[i%N][0], &r[i%N][0]+DIM, &point[0]);
            for (int j = 0; j < DIM; ++j) {
                point[j] = point[j] - r[ipivot[1]][j];
            }
            MatrixMulVector(matrix, point, DIM);
            for (int j = 0; j < DIM; ++j) {
                point[j] = point[j] + r[ipivot[1]][j];
            }
            std::copy(&point[0], &point[0]+DIM, &r[i%N][0]);
        }
        for (int i = N-1; i >= 0; --i) {
            for (int j = 0; j < DIM; ++j) {
                r[i][j] = r[i][j] - r[0][j];
            }
        }
    } else {
        for (int i = ipivot[0]+1; i < ipivot[1]; ++i) {
            double point[DIM];
            std::copy(&r[i][0], &r[i][0]+DIM, &point[0]);
            for (int j = 0; j < DIM; ++j) {
                point[j] = point[j] - r[ipivot[0]][j];
            }
            MatrixMulVector(matrix, point, DIM);
            for (int j = 0; j < DIM; ++j) {
                point[j] = point[j] + r[ipivot[0]][j];
            }
            std::copy(&point[0], &point[0]+DIM, &r[i][0]);
        }
    }

    delete2DArray(matrix);
}

void Montecarlo::moveChain(int N, double** r)
{
    // pre-move, randomly select pivot, calculate rotate axis and rotate matrix
    int ipivot[2];
    pivot(N, ipivot);
    
    double axis[DIM];
    rotateAxis(ipivot, r, axis);

    double theta;
    double phi = PI;
    theta = (2*Ran(seed)-1)*phi;
    double** matrix = create2DArray<double>(DIM, DIM);
    rotateMatrix(theta, axis, matrix);

    for (int i = ipivot[0]+1; i < ipivot[1]; ++i) {
        double point[DIM];
        std::copy(&r[i][0], &r[i][0]+DIM, &point[0]);
        for (int j = 0; j < DIM; ++j) {
            point[j] = point[j] - r[ipivot[0]][j];
        }
        MatrixMulVector(matrix, point, DIM);
        for (int j = 0; j < DIM; ++j) {
            point[j] = point[j] + r[ipivot[0]][j];
        }
        std::copy(&point[0], &point[0]+DIM, &r[i][0]);
    }

    delete2DArray(matrix);
}

void Montecarlo::moveRingPair(int N, double** r)
{
    int ringSize = (N+1)/2;
    double** ring1 = create2DArray<double>(ringSize,DIM);
    double** ring2 = create2DArray<double>(ringSize,DIM);

    std::copy(&r[0][0], &r[0][0]+ringSize*DIM, &ring1[0][0]);
    std::copy(&r[0][0], &r[0][0]+DIM, &ring2[0][0]);
    std::copy(&r[ringSize][0], &r[ringSize][0]+(ringSize-1)*DIM, &ring2[1][0]);

    moveRing(ringSize, ring1);
    moveRing(ringSize, ring2);

    std::copy(&ring1[0][0], &ring1[0][0]+ringSize*DIM, 
            &r[0][0]);
    std::copy(&ring2[1][0], &ring2[1][0]+(ringSize-1)*DIM, 
            &r[ringSize][0]);

    delete2DArray(ring1);
    delete2DArray(ring2);
}

void Montecarlo::moveThreeRingPair(int N, double** r)
{
    int m[3];
    int monomer = (N+5)/2;
    m[1] = monomer * 245/1257;
    m[2] = monomer * 454/1257;
    m[0] = monomer - m[1] -m[2];

    int pairSize1 = 2*m[0]-1;
    int pairSize2 = 2*m[1]-1;
    int pairSize3 = 2*m[2]-1;
    double** pair1 = create2DArray<double>(pairSize1,DIM);
    double** pair2 = create2DArray<double>(pairSize2,DIM);
    double** pair3 = create2DArray<double>(pairSize3,DIM);

    std::copy(&r[0][0], &r[0][0]+pairSize1*DIM, &pair1[0][0]);
    std::copy(&r[0][0], &r[0][0]+DIM, &pair2[0][0]);
    std::copy(&r[0][0], &r[0][0]+DIM, &pair3[0][0]);
    std::copy(&r[pairSize1][0], 
            &r[pairSize1][0] + (pairSize2-1)*DIM, 
            &pair2[1][0]);
    std::copy(&r[pairSize1+pairSize2-1][0], 
            &r[pairSize1+pairSize2-1][0] + (pairSize3-1)*DIM,
            &pair3[1][0]);

    moveRingPair(pairSize1, pair1);
    moveRingPair(pairSize2, pair2);
    moveRingPair(pairSize3, pair3);

    std::copy(&pair1[0][0], &pair1[0][0]+pairSize1*DIM, 
            &r[0][0]);
    std::copy(&pair2[1][0], 
            &pair2[1][0] + (pairSize2-1)*DIM,
            &r[pairSize1][0]);
    std::copy(&pair3[1][0], 
            &pair3[1][0] + (pairSize3-1)*DIM,
            &r[pairSize1+pairSize2-1][0]);

    delete2DArray(pair1);
    delete2DArray(pair2);
    delete2DArray(pair3);
}

void Montecarlo::moveRingPairWithCentromere(int N, double** r)
{
    int cm = 10;        // To change accordingly
    int ringSize = (N+2)/2;
    double** ring = create2DArray<double>(ringSize, DIM);
    std::copy(&r[0][0], &r[0][0]+ringSize*DIM, &ring[0][0]);
    moveRing(ringSize, ring);
    std::copy(&ring[0][0], &ring[0][0]+ringSize*DIM, 
            &r[0][0]);
    delete2DArray(ring);

    Config *config = new Config();

    double** chain1 = create2DArray<double>(cm+1, DIM);
    std::fill(&chain1[0][0], &chain1[0][0]+DIM, 0);
    std::copy(&r[cm][0], &r[cm][0]+DIM, &chain1[cm][0]);
    config->fixedChain(cm+1, chain1);
    moveChain(cm+1, chain1);
    std::copy(&chain1[1][0], &chain1[1][0]+(cm-1)*DIM,
            &r[ringSize][0]);
    delete2DArray(chain1);

    double** chain2 = create2DArray<double>(ringSize-cm+1,
            DIM);
    std::fill(&chain2[0][0], &chain2[0][0]+DIM, 0);
    std::copy(&r[cm][0], &r[cm][0]+DIM, 
            &chain2[ringSize-cm][0]);
    config->fixedChain(ringSize-cm+1, chain2);
    moveChain(ringSize-cm+1, chain2);
    for (int i = 1; i < ringSize-cm; ++i) {
        std::copy(&chain2[i][0], &chain2[i][0]+DIM,
                &r[nBead-i][0]);
    }
    delete2DArray(chain2);

    delete(config);
}

void Montecarlo::moveTry(int N, double** r)
{
    switch (topoType) {
        case 0:
            moveRing(N, r);
            break;
        case 1:
            moveChain(N, r);
            break;
        case 2:
            moveRingPair(N, r);
            break;
        case 3:
            moveRingPairWithCentromere(N, r);
            break;
        case 4:
            moveThreeRingPair(N, r);
            break;
        default:
            throw "Invalid topoType!";
    }
}

double Montecarlo::energy(int N, double **r) 
{
    double eTotal = 0;
    for (int i = 0; i < N; ++i) {
        eTotal -= r[i][0];
    }
    // eTotal += potential->LennardJones(nBead, r);
    return eTotal;
}


int Montecarlo::move(double **r)
{
    double** rTry = create2DArray<double>(nBead, DIM);
    double dE;
    std::copy(&r[0][0], &r[0][0]+nBead*DIM, &rTry[0][0]);
    moveTry(nBead, rTry);
    dE = energy(nBead, rTry) - energy(nBead, r);
    if (dE <= 0) {
        std::copy(&rTry[0][0], &rTry[0][0]+nBead*DIM, 
                &r[0][0]);
        delete2DArray(rTry);
        return 1;
    } else if (Ran(seed) < exp(-dE/tempEff)) {
        std::copy(&rTry[0][0], &rTry[0][0]+nBead*DIM, 
                &r[0][0]);
        delete2DArray(rTry);
        return 1;
    }

    delete2DArray(rTry);
    return 0;
}

void Montecarlo::randomize(double **r)
{
    int moveStep = 1e6*Ran(seed) + 1e3;
    for (int i = 0; i < moveStep; ++i) {
        moveTry(nBead, r);
    }
}

void Montecarlo::equilibrate(double **r)
{
    int count = 0;
    do {
        count += move(r);
    } while (count < nBead*nBead);
}
