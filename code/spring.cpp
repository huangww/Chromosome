#include "spring.hpp"
#include "input.hpp"
#include "bead.hpp"
#include "ultilities.hpp"
#include "topo.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

Spring::Spring(Bead *beadPointer)
{ 
    bead = beadPointer;
    link = nullptr;
    g = nullptr;
    linkPair = nullptr;
    nPair = nullptr;
    u = nullptr;
    topo = nullptr;
}
Spring::~Spring() 
{ 
    delete2DArray(u);
    delete topo;
}

void Spring::setParameter(Input *input)
{
    nBead =int(input->parameter["nBead"]);
    // nLink = int(input->parameter["nLink"]);

    topo = new Topo();
    topo->setParameter(input);
    topo->init();
    nLink = topo->getNumLink();
    link = topo->getLink();
    g = topo->getMetricTensor();
    nPair = topo->getNumPair();
    linkPair = topo->getLinkPair();

    u = create2DArray<double>(nLink, DIM);
}


double** Spring::linkVector(double** r)
{
    for (int i = 0; i < nLink; i++) {
        double uLength = 0;
        for (int j = 0; j < DIM; j++) {
            u[i][j] = r[link[i][1]][j] - r[link[i][0]][j]; 
            uLength = uLength + u[i][j]*u[i][j];
        }
        uLength = sqrt(uLength);
        for (int j = 0; j < DIM; j++) {
            u[i][j] = u[i][j]/uLength;
        }
    }

    return u;
}


double** Spring::bending(double** r, double** f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    double kappa = 10.0;
    u = linkVector(r);
    int count = 0;
    for (int k = 0; k < nBead; ++k) {
        for (int l = 0; l < nPair[k]; ++l) {
            int i,i0,i1,j,j0,j1;
            i = linkPair[count][0];
            i0 = linkPair[count][1];
            i1 = linkPair[count][2];
            j = linkPair[count][3];
            j0 = linkPair[count][4];
            j1 = linkPair[count][5];
            count++;
            double uij;
            uij = Dot(&u[i][0], &u[j][0], DIM); 
            double pgr[DIM];
            for (int m = 0; m < DIM; ++m) {
                pgr[m] = (Delta(k,i1) - Delta(k,i0)) *
                    (u[j][m] - uij * u[i][m]) +
                    (Delta(k,j1) - Delta(k, j0)) *
                    (u[i][m] - uij * u[j][m]);			
                f[k][m] = f[k][m] + g[j][i] * kappa * pgr[m];
            }
        }
    }

    return f;
}


double** Spring::harmonic(double** r, double** f) 
{
    std::fill(&f[0][0], &f[0][0] + nBead*DIM, 0);

    double k = 30.0;
    for (int i = 0; i < nLink; ++i) {
        int i0 = link[i][0];
        int i1 = link[i][1];
        for (int j = 0; j < DIM; ++j) {
            double df;
            df = -k*(r[i0][j]-r[i1][j]); 
            f[i0][j] += df;
            f[i1][j] -= df;
        }
    }
    
    return f;
}
    
double** Spring::fene(double** r, double** f) 
{
    std::fill(&f[0][0], &f[0][0] + nBead*DIM, 0);

    double k = 30.0;
    double R0 = 1.5;
    for (int i = 0; i < nLink; ++i) {
        int i0 = link[i][0];
        int i1 = link[i][1];
        double d = Distance(&r[i0][0], &r[i1][0], DIM);
        for (int j = 0; j < DIM; ++j) {
            double rlogarg = 1 - (d*d)/(R0*R0);
            if (rlogarg < 0.01) rlogarg = 0.01;
            double df;
            df = -k*(r[i0][j]-r[i1][j])/rlogarg; 
            f[i0][j] += df;
            f[i1][j] -= df;
        }
    }
    return f;
}
