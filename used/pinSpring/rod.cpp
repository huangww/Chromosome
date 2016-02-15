#include "rod.hpp"
#include "simulation.hpp"
#include "ultilities.hpp"
#include "topo.hpp"
#include "bead.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <iomanip>

Rod::Rod(Simulation *simu) : Force(simu)
{
    link = create2DArray<int>(nRod, 2);
    g = create2DArray<int>(nRod, nRod);
    u = create2DArray<double>(nRod, DIM);
    b = create2DArray<double>(nRod, DIM);
    linkTable.nLinks = new int[nBead];
    linkTable.table = create2DArray<int>(3*nRod, 6);

    init();
}
Rod::~Rod() 
{
    delete2DArray(link);
    delete2DArray(g);
    delete2DArray(u);
    delete2DArray(b);
}

void Rod::init() 
{
    // init the link topology
    Topo *topo = new Topo(simulation);
    link = topo->init(link);
    delete topo;
    outputLinks();
    // printLinks();

    // init the metric matrix
    g = metricTensor();
    // printMetric();
    setLinkTable();
    
}

void Rod::printLinks() 
{
    std::cout << " Rod Topology: " << std::endl;
    for (int i = 0; i < nRod; ++i) {
        std::cout << "Rod " << i << ":  " << 
            link[i][0] << '\t' << link[i][1] 
            << std::endl;
    }
}

void Rod::outputLinks() 
{
    std::ofstream output("data/topo.dat");

    for (int i = 0; i < nRod; ++i) {
            output << std::setw(9) << link[i][0] << '\t' 
                << std::setw(9) << link[i][1] << std::endl;
    }

    output.close();
}

void Rod::printMetric()
{
    std::cout << "Metric Tensor: " << std::endl;
    for (int i = 0; i < nRod; ++i) {
        for (int j = 0; j < nRod; ++j) {
            std::cout << std::setw(3) << g[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}


double** Rod::linkVectorU() 
{
    double **r = bead->r;
    for (int i = 0; i < nRod; i++) {
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

double** Rod::linkVectorB()
{
    double **rs = bead->rs;
    for (int i = 0; i < nRod; i++) {
        for (int j = 0; j < DIM; j++) {
            b[i][j] = rs[link[i][1]][j] - rs[link[i][0]][j]; 
        }
    }
    
    return b;
}

int** Rod::metricTensor() 
{
    std::fill(&g[0][0], &g[0][0] + nRod * nRod, 0);

    for (int i = 0; i < nRod; i++) {
        g[i][i] = -2;
        for (int j = 0; j < nRod; j++) {
            if (i==j) continue;
            int condition1 = (link[i][0]-link[j][1]) *
                (link[i][1]-link[j][0]);
            int condition2 = (link[i][0]-link[j][0]) *
                (link[i][1]-link[j][1]);
            if (condition1 == 0) {
                g[i][j] = 1;
            }
            else if (condition2 == 0) {
                g[i][j] = -1;
            }
        }
    }

    return g;
}

void Rod::setLinkTable()
{
    int count = 0;
    for (int l = 0; l < nBead; ++l) {
        int countLink = 0;
        for (int i = 0; i < nRod; ++i) {
            for (int j = i+1; j < nRod; ++j) {
                if (g[i][j]!=0) {
                    int i0,i1,j0,j1;
                    i0 = link[i][0];
                    i1 = link[i][1];
                    j0 = link[j][0];
                    j1 = link[j][1];
                    if ((l-i0)*(l-i1)*(l-j0)*(l-j1)==0) {
                        linkTable.table[count][0] = i;
                        linkTable.table[count][1] = i0;
                        linkTable.table[count][2] = i1;
                        linkTable.table[count][3] = j;
                        linkTable.table[count][4] = j0;
                        linkTable.table[count][5] = j1;
                        count++;
                        countLink++;
                    }
                }
            }
        }
        linkTable.nLinks[l] = countLink;
    }
}

void Rod::matrixA(double *A) 
{
    std::fill(&A[0], &A[0] + nRod * nRod, 0);

    for (int i = 0; i < nRod; i++) {
        for (int j = 0; j < nRod; j++) {
            if (g[i][j] != 0) {
                // Indices transfer: A[j*Rodnumber+i] = A[i][j]
                A[j*nRod+i] = g[i][j]*Dot(&b[i][0],&u[j][0], DIM);
            }
        }
    }
}

void Rod::vectorB(double *x, double *B)
{
    for (int i = 0; i < nRod; i++) {
        double tmp[DIM];
        std::fill(&tmp[0], &tmp[0] + DIM, 0);
        for (int j = 0; j < nRod; j++) {
            if (g[i][j] != 0) {
                for (int k = 0; k < DIM; k++) {
                    tmp[k] += g[i][j]*x[j]*u[j][k];
                }
            }
        }
        
        double dotBiBi = 0.0;
        double dotTiTi = 0.0;
        for (int k = 0; k < DIM; k++) {
            dotBiBi += b[i][k]*b[i][k];
            dotTiTi += tmp[k]*tmp[k];
        }
        B[i] = (1.0 - dotBiBi)/(2.0*dt) -
            dt * dotTiTi/2.0;
    }
}

void Rod::solverPicard(double *x) 
{
    double A[nRod*nRod];
    matrixA(A);

    int n = nRod;
    int lda = nRod;
    int iPIv[nRod];
    int info;
    dgetrf_(&n, &n, A, &lda, iPIv, &info);

    double xold[nRod];
    for (int step = 0; step < 500; step++)
    {
        std::copy(&x[0], &x[0] + nRod, &xold[0]);
        double B[nRod];
        vectorB(x, B);

        char s = 'N';
        int nrhs = 1;
        dgetrs_(&s, &n, &nrhs, A, &n, iPIv, B, &n, &info);

        std::copy(&B[0], &B[0] + nRod, &x[0]);
        double maxDiff = fabs(xold[0] - x[0]);
        for (int i = 1; i < nRod; i++) {
            if (fabs(xold[i] - x[i]) > maxDiff) {
                maxDiff = fabs(xold[i] - x[i]);
            }
        }
        if (maxDiff < 1e-6) return;
    }

    bead->print();
    std::cout << "MaxStep exceeded in Picard!" << std::endl;
    exit(EXIT_FAILURE);
}

double** Rod::constraint(double** f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    u = linkVectorU();
    b = linkVectorB();
    double tension[nRod];
    std::fill(&tension[0], &tension[0] + nRod, 0);
    solverPicard(tension);

    for (int i = 0; i < nBead; i++) {
        for (int j = 0; j < nRod; j++) {
            if (link[j][0] == i) {
                for (int k = 0; k < DIM; k++) {
                    f[i][k] = f[i][k] + tension[j] * u[j][k];
                }
            }
            if (link[j][1] == i) {
                for (int k = 0; k < DIM; k++) {
                    f[i][k] = f[i][k] - tension[j] * u[j][k];
                }
            }
        }
    }

    return f;
}

double** Rod::pseudo(double **f) 
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    u = linkVectorU();
    double metric[nRod*nRod];
    std::fill(&metric[0], &metric[0] + nRod*nRod, 0);
    for (int i = 0; i < nRod; i++) {
        for (int j = 0; j < nRod; j++) {
            if (g[i][j] != 0) {
                metric[j*nRod+i] = -g[i][j]*Dot(&u[i][0],&u[j][0], DIM);
            }
        }
    }

    // calculate the inverse of the metric matrix
    int n = nRod;
    int iPIv[nRod];
    int info;
    double work[nRod];
    dgetrf_(&n, &n, metric, &n, iPIv, &info);
    dgetri_(&n, metric, &n, iPIv, work, &n, &info);

    int count = 0;
    for (int k = 0; k < nBead; ++k) {
        for (int l = 0; l < linkTable.nLinks[k]; ++l) {
            int i,i0,i1,j,j0,j1;
            i = linkTable.table[count][0];
            i0 = linkTable.table[count][1];
            i1 = linkTable.table[count][2];
            j = linkTable.table[count][3];
            j0 = linkTable.table[count][4];
            j1 = linkTable.table[count][5];
            count++;
            double uij;
            uij = Dot(&u[i][0], &u[j][0], DIM); 
            double pgr[DIM];
            for (int m = 0; m < DIM; ++m) {
                pgr[m] = (Delta(k,i1) - Delta(k,i0)) *
                    (u[j][m] - uij * u[i][m]) +
                    (Delta(k,j1) - Delta(k, j0)) *
                    (u[i][m] - uij * u[j][m]);			
                f[k][m] = f[k][m] + g[j][i] * metric[i*nRod+j] * pgr[m];
            }
        }
    }

    return f;
}

double Rod::detBandMetric(int n, double *coeff)
{
    double det[n+1];
    det[0] = 1; det[1] = 2;
    for (int i = 1; i < n; ++i) {
        det[i+1] = 2 * det[i] - coeff[i-1] * coeff[i-1] * det[i-1];
    }
    return det[n];
}

double **Rod::pseudoRing(double **f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    u = linkVectorU();
    double coeff[nRod];
    coeff[0] = -Dot(&u[0][0], &u[nRod-1][0], 3);
    double prodCoeff = coeff[0];
    for (int i = 1; i < nRod; i++) {
        coeff[i] = -Dot(&u[i][0], &u[i-1][0], 3);
        prodCoeff *= coeff[i];  
    }
    double det1 = detBandMetric(nRod-2, &coeff[2]);
    double det2 = detBandMetric(nRod-2, &coeff[1]);
    double det3 = detBandMetric(nRod-1, &coeff[1]);
    int sign = 1 - 2*(nRod%2);
    double detG = -coeff[0]*coeff[0]*det1 
        - coeff[nRod-1]*coeff[nRod-1]*det2
        + 2*det3 - 2*sign*prodCoeff; // check
    
    double metricInv[nRod][nRod];
    std::fill(&metricInv[0][0], &metricInv[0][0] + nRod*nRod, 0);
    double detTop, detBottom, cofG;

    detTop = 1; detBottom = detBandMetric(nRod-2, &coeff[3]);
    cofG = -coeff[1]*detTop*detBottom - sign*prodCoeff/coeff[1];
    metricInv[1][0] = cofG / detG;
    metricInv[0][1] = metricInv[1][0];

    detTop = 1; detBottom = detBandMetric(nRod-4, &coeff[4]);
    det1 = coeff[2]*detTop*detBottom;
    detTop = 1; detBottom = detBandMetric(nRod-3, &coeff[4]);
    det3 = coeff[2]*detTop*detBottom;
    cofG = coeff[0]*coeff[0]*det1 - 2*det3
        - sign*prodCoeff/coeff[2];
    metricInv[2][1] = cofG / detG;
    metricInv[1][2] = metricInv[2][1];

    for (int i = 3; i < nRod-2; ++i) {
        detTop = detBandMetric(i-2, &coeff[2]);
        detBottom = detBandMetric(nRod-i-2, &coeff[i+2]);
        det1 = coeff[i]*detTop*detBottom;
        detTop = detBandMetric(i-1, &coeff[1]);
        detBottom = detBandMetric(nRod-i-3, &coeff[i+2]);
        det2 = coeff[i]*detTop*detBottom;
        detTop = detBandMetric(i-1, &coeff[1]);
        detBottom = detBandMetric(nRod-i-2, &coeff[i+2]); //check
        det3 = coeff[i]*detTop*detBottom;
        cofG = coeff[0]*coeff[0]*det1 
        + coeff[nRod-1]*coeff[nRod-1]*det2 -
        2*det3 - sign*prodCoeff/coeff[i]; // check
        metricInv[i][i-1] = cofG / detG;
        metricInv[i-1][i] = metricInv[i][i-1];
    }

    detTop = detBandMetric(nRod-4, &coeff[2]); detBottom = 1;
    det1 = coeff[nRod-2]*detTop*detBottom;
    detTop = detBandMetric(nRod-3, &coeff[1]); detBottom = 1;
    det3 = coeff[nRod-2]*detTop*detBottom;
    cofG = coeff[0]*coeff[0]*det1 - 2*det3
        - sign*prodCoeff/coeff[nRod-2];
    metricInv[nRod-2][nRod-3] = cofG / detG;
    metricInv[nRod-3][nRod-2] = metricInv[nRod-2][nRod-3];

    detTop = detBandMetric(nRod-2, &coeff[1]); detBottom = 1;
    cofG = -coeff[nRod-1]*detTop*detBottom
        - sign*prodCoeff/coeff[nRod-1];
    metricInv[nRod-1][nRod-2] = cofG / detG;
    metricInv[nRod-2][nRod-1] = metricInv[nRod-1][nRod-2];
    
    cofG = -coeff[0]*detBandMetric(nRod-2, &coeff[2]) 
        - sign*prodCoeff / coeff[0];
    metricInv[0][nRod-1] = cofG / detG;
    metricInv[nRod-1][0] = metricInv[0][nRod-1];

    int count = 0;
    for (int k = 0; k < nBead; ++k) {
        for (int l = 0; l < linkTable.nLinks[k]; ++l) {
            int i,i0,i1,j,j0,j1;
            i = linkTable.table[count][0];
            i0 = linkTable.table[count][1];
            i1 = linkTable.table[count][2];
            j = linkTable.table[count][3];
            j0 = linkTable.table[count][4];
            j1 = linkTable.table[count][5];
            count++;
            double uij;
            uij = Dot(&u[i][0], &u[j][0], DIM); 
            double pgr[DIM];
            for (int m = 0; m < DIM; ++m) {
                pgr[m] = (Delta(k,i1) - Delta(k,i0)) *
                    (u[j][m] - uij * u[i][m]) +
                    (Delta(k,j1) - Delta(k, j0)) *
                    (u[i][m] - uij * u[j][m]);			
                f[k][m] = f[k][m] + g[j][i]*metricInv[j][i]*pgr[m];
            }
        }
    }

    return f;
}
