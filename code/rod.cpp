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
#include "Eigen/Sparse"
#include "Eigen/SparseLU"

Rod::Rod(Simulation *simu) : Force(simu)
{
    link = create2DArray<int>(nRod, 2);
    u = create2DArray<double>(nRod+DIM, DIM);
    b = create2DArray<double>(nRod+DIM, DIM);
    g = create2DArray<int>(nRod, nRod);
    gSparse = SpMatD(nRod, nRod);
    linkTable.nLinks = new int[nBead];
    linkTable.table = create2DArray<int>(3*nRod, 6);

    init();
}
Rod::~Rod() 
{
    delete2DArray(link);
    delete2DArray(u);
    delete2DArray(b);
    delete2DArray(g);
    delete linkTable.nLinks;
    delete2DArray(linkTable.table);
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
    metricTensorSparse();
    
    // init linkTable
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

    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
            u[nRod+i][j] = 0;
        }
        u[nRod+i][i] = 1;
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
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
            b[nRod+i][j] = 0;
        }
        b[nRod+i][i] = 1;
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

void Rod::metricTensorSparse() 
{
    std::vector<T> tripletList;
    tripletList.reserve(3*nRod);
    for (int i = 0; i < nRod; i++) {
        tripletList.push_back(T(i,i,-2));
        for (int j = 0; j < nRod; j++) {
            if (i==j) continue;
            int condition1 = (link[i][0]-link[j][1]) *
                (link[i][1]-link[j][0]);
            int condition2 = (link[i][0]-link[j][0]) *
                (link[i][1]-link[j][1]);
            if (condition1 == 0) {
                tripletList.push_back(T(i,j,1));
            }
            else if (condition2 == 0) {
                tripletList.push_back(T(i,j,-1));
            }
        }
    }
    gSparse.setFromTriplets(tripletList.begin(), tripletList.end());
    gSparse.makeCompressed();
}

void Rod::setLinkTable()
{
    int count = 0;
    for (int l = 0; l < nBead; ++l) {
        int countLink = 0;
        for (int k = 0; k < gSparse.outerSize(); ++k) {
            for (SpMatD::InnerIterator it(gSparse, k); it; ++it) {
                int i = it.row();
                int j = it.col();
                if (j>i) {
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
    std::fill(&A[0], &A[0] + (nRod+DIM) * (nRod+DIM), 0);

    for (int i = 0; i < nRod; i++) {
        for (int j = 0; j < nRod; j++) {
            if (g[i][j] != 0) {
                // Indices transfer: A[j*Rodnumber+i] = A[i][j]
                A[j*(nRod+DIM)+i] = g[i][j]*Dot(&b[i][0],&u[j][0], DIM);
            }
        }
    }

    for (int i = 0; i < DIM; ++i) {
        A[(nRod+i)*(nRod+DIM)] = 
            -Dot(&b[0][0],&u[nRod+i][0], DIM);
        A[(nRod+i)*(nRod+DIM)+(nRod-1)] = 
            Dot(&b[nRod-1][0],&u[nRod+i][0], DIM);
        A[nRod+i] = -Dot(&b[nRod+i][0],&u[0][0], DIM);
        A[(nRod-1)*(nRod+DIM)+(nRod+i)] = 
            Dot(&b[nRod+i][0],&u[nRod-1][0], DIM);
        A[(nRod+i)*(nRod+DIM)+(nRod+i)] = 
            -Dot(&b[nRod+i][0],&u[nRod+i][0], DIM);
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
        if (i==0) {
            for (int k = 0; k < DIM; ++k) {
                tmp[k] -= x[nRod+k]*u[nRod+k][k];
            }
        }
        if (i==(nRod-1)) {
            for (int k = 0; k < DIM; ++k) {
                tmp[k] += x[nRod+k]*u[nRod+k][k];
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
    for (int i = 0; i < DIM; ++i) {
        B[nRod+i] = bead->rs[0][i]/dt;
        // prescribe moving term
        // B[nRod] = (bead->rs[0][0]-x(t))/dt;
    }
}

void Rod::solverPicard(double *x) 
{
    double A[(nRod+DIM)*(nRod+DIM)];
    matrixA(A);

    int n = nRod+DIM;
    int lda = nRod+DIM;
    int iPIv[nRod+DIM];
    int info;
    dgetrf_(&n, &n, A, &lda, iPIv, &info);

    double xold[nRod+DIM];
    for (int step = 0; step < 500; step++)
    {
        std::copy(&x[0], &x[0] + nRod + DIM, &xold[0]);
        double B[nRod+DIM];
        vectorB(x, B);

        //tmp output
        // if (fabs(bead->t - 0.0050) < 1e-4) {
        //     for (int i = 0; i < nRod+DIM; ++i) {
        //         std::cout << x[i] << ' ' << B[i] << std::endl;
        //     }
        //     std::cin.get();
        // }
        //tmp output

        char s = 'N';
        int nrhs = 1;
        dgetrs_(&s, &n, &nrhs, A, &n, iPIv, B, &n, &info);

        std::copy(&B[0], &B[0] + nRod + DIM, &x[0]);
        double maxDiff = fabs(xold[0] - x[0]);
        for (int i = 1; i < nRod+DIM; i++) {
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
    double tension[nRod+DIM];
    std::fill(&tension[0], &tension[0] + nRod + DIM, 0);
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

    for (int k = 0; k < DIM; ++k) {
        f[0][k] = f[0][k] + tension[nRod+k]*u[nRod+k][k];
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

double** Rod::pseudoSparse(double **f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    u = linkVectorU();
    SpMatD metric(gSparse);
    SpMatD uij(gSparse);
    for (int k = 0; k < uij.outerSize(); ++k) {
        for (SpMatD::InnerIterator it(uij, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            it.valueRef() = Dot(&u[i][0], &u[j][0], DIM);
        }
    }
    metric = -uij.cwiseProduct(gSparse);

    // calculate the inverse of the metric matrix
    SpMatD I(nRod,nRod);
    I.setIdentity();
    // Eigen::SparseLU<SpMatD> solver;
    Eigen::SimplicialLDLT<SpMatD> solver;
    solver.compute(metric);
    metric = solver.solve(I);
    metric = metric.cwiseProduct(gSparse);

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
            double pgr[DIM];
            for (int m = 0; m < DIM; ++m) {
                pgr[m] = (Delta(k,i1) - Delta(k,i0)) *
                    (u[j][m] - uij.coeff(i,j) * u[i][m]) +
                    (Delta(k,j1) - Delta(k, j0)) *
                    (u[i][m] - uij.coeff(i, j) * u[j][m]);			
                f[k][m] = f[k][m] + metric.coeff(i,j)*pgr[m];
            }
        }
    }


    return f;
}
