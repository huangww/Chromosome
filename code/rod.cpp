#include "rod.hpp"
#include "input.hpp"
#include "project.hpp"
#include "bead.hpp"
#include "constant.hpp"
#include "ultilities.hpp"
#include "topo.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Eigen/Sparse"
#include "Eigen/SparseLU"


Rod::Rod(Bead *beadPointer)
{
    bead = beadPointer;
    link = NULL;
    u = NULL;
    b = NULL;
    g = NULL;
    gSparse = NULL;
    nPair = NULL;
    linkPair = NULL;
    topo = NULL;
}
Rod::~Rod() 
{
    delete2DArray(u);
    delete2DArray(b);
    delete topo;
}

void Rod::setParameter(Input *input) 
{

    nBead = int(input->parameter["nBead"]);
    // nLink = int(input->parameter["nLink"]);
    dt = input->parameter["dt"];

    // setup the link topology
    topo = new Topo();
    topo->setParameter(input);
    topo->init();
    // topo->initSparse();
    nLink = topo->getNumLink();
    link = topo->getLink();
    g = topo->getMetricTensor();
    // gSparse = topo->getMetricTensorSparse();
    linkPair = topo->getLinkPair();
    nPair = topo->getNumPair();

    u = create2DArray<double>(nLink+DIM, DIM);
    b = create2DArray<double>(nLink+DIM, DIM);
    // init();
}

void Rod::init() 
{

    // metric in sparse format using Egen lib
    // metricTensorSparse();

}


void Rod::printMetric()
{
    std::cout << "Metric Tensor: " << std::endl;
    for (int i = 0; i < nLink; ++i) {
        for (int j = 0; j < nLink; ++j) {
            std::cout << std::setw(3) << g[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}

double** Rod::linkVectorU() 
{
    double **r = bead->r;
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

    // additional unit vector for first bead costraints
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
            u[nLink+i][j] = 0;
        }
        u[nLink+i][i] = 1;
    }

    return u;
}

double** Rod::linkVectorB()
{
    double **rs = bead->rs;
    for (int i = 0; i < nLink; i++) {
        for (int j = 0; j < DIM; j++) {
            b[i][j] = rs[link[i][1]][j] - rs[link[i][0]][j]; 
        }
    }
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
            b[nLink+i][j] = 0;
        }
        b[nLink+i][i] = 1;
    }
    return b;
}

void Rod::matrixA(double *A) 
{
    std::fill(&A[0], &A[0] + (nLink+DIM) * (nLink+DIM), 0);

    for (int i = 0; i < nLink; i++) {
        // normal A entities 
        for (int j = 0; j < nLink; j++) {
            if (g[i][j] != 0) {
                // Indices transfer: A[j*nLink+i] = A[i][j]
                A[j*(nLink+DIM)+i] = g[i][j]*Dot(&b[i][0],&u[j][0], DIM);
            }
        }
        // the one for pinned bead Fc[1]*b[k], symetric matrix
        // for a ring, A[N+i, 1] = u[1] . b[N+i]
        if (link[i][0]==0) {
            for (int j = 0; j < DIM; ++j) {
                A[i*(nLink+DIM)+nLink+j] = 
                    Dot(&b[nLink+j][0],&u[i][0], DIM);
                A[(nLink+j)*(nLink+DIM)+i] = 
                    Dot(&b[i][0],&u[nLink+j][0], DIM);
            }
        } 
        // the one for pinned bead Fc[1]*b[k]
        // for a ring, A[N+i, N] = u[N] . b[N+i]
        if (link[i][1]==0) {
            for (int j = 0; j < DIM; ++j) {
                A[(nLink+j)*(nLink+DIM)+i] = 
                    -Dot(&b[i][0],&u[nLink+j][0], DIM);
                A[i*(nLink+DIM)+(nLink+j)] = 
                    -Dot(&b[nLink+j][0],&u[i][0], DIM);
            }
        }
    }

    for (int i = 0; i < DIM; ++i) {
        A[(nLink+i)*(nLink+DIM)+(nLink+i)] = 
            Dot(&b[nLink+i][0],&u[nLink+i][0], DIM);
    }

}

void Rod::vectorB(double *x, double *B)
{
    // B =  (1 - b^2[i] - nonlinear term)/(2dt)
    // nonlinear term = ((Fc[i+1]-Fc[i])*dt)^2
    for (int i = 0; i < nLink; i++) {
        double tmp[DIM];
        std::fill(&tmp[0], &tmp[0] + DIM, 0);
        for (int j = 0; j < nLink; j++) {
            if (g[i][j] != 0) {
                for (int k = 0; k < DIM; k++) {
                    tmp[k] += g[i][j]*x[j]*u[j][k];
                }
            }
        }
        if (i==0) {
            for (int k = 0; k < DIM; ++k) {
                tmp[k] -= x[nLink+k]*u[nLink+k][k];
            }
        }
        if (i==(nLink-1)) {
            for (int k = 0; k < DIM; ++k) {
                tmp[k] += x[nLink+k]*u[nLink+k][k];
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

    // Additional constraint for the first bead
    for (int i = 0; i < DIM; ++i) {
        B[nLink+i] = -bead->rs[0][i]/dt;
        // prescribe moving term
        // B[nLink] = -(bead->rs[0][0]-x(t))/dt;
    }
}

void Rod::solverPicard(double *x) 
{
    double A[(nLink+DIM)*(nLink+DIM)];
    matrixA(A);

    int n = nLink+DIM;
    int lda = nLink+DIM;
    int iPIv[nLink+DIM];
    int info;
    dgetrf_(&n, &n, A, &lda, iPIv, &info);

    double xold[nLink+DIM];
    for (int step = 0; step < 1000; step++)
    {
        std::copy(&x[0], &x[0] + nLink + DIM, &xold[0]);
        double B[nLink+DIM];
        vectorB(x, B);

        //tmp output
        // if (fabs(bead->t - 0.0050) < 1e-4) {
        // for (int i = 0; i < nLink+DIM; ++i) {
        //     std::cout << x[i] << ' ' << B[i] << std::endl;
        // }
        // std::cin.get();
        // }
        //tmp output

        char s = 'N';
        int nrhs = 1;
        dgetrs_(&s, &n, &nrhs, A, &n, iPIv, B, &n, &info);

        std::copy(&B[0], &B[0] + nLink + DIM, &x[0]);
        double maxDiff = fabs(xold[0] - x[0]);
        for (int i = 1; i < nLink+DIM; i++) {
            if (fabs(xold[i] - x[i]) > maxDiff) {
                maxDiff = fabs(xold[i] - x[i]);
            }
        }
        if (maxDiff < 1e-6) return;
    }

    // bead->print();
    throw "MaxStep exceeded in Picard!";
}

void Rod::jacobian(double *x, double *A, double *B)
{
    // F[i] = \sum_j g[i,j] * b[i] . u[j] * x[j]
    // + (\sum_j g[i,j] * b[i] . u[j] * x[j])^2 * dt / 2
    //  + (b[i]^2 - 1)/(2dt)
    // J[i,j] = \partial F[i] / \partial x[j]
    // J[i,j] = g[i,j] * b[i].u[j] + dt * g[i,j]
    // * b[i].u[j] * (\sum_j g[i,j] * b[i].u[j] * x[j])

    std::fill(&A[0], &A[0] + (nLink+DIM) * (nLink+DIM), 0);
    std::fill(&B[0], &B[0] + (nLink+DIM), 0);

    double tmp[nLink+DIM]; //  sum_j g[i,j] * b[i].u[j] * x[j]
    std::fill(&tmp[0], &tmp[0] + (nLink+DIM), 0);
    for (int i = 0; i < nLink; i++) {
        // normal A entities 
        for (int j = 0; j < nLink; j++) {
            if (g[i][j] != 0) {
                // Indices transfer: A[j*nLink+i] = A[i][j]
                A[j*(nLink+DIM)+i] = g[i][j]*Dot(&b[i][0],&u[j][0], DIM);
                tmp[i] += A[j*(nLink+DIM)+i]*x[j];
            }
        }

        // the one for pinned bead Fc[1]*b[k], symetric matrix
        // for a ring, A[N+j, 1] = u[1] . b[N+j]
        if (link[i][0]==0) {
            for (int j = 0; j < DIM; ++j) {
                A[i*(nLink+DIM)+nLink+j] = 
                    Dot(&b[nLink+j][0],&u[i][0], DIM);
                tmp[nLink+j] += A[i*(nLink+DIM)+nLink+j]*x[i];
                A[(nLink+j)*(nLink+DIM)+i] = 
                    Dot(&b[i][0],&u[nLink+j][0], DIM);
                tmp[i] += A[i*(nLink+DIM)+nLink+j]*x[nLink+j];
            }
        } 
        // the one for pinned bead Fc[1]*b[k]
        // for a ring, A[N+i, N] = u[N] . b[N+i]
        if (link[i][1]==0) {
            for (int j = 0; j < DIM; ++j) {
                A[(nLink+j)*(nLink+DIM)+i] = 
                    -Dot(&b[i][0],&u[nLink+j][0], DIM);
                tmp[i] += A[(nLink+j)*(nLink+DIM)+i]*x[nLink+j];
                A[i*(nLink+DIM)+(nLink+j)] = 
                    -Dot(&b[nLink+j][0],&u[i][0], DIM);
                tmp[nLink+j] += A[i*(nLink+DIM)+(nLink+j)]*x[i];
            }
        }
    }

    for (int i = 0; i < DIM; ++i) {
        A[(nLink+i)*(nLink+DIM)+(nLink+i)] = 
            Dot(&b[nLink+i][0],&u[nLink+i][0], DIM);
        tmp[nLink+i] += 
        A[(nLink+i)*(nLink+DIM)+(nLink+i)] * x[nLink+i];
    }

    for (int i = 0; i < nLink+DIM; ++i) {
        for (int j = 0; j < nLink+DIM; ++j) {
            if (A[j*(nLink+DIM)+i]!=0) {
                A[j*(nLink+DIM)+i] *= (1. + tmp[i]*dt);
            }
        }
        if (i<nLink) {
            // if inverse is used, must use B[i]
            B[i] = -(tmp[i] + tmp[i]*tmp[i]*dt/2 +
                (Dot(&b[i][0], &b[i][0], DIM) - 1) / (2.*dt));
        } else {
            B[i] = -(tmp[i] + bead->rs[0][i-nLink]/dt);
        }
    }

}

void Rod::inverse(double *A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete [] IPIV;
    delete [] WORK;
}

void Rod::solverNewton(double *x) 
{
    double xold[nLink+DIM];
    for (int step = 0; step < 1000; step++)
    {
        std::copy(&x[0], &x[0] + nLink + DIM, &xold[0]);
        double A[(nLink+DIM)*(nLink+DIM)];
        double B[nLink+DIM];
        jacobian(xold, A, B);

        int n = nLink+DIM;
        int lda = nLink+DIM;
        int iPIv[nLink+DIM+1];
        int info;
        char s = 'N';
        int nrhs = 1;
        dgetrf_(&n, &n, A, &lda, iPIv, &info);
        dgetrs_(&s, &n, &nrhs, A, &n, iPIv, B, &n, &info);
        for (int i = 0; i < nLink+DIM; ++i) {
            x[i] = xold[i] + B[i];
        }

        // inverse(A, nLink+DIM);
        // for (int i = 0; i < nLink+DIM; ++i) {
        //     double sumAB = 0;
        //     for (int j = 0; j < nLink+DIM; ++j) {
        //         sumAB+=A[j*(nLink+DIM)+i] * B[j];
        //     }
        //     x[i] = xold[i] - sumAB;
        // }

        double maxDiff = fabs(xold[0] - x[0]);
        for (int i = 1; i < nLink+DIM; i++) {
            if (fabs(xold[i] - x[i]) > maxDiff) {
                maxDiff = fabs(xold[i] - x[i]);
            }
        }
        if (maxDiff < 1e-6) return;
    }

    // bead->print();
    throw "MaxStep exceeded in Newton!";
}

double** Rod::constraint(double** f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    u = linkVectorU();
    b = linkVectorB();
    double tension[nLink+DIM];
    std::fill(&tension[0], &tension[0] + nLink + DIM, 0);
    // solverPicard(tension);
    solverNewton(tension);

    for (int i = 0; i < nBead; i++) {
        for (int j = 0; j < nLink; j++) {
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
        f[0][k] = f[0][k] + tension[nLink+k]*u[nLink+k][k];
    }
    return f;
}

double** Rod::pseudo(double **f) 
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    u = linkVectorU();
    double metric[nLink*nLink];
    std::fill(&metric[0], &metric[0] + nLink*nLink, 0);
    for (int i = 0; i < nLink; i++) {
        for (int j = 0; j < nLink; j++) {
            if (g[i][j] != 0) {
                metric[j*nLink+i] = -g[i][j]*Dot(&u[i][0],&u[j][0], DIM);
            }
        }
    }

    // calculate the inverse of the metric matrix
    int n = nLink;
    int iPIv[nLink];
    int info;
    double work[nLink];
    dgetrf_(&n, &n, metric, &n, iPIv, &info);
    dgetri_(&n, metric, &n, iPIv, work, &n, &info);

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
                f[k][m] = f[k][m] + g[j][i] * metric[i*nLink+j] * pgr[m];
            }
        }
    }

    return f;
}

double** Rod::pseudoSparse(double **f)
{
    std::fill(&f[0][0], &f[0][0] + nBead * DIM, 0);

    u = linkVectorU();
    SpMatD metric(*gSparse);
    SpMatD uij(*gSparse);
    for (int k = 0; k < uij.outerSize(); ++k) {
        for (SpMatD::InnerIterator it(uij, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            it.valueRef() = Dot(&u[i][0], &u[j][0], DIM);
        }
    }
    metric = -uij.cwiseProduct(*gSparse);

    // calculate the inverse of the metric matrix
    SpMatD I(nLink,nLink);
    I.setIdentity();
    // Eigen::SparseLU<SpMatD> solver;
    Eigen::SimplicialLDLT<SpMatD> solver;
    solver.compute(metric);
    // std::cout << solver.determinant() << std::endl;
    metric = solver.solve(I);
    metric = metric.cwiseProduct(*gSparse);

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
    double coeff[nLink];
    coeff[0] = -Dot(&u[0][0], &u[nLink-1][0], 3);
    double prodCoeff = coeff[0];
    for (int i = 1; i < nLink; i++) {
        coeff[i] = -Dot(&u[i][0], &u[i-1][0], 3);
        prodCoeff *= coeff[i];  
    }
    double det1 = detBandMetric(nLink-2, &coeff[2]);
    double det2 = detBandMetric(nLink-2, &coeff[1]);
    double det3 = detBandMetric(nLink-1, &coeff[1]);
    int sign = 1 - 2*(nLink%2);
    double detG = -coeff[0]*coeff[0]*det1 
        - coeff[nLink-1]*coeff[nLink-1]*det2
        + 2*det3 - 2*sign*prodCoeff; // check

    double metricInv[nLink][nLink];
    std::fill(&metricInv[0][0], &metricInv[0][0] + nLink*nLink, 0);
    double detTop, detBottom, cofG;

    detTop = 1; detBottom = detBandMetric(nLink-2, &coeff[3]);
    cofG = -coeff[1]*detTop*detBottom - sign*prodCoeff/coeff[1];
    metricInv[1][0] = cofG / detG;
    metricInv[0][1] = metricInv[1][0];

    detTop = 1; detBottom = detBandMetric(nLink-4, &coeff[4]);
    det1 = coeff[2]*detTop*detBottom;
    detTop = 1; detBottom = detBandMetric(nLink-3, &coeff[4]);
    det3 = coeff[2]*detTop*detBottom;
    cofG = coeff[0]*coeff[0]*det1 - 2*det3
        - sign*prodCoeff/coeff[2];
    metricInv[2][1] = cofG / detG;
    metricInv[1][2] = metricInv[2][1];

    for (int i = 3; i < nLink-2; ++i) {
        detTop = detBandMetric(i-2, &coeff[2]);
        detBottom = detBandMetric(nLink-i-2, &coeff[i+2]);
        det1 = coeff[i]*detTop*detBottom;
        detTop = detBandMetric(i-1, &coeff[1]);
        detBottom = detBandMetric(nLink-i-3, &coeff[i+2]);
        det2 = coeff[i]*detTop*detBottom;
        detTop = detBandMetric(i-1, &coeff[1]);
        detBottom = detBandMetric(nLink-i-2, &coeff[i+2]); //check
        det3 = coeff[i]*detTop*detBottom;
        cofG = coeff[0]*coeff[0]*det1 
            + coeff[nLink-1]*coeff[nLink-1]*det2 -
            2*det3 - sign*prodCoeff/coeff[i]; // check
        metricInv[i][i-1] = cofG / detG;
        metricInv[i-1][i] = metricInv[i][i-1];
    }

    detTop = detBandMetric(nLink-4, &coeff[2]); detBottom = 1;
    det1 = coeff[nLink-2]*detTop*detBottom;
    detTop = detBandMetric(nLink-3, &coeff[1]); detBottom = 1;
    det3 = coeff[nLink-2]*detTop*detBottom;
    cofG = coeff[0]*coeff[0]*det1 - 2*det3
        - sign*prodCoeff/coeff[nLink-2];
    metricInv[nLink-2][nLink-3] = cofG / detG;
    metricInv[nLink-3][nLink-2] = metricInv[nLink-2][nLink-3];

    detTop = detBandMetric(nLink-2, &coeff[1]); detBottom = 1;
    cofG = -coeff[nLink-1]*detTop*detBottom
        - sign*prodCoeff/coeff[nLink-1];
    metricInv[nLink-1][nLink-2] = cofG / detG;
    metricInv[nLink-2][nLink-1] = metricInv[nLink-1][nLink-2];

    cofG = -coeff[0]*detBandMetric(nLink-2, &coeff[2]) 
        - sign*prodCoeff / coeff[0];
    metricInv[0][nLink-1] = cofG / detG;
    metricInv[nLink-1][0] = metricInv[0][nLink-1];

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
                f[k][m] = f[k][m] + g[j][i]*metricInv[j][i]*pgr[m];
            }
        }
    }

    return f;
}
