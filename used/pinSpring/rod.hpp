#ifndef ROD_HPP_FFKADSOB
#define ROD_HPP_FFKADSOB

#include "force.hpp"

class Rod: protected Force
{
public:
    Rod (Simulation *simu);
    virtual ~Rod ();

    double** constraint(double** f);
    double** pseudo(double** f);

private:
    int **link;         // index of pair of beads linked
    int **g;            // metric matrix
    double **u;         // unit vector of rod
    double **b;         // vector of rod unnormalized

    void init();
    void printLinks();
    void outputLinks();
    void printMetric();
    double** linkVectorU(); 
    double** linkVectorB();
    int** metricTensor();
    void matrixA(double *A);
    void vectorB(double *x, double *B);
    void solverPicard(double *x);
};

extern "C" int dgetrf_(int* m, int* n, double* A, 
        int* lda, int* iPIv, int* info);
extern "C" int dgetri_(int* m, double* a, int* lda,
        int* iPIv, double* work, int* lwork, int* info);
extern "C" int dgetrs_(char* s, int* n, int* nrhs, 
        double* A, int* lda, int* iPIv, 
        double* B, int* ldb, int* info);

#endif /* end of include guard: ROD_HPP_FFKADSOB */
