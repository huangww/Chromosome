#ifndef ROD_HPP_FFKADSOB
#define ROD_HPP_FFKADSOB

#include "force.hpp"
#include "Eigen/Sparse"

class Rod: protected Force
{
public:
    Rod (Simulation *simu);
    virtual ~Rod ();

    double** constraint(double** f);
    double** pseudo(double** f);
    double** pseudoRing(double** f);
    double** pseudoSparse(double **f);

private:
    typedef Eigen::SparseMatrix<double> SpMatD;
    typedef Eigen::Triplet<int> T;

    int **link;         // index of pair of beads linked
    double **u;         // unit vector of rod
    double **b;         // vector of rod unnormalized
    int **g;            // metric matrix
    SpMatD gSparse;     // metric matrix sparse representation

    struct LinkTable {
        int *nLinks;
        int **table;
    } linkTable;

    void init();
    void printLinks();
    void outputLinks();
    void printMetric();
    void setLinkTable();
    double** linkVectorU(); 
    double** linkVectorB();
    int** metricTensor();
    void metricTensorSparse();
    void matrixA(double *A);
    void vectorB(double *x, double *B);
    void solverPicard(double *x);
    double detBandMetric(int n, double *coeff);
};


extern "C" int dgetrf_(int* m, int* n, double* A, 
        int* lda, int* iPIv, int* info);
extern "C" int dgetri_(int* m, double* a, int* lda,
        int* iPIv, double* work, int* lwork, int* info);
extern "C" int dgetrs_(char* s, int* n, int* nrhs, 
        double* A, int* lda, int* iPIv, 
        double* B, int* ldb, int* info);

#endif /* end of include guard: ROD_HPP_FFKADSOB */
