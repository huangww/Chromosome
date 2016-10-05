#ifndef ROD_HPP_FFKADSOB
#define ROD_HPP_FFKADSOB

#include "input.hpp"
#include "bead.hpp"


class Rod
{
public:
    Rod (Bead *beadPointer);
    virtual ~Rod ();

    void setParameter(Input *input);
    double** constraint(double** f);
    double** pseudo(double** f);
    double** pseudoRing(double** f);
    double** pseudoSparse(double **f);

private:

    int **link;         // index of pair of beads linked
    double **u;         // unit vector of rod
    double **b;         // vector of rod unnormalized
    int **g;            // metric matrix, first bead pinned
    SpMatD* gSparse;    // metric matrix sparse representation

    int *nPair;         // # of link pairs related to each bead
    int **linkPair;     // table of link pairs

    class Bead *bead;
    class Topo *topo;

    void init();
    void printMetric();
    double** linkVectorU(); 
    double** linkVectorB();
    void matrixA(double *A);
    void vectorB(double *x, double *B);
    void solverPicard(double *x);
    void jacobian(double *x, double *A, double *B);
    void inverse(double *A, int N);
    void solverNewton(double *x);
    void solverHydrj(double *x);
    double detBandMetric(int n, double *coeff);

    // parameters
    int nLink;
    int nBead;
    double dt;

};


extern "C" int dgetrf_(int* m, int* n, double* A, 
        int* lda, int* iPIv, int* info);
extern "C" int dgetri_(int* m, double* a, int* lda,
        int* iPIv, double* work, int* lwork, int* info);
extern "C" int dgetrs_(char* s, int* n, int* nrhs, 
        double* A, int* lda, int* iPIv, 
        double* B, int* ldb, int* info);

#endif /* end of include guard: ROD_HPP_FFKADSOB */
