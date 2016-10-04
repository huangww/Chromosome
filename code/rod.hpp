#ifndef ROD_HPP_FFKADSOB
#define ROD_HPP_FFKADSOB

#include "Eigen/Sparse"
#include "input.hpp"
#include "bead.hpp"
#include "Eigen/unsupported/Eigen/NonLinearOptimization"

typedef Eigen::SparseMatrix<double> SpMatD;

template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  const int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  // you should define that in the subclass :
//  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};


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
