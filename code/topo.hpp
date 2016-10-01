#ifndef TOPO_HPP_FJPVDXQX
#define TOPO_HPP_FJPVDXQX


#include "input.hpp"
#include "Eigen/Sparse"

typedef Eigen::SparseMatrix<double> SpMatD;
typedef Eigen::Triplet<int> T;

class Topo
{
public:
    Topo ();
    virtual ~Topo ();

    void setParameter(Input *input);
    void init();
    void initSparse();
    int getNumLink();
    int** getLink();
    int* getNumPair();
    int** getLinkPair();
    int **getMetricTensor();
    SpMatD* getMetricTensorSparse();
    
private:
    // parameters
    int topoType;
    int nBead;
    int nLink;
    int** link;
    int* nPair;         // # of link pair for each bead
    int **linkPair;
    int **g;            // metric tensor
    SpMatD* gSparse;    // metric matrix sparse representation
    int initFlag;

    void setLink();
    void setLinkPair();
    void setMetricTensor();
    void setMetricTensorSparse();
    void printLinks();
    void outputLinks();
    int setNumLink();
    int getTotalLinkPair();

    void ring(int N, int **link);
    void chain(int N, int **link);
    void ringPair(int N, int **link);
    void ringPairWithCentromere(int N, int **link);
    void threeRingPair(int N, int **link);
    void listTopo();
};

#endif /* end of include guard: TOPO_HPP_FJPVDXQX */
