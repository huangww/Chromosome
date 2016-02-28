#ifndef TOPO_HPP_FJPVDXQX
#define TOPO_HPP_FJPVDXQX


#include "input.hpp"

class Topo
{
public:
    Topo ();
    virtual ~Topo ();

    void setParameter(Input *input);
    int** init(int **link);

private:
    // parameters
    int topoType;
    int nLink;
    int nBead;

    void ring(int N, int **link);
    void chain(int N, int **link);
    void ringPair(int N, int **link);
    void ringPairWithCentromere(int N, int **link);
    void threeRingPair(int N, int **link);
};

#endif /* end of include guard: TOPO_HPP_FJPVDXQX */
