#ifndef TOPO_HPP_FJPVDXQX
#define TOPO_HPP_FJPVDXQX

#include "simulation.hpp"
#include "parameter.hpp"

class Topo: protected Parameter
{
public:
    Topo (Simulation *simu);
    virtual ~Topo ();

    int** init(int **link);

private:
    void ring(int N, int **link);
    void chain(int N, int **link);
    void ringPair(int N, int **link);
    void ringPairWithCentromere(int N, int **link);
    void threeRingPair(int N, int **link);
};

#endif /* end of include guard: TOPO_HPP_FJPVDXQX */
