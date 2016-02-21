#ifndef SPRING_HPP_QBXE97YK
#define SPRING_HPP_QBXE97YK

#include "simulation.hpp"
#include "force.hpp"

class Spring: protected Force
{
public:
    Spring (Simulation *simu);
    virtual ~Spring ();

    double** harmonic(double** f);
    double** fene(double** f);

private:
    int **link;         // index of pair of beads linked

    void init();
    void outputLinks();
};

#endif /* end of include guard: SPRING_HPP_QBXE97YK */
