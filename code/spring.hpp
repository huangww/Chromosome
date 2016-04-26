#ifndef SPRING_HPP_QBXE97YK
#define SPRING_HPP_QBXE97YK

#include "input.hpp"

class Spring
{
public:
    Spring ();
    virtual ~Spring ();

    void setParameter(Input* input);
    double** harmonic(double** r, double** f);
    double** fene(double** r, double** f);

private:
    // parameters
    int nBead;
    int nLink;

    int** link;         // index of pair of beads linked

    void outputLinks();
};

#endif /* end of include guard: SPRING_HPP_QBXE97YK */
