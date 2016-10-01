#ifndef SPRING_HPP_QBXE97YK
#define SPRING_HPP_QBXE97YK

#include "input.hpp"
#include "bead.hpp"

class Spring
{
public:
    Spring (Bead *beadPointer);
    virtual ~Spring ();

    void setParameter(Input* input);
    double** bending(double** r, double** f);
    double** harmonic(double** r, double** f);
    double** fene(double** r, double** f);

private:
    // parameters
    int nBead;
    int nLink;

    int** link;         // index of pair of beads linked
    int* nPair;
    int** linkPair;
    int** g;            // metric tensor, linking matrix
    double** u;         // unit vector of links

    class Bead *bead;
    class Topo *topo;

    // functions
    double** linkVector(double** r);
};

#endif /* end of include guard: SPRING_HPP_QBXE97YK */
