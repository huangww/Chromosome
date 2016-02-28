#ifndef POTENTIAL_HPP_RJ98ISVS
#define POTENTIAL_HPP_RJ98ISVS

#include "input.hpp"

class Potential
{
public:
    Potential ();
    virtual ~Potential ();

    double LennardJones(int N, double **r);


};

#endif /* end of include guard: POTENTIAL_HPP_RJ98ISVS */
