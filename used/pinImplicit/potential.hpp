#ifndef POTENTIAL_HPP_RJ98ISVS
#define POTENTIAL_HPP_RJ98ISVS

#include "parameter.hpp"
#include "simulation.hpp"

class Potential: protected Parameter
{
public:
    Potential (Simulation *simu);
    virtual ~Potential ();
    
    double LennardJones();

};

#endif /* end of include guard: POTENTIAL_HPP_RJ98ISVS */
