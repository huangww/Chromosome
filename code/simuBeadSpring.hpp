#ifndef SIMUBEADSPRING_HPP_CPGHOHFE
#define SIMUBEADSPRING_HPP_CPGHOHFE

#include "simulation.hpp"

class SimuBeadSpring: public Simulation
{
public:
    SimuBeadSpring ();
    virtual ~SimuBeadSpring ();

    void print();
    void run();
};

#endif /* end of include guard: SIMUBEADSPRING_HPP_CPGHOHFE */
