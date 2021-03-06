#ifndef SIMUASEP_HPP_AUILPX1S
#define SIMUASEP_HPP_AUILPX1S

#include "simulation.hpp"

class SimuASEP: public Simulation
{
public:
    SimuASEP ();
    virtual ~SimuASEP ();

    void print();
    void run();
private:
    class State *state;
};

#endif /* end of include guard: SIMUASEP_HPP_AUILPX1S */
