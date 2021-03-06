#ifndef SIMUBEADROD_HPP_TGEVUJIM
#define SIMUBEADROD_HPP_TGEVUJIM

#include "simulation.hpp"

class SimuBeadRod: public Simulation
{
public:
    SimuBeadRod ();
    virtual ~SimuBeadRod ();

    void print();
    void run();
private:
    void runMD();
    void runMonteCarlo();
};

#endif /* end of include guard: SIMUBEADROD_HPP_TGEVUJIM */
