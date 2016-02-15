#ifndef SIMUSINGLEFILE_HPP_RUAETXWZ
#define SIMUSINGLEFILE_HPP_RUAETXWZ

#include "simulation.hpp"
class SimuSingleFile: public Simulation
{
public:
    SimuSingleFile ();
    virtual ~SimuSingleFile ();

    void print();
    void run();
};

#endif /* end of include guard: SIMUSINGLEFILE_HPP_RUAETXWZ */
