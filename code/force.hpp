#ifndef FORCE_HPP_WU0ZVG8B
#define FORCE_HPP_WU0ZVG8B

#include "parameter.hpp"
#include "simulation.hpp"

class Force: protected Parameter
{
public:
    Force (Simulation *simu);
    virtual ~Force ();

    void print(double* f);
    void print(double** f);
    double* repulsive(double* f);
    double** repulsive(double** f);
    double* brownian(double* f);
    double** brownian(double** f);
    double* external(double* f);
    double** external(double** f);
    double* boundary(double* f);
protected:
    /* data */
private:
    /* data */
};

#endif /* end of include guard: FORCE_HPP_WU0ZVG8B */
