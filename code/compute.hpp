#ifndef COMPUTE_HPP_6KHWKDTL
#define COMPUTE_HPP_6KHWKDTL

#include "parameter.hpp"
#include "simulation.hpp"

class Compute: protected Parameter
{
public:
    Compute (Simulation *simu);
    virtual ~Compute ();
    double gyrationRadius(double **r);

protected:
    /* data */
private:
    /* data */
};


#endif /* end of include guard: COMPUTE_HPP_6KHWKDTL */
