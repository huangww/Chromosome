#ifndef FORCE_HPP_WU0ZVG8B
#define FORCE_HPP_WU0ZVG8B


#include "input.hpp"

class Force
{
public:
    Force ();
    virtual ~Force ();

    void setParameter(Input *input);
    void print(double* f);
    void print(double** f);
    double* repulsive(double* x, double* f);
    double** repulsive(double** r, double** f);
    double* brownian(double* f);
    double** brownian(double** f);
    double* external(double* f);
    double** external(double** f);
    double* boundary(double* x, double* f);

private:
    int nBead;
    int nPar;
    int nSite;
    unsigned long seed;
    double tempEff;
    double dt;

    /* data */
};

#endif /* end of include guard: FORCE_HPP_WU0ZVG8B */
