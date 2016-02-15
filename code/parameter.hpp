#ifndef PARAMETER_HPP_H4YIFEO2
#define PARAMETER_HPP_H4YIFEO2

#include "simulation.hpp"


class Parameter
{
public:
    Parameter (Simulation *simu);
    virtual ~Parameter ();

    int taskID;

    int nSite;
    int nPar;
    int nSample;

    double rateToLeft;
    double rateToRight;

    int topoType;
    int nBead;
    int nRod;

    double dt;
    double tEnd;
    int outputStep;
    double tempEff;
    unsigned long seed;

protected:
    // some constant
    const double PI;
    const int DIM;

    Simulation *simulation;
    Particle *&particle;
    Bead *&bead;

private:
    void setup();
    void setupASEP();
    void setupSingleFile();
    void setupBeadRod();
    void setupBeadSpring();
};

#endif /* end of include guard: PARAMETER_HPP_H4YIFEO2 */
