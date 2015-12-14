#include "parameter.hpp"
#include "simulation.hpp"
#include <iostream>
#include <random>

Parameter::Parameter(Simulation *simu) :
    PI (3.141592653589793238463),
    DIM (3),
    simulation(simu),
    particle(simu->particle),
    bead(simu->bead)
{
    setup();
}
Parameter::~Parameter() 
{ }

void Parameter::setup() 
{
    // setupASEP();
    setupBeadRod();
    // setupSingleFile();
}

void Parameter::setupASEP() 
{
    nSite = 100;
    nPar = 50;
    nSample = 1;
    tempEff = 200;

    dt = 5;
    tEnd = 1.e+7;
    outputStep = 1;
    double jumpRate = 1.0;
    double factor = exp(-1.0/tempEff);
    // factor = 1.0;
    rateToLeft = jumpRate*1.0/(1+factor);
    rateToRight = jumpRate*factor/(1.0+factor);

    std::random_device rd;
    seed = rd();
    // seed = 5489;
}

void Parameter::setupSingleFile()
{
    nSite = 100;
    nPar = 50;
    nSample = 1000;
    tempEff = 1.0;

    dt = 1e-4;
    tEnd = 1.e+2;
    outputStep = 1e3;
    std::random_device rd;
    // seed = rd();
    seed = 5489;
}

void Parameter::setupBeadRod() 
{
    topoType = 0;
    nBead = 100;
    nRod = 100;
    tempEff = 1;

    dt = 1e-4;
    tEnd = 1e+1;
    outputStep = 1e3;

    std::random_device rd;
    // seed = rd();
    seed = 5489;
}
