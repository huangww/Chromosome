#include "simulation.hpp"
#include "input.hpp"

Simulation::Simulation() 
{
    input = new Input();
    simuASEP = NULL;
    simuBeadRod = NULL;
    simuBeadSpring = NULL;
    simuSingleFile = NULL;
    
}
Simulation::~Simulation() 
{
    delete input;
    delete simuASEP;
    delete simuBeadRod;
    delete simuBeadSpring;
    delete simuSingleFile;
}


