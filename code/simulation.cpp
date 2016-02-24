#include "simulation.hpp"
#include "input.hpp"
#include "project.hpp"

Simulation::Simulation() 
{
    input = new Input();
    project = NULL;
    
}
Simulation::~Simulation() 
{
    delete input;
    delete project;
}

void Simulation::init() 
{
    // input->parameter->paraName[0] ==
    
}


