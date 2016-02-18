#include "simulation.hpp"
#include "input.hpp"
#include "output.hpp"

Simulation::Simulation() 
{
    input = new Input();
    output = new Output();
}
Simulation::~Simulation() 
{
    delete input;
    delete output;
}

