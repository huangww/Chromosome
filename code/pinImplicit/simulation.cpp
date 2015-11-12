#include "simulation.hpp"
#include "particle.hpp"
#include "bead.hpp"

Simulation::Simulation() 
{
    particle = NULL;
    bead = NULL;
}
Simulation::~Simulation() 
{
    delete particle;
    delete bead;
}

