#include "simulation.hpp"
#include "particle.hpp"
#include "bead.hpp"
#include "state.hpp"

Simulation::Simulation() 
{
    particle = NULL;
    bead = NULL;
    state = NULL;
}
Simulation::~Simulation() 
{
    delete particle;
    delete bead;
    delete state;
}

