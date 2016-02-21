#include "parameter.hpp"
#include "simulation.hpp"
#include <iostream>
#include <sstream>
#include <random>

Parameter::Parameter(Simulation *simu) :
    PI (3.141592653589793238463),
    DIM (3),
    simulation(simu),
{ }
Parameter::~Parameter() 
{ }


void Parameter::set(std::string key, std::string value)
{

}
