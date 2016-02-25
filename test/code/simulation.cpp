#include "simulation.hpp"
#include "input.hpp"
#include "project.hpp"
#include "asep.hpp"

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
    if (input->projectName=="Asep") {
        project = new Asep;
    } else {
        project = new Asep;
    }
    
    input->print();
}

void Simulation::run() 
{

    project->setup(input);
    project->run();
}


