#include "simulation.hpp"
#include "input.hpp"
#include "project.hpp"
#include "asep.hpp"
#include "beadRod.hpp"
#include "beadSpring.hpp"
#include "singleFile.hpp"

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
    if (input->projectName=="Asep") 
        project = new Asep;

    if (input->projectName=="BeadRod") 
        project = new BeadRod;
    
    if (input->projectName=="BeadSpring") 
        project = new BeadSpring;
    
    if (input->projectName=="SingleFile") 
        project = new SingleFile;

    input->print();

    if (project==NULL) {
        throw "Wrong projet Name!";
     }
}

void Simulation::run() 
{
    project->setup(input);
    project->run();
}


