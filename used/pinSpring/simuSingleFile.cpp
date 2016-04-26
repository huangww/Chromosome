#include "simuSingleFile.hpp"
#include "particle.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

SimuSingleFile::SimuSingleFile() : Simulation()
{
    particle = new Particle(this);
}
SimuSingleFile::~SimuSingleFile() { }


void SimuSingleFile::print() 
{
    std::cout << "================================="
        << std::endl;
    std::cout << "Parameters of the Simulation System" 
        << std::endl;
    std::cout << "=> Number of Sites: "
       << particle->nSite << std::endl;
    std::cout << "=> Number of Particles: "
       << particle->nPar << std::endl;
    std::cout << "=> Time step: "
       << particle->dt << std::endl;
    std::cout << "=> Total evolving time: "
       << particle->tEnd << std::endl;
    std::cout << "=> Effective temperature: "
       << particle->tempEff << std::endl;
    std::cout << "================================="
        << std::endl;
}

void SimuSingleFile::run() 
{
    std::stringstream fname;
    fname << "data/par_N" << particle->nSite << "_T"
        << particle->tempEff << ".dat";
    std::cout << fname.str() << std::endl;
    std::ofstream output(fname.str());

    particle->init();

    int maxStep = int(particle->tEnd / particle->dt);
    for (int step = 0; step < maxStep; ++step) {
        // output to data file
        // if (step % int(1.0/particle->dt) == 0) {
        if (step % particle->outputStep == 0) {
            particle->output(output);
        }

        // output progressing to screen
        if (step % (maxStep/100) == 0) {
           std::cout << step / (maxStep/100) <<" % done!" 
               << std::endl; 
        }

        particle->updateBD();
    }

    output.close();
}
