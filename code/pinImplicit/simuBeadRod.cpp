#include "simuBeadRod.hpp"
#include "bead.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

SimuBeadRod::SimuBeadRod() : Simulation()
{
    bead = new Bead(this);
}
SimuBeadRod::~SimuBeadRod() { }

void SimuBeadRod::print() 
{
    std::cout << "================================="
        << std::endl;
    std::cout << "Parameters of the Simulation System" 
        << std::endl;
    std::cout << "=> Number of Bead: "
       << bead->nBead << std::endl;
    std::cout << "=> Number of Rod: "
       << bead->nRod << std::endl;
    std::cout << "=> Topological Type: "
       << bead->topoType << std::endl;
    std::cout << "=> Time Step: "
       << bead->dt << std::endl;
    std::cout << "=> Total Evolving Time: "
       << bead->tEnd << std::endl;
    std::cout << "=> Output Time Interval: "
       << bead->dt*bead->outputStep << std::endl;
    std::cout << "=> Effective Temperature: "
       << bead->tempEff << std::endl;
    std::cout << "=> Random Seed: "
       << bead->seed << std::endl;
    std::cout << "=> Task ID: "
       << bead->taskID << std::endl;
    std::cout << "================================="
        << std::endl;
}

void SimuBeadRod::run() 
{
    std::stringstream fname;
    fname << "data/r_N" << bead->nBead 
        << "_T" << bead->tempEff 
        << "_" << bead->taskID << ".dat";
    std::cout << fname.str() << std::endl;
    std::ofstream output(fname.str());

    bead->init();

    double t = 0;
    int maxStep = int(bead->tEnd / bead->dt);
    for (int step = 0; step < maxStep; ++step) {
        // output bead position to data file
        // if (step % int(1.0/bead->dt) == 0) {
        if (step % bead->outputStep == 0) {
            bead->output(output);
        }

        // output progressing to screen
        if (step % (maxStep/100) == 0) {
           std::cout << step / (maxStep/100) <<" % done!" 
               << std::endl; 
        }

        // Monte-Carlo Move
        // bead->montecarloUpdate();

        // MD Move
        bead->predict();
        bead->correct();
        
        t += bead->dt;
    }

    output.close();
}
