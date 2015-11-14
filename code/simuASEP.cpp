#include "simuASEP.hpp"
#include "simulation.hpp"
#include "state.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

SimuASEP::SimuASEP() : Simulation()
{
   state = new State(this);
}
SimuASEP::~SimuASEP() 
{
    delete state;
}

void SimuASEP::print() 
{
    std::cout << "================================="
        << std::endl;
    std::cout << "Parameters of the Simulation System" 
        << std::endl;
    std::cout << "=> Number of Sites: "
       << state->nSite << std::endl;
    std::cout << "=> Number of states: "
       << state->nPar << std::endl;
    std::cout << "=> Time step: "
       << state->dt << std::endl;
    std::cout << "=> Total evolving time: "
       << state->tEnd << std::endl;
    std::cout << "=> Effective temperature: "
       << state->tempEff << std::endl;
    std::cout << "================================="
        << std::endl;
}

void SimuASEP::run() 
{
    std::stringstream fname;
    std::ofstream *output = new std::ofstream[2];
    fname << "data/par_N" << state->nSite << "_T"
        << state->tempEff << ".dat";
    std::cout << fname.str() << std::endl;
    output[0].open(fname.str());
    fname.str("");
    fname << "data/rg1D_N" << state->nSite << "_T"
        << state->tempEff << ".dat";
    std::cout << fname.str() << std::endl;
    output[1].open(fname.str());
    
    double jumpTime = 0;
    for (int i = 0; i < state->nSample; ++i) {

        state->init();

        int step = 0;
        while (state->t < state->tEnd) {
            // output to data file
            // if (step % int(1.0/state->dt) == 0) {
            if (step % state->outputStep == 0) {
            state->output(output);
            }

            state->update();
            step++;
        }

        // output progressing to screen
        jumpTime += state->tEnd/step;
        std::cout << i <<  std::endl;
    }
    std::cout << jumpTime/state->nSample << std::endl;

    output[0].close();
    output[1].close();
}
