#include "simuASEP.hpp"
#include "simulation.hpp"
#include "state.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

SimuASEP::SimuASEP() : Simulation()
{
   state = new State(this);
}
SimuASEP::~SimuASEP() { }

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
        << state->tempEff << "_eq.dat";
    std::cout << fname.str() << std::endl;
    output[1].open(fname.str());
    
    int nRg = int(state->tEnd / state->dt);
    double *rgMean = new double[nRg];
    std::fill(&rgMean[0], &rgMean[0] + nRg, 0);
   
    for (int i = 0; i < state->nSample; ++i) {

        state->init();

        int step = 0;
        while (state->t < state->tEnd) {
            // output to data file
            // if (step % state->outputStep == 0) {
            //     state->output(output);
            // }

            rgMean[step] += state->rg;
            state->update();
            step++;
        }

        // output progressing to screen
        std::cout << i <<  std::endl;
    }

    for (int i = 0; i < nRg; ++i) {
        output[1] << std::setprecision(9) << std::setw(10);
        output[1] << i*state->dt << '\t';
        output[1] << rgMean[i] / state->nSample 
            << std::endl;
    }
    delete [] rgMean;
    output[0].close();
    output[1].close();
}

