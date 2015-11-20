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
        << state->tempEff << "_eq.dat";
    std::cout << fname.str() << std::endl;
    output[1].open(fname.str());
    
    // const int nbins = int(state->tEnd);
    // int *binCount = new int[nbins];
    // double *rgMean = new double[nbins];
    // std::fill(&binCount[0], &binCount[0]+nbins, 0);
    // std::fill(&rgMean[0], &rgMean[0]+nbins, 0);

    for (int i = 0; i < state->nSample; ++i) {

        state->init();

        int step = 0;
        // int i0 = 0;
        while (state->t < state->tEnd) {
            // output to data file
            if (step % state->outputStep == 0) {
                state->output(output);
            }

            state->update();
            // i0 = binRg(i0, nbins, binCount, rgMean);
            step++;
        }

        // output progressing to screen
        std::cout << i <<  std::endl;
    }

    // for (int i = 0; i < nbins; ++i) {
    //     output[1] << float(i)*state->tEnd/nbins << '\t';
    //     output[1] << rgMean[i] / binCount[i] << std::endl;
    // }

    // delete [] binCount;
    // delete [] rgMean;
    output[0].close();
    output[1].close();
}

int SimuASEP::binRg(int i0, int nbins, int* binCount, double* rgMean) 
{
    for (int i = i0; i < nbins; ++i) {
        bool inBin;
        inBin = (state->t > i*state->tEnd/nbins) && 
            (state->t < (i+1)*state->tEnd/nbins); 
        if (inBin) {
            binCount[i]++;
            rgMean[i] += state->rg;
            return i;
        }
    }
    return nbins;
}
