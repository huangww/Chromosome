#include "asep.hpp"
#include "project.hpp"
#include "state.hpp"
#include "input.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

Asep::Asep() : Project()
{
    state = new State();
    outFile = new std::ofstream[2];  
}
Asep::~Asep() 
{ 
    delete state;
    outFile[0].close();
    outFile[1].close();
    delete[] outFile;
}

void Asep::setup(Input *input) 
{
    if (input->parameter.count("nSample") == 0) {
        nSample = 1;
    }
    nSample = int(input->parameter["nSample"]);
    if (input->parameter.count("outputStep") == 0) {
        outputStep = 1;
    }
    outputStep = int(input->parameter["outputStep"]);
    if (input->parameter.count("tEnd") == 0) {
        throw "Parameter \"tEnd\" is not specified!";
    }
    tEnd = input->parameter["tEnd"];
    state->setParameter(input);

    // Open files for output
    std::ostringstream fname;
    fname << "data/" << input->projectName <<
        "_pos_N" << input->parameter["nSite"] 
        << "_T" << input->parameter["tempEff"] 
        << "_" << input->parameter["taskID"] << ".dat";
    outFile[0].open(fname.str());
    std::cout << fname.str() << std::endl;
    fname.str("");
    fname << "data/" << input->projectName <<
        "_rg_N" << input->parameter["nSite"] 
        << "_T" << input->parameter["tempEff"] 
        << "_" << input->parameter["taskID"] << ".dat";
    std::cout << fname.str() << std::endl;
    outFile[1].open(fname.str());
    outFile[1] << "#\t" << "t\t\t" << "Rg\t\t" 
        << "r[N/2]" << std::endl;
}

void Asep::run() 
{
    // int nRg = int(tEnd / state->dt);
    // double *rgMean = new double[nRg];
    // std::fill(&rgMean[0], &rgMean[0] + nRg, 0);
   
    for (int i = 0; i < nSample; ++i) {

        state->init();

        int step = 0;
        while (state->t < tEnd) {
            // output to data file
            if (step % outputStep == 0) {
                state->output(outFile);
                // this->output();
            }

            // rgMean[step] += state->rg;
            state->update();
            step++;
        }

        // output progressing to screen
        std::cout << i <<  std::endl;
    }

    // for (int i = 0; i < nRg; ++i) {
    //     output[1] << std::setprecision(9) << std::setw(10);
    //     output[1] << i*state->dt << '\t';
    //     output[1] << rgMean[i] / state->nSample 
    //         << std::endl;
    // }
    // delete [] rgMean;
}

