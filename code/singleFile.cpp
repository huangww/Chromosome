#include "singleFile.hpp"
#include "project.hpp"
#include "input.hpp"
#include "particle.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

SingleFile::SingleFile() : Project()
{
    particle = NULL;
}
SingleFile::~SingleFile() 
{
    delete particle;
    outFile.close();
}

void SingleFile::setup(Input* input) 
{
    if (input->parameter.count("tEnd") == 0) {
        throw "Parameter \"tEnd\" is not specified!";
    }
    tEnd = input->parameter["tEnd"];
    if (input->parameter.count("tEnd") == 0) {
        outputStep = 1;
    }
    outputStep = int(input->parameter["outputStep"]);

    particle = new Particle();
    particle->setParameter(input);

    // open files for output
    std::stringstream fname;
    fname << "data/" << input->projectName <<
        "_pos_N" << input->parameter["nSite"] 
        << "_T" << input->parameter["tempEff"] 
        << "_" << input->parameter["taskID"] << ".dat";
    std::cout << fname.str() << std::endl;
    outFile.open(fname.str());
}

void SingleFile::run() 
{
    
    particle->init();

    int maxStep = int(tEnd / particle->dt);
    for (int step = 0; step < maxStep; ++step) {
        // output to data file
        // if (step % int(1.0/particle->dt) == 0) {
        if (step % outputStep == 0) {
            particle->output(outFile);
        }

        // output progressing to screen
        if (step % (maxStep/100) == 0) {
           std::cout << step / (maxStep/100) <<" % done!" 
               << std::endl; 
        }

        particle->updateBD();
    }
}
