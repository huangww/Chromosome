#include "beadSpring.hpp"
#include "project.hpp"
#include "input.hpp"
#include "bead.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

BeadSpring::BeadSpring() : Project()
{
    bead = NULL;
}
BeadSpring::~BeadSpring() 
{
    delete bead;
    outFile[0].close();
    outFile[1].close();
    outFile[2].close();
    delete[] outFile;
}


void BeadSpring::setup(Input *input) 
{
    outputStep = int(input->parameter["outputStep"]);
    tEnd = input->parameter["tEnd"];

    bead = new Bead();
    bead->setParameter(input);

    // open files for output
    std::stringstream fname;
    fname << "data/r_N" << input->parameter["nBead"]
        << "_T" << input->parameter["tempEff"]
        << "_" << input->parameter["taskID"] << ".dat";
    std::cout << fname.str() << std::endl;
    outFile[0].open(fname.str());
    fname.str("");
    fname << "data/rg_N" << input->parameter["nBead"]
        << "_T" << input->parameter["tempEff"]
        << "_" << input->parameter["taskID"] << ".dat";
    std::cout << fname.str() << std::endl;
    outFile[1].open(fname.str());
    fname.str("");
    fname << "data/rd_N" << input->parameter["nBead"]
        << "_T" << input->parameter["tempEff"]
        << "_" << input->parameter["taskID"] << ".dat";
    std::cout << fname.str() << std::endl;
    outFile[2].open(fname.str());
}
    

void BeadSpring::run() 
{
    bead->init();

    double t = 0;
    int maxStep = int(tEnd / bead->dt);
    for (int step = 0; step < maxStep; ++step) {
        // output bead position to data file
        // if (step % int(1.0/bead->dt) == 0) {
        if (step % outputStep == 0) {
            bead->output(outFile);
        }

        // output progressing to screen
        if (step % (maxStep/100) == 0) {
           std::cout << step / (maxStep/100) <<" % done!" 
               << std::endl; 
        }

        // MD Move
        bead->update();
        
        t += bead->dt;
    }

}
