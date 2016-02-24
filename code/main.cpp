#include <iostream>
#include <cstdlib>
#include "simulation.hpp"
#include "input.hpp"

int main(int argc, char *argv[])
{
    int startTime = time(NULL);
   
    Simulation *simu = new Simulation;

    simu->input->getInput(argc, argv);
    simu->input->file();

    simu->run();

    delete simu;

    int endTime = time(NULL);
    int elapsedTime = endTime - startTime;
    std::cout << "Runing Time: " << elapsedTime 
    << "seconds" << std::endl;
    return 0;
}
