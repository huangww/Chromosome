#include <iostream>
#include <cstdlib>
#include "simulation.hpp"
#include "input.hpp"

int main(int argc, char *argv[])
{
    int startTime = time(NULL);

    try{  
        Simulation *simu = new Simulation;

        simu->input->getInput(argc, argv);
        simu->input->file();

        std::cout << "Start The Simulation: " << std::endl;
        simu->init();
        simu->run();

        delete simu;

    } catch (const char* error) {

        std::cout << error << std::endl;
        return EXIT_FAILURE;

    }

    int endTime = time(NULL);
    int elapsedTime = endTime - startTime;
    std::cout << "Runing Time: " << elapsedTime 
    << "seconds" << std::endl;
    return EXIT_SUCCESS;
}
