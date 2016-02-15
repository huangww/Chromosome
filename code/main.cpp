#include <iostream>
#include <cstdlib>
#include "test.hpp"
#include "simulation.hpp"

int main(int argc, char *argv[])
{
    int startTime = time(NULL);

    Simulation *simu = new Simulation();

    simu->input->getArg(argc, argv);
    simu->input->file();

    // test();
    simu->print();
    simu->run();

    delete simu;

    int endTime = time(NULL);
    int elapsedTime = endTime - startTime;
    std::cout << "Runing Time: " << elapsedTime 
    << "seconds" << std::endl;
    return 0;
}
