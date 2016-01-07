#include <iostream>
#include <cstdlib>
#include "simuBeadRod.hpp"
#include "simuSingleFile.hpp"
#include "simuASEP.hpp"
#include "bead.hpp"
#include "test.hpp"

int main(int argc, char *argv[])
{
    int startTime = time(NULL);

    Simulation *simu = new SimuASEP();

    if (argc > 1 && atoi(argv[1])>0) {
        // simu->bead->taskID = atoi(argv[1]);
    }

    simu->print();
    // test();
    simu->run();

    delete simu;

    int endTime = time(NULL);
    int elapsedTime = endTime - startTime;
    std::cout << "Runing Time: " << elapsedTime 
    << "seconds" << std::endl;
    return 0;
}
