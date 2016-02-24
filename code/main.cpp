#include <iostream>
#include <cstdlib>
#include "input.hpp"

int main(int argc, char *argv[])
{
    // int startTime = time(NULL);
   
    Input *input = new Input;
    input->getInput(argc, argv);
    input->file();

    delete input;

    // int endTime = time(NULL);
    // int elapsedTime = endTime - startTime;
    // std::cout << "Runing Time: " << elapsedTime 
    // << "seconds" << std::endl;
    return 0;
}
