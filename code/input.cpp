#include "input.hpp"
#include "simulation.hpp"
#include "parameter.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

Input::Input(Simulation *simu) : Parameter(simu) 
{
        

}
Input::~Input() 
{

}


void Input::getArg(int argc, char* argv[]) 
{
    if (argc > 1) {
        infile.open(argv[1]);
    } 

    if (argc >2 && atoi(argv[2])>0) {
        parameter->taskID = atoi(argv[2]);
    }
    
}

void Input::file() 
{
    std::string line;
    while (1) {
        // read a line 
        std::getline(infile, line);
        
        
    }
}

