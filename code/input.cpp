#include "input.hpp"
#include "simulation.hpp"
#include "parameter.hpp"
#include <string>
#include <iostream>

Input::Input(Simulation *simu) : Parameter(simu) 
{
        

}
Input::~Input() 
{

}



void Input::file() 
{
    std::string line;
    while (1) {
        // read a line 
        std::getline(infile, line);
        
        
    }
}
