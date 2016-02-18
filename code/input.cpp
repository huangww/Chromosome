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

    // if (argc >2 && atoi(argv[2])>0) {
    //     parameter->taskID = atoi(argv[2]);
    // }
    
}

void Input::file() 
{
    while (std::getline(infile, line)) {
        if (line.compare(0,1,"#")) {
            parse();
            excute();
        }
    }
}

void Input::parse() 
{
    std::string key;
    std::istringstream lineStream(line);
    if (std::getline(lineStream, key, '=')) {
        std::string value;
        if (std::getline(lineStream, value)) {
            store_line(key, value);
        }
    }
}

void Input::excute() 
{
    if (nline == 1 && key.compare("setup")) {
        simu->init(value);
    } else {
        simu->parameter->setPara(key, value);
}

