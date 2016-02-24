#include "input.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

Input::Input()
{
    parameter = new Parameter;
}
Input::~Input() 
{
    delete [] parameter->paraName;
    delete [] parameter->paraValue;
    delete parameter;
}


void Input::getInput(int argc, char* argv[]) 
{
    std::ostringstream fname;
    if (argc > 1) {
        fname << argv[1];
    } else {
        fname << "input.in";
    }
    inputfname = fname.str();
}

int Input::getParaNumber()
{
    std::ifstream infile(inputfname);
    nPara = 1;
    std::string line;
    while (std::getline(infile, line)) {
        if (line.compare(0,1,"#")) {
            nPara++;
        }
    }
    return nPara;
}

void Input::init()
{
    nPara = getParaNumber();
    parameter->paraName = new std::string[nPara];
    parameter->paraValue = new double[nPara];
}

void Input::file() 
{
    init();
    std::ifstream infile(inputfname);
    std::string line;
    int index = 0;
    while (std::getline(infile, line)) {
        if (line.compare(0,1,"#")) {
            parse(index, line);
            index++;
        }
    }
    print();
}

void Input::parse(int index, std::string line) 
{
    std::istringstream lineStream(line);
    std::string key;
    if (std::getline(lineStream, key, '=')) {
        parameter->paraName[index] = key;
        lineStream >> parameter->paraValue[index];
    }
}

void Input::print()
{

    std::cout << "Input File: " << inputfname << std::endl;
    std::cout << "Total Number of parameters: "<< nPara << std::endl;
    std::cout << "    Check Parameters    " << std::endl;
    for (int i = 0; i < nPara; ++i) {
        std::cout << parameter->paraName[i] << "=" 
           << parameter->paraValue[i] << std::endl;
    }
}
