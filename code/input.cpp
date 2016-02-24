#include "input.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

Input::Input() { }
Input::~Input() { }


void Input::getInput(int argc, char* argv[]) 
{
    std::ostringstream fname;
    if (argc > 1) {
        fname << argv[1];
    } else {
        fname << "input.in";
    }
    infname = fname.str();
}

void Input::file() 
{
    std::ifstream infile(infname);
    std::string line;
    while (std::getline(infile, line)) {
        if (line.compare(0,1,"#")) {
            parse(line);
        }
    }
    print();
    infile.close();
}

void Input::parse(std::string line) 
{
    std::istringstream lineStream(line);
    std::string key;
    if (std::getline(lineStream, key, '=')) {
        if (key.compare("projectName")) {
            double value;
            lineStream >> value;
            parameter[key] = value;
        } else {
            lineStream >> projectName;
        } 
    }
}

void Input::print()
{
    std::cout << "Input File: " << infname << std::endl;
    std::cout << "Total Number of parameters: "<< parameter.size() << std::endl;
    std::cout << "    Check Parameters    " << std::endl;
    for (std::map<std::string, double>::iterator
        it=parameter.begin(); it!=parameter.end(); ++it) {
        std::cout << it->first << "=" 
           << it->second << std::endl;
    }
}
