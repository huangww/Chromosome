#include "input.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>

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
        line.erase(std::remove_if(line.begin(), line.end(),
                    isspace), line.end());
        if (line.compare(0,1,"#")) {
            parse(line);
        }
    }
    infile.close();
}

void Input::parse(std::string line) 
{
    std::istringstream lineStream(line);
    std::string key;
    if (std::getline(lineStream, key, '=')) {
        if (key=="projectName") {
            lineStream >> projectName;
        } else {
            double value;
            lineStream >> value;
            parameter[key] = value;
        } 
    }
}

void Input::print()
{
    std::cout << "================================="
        << std::endl;
    std::cout << '\t' << projectName << " Simulation" << std::endl;
    std::cout << "Input File: " << infname << std::endl;
    std::cout << "Total Number of Parameters: "<< parameter.size() << std::endl;
    std::cout << "    Check   Parameters    " << std::endl;
    int count = 0;
    for (std::map<std::string, double>::iterator
        it=parameter.begin(); it!=parameter.end(); ++it)
    {
        count++;
        std::cout << count << ", ";
        std::cout << it->first << " = " 
           << it->second << std::endl;
    }
    std::cout << "================================="
        << std::endl;
}
