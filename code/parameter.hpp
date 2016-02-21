#ifndef PARAMETER_HPP_H4YIFEO2
#define PARAMETER_HPP_H4YIFEO2

#include "simulation.hpp"
#include <iostream>
#include <sstream>


class Parameter
{
public:
    Parameter (Simulation *simu);
    virtual ~Parameter ();

    int taskID;
    
    void set(std::string key, std::string value);

protected:
    // some constant
    const double PI;
    const int DIM;

    Simulation *&simulation;

private:

};

#endif /* end of include guard: PARAMETER_HPP_H4YIFEO2 */
