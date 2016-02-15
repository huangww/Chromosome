#ifndef INPUT_HPP_DRBPSRM7
#define INPUT_HPP_DRBPSRM7

#include "parameter.hpp"
#include "simulation.hpp"
#include <iostream>
#include <fstream>

class Input: protected Parameter
{
public:
    Input (Simulation *simu);
    virtual ~Input ();

    void getArg(int argc, char *argv[]);
    void file();

private:
    std::fstream infile;
    void parse();
    void excute();
};

#endif /* end of include guard: INPUT_HPP_DRBPSRM7 */
