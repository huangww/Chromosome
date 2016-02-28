#ifndef SINGLEFILE_HPP_ILK3AJTE
#define SINGLEFILE_HPP_ILK3AJTE

#include "project.hpp"
#include "input.hpp"
#include <fstream>

class SingleFile: public Project
{
public:
    SingleFile ();
    virtual ~SingleFile ();

    void setup(Input* input);
    void run();

private:
    class Particle *particle;
    std::ofstream outFile;

    // parameters
    int outputStep;
    double tEnd;

};

#endif /* end of include guard: SINGLEFILE_HPP_ILK3AJTE */
