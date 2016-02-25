#ifndef BEADROD_HPP_2YIZDR6W
#define BEADROD_HPP_2YIZDR6W

#include "project.hpp"
#include "input.hpp"

class BeadRod: public Project
{
public:
    BeadRod ();
    virtual ~BeadRod ();

    void setup(Input *input);
    void run();

private:
    class State *bead;
    std::ofstream *outFile;
    
    // parameters
    int outputStep;
    double tEnd;
};

#endif /* end of include guard: BEADROD_HPP_2YIZDR6W */
