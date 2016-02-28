#ifndef BEADSPRING_HPP_BTRDNUMP
#define BEADSPRING_HPP_BTRDNUMP

#include "project.hpp"
#include "input.hpp"
#include <fstream>

class BeadSpring: public Project
{
public:
    BeadSpring ();
    virtual ~BeadSpring ();

    void setup(Input *input);
    void run();

private:
    class Bead *bead;
    std::ofstream *outFile;

    // parameters
    int outputStep;
    double tEnd;
};

#endif /* end of include guard: BEADSPRING_HPP_BTRDNUMP */
