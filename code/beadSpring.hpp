#ifndef BEADSPRING_HPP_BTRDNUMP
#define BEADSPRING_HPP_BTRDNUMP

#include "project.hpp"
#include "input.hpp"

class BeadSpring: public Project
{
public:
    BeadSpring ();
    virtual ~BeadSpring ();

    void setup(Input *input);
    void run();

private:
    // parameters
    int outputStep;
    double tEnd;

    class Bead *bead;
    std::ofstream *outFile;

};

#endif /* end of include guard: BEADSPRING_HPP_BTRDNUMP */
