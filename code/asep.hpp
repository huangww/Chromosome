#ifndef ASEP_HPP_MB6PNS13
#define ASEP_HPP_MB6PNS13

#include "project.hpp"
#include "input.hpp"

class Asep: public Project
{
public:
    Asep ();
    virtual ~Asep ();

    void setup(Input *input);
    void run();

private:
    class State *state;
    std::ofstream *outFile;
    
    // parameters
    int nSample;
    int outputStep;
    double tEnd;

};

#endif /* end of include guard: ASEP_HPP_MB6PNS13 */
