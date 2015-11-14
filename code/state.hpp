#ifndef STATE_HPP_3TUGIGTC
#define STATE_HPP_3TUGIGTC

#include <fstream>
#include "simulation.hpp"
#include "parameter.hpp"

class State: public Parameter
{
public:
    State (Simulation *simu);
    virtual ~State ();

    void init();
    double update();
    void output(std::ofstream* output);

    double t;

private:
    bool *site;
    int *pos;
    double *beadPos;
    double *rate;
    double totalRate;
    class Compute *compute;

    void initStretch();
    void initRandom();
    void print();
    void par2bead();
    void outputSite(std::ofstream& output);
    void outputPar(std::ofstream& output);
    void outputPos(std::ofstream& output);
    void outputRg(std::ofstream& output);
};

#endif /* end of include guard: STATE_HPP_3TUGIGTC */
