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
    void output(std::ofstream& output);

    double t;

private:
    bool *site;
    int *pos;
    double *rate;
    double totalRate;

    void initStretch();
    void initRandom();
    void print();
};

#endif /* end of include guard: STATE_HPP_3TUGIGTC */
