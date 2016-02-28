#ifndef STATE_HPP_3TUGIGTC
#define STATE_HPP_3TUGIGTC

#include <fstream>
#include "input.hpp"

class State 
{
public:
    State ();
    virtual ~State ();

    void setParameter(Input *input);
    void init();
    void update();
    void output(std::ofstream* output);

    double t;       // real time
    double tGrid;   // grid time for output
    double rg;      // radius of gyration
    double dt;      // time step of the update

private:
    // parameters
    int nSite;
    int nPar;

    double tempEff;
    double rateToLeft;
    double rateToRight;
    unsigned long seed;

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
    void monteCarloMove();
    void equilibrate();
    void outputSite(std::ofstream& output);
    void outputPar(std::ofstream& output);
    void outputPos(std::ofstream& output);
    void outputRg(std::ofstream& output);
};

#endif /* end of include guard: STATE_HPP_3TUGIGTC */
