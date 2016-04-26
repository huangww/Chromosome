#ifndef BEAD_HPP_JGQK9V2L
#define BEAD_HPP_JGQK9V2L

#include "simulation.hpp"
#include "parameter.hpp"
#include <fstream>

class Bead: public Parameter
{
public:
    Bead (Simulation *simu);
    virtual ~Bead ();

    void print();
    void init();
    void predict();
    void correct();
    void output(std::ofstream* output);
    void montecarloUpdate();

    double t;
// protected:
    double **r;     // postion of beads
    double **rs;    // predicted postion of beads
    double **f;     // forces applied on beads
    // double **v;     // velocity of beads

private:
    double **ftotal;    // total force

    class Force *force;
    class Rod *rod;
    class Montecarlo *montecarlo; 
    class Compute *compute;

    void addForce(double **);
    void pinSPB();
    void create();
    void destroy();
    void outputPos(std::ofstream& output);
    void outputRg(std::ofstream& output);
    void outputRd(std::ofstream& output);
};

#endif /* end of include guard: BEAD_HPP_JGQK9V2L */
