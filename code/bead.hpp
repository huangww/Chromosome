#ifndef BEAD_HPP_JGQK9V2L
#define BEAD_HPP_JGQK9V2L

#include "input.hpp"
#include <fstream>

class Bead
{
public:
    Bead ();
    virtual ~Bead ();

    void setParameter(Input *input);
    void print();
    void init();
    void init(const char* mode);
    void predict();
    void correct();
    void eulerUpdate();
    void rungerKuttaUpdate();
    void montecarloUpdate();
    void output(std::ofstream* outFile);
    void outputPos(std::ofstream& output);
    void outputRd(std::ofstream& output);
    void outputRg(std::ofstream& output);

    double t;
    double dt;
// protected:
    double **r;     // postion of beads
    double **rs;    // predicted postion of beads
    double **f;     // forces applied on beads
    // double **v;     // velocity of beads

private:
    double **ftotal;    // total force

    class Force *force;
    class Rod *rod;
    class Spring *spring;
    class Config *config;
    class Montecarlo *montecarlo; 
    class Compute *compute;

    void addForce(double **);
    void pinSPB();
    void addDrivenSPB();
    void create();
    void destroy();
    
    // parameters
    int nBead;
};

#endif /* end of include guard: BEAD_HPP_JGQK9V2L */
