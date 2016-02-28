#ifndef PARTICLE_HPP_VUOZ3OXC
#define PARTICLE_HPP_VUOZ3OXC

#include "input.hpp"
#include <fstream>

class Particle
{
public:
    Particle ();
    virtual ~Particle ();
    
    double t;       // time
    double dt;      // time step
    double *x;      // position
    double *v;      // velocity
    double *f;      // force

    void setParameter(Input* input);
    void init();
    void print();
    void update();
    void updateBD();
    void output(std::ofstream &outFile);

private:
    double *ftotal; // total force
    bool *site;     // 0 empty, 1 occupied

    class Force *force;

    void initRandom();
    void addForce(double *);
    /* data */
    // parameters
    int nPar;
    int nSite;
    unsigned long seed;

};

#endif /* end of include guard: PARTICLE_HPP_VUOZ3OXC */
