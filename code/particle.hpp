#ifndef PARTICLE_HPP_VUOZ3OXC
#define PARTICLE_HPP_VUOZ3OXC

#include "parameter.hpp"
#include "simulation.hpp"
#include <fstream>

class Particle: public Parameter
{
public:
    Particle (Simulation *simu);
    virtual ~Particle ();
    
    double t;       // time
    double *x;      // position
    double *v;      // velocity
    double *f;      // force

    void init();
    void print();
    void update();
    void updateBD();
    void output(std::ofstream &output);

private:
    double *ftotal; // total force
    bool *site;     // 0 empty, 1 occupied

    class Force *force;

    void initRandom();
    void addForce(double *);
    /* data */
};

#endif /* end of include guard: PARTICLE_HPP_VUOZ3OXC */
