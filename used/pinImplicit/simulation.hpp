#ifndef SIMULATION_HPP_ZW9DJGS5
#define SIMULATION_HPP_ZW9DJGS5

class Simulation
{
public:
    Simulation ();
    virtual ~Simulation ();

    class Particle *particle;       // unit of 1D
    class Bead *bead;               // unit of 3D

    virtual void print() = 0;
    virtual void run() = 0;
};

#endif /* end of include guard: SIMULATION_HPP_ZW9DJGS5 */
