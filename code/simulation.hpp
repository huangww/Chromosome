#ifndef SIMULATION_HPP_ZW9DJGS5
#define SIMULATION_HPP_ZW9DJGS5

class Simulation
{
public:
    Simulation ();
    virtual ~Simulation ();

    class Input *input;             // input scripting
    class Output *output;           // output streams

    virtual void print() = 0;
    virtual void run() = 0;
};

#endif /* end of include guard: SIMULATION_HPP_ZW9DJGS5 */
