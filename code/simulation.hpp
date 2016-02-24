#ifndef SIMULATION_HPP_ZW9DJGS5
#define SIMULATION_HPP_ZW9DJGS5

#include <sstream>
class Simulation
{
public:
    Simulation ();
    virtual ~Simulation ();

    class Input *input;     // input scripting
    class Box *box;         // box to do the simulation
    class Output *output;   // output streams
    
    void init(std::string simuName);
    virtual void print() = 0;
    virtual void run() = 0;
};

#endif /* end of include guard: SIMULATION_HPP_ZW9DJGS5 */
