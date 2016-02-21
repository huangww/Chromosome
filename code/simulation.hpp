#ifndef SIMULATION_HPP_ZW9DJGS5
#define SIMULATION_HPP_ZW9DJGS5

#include <sstream>
class Simulation
{
public:
    Simulation ();
    virtual ~Simulation ();

<<<<<<< HEAD
    class Input *input;         // input scripting
    class Box *box;             // box to do simulation
    class Output *output;       // output streams

=======
    class Input *input;     // input scripting
    class Box *box;         // box to do the simulation
    class Output *output;   // output streams
    
    void init(std::string simuName);
>>>>>>> 06d6606705b9676fb9bda9404137c07e1b7d806c
    virtual void print() = 0;
    virtual void run() = 0;
};

#endif /* end of include guard: SIMULATION_HPP_ZW9DJGS5 */
