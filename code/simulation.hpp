#ifndef SIMULATION_HPP_ZW9DJGS5
#define SIMULATION_HPP_ZW9DJGS5

class Simulation
{
public:
    Simulation ();
    virtual ~Simulation ();

    class Input *input;     // input scripting
    class Project *project; // projects of all types

    void init();
    void print();
};

#endif /* end of include guard: SIMULATION_HPP_ZW9DJGS5 */
