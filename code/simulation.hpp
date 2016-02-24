#ifndef SIMULATION_HPP_ZW9DJGS5
#define SIMULATION_HPP_ZW9DJGS5

class Simulation
{
public:
    Simulation ();
    virtual ~Simulation ();

    class Input *input;     // input scripting

    void init();
    void run();
    void print();


private:
    class SimuASEP *simuASEP;
    class SimuSingleFile *simuSingleFile;
    class SimuBeadRod *simuBeadRod;
    class SimuBeadSpring *simuBeadSpring;

};

#endif /* end of include guard: SIMULATION_HPP_ZW9DJGS5 */
