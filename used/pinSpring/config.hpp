#ifndef CONFIG_HPP_4IBHETEB
#define CONFIG_HPP_4IBHETEB

#include "simulation.hpp"
#include "parameter.hpp"

class Config: protected Parameter
{
public:
    Config (Simulation *simu);
    virtual ~Config ();

    double** init();
    void fixedChain(int N, double **pos);

private:
    void ring(int N, double **pos);
    void chain(int N, double **pos);
    // void fixedChain(int N, double **pos);
    void straightRing(int N, double **pos); 
    void quenchedRing(int N, double **pos); 
    void ringPair(int N, double **pos); 
    void ringPairWithCentromere(int N, double **pos);
    void threeRingPair(int N, double **pos); 

};

#endif /* end of include guard: CONFIG_HPP_4IBHETEB */
