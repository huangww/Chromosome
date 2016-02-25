#ifndef CONFIG_HPP_4IBHETEB
#define CONFIG_HPP_4IBHETEB

#include "input.hpp"

class Config
{
public:
    Config ();
    virtual ~Config ();

    double** init(double **r);
    void setParameter(Input *input);
    void fixedChain(int N, double **pos);

private:
    // parameters
    int topoType;
    int nBead;

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
