#ifndef CONFIG_HPP_4IBHETEB
#define CONFIG_HPP_4IBHETEB

#include "input.hpp"

class Config
{
public:
    Config ();
    virtual ~Config ();

    void setParameter(Input *input);

    Vec init(Vec r);
    double** init(double **r);
    void fixedChain(int N, double **pos);

private:
    // parameters
    int topoType;
    int nBead;

    void ring(int N, double **pos);
    void chain(int N, double **pos);
    void straightRing(int N, double **pos); 
    void quenchedRing(int N, double **pos); 
    void ringPair(int N, double **pos); 
    void ringPairWithCentromere(int N, double **pos);
    void threeRingPair(int N, double **pos); 
    void ring(int N, Vec pos);
    void chain(int N, Vec pos);
    void straightRing(int N, Vec pos); 
    void quenchedRing(int N, Vec pos); 
    void ringPair(int N, Vec pos); 
    void ringPairWithCentromere(int N, Vec pos);
    void threeRingPair(int N, Vec pos); 


};

#endif /* end of include guard: CONFIG_HPP_4IBHETEB */
