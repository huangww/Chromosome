#ifndef MONTECARLO_HPP_9WRC7ZF3
#define MONTECARLO_HPP_9WRC7ZF3


#include "input.hpp"

class Montecarlo
{
public:
    Montecarlo ();
    virtual ~Montecarlo ();

    void setParameter(Input *input);
    int move(double **r);
    void randomize(double **r);
    void equilibrate(double **r);

private:
    // parameters
    int nBead;
    int topoType;
    unsigned long seed;
    double tempEff;

    void pivot(int N, int *ipivot);
    void rotateAxis(int* ipivot, double** r, double* axis);
    void rotateMatrix( double theta, double* axis, 
            double** matrix);
    void moveRing(int N, double** r); 
    void moveChain(int N, double** r);
    void moveRingPair(int N, double** r);
    void moveThreeRingPair(int N, double** r);
    void moveRingPairWithCentromere(int N, double** r);
    void moveTry(int N, double** r);
    double energy(int N, double** r);

};

#endif /* end of include guard: MONTECARLO_HPP_9WRC7ZF3 */
