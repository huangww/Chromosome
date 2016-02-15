#ifndef MONTECARLO_HPP_9WRC7ZF3
#define MONTECARLO_HPP_9WRC7ZF3

#include "parameter.hpp"
#include "simulation.hpp"

class Montecarlo: protected Parameter
{
public:
    Montecarlo (Simulation *simu);
    virtual ~Montecarlo ();

    int move();
    void randomize();
    void equilibrate();

private:
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

    // class Potential *potential;
};

#endif /* end of include guard: MONTECARLO_HPP_9WRC7ZF3 */
