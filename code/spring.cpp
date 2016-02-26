#include "spring.hpp"
#include "input.hpp"
#include "bead.hpp"
#include "ultilities.hpp"
#include "constant.hpp"
#include "topo.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

Spring::Spring()
{ 
    link = NULL;
}
Spring::~Spring() 
{ 
    delete2DArray(link);
}

void Spring::setParameter(Input *input)
{

    nBead =int(input->parameter["nRod"]);
    nRod = int(input->parameter["nRod"]);

    link = create2DArray<int>(nRod, 2);
    Topo *topo = new Topo();
    topo->setParameter(input);
    link = topo->init(link);
    delete topo;

    outputLinks();
}


void Spring::outputLinks() 
{
    std::ofstream output("data/topo.dat");

    for (int i = 0; i < nRod; ++i) {
            output << std::setw(9) << link[i][0] << '\t' 
                << std::setw(9) << link[i][1] << std::endl;
    }

    output.close();
}


double** Spring::harmonic(double** r, double** f) 
{
    std::fill(&f[0][0], &f[0][0] + nBead*DIM, 0);

    double k = 3.0;
    for (int i = 0; i < nRod; ++i) {
        int i0 = link[i][0];
        int i1 = link[i][1];
        for (int j = 0; j < DIM; ++j) {
            double df;
            df = -k*(r[i0][j]-r[i1][j]); 
            f[i0][j] += df;
            f[i1][j] -= df;
        }
    }
    
    return f;
}
    
double** Spring::fene(double** r, double** f) 
{
    std::fill(&f[0][0], &f[0][0] + nBead*DIM, 0);

    double k = 100.0;
    double R0 = 1.1;
    for (int i = 0; i < nRod; ++i) {
        int i0 = link[i][0];
        int i1 = link[i][1];
        double d = Distance(&r[i0][0], &r[i1][0], DIM);
        for (int j = 0; j < DIM; ++j) {
            double rlogarg = 1 - (d*d)/(R0*R0);
            if (rlogarg < 0.01) rlogarg = 0.01;
            double df;
            df = -k*(r[i0][j]-r[i1][j])/rlogarg; 
            f[i0][j] += df;
            f[i1][j] -= df;
        }
    }
    return f;
}
