// #include "bead.hpp"
#include "ultilities.hpp"
#include "config.hpp"
#include "force.hpp"
#include "rod.hpp"
#include "spring.hpp"
#include "montecarlo.hpp"
#include "constant.hpp"
#include "compute.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>


Bead::Bead()
{
    r = NULL;
    rs = NULL;
    // v = NULL;
    f = NULL;
    ftotal = NULL;

    force = NULL;
    rod = NULL;
    spring = NULL;
    config = NULL;
    montecarlo = NULL;
    compute = NULL;
}
Bead::~Bead() 
{
    delete2DArray(r);
    delete2DArray(rs);
    // delete2DArray(v);
    delete2DArray(f);
    delete2DArray(ftotal);

    delete force;
    delete rod;
    delete spring;
    delete config;
    delete montecarlo;
    delete compute;
}

void Bead::setParameter(Input *input)
{
    if (input->parameter.count("nBead") == 0) {
        throw "Parameter \"nBead\" is not specified!";
    }
    nBead = int(input->parameter["nBead"]);
    if (input->parameter.count("dt") == 0) {
        throw "Parameter \"dt\" is not specified!";
    }
    dt = input->parameter["dt"];
    
    // set up the state class
    r = create2DArray<double>(nBead, DIM);
    rs = create2DArray<double>(nBead, DIM);
    // v = create2DArray<double>(nBead, DIM);
    f = create2DArray<double>(nBead, DIM);
    ftotal = create2DArray<double>(nBead, DIM);

    force = new Force();
    force->setParameter(input);
    if (input->projectName=="BeadRod") {
        rod = new Rod(this);
        rod->setParameter(input);
    }
    if (input->projectName=="BeadSpring") {
        spring = new Spring(this);
        spring->setParameter(input);
    }
    config = new Config();
    config->setParameter(input);
    montecarlo = new Montecarlo();
    montecarlo->setParameter(input);
    compute = new Compute();
}

void Bead::init() 
{
    t = 0.0;

    r = config->init(r);
    
    std::copy(&r[0][0], &r[0][0] + nBead * DIM, &rs[0][0]);
    // std::fill(&v[0][0], &v[0][0]+nBead*DIM, 0);
    std::fill(&f[0][0], &f[0][0]+nBead*DIM, 0);
    std::fill(&ftotal[0][0], &ftotal[0][0]+nBead*DIM, 0);
}

void Bead::init(const char* mode) 
{
    t = 0.0;

    r = config->init(r);
    if (strcmp(mode,"random")==0) {
        montecarlo->randomize(r);
    } else if (strcmp(mode,"equilibrate")==0) {
        montecarlo->equilibrate(r);
    }
    
    std::copy(&r[0][0], &r[0][0] + nBead * DIM, &rs[0][0]);
    // std::fill(&v[0][0], &v[0][0]+nBead*DIM, 0);
    std::fill(&f[0][0], &f[0][0]+nBead*DIM, 0);
    std::fill(&ftotal[0][0], &ftotal[0][0]+nBead*DIM, 0);
}


void Bead::print() 
{
    std::cout << "t = " << t << std::endl;
    std::cout << "Position: " << std::endl;
    for (int i = 0; i < nBead; ++i) {
        std::cout << i << '\t';
        for (int j = 0; j < DIM; ++j) {
            std::cout << std::setw(9) << r[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << "Specific Force: " << std::endl;
    for (int i = 0; i < nBead; ++i) {
        std::cout << i << '\t';
        for (int j = 0; j < DIM; ++j) {
            std::cout << std::setw(9) << f[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << "Total Force: " << std::endl;
    for (int i = 0; i < nBead; ++i) {
        std::cout << i << '\t';
        for (int j = 0; j < DIM; ++j) {
            std::cout << std::setw(9) << ftotal[i][j] << '\t';
        }
        std::cout << std::endl;
    }
}

void Bead::pinSPB() 
{
    for (int i = 0; i < DIM; ++i) {
        ftotal[0][i] += -1000.0*r[0][i];
    }
}

void Bead::addDrivenSPB() 
{
    ftotal[0][0] += 300.0*sin(t/10.0);
}

void Bead::addForce(double **f)
{
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            ftotal[i][j] += f[i][j];
        }
    }
}

void Bead::predict()
{
    std::fill(&ftotal[0][0], &ftotal[0][0] + nBead * DIM, 0);
    // addForce(rod->pseudo(f));
    // addForce(rod->pseudoSparse(f));
    addForce(rod->pseudoRing(f));
    addForce(force->brownian(f));
    addForce(force->external(f));
    // addForce(force->repulsive(r, f));
    
    // predict the next step position as rs
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            rs[i][j] = r[i][j] + ftotal[i][j] * dt;
        }
    }
    
}

void Bead::correct()
{
    f = rod->constraint(f);

    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            r[i][j] = rs[i][j] + f[i][j] * dt;
        }
    }

    t += dt;
}

void Bead::eulerUpdate()
{
    std::fill(&ftotal[0][0], &ftotal[0][0] + nBead * DIM, 0);

    pinSPB();
    addForce(spring->harmonic(r, f));
    // addForce(spring->fene(r, f));
    // addForce(spring->bending(r, f));
    // addForce(force->repulsive(r, f));
    addForce(force->external(f));
    addForce(force->brownian(f));
    
    // predict the next step position
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            r[i][j] = r[i][j] + ftotal[i][j] * dt;
        }
    }

    t += dt;
}

void Bead::rungerKuttaUpdate()
{
    std::fill(&ftotal[0][0], &ftotal[0][0] + nBead * DIM, 0);

    addForce(spring->fene(r, f));
    addForce(spring->bending(r, f));
    addForce(force->repulsive(r, f));
    addDrivenSPB();
    f = force->brownian(f);
    
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            rs[i][j] = r[i][j] + ftotal[i][j]*dt + f[i][j]*dt;
        }
    }

    addForce(spring->fene(rs, f));
    addForce(spring->bending(rs, f));
    addForce(force->repulsive(rs, f));
    addDrivenSPB();
    f = force->brownian(f);

    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            r[i][j] = r[i][j] + ftotal[i][j]*dt/2.0 + f[i][j]*dt;
        }
    }

    t += dt;
}

void Bead::montecarloUpdate()
{
    montecarlo->move(r);
    t += dt;
}

void Bead::output(std::ofstream* outFile) 
{
    if (std::isnan(r[0][0])) {
       throw "NaN error in r!"; 
    }
    outputPos(outFile[0]);
    // outputRg(outFile[1]);
    outputRd(outFile[2]);
}

void Bead::outputPos(std::ofstream& outFile) 
{
    outFile << "# t = " << t << '\n';
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            outFile << std::setw(9) << r[i][j] << '\t';
        }
        outFile << '\n';
    } 
}

void Bead::outputRg(std::ofstream& outFile)
{
    double rg = compute->gyrationRadius(nBead, r);
    outFile << std::setw(9) << t << '\t'
        << std::setw(9) << rg << '\n';
}

void Bead::outputRd(std::ofstream& outFile)
{
    outFile << std::setw(9) << t << '\t';
    for (int i = 0; i < DIM; ++i) {
        outFile << std::setw(9) << 
            r[nBead/2][i] - r[0][i] << '\t';
    }
    outFile << '\n';
}
