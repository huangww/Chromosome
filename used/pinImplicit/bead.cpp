#include "bead.hpp"
#include "simulation.hpp"
#include "parameter.hpp"
#include "ultilities.hpp"
#include "config.hpp"
#include "force.hpp"
#include "rod.hpp"
#include "montecarlo.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

Bead::Bead(Simulation *simu) : Parameter(simu)
{
    r = create2DArray<double>(nBead, DIM);
    rs = create2DArray<double>(nBead, DIM);
    // v = create2DArray<double>(nBead, DIM);
    f = create2DArray<double>(nBead, DIM);
    ftotal = create2DArray<double>(nBead, DIM);

    force = new Force(simu);
    rod = new Rod(simu);
    montecarlo = new Montecarlo(simu);
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
    delete montecarlo;
}

void Bead::init() 
{
    t = 0.0;

    Config *config = new Config(simulation); 
    r = config->init();
    delete config;
    montecarlo->randomize();
    // montecarlo->equilibrate();
    
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
    addForce(rod->pseudo(f));
    addForce(force->brownian(f));
    addForce(force->external(f));
    // addForce(force->repulsive(f));

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

void Bead::montecarloUpdate()
{
        montecarlo->move();
        t += dt;
}

void Bead::output(std::ofstream& output) 
{
    output << "# t = " << t << std::endl;
    for (int i = 0; i < nBead; ++i) {
        for (int j = 0; j < DIM; ++j) {
            output << std::setw(9) << r[i][j] << '\t';
        }
        output << std::endl;
    } 
}
