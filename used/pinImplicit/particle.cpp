#include "particle.hpp"
#include "parameter.hpp"
#include "random.hpp"
#include "force.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

Particle::Particle(Simulation *simu) : Parameter(simu)
{
    x = new double[nPar];
    v = new double[nPar];
    f = new double[nPar];
    ftotal = new double[nPar];
    site = new bool[nSite];

    force = new Force(simu);
}
Particle::~Particle() 
{
    delete [] x;
    delete [] v;
    delete [] f;
    delete [] ftotal;
    delete [] site;

    delete force;
}

void Particle::init() 
{
    t = 0;
    initRandom();
    std::fill(&v[0], &v[0] + nPar, 0);
    std::fill(&f[0], &f[0] + nPar, 0);
    std::fill(&ftotal[0], &ftotal[0] + nPar, 0);
}

void Particle::initRandom() 
{
    std::fill(&x[0], &x[0]+nPar, 0);
    std::fill(&site[0], &site[0]+nSite, 0);

    int count = 0;
    while (count < nPar) {
        double r = Ran(seed);
        int i = round(r*(nSite-1));
        if (!site[i]) {
            site[i] = true;
            count++;
        }
    }

    count = 0;
    for (int i = 0; i < nSite; ++i) {
        if (site[i]) {
            x[count] = i + 0.5;
            count++;
        }
    }
    
    // for (int i = 0; i < nPar; ++i) {
    //     v[i] = GaussRan(seed);
    // }
}

void Particle::print() 
{
    std::cout << "t = " << t << std::endl;
    std::cout << "Position: " << std::endl;
    for (int i = 0; i < nPar; ++i) {
        std::cout << i << '\t' << x[i] << std::endl;
    }
    // std::cout << "Velocity: " << std::endl;
    // for (int i = 0; i < nPar; ++i) {
    //     std::cout << i << '\t' << v[i]  << std::endl;
    // }
    std::cout << "Specific Force: " << std::endl;
    for (int i = 0; i < nPar; ++i) {
        std::cout << i << '\t' << f[i]  << std::endl;
    }
    std::cout << "Total Force: " << std::endl;
    for (int i = 0; i < nPar; ++i) {
        std::cout << i << '\t' << ftotal[i]  << std::endl;
    }
}

void Particle::addForce(double *f)
{
    for (int i = 0; i < nPar; ++i) {
        ftotal[i] += f[i];
    }
}

void Particle::update() 
{
    std::fill(&site[0], &site[0] + nSite, 0);
    for (int i = 0; i < nPar; ++i) {
        v[i] += 0.5*ftotal[i]*dt;
        x[i] += v[i]*dt;
        // site[int(x[i])] = true;
    }

    std::fill(&ftotal[0], &ftotal[0] + nPar, 0);
    addForce(force->brownian(f));
    addForce(force->repulsive(f));
    addForce(force->external(f));
    if (x[0] < 1.0 || x[nPar-1] > nSite-1.0) {
        addForce(force->boundary(f));
    }

    for (int i = 0; i < nPar; ++i) {
        v[i] += 0.5*ftotal[i]*dt;
    }
}

void Particle::updateBD() 
{
    std::fill(&ftotal[0], &ftotal[0] + nPar, 0);
    addForce(force->brownian(f));
    addForce(force->repulsive(f));
    addForce(force->external(f));
    if (x[0] < 1.0 || x[nPar-1] > nSite-1.0) {
        addForce(force->boundary(f));
    }

    // output for debug
    // double xmax, xmin;
    // xmin = *std::min_element(&x[0], &x[0]+nPar);
    // xmax = *std::max_element(&x[0], &x[0]+nPar); 
    // if (xmin < 0 or xmax > nSite) {
    //    print(); 
    //    std::cin.get();
    // }

    for (int i = 0; i < nPar; ++i) {
        x[i] += ftotal[i] * dt;
    }

    t += dt;
}

void Particle::output(std::ofstream& output) 
{
    output << t << '\t';
    // double beadPos = 0;
    // for (int i = 0; i < nSite; ++i) {
    //     output << beadPos << '\t';
    //     beadPos += 2*site[i] - 1;
    // }
    for (int i = 0; i < nPar; ++i) {
        output << std::setw(9) << x[i] << '\t';
    }
    output << std::endl;
}
