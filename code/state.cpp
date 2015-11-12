#include "state.hpp"
#include "simulation.hpp"
#include "random.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>


State::State(Simulation *simu) : Parameter(simu)
{
    site = new bool[nSite];
    pos = new int[nPar];
    rate = new double[nPar*2];
}

State::~State() 
{
    delete[] site;
    delete[] pos;
    delete[] rate;
}

void State::init() 
{
    initRandom();
    // initStretch();
}

void State::initRandom() 
{
    std::fill(&pos[0], &pos[0]+nPar, 0);
    std::fill(&site[0], &site[0]+nSite, 0);
    std::fill(&rate[0], &rate[0]+2*nPar, 0);

    t = 0;

    int count = 0;
    while(count < nPar) {
        double r = Ran(seed);
        int i = round(r*(nSite-1));
        if (!site[i])
        {
            site[i] = true;
            count++;
        }
    }

    count = 0;
    totalRate = 0;
    if (site[0])
    {
        pos[count] = 0;
        rate[0] = 0;
        if (site[1])
            rate[1] = 0;
        else
            rate[1] = rateToRight;
        totalRate += (rate[0]+rate[1]);
        count++;
    }
    for (int i = 1; i < nSite-1; ++i)
    {
        if (site[i])
        {
            pos[count] = i;
            if (!site[i-1]) 
                rate[2*count] = rateToLeft;
            if (!site[i+1])
                rate[2*count+1] = rateToRight;
            totalRate += (rate[2*count]+rate[2*count+1]);
            count++; 
        }
    }
    if (site[nSite-1])
    {
        pos[count] = nSite - 1;
        rate[2*count+1] = 0;
        if (site[nSite-2])
            rate[2*count] = 0;
        else
            rate[2*count] = rateToLeft;
        totalRate += (rate[2*count]+rate[2*count+1]);
        count++;
    }

}

void State::initStretch() 
{
    t = 0;

    for (int i = 0; i < nSite; ++i)
    {
        if (i < nPar)
        {
            site[i] = true;
            pos[i] = i;
        }
        else
        {
            site[i] = false;
        }
    }

    for (int i = 0; i < 2*nPar; ++i)
    {
        rate[i] = 0;
        if (i == 2*nPar-1) rate[i] = rateToRight;
    }

    totalRate = rateToRight;
}

void State::print() 
{
    std::cout << "================================="
        << std::endl;
    std::cout << "Time:" << t
        << std::endl;

    std::cout << "Lattice Site Occupation State:"
        << std::endl;
    for (int i = 0; i < nSite; ++i)
    {
        std::cout << i << "\t"
            << site[i] << std::endl;
    }

    std::cout << "Particle Position:"
        << std::endl;
    for (int i = 0; i < nPar; ++i)
    {
        std::cout << i << "\t"
            << pos[i] << std::endl;
    }

    std::cout << "Jumping Rate:"
        << std::endl;
    for (int i = 0; i < nPar; ++i)
    {
        std::cout << i << "\t"
            << rate[2*i] << "\t" 
            << rate[2*i+1] << std::endl;
    }

    std::cout << "Total Jumping Rate:" << totalRate
        << std::endl;

}

double State::update() 
{
    double ran1 = Ran(seed);
    double ran2 = Ran(seed);

    double dt = -1.0 / totalRate * log(ran1);

    // pick which particle to jump
    double cumulateRate = 0.0;
    int index;      // index of particle to move
    int direction;  // 0 move left, 1 move right
    int i = 0;
    while(cumulateRate <= ran2*totalRate ) {
        cumulateRate += rate[i]; 
        index = i / 2;
        direction = i % 2;
        i++;
    }

    // update state, carefully
    bool checkSite;
    totalRate -= (rate[2*index-1]+rate[2*index]
            +rate[2*index+1]+rate[2*index+2]);
    if (direction)
    {
        // move particle right
        site[pos[index]] = false;
        site[pos[index]+1] = true;
        pos[index]++;

        // update the rate of the ith particle
        // no need to check the left side
        rate[2*index] = rateToLeft;
        // check the right side
        checkSite = site[pos[index]+1];
        if (checkSite || pos[index]==nSite-1) 
            rate[2*index+1] = 0;
    }
    else
    {
        // move particle to left
        site[pos[index]] = false;
        site[pos[index]-1] = true;
        pos[index]--;

        // update the rate of the ith particle
        // check the left side
        checkSite = site[pos[index]-1];
        if (checkSite || pos[index]==0) 
            rate[2*index] = 0;
        // no need to check the right side
        rate[2*index+1] = rateToRight;
    }

    // update rate of (i-1)th and (i+1)th particle
    if (index != 0)
    {
        // check the right side of (i-1)th particle
        checkSite = site[pos[index-1]+1];
        if (checkSite) 
            rate[2*index-1] = 0;
        else
            rate[2*index-1] = rateToRight;
    }
    if (index != nPar-1)
    {
        // check the left side of (i+1)th particle
        checkSite = site[pos[index+1]-1];
        if (checkSite) 
            rate[2*index+2] = 0;
        else
            rate[2*index+2] = rateToLeft;
    }
    totalRate += (rate[2*index-1]+rate[2*index]
            +rate[2*index+1]+rate[2*index+2]);

    // update done

    t += dt;
    return dt;
}

void State::output(std::ofstream& output) 
{
    // output site state
    // output << t << '\t';
    // for (int i = 0; i < nSite; ++i) {
    //     output << std::setw(6) << site[i] << ' ';
    // } 
    // output << std::endl;

    // output particle position
    // output << t << '\t';
    // for (int i = 0; i < nPar; ++i) {
    //     output << std::setw(6) << pos[i]; 
    // }
    // output << std::endl;

    // output corresponding polymer bead position
    output << t << '\t';
    double beadPos = 0;
    for (int i = 0; i < nSite; ++i)
    {
        output << std::setw(6) << beadPos << ' ';
        beadPos += 2*site[i] - 1;
    } 
    output << std::endl;
    
}
