#include "state.hpp"
#include "input.hpp"
#include "random.hpp"
#include "compute.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

State::State()
{
    site = NULL;
    pos = NULL;
    beadPos = NULL;
    rate =  NULL;
    compute = NULL;
}

State::~State() 
{
    delete[] site;
    delete[] pos;
    delete[] beadPos;
    delete[] rate;
    delete compute;
}

void State::setParameter(Input *input)
{
// Todo: use try catch here
    nSite = int(input->parameter["nSite"]);
    nPar = int(input->parameter["nPar"]);

    tempEff = input->parameter["tempEff"];
    dt = input->parameter["dt"];
    double jumpRate = 1.0;
    double factor = exp(-1.0/tempEff);
    rateToLeft = jumpRate/(1+factor);
    rateToRight = jumpRate*factor/(1+factor);
    
    std::random_device rd;
    // seed = rd();
    seed = long(input->parameter["seed"]);

    // set up the state class
    site = new bool[nSite];
    pos = new int[nPar];
    beadPos = new double[nSite];
    rate = new double[nPar*2];
    compute = new Compute();
}

void State::init() 
{

    std::fill(&pos[0], &pos[0]+nPar, 0);
    std::fill(&site[0], &site[0]+nSite, 0);
    std::fill(&rate[0], &rate[0]+2*nPar, 0);

    initRandom();
    // initStretch();
    equilibrate();

    int count = 0;
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

    t = 0;
    tGrid = 0;
    par2bead();
    rg = compute->gyrationRadius(nSite, beadPos);
}

void State::initRandom() 
{
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
}

void State::initStretch() 
{
    for (int i = 0; i < nSite; ++i) {
        if (i < nPar) {
            site[i] = true;
        } else {
            site[i] = false;
        }
    }
}

void State::monteCarloMove()
{
    int i1 = round(Ran(seed)*(nSite-1));
    int i2 = round(Ran(seed)*(nSite-1));
    do {
        i1 = round(Ran(seed)*(nSite-1));
        i2 = round(Ran(seed)*(nSite-1));
    } while (site[i1] == site[i2]);

    double dE = 0;
    if (site[i1]) {
        dE = i2 - i1;
    } else {
        dE = i1 - i2;
    }
    if (dE < 0) {
        bool tmp = site[i1];
        site[i1] = site[i2];
        site[i2] = tmp;
    } else {
        if (Ran(seed) < exp(-dE/tempEff)) {
            bool tmp = site[i1];
            site[i1] = site[i2];
            site[i2] = tmp;
        } 
    }
}

void State::equilibrate()
{
    for (int i = 0; i < int(nSite*tempEff*100); ++i) {
        monteCarloMove();
    }
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

void State::par2bead()
{
    std::fill(&beadPos[0], &beadPos[0] + nSite, 0);
    double tmp = 0;
    for (int i = 0; i < nSite; ++i) {
        beadPos[i] = tmp;
        tmp += 2* site[i] - 1;
    }
}


void State::update() 
{
    tGrid += dt;
    while (t < tGrid) {
        double ran1 = Ran(seed);
        double dtJump = -1.0 / totalRate * log(ran1);
        if (t + dtJump < tGrid) {
            // pick which particle to jump
            double cumulateRate = 0.0;
            int index=0;      // index of particle to move
            int direction=0;  // 0 move left, 1 move right
            int i = 0;
            double ran2 = Ran(seed);
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
            if (direction) {
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
            } else {
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
            if (index != 0) {
                // check the right side of (i-1)th particle
                checkSite = site[pos[index-1]+1];
                if (checkSite) 
                    rate[2*index-1] = 0;
                else
                    rate[2*index-1] = rateToRight;
            }
            if (index != nPar-1) {
                // check the left side of (i+1)th particle
                checkSite = site[pos[index+1]-1];
                if (checkSite) 
                    rate[2*index+2] = 0;
                else
                    rate[2*index+2] = rateToLeft;
            }
            totalRate += (rate[2*index-1]+rate[2*index]
                    +rate[2*index+1]+rate[2*index+2]);
            t += dtJump;
        } else {
            t = tGrid;
        }
    }

    // update corresponding bead cofiguration
    par2bead();
    rg = compute->gyrationRadius(nSite, beadPos);
    // update done
}


void State::output(std::ofstream* output) 
{
    // outputPos(output[0]);
    outputRg(output[1]);
}

void State::outputSite(std::ofstream& output) 
{
    // output site state
    output << tGrid << '\t';
    for (int i = 0; i < nSite; ++i) {
        output << std::setw(6) << site[i] << ' ';
    } 
    output << std::endl;
}

void State::outputPar(std::ofstream& output) 
{
    // output particle position
    output << tGrid << '\t';
    for (int i = 0; i < nPar; ++i) {
        output << std::setw(6) << pos[i]; 
    }
    output << std::endl;
}

void State::outputPos(std::ofstream& output) 
{
     // output corresponding polymer bead position
    output << tGrid << '\t';
    for (int i = 0; i < nSite; ++i) {
        output << std::setw(6) << beadPos[i] << ' ';
    } 
    output << std::endl;
}

void State::outputRg(std::ofstream& output) 
{
    // output gyration radius
    output << std::setprecision(12) << tGrid << '\t';
    output << std::setprecision(12) << rg << '\t';
    output << std::setprecision(12) << beadPos[1] << '\t';
    output << std::setprecision(12) << beadPos[nSite/2] << '\t';
    output << std::endl;
}
