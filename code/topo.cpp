#include "topo.hpp"
#include "input.hpp"
#include "ultilities.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>


Topo::Topo() 
{
    initFlag = 0;
    link = nullptr;
    nPair = nullptr;
    linkPair = nullptr;
    g = nullptr;
    gSparse = nullptr;
}
Topo::~Topo() 
{
    delete2DArray(link);
    delete2DArray(linkPair);
    delete2DArray(g);
    delete gSparse;
    delete [] nPair;
}

void Topo::setParameter(Input *input) 
{
    if (input->parameter.count("topoType") == 0) {
        throw "Parameter \"topoType\" is not specified!";
    }

    topoType = int(input->parameter["topoType"]);
    nBead = int(input->parameter["nBead"]);
    // nLink = int(input->parameter["nLink"]);
    nLink = setNumLink();

    link = create2DArray<int>(nLink, 2);
}

void Topo::init()
{
    if (!initFlag) {
        g = create2DArray<int>(nLink, nLink);
        setLink();
        setMetricTensor();
        int totalLinkPair = getTotalLinkPair();
        linkPair = create2DArray<int>(totalLinkPair, 6);
        nPair = new int[nBead];
        setLinkPair();
        initFlag = 1;
    }
}

void Topo::initSparse()
{
    if (!initFlag) {
        g = create2DArray<int>(nLink, nLink);
        setLink();
        setMetricTensor();
        int totalLinkPair = getTotalLinkPair();
        linkPair = create2DArray<int>(totalLinkPair, 6);
        nPair = new int[nBead];
        gSparse = new SpMatD(nLink, nLink);
        setMetricTensorSparse();
        setLinkPair();
        initFlag = 1;
    }
}

void Topo::setMetricTensor() 
{
    std::fill(&g[0][0], &g[0][0] + nLink * nLink, 0);

    for (int i = 0; i < nLink; i++) {
        g[i][i] = -2;
        for (int j = 0; j < nLink; j++) {
            if (i==j) continue;
            int condition1 = (link[i][0]-link[j][1]) *
                (link[i][1]-link[j][0]);
            int condition2 = (link[i][0]-link[j][0]) *
                (link[i][1]-link[j][1]);
            if (condition1 == 0) {
                g[i][j] = 1;
            }
            else if (condition2 == 0) {
                g[i][j] = -1;
            }
        }
    }

    // check \sum_j g[i][j] = 0
    // for (int i = 0; i < nLink; ++i) {
    //     int count = 0;
    //     for (int j = 0; j < nLink; ++j) {
    //         count += g[i][j];
    //     }
    //     if (count!=0) {
    //         std::cout << i << '\t' << count << std::endl;
    //     }
    // }
    // std::cout << "check is done!" << std::endl;
    // std::cin.get();

}

void Topo::setMetricTensorSparse() 
{
    std::vector<T> tripletList;
    tripletList.reserve(3*nLink);
    for (int i = 0; i < nLink; i++) {
        tripletList.push_back(T(i,i,-2));
        for (int j = 0; j < nLink; j++) {
            if (i==j) continue;
            int condition1 = (link[i][0]-link[j][1]) *
                (link[i][1]-link[j][0]);
            int condition2 = (link[i][0]-link[j][0]) *
                (link[i][1]-link[j][1]);
            if (condition1 == 0) {
                tripletList.push_back(T(i,j,1));
            }
            else if (condition2 == 0) {
                tripletList.push_back(T(i,j,-1));
            }
        }
    }
    gSparse->setFromTriplets(tripletList.begin(), tripletList.end());
    gSparse->makeCompressed();
}

int Topo::getTotalLinkPair()
{
    int count = 0;
    for (int l = 0; l < nBead; ++l) {
        for (int i = 0; i < nLink; ++i) {
            for (int j = i+1; j < nLink; ++j) {
                if (g[i][j]!=0) {
                    int i0,i1,j0,j1;
                    i0 = link[i][0];
                    i1 = link[i][1];
                    j0 = link[j][0];
                    j1 = link[j][1];
                    if ((l-i0)*(l-i1)*(l-j0)*(l-j1)==0) {
                        count++;
                    }
                }
            }
        }
    }
    return count;
}

void Topo::setLinkPair()
{
    int count = 0;
    for (int l = 0; l < nBead; ++l) {
        int countLink = 0;
        for (int i = 0; i < nLink; ++i) {
            for (int j = i+1; j < nLink; ++j) {
                if (g[i][j]!=0) {
                    int i0,i1,j0,j1;
                    i0 = link[i][0];
                    i1 = link[i][1];
                    j0 = link[j][0];
                    j1 = link[j][1];
                    if ((l-i0)*(l-i1)*(l-j0)*(l-j1)==0) {
                        linkPair[count][0] = i;
                        linkPair[count][1] = i0;
                        linkPair[count][2] = i1;
                        linkPair[count][3] = j;
                        linkPair[count][4] = j0;
                        linkPair[count][5] = j1;
                        count++;
                        countLink++;
                    }
                }
            }
        }
        nPair[l] = countLink;
    }
}

int Topo::getNumLink() 
{ 
    if (!nLink) setNumLink();
    return nLink; 
}

int** Topo::getLink()
{ 
    if (!initFlag) setLink();
    outputLinks();
    return link;
}

int** Topo::getLinkPair() 
{
    if (!initFlag) setLinkPair();
    return linkPair; 
}

int* Topo::getNumPair() 
{ 
    if (!initFlag) setLinkPair();
    return nPair; 
}

int** Topo::getMetricTensor() 
{ 
    if (!initFlag) setMetricTensor();
    return g; 
}

SpMatD* Topo::getMetricTensorSparse() 
{ 
    if (!initFlag) setMetricTensorSparse();
    return gSparse; 
}

void Topo::outputLinks() 
{
    std::ofstream output("data/topo.dat");

    for (int i = 0; i < nLink; ++i) {
            output << std::setw(9) << link[i][0] << '\t' 
                << std::setw(9) << link[i][1] << std::endl;
    }

    output.close();
}

void Topo::printLinks() 
{
    for (int i = 0; i < nLink; ++i) {
        std::cout << std::setw(9) << link[i][0] << '\t' 
                << std::setw(9) << link[i][1] << std::endl;
    }
}

void Topo::ring(int N, int **link)
{
    std::fill(&link[0][0], &link[0][0] + N * 2, 0);
    for (int i = 0; i < N; i++) 
    {
        link[i][0] = i;
        link[i][1] = (i+1) % N;
    }
}

void Topo::chain(int N, int **link)
{
    std::fill(&link[0][0], &link[0][0] + N * 2, 0);
    for (int i = 0; i < N; i++) 
    {
        link[i][0] = i;
        link[i][1] = i+1;
    }
}

void Topo::ringPair(int N, int **link)
{
    std::fill(&link[0][0], &link[0][0] + N * 2, 0);

    if (N%2!=0) {
        throw "Odd number of beads required for this topo type!";
    }

    for (int i = 0; i < N/2-1; ++i)
    {
        link[i][0] = i;
        link[i][1] = i + 1;
    }
    link[N/2-1][0] = N/2 - 1;
    link[N/2-1][1] = 0;
    link[N/2][0] = 0;
    link[N/2][1] = N/2;
    for (int i = N/2+1; i < nBead; ++i)
    {
        link[i][0] = i - 1;
        link[i][1] = i;
    }
    link[N-1][0] = nBead - 1;
    link[N-1][1] = 0;
}

void Topo::ringPairWithCentromere(int N, int **link)
{
    std::fill(&link[0][0], &link[0][0] + N * 2, 0);

    if (N%2!=0) {
        throw "Even number of beads required for this topo type!";
    }
    int index = 0;
    int ringSize = N/2;
    for (int i = 0; i < ringSize; ++i)
    {
        link[index][0] = i;
        link[index][1] = (i+1)%ringSize;
        index = index + 1;
    }
    link[index][0] = 0;
    link[index][1] = ringSize;
    index = index + 1;
    int cm = ringSize/2;
    for (int i = ringSize; i < ringSize+cm-2; ++i)
    {
        link[index][0] = i;
        link[index][1] = i+1;
        index = index + 1;
    }
    link[index][0] = ringSize + cm - 2;
    link[index][1] = cm;
    index = index + 1;
    link[index][0] = cm;
    link[index][1] = ringSize + cm - 1;
    index = index + 1;
    for (int i = ringSize+cm-1; i < nBead; ++i)
    {
        link[index][0] = i;
        link[index][1] = (i+1)%nBead;
        index = index + 1;
    }

}

void Topo::threeRingPair(int N, int **link)
{
    std::fill(&link[0][0], &link[0][0] + N * 2, 0);
    if (N%2!=0) {
        throw "Odd number of beads required for this topo type!";
    }
    // int m[3] = {558, 454, 245};
    int m[3];
    int monomer = N/2;
    m[1] = monomer * 245/1257;
    m[2] = monomer * 454/1257;
    m[0] = monomer - m[1] -m[2];
    int k = 0;

    /* the first pair */
    link[k][0] = 0;
    link[k][1] = 1;
    k = k + 1;
    for (int i = 0; i < m[0] - 1; ++i)
    {
        link[k][0] = k;
        link[k][1] = (k+1) % m[0];
        k = k + 1;
    }
    link[k][0] = 0;
    link[k][1] = m[0];
    k = k + 1;
    for (int i = 0; i < m[0] - 1; ++i)
    {
        link[k][0] = k - 1;
        link[k][1] = k % (2*m[0]-1);
        k = k + 1;
    }

    /* the second pair */
    link[k][0] = 0;
    link[k][1] = 2*m[0] - 1;
    k = k + 1;
    for (int i = 0; i < m[1] - 1; ++i)
    {
        link[k][0] = k - 2;
        link[k][1] = (k-1) % (2*m[0]+m[1]-2);
        k = k + 1;
    }
    link[k][0] = 0;
    link[k][1] = 2*m[0] + m[1] - 2;
    k = k + 1;
    for (int i = 0; i < m[1] - 1; ++i)
    {
        link[k][0] = k - 3;
        link[k][1] = (k-2) % (2*m[0]+2*m[1]-3);
        k = k + 1;
    }

    /* the third pair */
    link[k][0] = 0;
    link[k][1] = 2*m[0] + 2*m[1] - 3;
    k = k + 1;
    for (int i = 0; i < m[2] - 1; ++i)
    {
        link[k][0] = k - 4;
        link[k][1] = (k-3) % (2*m[0]+2*m[1]+m[2]-4);
        k = k + 1;
    }
    link[k][0] = 0;
    link[k][1] = 2*m[0] + 2*m[1] + m[2] - 4;
    k = k + 1;
    for (int i = 0; i < m[2] - 1; ++i)
    {
        link[k][0] = k - 5;
        link[k][1] = (k-4) % (2*(m[0]+m[1]+m[2])-5);
        k = k + 1;
    }
}

int Topo::setNumLink()
{
    switch (topoType) {
        case 0:
            nLink = nBead;
            break;
        case 1:
            nLink = nBead - 1;
            break;
        case 2:
            nLink = nBead + 1;
            break;
        case 3:
            nLink = nBead + 2;
            break;
        case 4:
            nLink = nBead + 5;
            break;
        default:
            nLink = nBead;
    }
    return nLink;
}

void Topo::listTopo()
{
    std::cout << " 0 : Ring " << std::endl;
    std::cout << " 1 : Chain " << std::endl;
    std::cout << " 2 : A pair of rings " << std::endl;
    std::cout << " 3 : A pair of rings with centromere" << std::endl;
    std::cout << " 4 : Three pairs of rings" << std::endl;
}
void Topo::setLink()
{
    switch (topoType) {
        case 0:
            ring(nLink, link);
            break;
        case 1:
            chain(nLink, link);
            break;
        case 2:
            ringPair(nLink, link);
            break;
        case 3:
            ringPairWithCentromere(nLink, link);
            break;
        case 4:
            threeRingPair(nLink, link);
            break;
        default:
            listTopo();
            std::cout << topoType << 
                " is not valid topoType!" << std::endl;
            throw "Wrong topoTpye specified!";
    }

}
