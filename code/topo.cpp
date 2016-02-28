#include "topo.hpp"
#include "input.hpp"
#include <algorithm>


Topo::Topo()  { }
Topo::~Topo() { }

void Topo::setParameter(Input *input) 
{
    if (input->parameter.count("topoType") == 0) {
        throw "Parameter \"topoType\" is not specified!";
    }
    topoType = int(input->parameter["topoType"]);
    if (input->parameter.count("nLink") == 0) {
        throw "Parameter \"nLink\" is not specified!";
    }
    nLink = int(input->parameter["nLink"]);
    nBead = int(input->parameter["nBead"]);
}

int** Topo::init(int** link)
{
    switch (topoType) {
        case 0:
            ring(nLink, link);
            break;
        case 1:
            chain(nLink, link);
            break;
        default:
            ring(nLink, link);
    }

    return link;
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
    int cm = 5;
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
    int m[3] = {10, 10, 10};
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
