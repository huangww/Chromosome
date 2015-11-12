#include "main.h"
#include "utilities.h"
#include <stdio.h>
#include <string.h>


int TopolSingleRing(int N, int link[N][2])
{
    memset(link, 0, sizeof(link[0][0]) * N * 2);
    for (int i = 0; i < N; i++) 
    {
        link[i][0] = i;
        link[i][1] = (i+1) % N;
    }
    return 0;
}

int TopolSingleChain(int N, int link[N][2])
{
    memset(link, 0, sizeof(link[0][0]) * N * 2);
    for (int i = 0; i < N; i++) 
    {
        link[i][0] = i;
        link[i][1] = i+1;
    }
    return 0;
}

int TopolRingPair(int N, int link[N][2])
{
    memset(link, 0, sizeof(link[0][0]) * N * 2);
    for (int i = 0; i < N/2-1; ++i)
    {
        link[i][0] = i;
        link[i][1] = i + 1;
    }
    link[N/2-1][0] = N/2 - 1;
    link[N/2-1][1] = 0;
    link[N/2][0] = 0;
    link[N/2][1] = N/2;
    for (int i = N/2+1; i < numUnit; ++i)
    {
        link[i][0] = i - 1;
        link[i][1] = i;
    }
    link[N-1][0] = numUnit - 1;
    link[N-1][1] = 0;
    return 0;
}

int TopolCentromerePair(int N, int link[N][2])
{
    memset(link, 0, sizeof(link[0][0]) * N * 2);
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
    for (int i = ringSize+cm-1; i < numUnit; ++i)
    {
        link[index][0] = i;
        link[index][1] = (i+1)%numUnit;
        index = index + 1;
    }

    return 0;
}

int TopolThreePair(int N, int link[N][2])
{
    memset(link, 0, sizeof(link[0][0]) * N * 2);
    int *m = monomer;
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

    return 0;
}

void BondList(int topologyType,
        int link[][2])
{
    memset(link, 0, sizeof(link[0][0]) * numLink * 2);

    /* Define topological constriants (rods)*/
    switch (topologyType) {
        case 0:
            TopolSingleRing(numLink, link);
            break;
        case 1:
            TopolSingleChain(numLink, link);
            break;
        case 2:
            TopolRingPair(numLink, link);
            break;
        case 3:
            TopolCentromerePair(numLink, link);
            break;
        case 4:
            TopolThreePair(numLink, link);
            break;
        default:
            TopolSingleRing(numLink, link);

    }

}

int AngleList(int link[][2], int angleList[][3])
{
    int count = 0;
    for (int i = 0; i < numLink; ++i)
    {
        for (int j = i+1; j < numLink; ++j)
        {
            if (link[i][1]==link[j][0] && link[i][1]!=0)
            {
                angleList[count][0] = link[i][0];
                angleList[count][1] = link[i][1];
                angleList[count][2] = link[j][1];
                count++;
            }
        }
    }
    return count;
}
