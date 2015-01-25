#include "main.h"
#include <stdio.h>
#include <string.h>


int TopolSingleRing(int N, int link[N][2])
{
	int topolType = 1;
	memset(link, 0, sizeof(link[0][0]) * N * 2);
	for (int i = 0; i < N; i++) 
	{
		link[i][0] = i;
		link[i][1] = (i+1) % N;
	}
	return topolType;
}

int TopolSingleChain(int N, int link[N][2])
{
	int topolType = 1;
	memset(link, 0, sizeof(link[0][0]) * N * 2);
	for (int i = 0; i < N; i++) 
	{
		link[i][0] = i;
		link[i][1] = i+1;
	}
	return topolType;
}

int TopolRingPair(int N, int link[N][2])
{
	int topolType = 2;
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
	for (int i = N/2+1; i < beadNumber; ++i)
	{
		link[i][0] = i - 1;
		link[i][1] = i;
	}
	link[N-1][0] = beadNumber - 1;
	link[N-1][1] = 0;
	return topolType;
}

int TopolCentromerePair(int N, int link[N][2])
{
	int topolType = 3;
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
	for (int i = ringSize+cm-1; i < beadNumber; ++i)
	{
		link[index][0] = i;
		link[index][1] = (i+1)%beadNumber;
		index = index + 1;
	}

	return topolType;
}

void Topology(int link[rodNumber][2], 
		int g[rodNumber][rodNumber])
{
	memset(link, 0, sizeof(link[0][0]) * rodNumber * 2);
	/* Define topological constriants (rods)*/
	/* TopolSingleChain(rodNumber, link); */
	TopolCentromerePair(rodNumber, link);
	/* TopolRingPair(rodNumber, link); */

	
	/*Calculate metric matrix*/
	memset(g, 0, sizeof(g[0][0]) * rodNumber * rodNumber);
	for (int i = 0; i < rodNumber; i++) 
	{
		g[i][i] = -2;
		for (int j = 0; j < rodNumber; j++) {
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
}
