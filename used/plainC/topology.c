#include "main.h"
#include <stdio.h>
#include <string.h>

void TopolSingleChain(int link[rodNumber][2])
{
	memset(link, 0, sizeof(link[0][0]) * rodNumber * 2);
	for (int i = 0; i < beadNumber-1; i++) 
	{
		link[i][0] = i;
		link[i][1] = i+1;
	}
}

void TopolSingleRing(int link[rodNumber][2])
{
	memset(link, 0, sizeof(link[0][0]) * rodNumber * 2);
	for (int i = 0; i < beadNumber; i++) 
	{
		link[i][0] = i;
		link[i][1] = (i+1) % beadNumber;
	}
}

void TopolRingPair(int link[rodNumber][2])
{
	memset(link, 0, sizeof(link[0][0]) * rodNumber * 2);
	for (int i = 0; i < rodNumber/2-1; ++i)
	{
		link[i][0] = i;
		link[i][1] = i + 1;
	}
	link[rodNumber/2-1][0] = rodNumber/2 - 1;
	link[rodNumber/2-1][1] = 0;
	link[rodNumber/2][0] = 0;
	link[rodNumber/2][1] = rodNumber/2;
	for (int i = rodNumber/2+1; i < beadNumber; ++i)
	{
		link[i][0] = i - 1;
		link[i][1] = i;
	}
	link[rodNumber-1][0] = beadNumber - 1;
	link[rodNumber-1][1] = 0;

}

void Topology(int link[rodNumber][2], 
		int g[rodNumber][rodNumber])
{
	memset(link, 0, sizeof(link[0][0]) * rodNumber * 2);
	/* Define topological constriants (rods)*/
	/* TopolSingleChain(link); */
	TopolSingleRing(link);
	/* TopolRingPair(link); */
	
	/*Calculate metric matrix*/
	memset(g, 0, sizeof(g[0][0]) * rodNumber * rodNumber);
	for (int i = 0; i < rodNumber; i++) {
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
