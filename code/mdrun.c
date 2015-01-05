#include "main.h"
#include "forces.h"

void MDRun(double r[beadNumber][dimension],
		int link[rodNumber][2],
		int g[rodNumber][rodNumber],
		double v0,
		unsigned long seed)
{
	/* Generate random force */
	double fb[beadNumber][dimension] = {{0}};
	GenerateFb(fb, seed);

	/* calculate pseudo-force */
	double fa[beadNumber][dimension] = {{0}};
	/* PseudoForce(r, link, g, fa); */

	/* calculate Lennard-Jones  */
	/* double flj[beadNumber][dimension] = {{0}}; */
	/* LennardJones(r, flj); */

	/* Add external force */
	/* 1) pinned SPB */
	for (int i = 0; i < dimension; ++i)
	{
		fa[0][i] = fa[0][i] - 2.0e3*r[0][i];
	}

	/* 2) external force field */
	for (int i = 1; i < beadNumber; ++i)
	{
		fa[i][0] = fa[i][0] + 1.0 * v0;
	}
	
	/* Predictive step */
	double rs[beadNumber][dimension];
	for (int i = 0; i < beadNumber; i++) 
	{
		for (int j = 0; j < dimension; j++) 
		{
			/* rs[i][j] = r[i][j] + (fa[i][j]+fb[i][j]+flj[i][j])*dt; */
			rs[i][j] = r[i][j] + (fa[i][j]+fb[i][j])*dt;
		}
	}

	/* Calculate the constraint force on rods */
	double fc[beadNumber][dimension];
	CalculateFc(r, rs, link, g, fc);

	/* Corrective step */
	for (int i = 0; i < beadNumber; i++) 
	{
		for (int j = 0; j < dimension; j++)
		{
			r[i][j] = rs[i][j] + fc[i][j]*dt;
		}
	}

}
