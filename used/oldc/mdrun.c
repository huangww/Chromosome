#include "main.h"
#include "interaction.h"
#include "output.h"
#include <string.h>

void MDRun(double r[numUnit][DIM],
        int link[numLink][2],
        int g[numLink][numLink],
        double v0,
        unsigned long seed)
{
    /* Generate random force */
    double fb[numUnit][DIM];
    GenerateFb(fb, seed);

    /* calculate pseudo-force */
    double fa[numUnit][DIM];
    memset(fa, 0, sizeof(fa[0][0])*numUnit*DIM);
    PseudoForce(r, link, g, fa);

    /* calculate Lennard-Jones  */
    /* double flj[numUnit][DIM] = {{0}}; */
    /* LennardJonesForce(r, flj); */

    /* Add external force */
    /* 1) PInned SPB */
    for (int i = 0; i < DIM; ++i)
    {
    	fa[0][i] = fa[0][i] - 1.0e3*r[0][i];
    }

    /* bond centromere */
    /* for (int i = 0; i < DIM; ++i) */
    /* { */
    /* 	fa[cm][i] = fa[cm][i] - 1.0e3*(r[cm][i]-r[cm+numUnit/2][i]); */
    /* 	fa[cm+numUnit/2][i] = fa[cm+numUnit/2][i] + */
    /* 		1.0e3*(r[cm][i]-r[cm+numUnit/2][i]); */
    /* } */

    /* 2) external force field */
    for (int i = 1; i < numUnit; ++i)
    {
        fa[i][0] = fa[i][0] + 1.0 * v0;
    }


    /* Predictive step */
    double rs[numUnit][DIM];
    for (int i = 0; i < numUnit; i++) 
    {
        for (int j = 0; j < DIM; j++) 
        {
            /* rs[i][j] = r[i][j] + (fa[i][j]+fb[i][j]+flj[i][j])*dt; */
            rs[i][j] = r[i][j] + (fa[i][j]+fb[i][j])*dt;
        }
    }

    /* Calculate the constraint force on rods */
    double fc[numUnit][DIM];
    CalculateFc(r, rs, link, g, fc);

    /* Corrective step */
    for (int i = 0; i < numUnit; i++) 
    {
        for (int j = 0; j < DIM; j++)
        {
            r[i][j] = rs[i][j] + fc[i][j]*dt;
        }
    }

}
