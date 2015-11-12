#include "main.h"
#include "utilities.h"
#include "random.h"
#include "output.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

void ForceFENE(double r[][DIM], 
        int link[][2],
        double  f[][DIM])
{
    memset(f, 0, sizeof(f[0][0])*DIM*numUnit);

    double k = 100.0;
    double R0 = 1.1;
    for (int i = 0; i < numLink; ++i)
    {
        int i0 = link[i][0];
        int i1 = link[i][1];
        double d = Distance(&r[i0][0], &r[i1][0]);
        for (int j = 0; j < DIM; ++j)
        {
            double rlogarg = 1 - (d*d)/(R0*R0);
            if (rlogarg < 0.01) rlogarg = 0.01;
            double df;
            df = -k*(r[i0][j]-r[i1][j])/rlogarg; 
            f[i0][j] += df;
            f[i1][j] -= df;
        }
    }
}

void ForceLJ(double r[][DIM],
        double f[][DIM])
{
    memset(f, 0, sizeof(f[0][0])*numUnit*DIM);

    double r0 = 1.0;
    double eps = 1.0;
    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = i+1; j < numUnit; ++j)
        {
            double d = Distance(&r[i][0],&r[j][0]);
            if (d <= r0)
            {
                double r6;
                r6 = pow(r0/d, 6);
                for (int k = 0; k < DIM; ++k)
                {
                    double df;
                    df = 48 * eps *
                        (r6 - 0.5) * r6 *
                        (r[i][k]-r[j][k])/(d*d);
                    f[i][k] += df;
                    f[j][k] -= df;
                }
            }
        }
    }

}

void ForceBend(double r[][DIM],
        int angleList[][3],
        double f[][DIM])
{
    memset(f, 0, sizeof(f[0][0])*numUnit*DIM);

    /* U = 0.5*kappa*(theta-theta0)^2 */
    double kappa = 1.0;
    double theta0 = 0;

    for (int i = 0; i < numAngle; ++i)
    {
        int i0, i1, i2;
        i0 = angleList[i][0];
        i1 = angleList[i][1];
        i2 = angleList[i][2];
        double d1[DIM],d2[DIM];
        double d1s = 0, d2s = 0, d12 = 0;

        for (int k = 0; k < DIM; ++k)
        {
            d1[k] = r[i0][k] - r[i1][k];
            d2[k] = r[i2][k] - r[i1][k];
            d1s += d1[k]*d1[k];
            d2s += d2[k]*d2[k];
            d12 += d1[k]*d2[k];
        }
        double cosAngle = d12/sqrt(d1s*d2s);
        double s = sqrt(1.0 - cosAngle*cosAngle);
        if (s < 0.001) s= 0.001;
        s = 1.0/s;

        double dtheta = acos(cosAngle) - theta0;
        double a = -kappa*dtheta*s;
        double a11 = a*cosAngle/d1s;
        double a12 = -a/sqrt(d1s*d2s);
        double a22 = a*cosAngle/d2s;

        double f1[DIM], f3[DIM];
        for (int k = 0; k < DIM; ++k)
        {
            f1[k] = a11*d1[k] + a12*d2[k];
            f3[k] = a22*d1[k] + a12*d2[k];
            f[i0][k] += f1[k];
            f[i1][k] -= f1[k] + f3[k];
            f[i2][k] += f3[k];
        }

    }

}

void ForceFriction(double v[][DIM],
        double f[][DIM])
{
    double gamma = 1.0;
    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = 0; j < DIM; ++j)
        {
            f[i][j] = - gamma*v[i][j];
        }
    }
}

void ForceRandom(unsigned long seed,
        double f[][DIM])
{
    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = 0; j < DIM; ++j)
        {
            f[i][j] = sqrt(2.0/dt) * GaussRan(seed);
        }
    }
}

void ForceWall(double r[][DIM],
        double f[][DIM])
{
    memset(f, 0, sizeof(f[0][0])*numUnit*DIM);

    double kappa = 100.0;
    double radius = 100.0, altitute = 200.0;
    for (int i = 0; i < numUnit; ++i)
    {
        double d = sqrt(r[i][1]*r[i][1] + r[i][2]*r[i][2]);
        if (d > radius)
        {
            f[i][1] = -kappa*(1.0-radius/d)*r[i][1];
            f[i][2] = -kappa*(1.0-radius/d)*r[i][2];
        }
        if (r[i][0]>altitute/2)
           f[i][0] = -kappa*(r[i][0] - altitute/2); 
        if (r[i][0]<-altitute/2)
           f[i][0] = -kappa*(r[i][0] + altitute/2); 
    }

}



void ForceAll(double r[][DIM],
        double v[][DIM],
        int link[][2],
        int angleList[][3],
        unsigned long seed,
        int step,
        double f[][DIM])
{
    memset(f, 0, sizeof(f[0][0])*numUnit*DIM);

    double fRandom[numUnit][DIM];
    double fFENE[numUnit][DIM];
    double fLJ[numUnit][DIM];
    double fBend[numUnit][DIM];
    double fFriction[numUnit][DIM];
    double fWall[numUnit][DIM];

    ForceRandom(seed, fRandom);
    ForceFENE(r, link, fFENE);
    ForceLJ(r, fLJ);
    ForceBend(r, angleList, fBend); 
    ForceFriction(v, fFriction);
    ForceWall(r, fWall);
    
    for (int i = 0; i < numUnit; ++i)
    {
        for (int j = 0; j < DIM; ++j)
        {
            f[i][j] =  fRandom[i][j] + fFriction[i][j]
                + fFENE[i][j] + fLJ[i][j]
                + fBend[i][j] + fWall[i][j];
        }
    }

    /* add external force */
    f[0][0] -= 100.0*sin(0.1*step*dt);
    for (int i = 1; i < DIM; ++i)
    {
        f[0][i] -= 100.0*r[0][i];
    }

    /* PrintFrame(f); */
    /* getchar(); */
}
