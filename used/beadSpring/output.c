#include <stdio.h>
#include "main.h"

void Output4lammps(FILE *lammpsFile,
		int link[numLink][2], 
		double r[numUnit][DIM])
{
	fprintf(lammpsFile, "# Input file for LAMMPS\n" );
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "%d\t atoms\n", numUnit);
	fprintf(lammpsFile, "%d\t bonds\n", numLink);
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "1\t atom types\n");
	fprintf(lammpsFile, "1\t bond types\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "# simulation box\n");
	fprintf(lammpsFile, "%f %f xlo xhi\n", -numUnit/2.0, numUnit/2.0);
	fprintf(lammpsFile, "%f %f ylo yhi\n", -numUnit/2.0, numUnit/2.0);
	fprintf(lammpsFile, "%f %f zlo zhi\n", -numUnit/2.0, numUnit/2.0);
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "Masses\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "1 1\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "Atoms\n");
	fprintf(lammpsFile, "\n");
	for (int i = 0; i < numUnit; ++i)
	{
		fprintf(lammpsFile, "%d\t 1\t 1\t%lf\t%lf\t%lf\n",
				i+1, r[i][0], r[i][1], r[i][2] );
	}
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "Bonds\n");
	fprintf(lammpsFile, "\n");
	for (int i = 0; i < numLink; ++i)
	{
		fprintf(lammpsFile, "%d\t 1\t %d\t%d\n",
				i+1, link[i][0]+1, link[i][1]+1);
	}

}


void OutputBonds(FILE *topolFile,
		int link[numLink][2])
{
	for (int i = 0; i < numLink; i++) 
	{
		fprintf(topolFile, "%d\t%d\n", link[i][0], link[i][1]);
	}
}

void OutputAngles(FILE *angleFile,
		int angleList[numAngle][3])
{
	for (int i = 0; i < numAngle; i++) 
	{
		fprintf(angleFile, "%d\t%d\t%d\n", 
                        angleList[i][0], 
                        angleList[i][1], 
                        angleList[i][2]);
	}
}


void OutputFrame(FILE *outputFile,
		double r[][DIM])
{

	for (int i = 0; i < numUnit; i++) 
	{
                for (int j = 0; j < DIM; ++j)
                {
                        fprintf(outputFile, "%lf\t", r[i][j]);
                        
                }
		fprintf(outputFile, "\n");
	}

}

void PrintFrame(double r[][DIM])
{

	for (int i = 0; i < numUnit; i++) 
	{
                for (int j = 0; j < DIM; ++j)
                {
                        printf("%lf\t", r[i][j]);
                }
		printf("\n");
	}

}
