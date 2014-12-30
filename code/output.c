#include <stdio.h>
#include "main.h"

void Output4lammps(FILE *lammpsFile,
		int link[rodNumber][2], 
		double r[beadNumber][dimension])

{
	fprintf(lammpsFile, "# Input file for LAMMPS\n" );
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "%d\t atoms\n", beadNumber);
	fprintf(lammpsFile, "%d\t bonds\n", rodNumber);
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "1\t atom types\n");
	fprintf(lammpsFile, "1\t bond types\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "# simulation box\n");
	fprintf(lammpsFile, "%f %f xlo xhi\n", -beadNumber/2.0, beadNumber/2.0);
	fprintf(lammpsFile, "%f %f ylo yhi\n", -beadNumber/2.0, beadNumber/2.0);
	fprintf(lammpsFile, "%f %f zlo zhi\n", -beadNumber/2.0, beadNumber/2.0);
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "Masses\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "1 1\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "Atoms\n");
	fprintf(lammpsFile, "\n");
	for (int i = 0; i < beadNumber; ++i)
	{
		fprintf(lammpsFile, "%d\t 1\t 1\t%lf\t%lf\t%lf\n",
				i+1, r[i][0], r[i][1], r[i][2] );
	}
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "\n");
	fprintf(lammpsFile, "Bonds\n");
	fprintf(lammpsFile, "\n");
	for (int i = 0; i < rodNumber; ++i)
	{
		fprintf(lammpsFile, "%d\t 1\t %d\t%d\n",
				i+1, link[i][0]+1, link[i][1]+1);
	}

}


void OutputTopol(FILE *topolFile,
		int link[rodNumber][2])
{
	for (int i = 0; i < rodNumber; i++) 
	{
		fprintf(topolFile, "%d\t%d\n", link[i][0], link[i][1]);
	}
}

void OutputConfiguration(FILE *outputFile,
		double r[beadNumber][dimension])
{
	for (int i = 0; i < beadNumber; i++) 
	{
		fprintf(outputFile, "%lf\t%lf\t%lf\n",
				r[i][0], r[i][1], r[i][2]);
	}

}
