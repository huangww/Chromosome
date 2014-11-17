#include <stdio.h>
#include "main.h"

void Output4lammps(int link[rodNumber][2], 
		double r[beadNumber][dimension],
		char *outputDir)

{
	FILE *file;
	char fileName[80];
	sprintf(fileName, "%slammpsIn_N%d.dat", outputDir, beadNumber);
	file = fopen(fileName,"w");
	
	fprintf(file, "#%s\n", fileName);
	fprintf(file, "\n");
	fprintf(file, "\n");
	fprintf(file, "%d\t atoms\n", beadNumber);
	fprintf(file, "%d\t bonds\n", rodNumber);
	fprintf(file, "\n");
	fprintf(file, "\n");
	fprintf(file, "1\t atom types\n");
	fprintf(file, "1\t bond types\n");
	fprintf(file, "\n");
	fprintf(file, "\n");
	fprintf(file, "# simulation box\n");
	fprintf(file, "%f %f xlo xhi\n", -beadNumber/2.0, beadNumber/2.0);
	fprintf(file, "%f %f ylo yhi\n", -beadNumber/2.0, beadNumber/2.0);
	fprintf(file, "%f %f zlo zhi\n", -beadNumber/2.0, beadNumber/2.0);
	fprintf(file, "\n");
	fprintf(file, "\n");
	fprintf(file, "Masses\n");
	fprintf(file, "\n");
	fprintf(file, "1 1\n");
	fprintf(file, "\n");
	fprintf(file, "\n");
	fprintf(file, "Atoms\n");
	fprintf(file, "\n");
	for (int i = 0; i < beadNumber; ++i)
	{
		fprintf(file, "%d\t 1\t 1\t%lf\t%lf\t%lf\n",
				i+1, r[i][0], r[i][1], r[i][2] );
	}
	fprintf(file, "\n");
	fprintf(file, "\n");
	fprintf(file, "Bonds\n");
	fprintf(file, "\n");
	for (int i = 0; i < rodNumber; ++i)
	{
		fprintf(file, "%d\t 1\t %d\t%d\n",
				i+1, link[i][0]+1, link[i][1]+1);
	}

	fclose(file);
	
}


void OutputTopol(int link[rodNumber][2],
		char *outputDir)
{
	FILE  *topolFile;
	char fileName[80];
	sprintf(fileName, "%stopol.dat", outputDir);
	topolFile = fopen(fileName,"w");
	for (int i = 0; i < rodNumber; i++) 
	{
		fprintf(topolFile, "%d\t%d\n", link[i][0], link[i][1]);
	}
	fclose(topolFile);

	
}
