#ifndef OUTPUT_H_Y0A1TE9S
#define OUTPUT_H_Y0A1TE9S

#include "main.h"
#include <stdio.h>


void Output4lammps(FILE *lammpsFile,
		int link[rodNumber][2], 
		double r[beadNumber][dimension]);

void OutputTopol(FILE *topolFile,
		int link[rodNumber][2]);

void OutputConfiguration(FILE *outputFile,
		double r[beadNumber][dimension]);

#endif /* end of include guard: OUTPUT_H_Y0A1TE9S */
