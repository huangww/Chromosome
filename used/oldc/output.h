#ifndef OUTPUT_H_Y0A1TE9S
#define OUTPUT_H_Y0A1TE9S

#include "main.h"
#include <stdio.h>


void Output4lammps(FILE *lammpsFile,
		int link[numLink][2], 
		double r[numUnit][DIM]);

void OutputTopol(FILE *topolFile,
		int link[numLink][2]);

void OutputConfiguration(FILE *outputFile,
		double r[numUnit][DIM]);

void OutputTmp(double r[numUnit][DIM]);
#endif /* end of include guard: OUTPUT_H_Y0A1TE9S */
