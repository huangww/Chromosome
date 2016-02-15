#ifndef OUTPUT_H_Y0A1TE9S
#define OUTPUT_H_Y0A1TE9S

#include "main.h"
#include <stdio.h>


void Output4lammps(FILE *lammpsFile,
		int link[numLink][2], 
		double r[numUnit][DIM]);

void OutputBonds(FILE *topolFile,
		int link[numLink][2]);

void OutputAngles(FILE *angleFile,
		int angleList[numAngle][3]);

void OutputFrame(FILE *outputFile,
		double r[numUnit][DIM]);

void PrintFrame(double r[numUnit][DIM]);

#endif /* end of include guard: OUTPUT_H_Y0A1TE9S */
