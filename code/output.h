#ifndef OUTPUT_H_Y0A1TE9S
#define OUTPUT_H_Y0A1TE9S

#include "main.h"

void Output4lammps(int link[rodNumber][2], 
		double r[beadNumber][dimension],
		char *outputDir);

void OutputTopol(int link[rodNumber][2],
		char *outputDir);

#endif /* end of include guard: OUTPUT_H_Y0A1TE9S */
