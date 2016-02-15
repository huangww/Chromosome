#include "main.h"

int DetermineNumLink(int topologyType, int numUnit)
{
	int numLink;
	switch (topologyType) {
		case 0:
			numLink = numUnit;
			break;
		case 1:
			numLink = numUnit - 1;
			break;
		case 2:
			numLink = numUnit + 1;
			break;
		case 3:
			numLink = numUnit + 2;
			break;
		case 4:
			numLink = numUnit + 5;
			break;
		default:
			numLink = numUnit;

	}
	return numLink;
}


