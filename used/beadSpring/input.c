#include "main.h"

int DetermineNumLink(int topologyType, int numUnit)
{
    int count;
    switch (topologyType) {
        case 0:
            count = numUnit;
            break;
        case 1:
            count = numUnit - 1;
            break;
        case 2:
            count = numUnit + 1;
            break;
        case 3:
            count = numUnit + 2;
            break;
        case 4:
            count = numUnit + 5;
            break;
        default:
            count = numUnit;

    }
    return count;
}

int DetermineNumAngle(int topologyType, int numUnit)
{
    int count;
    switch (topologyType) {
        case 0:
            count = numUnit - 1;
            break;
        case 1:
            count = numUnit - 2;
            break;
        case 2:
            count = numUnit - 1;
            break;
        case 3:
            /* need to check later */
            count = numUnit;
            break;
        case 4:
            count = numUnit - 1;
            break;
        default:
            count = numUnit - 1;
    }
    return count;
}

