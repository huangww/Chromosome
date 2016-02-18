#ifndef MAIN_H

#define MAIN_H

// include opencl lib
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

// example for a single ring
#define chainNumber  1
#define beadNumber  100
#define rodNumber  100

// example for a ring pair
// #define chainNumber  2
// #define beadNumber  9
// #define rodNumber  10

// example for a ring pair with centromere
// #define chainNumber  2
// #define beadNumber  8
// #define rodNumber  10
// #define centromere  2

#define dimension  3
#define dt  1e-4
#define pi  3.14159365359
#define maxStep  1e5

#endif /* end of include guard: MAIN_H */
