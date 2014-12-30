/*
    SAW_pivot implements the pivot algorithm for Monte Carlo simulations
      of self-avoiding walks
    Copyright (C) 2003 Tom Kennedy

    This file is part of SAW_pivot (version 1.0).

    SAW_pivot is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    SAW_pivot is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SAW_pivot; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    E-mail contact: tgk@math.arizona.edu
    Ordinary mail: 
	Tom Kennedy
        Mathematics Department
	University of Arizona
	Tucson, ZA 85721, USA

*/



// file name : src/public/local.h
//

//////////////////////////////////////////////////////////////
// Random number generator
//////////////////////////////////////////////////////////////

// RANDOM_NUMBER_GENERATOR determines the random number generator to be used
// RANDOM_NUMBER_GENERATOR=1 uses drand48(). This is the easy one since it 
//  is probably already installed on your system.
// RANDOM_NUMBER_GENERATOR=2 uses sprng(). This is what I have been using.
//  sprng stands for Scalable Parallel Random Number Generators.
//  It is available at http://sprng.cs.fsu.edu/ 
//  Depending on where you install it, you may need to modify the file 
//  'Makefile'. Makefile presently contains paths to where I have 
//  sprng installed (/home/tgk/sprng/include and /home/tgk/sprng/lib).
//  If you install it somewhere other than the places where the compiler
//  usually searches you will need to modify these paths appropriately.
#define RANDOM_NUMBER_GENERATOR 1

// headers for SAW routines: these include statements are a mess.
// this needs to be cleaned up
//#include "stdio.h"
#include <cstdio> 
using namespace std; 
#include "math.h"
#include <iostream>
using namespace std; 
//#include "iostream.h"
#include <cstdlib>
#include "sys/time.h"
#include "time.h"
#include "limits.h"
#include <vector>
#include <algorithm>
#include <assert.h>

#if RANDOM_NUMBER_GENERATOR==2
#define SIMPLE_SPRNG 1 // used by sprng(); this line must precede following 
#include "sprng.h" // header file for sprng() random number generator
#endif

#include "../public/variable.h"
// variable.h contains constants that must be defined outside this file. 
// Usually variable.h will be created by a run_ script 
// The file should contain the following #define statements:
// #define INNER_LOOP 1000000
// #define MAX_NPIVOT 1000 // maximum number of pivots in class walk
// #define NO_SAW 0 // when equal to 1 we just do ordinary RW 
// #define ALGORITHM_NUMBER ?
//   ALGORITHM_NUMBER should be one of following to to specify the algorithm
//   7  square lattice, 7 nontrivial lattice symmetries 
//   2 square lattice, 2 nontrivial lattice symmetries 
//   5 hexagonal lattice, 5 nontrivial lattice symmetries 
//   11 triangular lattice, 11 nontrivial lattice symmetries 
//   1 Manhattan lattice
//   47 Cubic lattice, all 47 nontrivial lattice symmetries 
//   383 4d hypercubic lattice, all 383 nontrivial lattice symmetries 
//
//   The following two lines can be included in variable.h (either one or both)
//   Each of them will change slightly how the pivot algorithm is implemented.
//   They will not change the sequence of walks generated. Depending
//   on the parameters values they will make the program run slightly
//   faster or slower.
// #define USE_FIND_SEGMENT
// #define USE_THIRD

// following #define's all depend on ALGORITHM
// NUM_SYM= number of nontrivial lattice symmetries used, identity not counted 

#if ALGORITHM_NUMBER==1
#define NUM_SYM 7 
#define SQUARE_LATTICE
#define TWO_DIMENSIONS
#endif

#if ALGORITHM_NUMBER==2
#define NUM_SYM 7 
#define SQUARE_LATTICE
#define TWO_DIMENSIONS
#endif

#if ALGORITHM_NUMBER==7
#define NUM_SYM 7 
#define SQUARE_LATTICE
#define TWO_DIMENSIONS
#endif

#if ALGORITHM_NUMBER==5
#define NUM_SYM 5 
#define HEXAGONAL_TRIANGULAR_LATTICE
#define TWO_DIMENSIONS
#endif

#if ALGORITHM_NUMBER==11
#define NUM_SYM 11 
#define HEXAGONAL_TRIANGULAR_LATTICE
#define TWO_DIMENSIONS
#endif

#if ALGORITHM_NUMBER==47
#define NUM_SYM 47
#define CUBIC_LATTICE
#define THREE_DIMENSIONS
#endif

#if ALGORITHM_NUMBER==383
#define NUM_SYM 383
#define HYPER_LATTICE
#define FOUR_DIMENSIONS
#endif

#ifdef TWO_DIMENSIONS
#define NU NU_TWOD
#endif
#ifdef THREE_DIMENSIONS
#define NU NU_THREED
#endif
#ifdef FOUR_DIMENSIONS
#define NU NU_FOURD
#endif

//////////////////////////////////////////////////////////////////////
//     simple routines  that don't involve any classes              //
//////////////////////////////////////////////////////////////////////

double mytime();

//////////////////////////////////////////////////////////////////////
//                          classes                                 //
//////////////////////////////////////////////////////////////////////

class excluded_class;
class parameters;

// Points are represented by two classes : point and rpoint.
// rpoint contains real coordinates and represents the actual coordinates of
// the point in R^d.
// point has integer coordinates. These are the coordinates of the point
// with respect to a basis. The bases used are 
// square :        e1=(1,0)  e2=(0,1)
// triangular :    e1=(1,0)  e2=(-1/2,root(3)/2)
// hexagonal :     e1=(1,0)  e2=(-1/2,root(3)/2)
// cubic :         e1=(1,0,0) e2=(0,1,0)  e3=(0,0,1)
// hypercubic :    e1=(1,0,0,0) e2=(0,1,0,0) e3=(0,0,1,0) e4=(0,0,0,1)  

class point {
 private:
  // these are coordinates wrt the two vectors that define the lattice
 #ifdef TWO_DIMENSIONS
  long x,y;
 #endif
 #ifdef THREE_DIMENSIONS
  long x,y,z;
 #endif
 #ifdef FOUR_DIMENSIONS
  long x,y,z,w;
 #endif
 public:
  inline point& operator+=(point p);
  inline point& operator-=(point p);
  point& operator*=(long c);
  void print(FILE *fptr);
  void scan(FILE *fptr);
  double distance(); // distance to the origin
  inline void euclidean_op(point* q, point* p,long isym);
  point& operator=(const point& p); 
  void zero(); // sets coords to zero 
  inline long lattice_distance();
  inline long excluded_distance(excluded_class* excluded);

  inline long coord_x();
  inline long coord_y();
  long coord_x_slow();
  long coord_y_slow();
  double real_coord_x();
  double real_coord_y();
 #ifdef TWO_DIMENSIONS
  void assign(long xx,long yy);
 #endif
 #ifdef THREE_DIMENSIONS
  inline long coord_z();
  long coord_z_slow();
  double real_coord_z();
  void assign(long xx,long yy,long zz);
 #endif
 #ifdef FOUR_DIMENSIONS
  inline long coord_z();
  long coord_z_slow();
  double real_coord_z();
  inline long coord_w();
  long coord_w_slow();
  double real_coord_w();
  void assign(long xx,long yy,long zz,long ww);
 #endif

}; 

point operator+(point p,point q);
point operator-(point p,point q);
int operator==(point p,point q);
point operator*(long c,point p);
long operator*(point p,point q);
void htab_clear(point *htab_point);

class rpoint {
 private:
  // these are actual x,y coordinates 
 #ifdef TWO_DIMENSIONS
  double x,y;
 #endif
 #ifdef THREE_DIMENSIONS
  double x,y,z;
 #endif
 #ifdef FOUR_DIMENSIONS
  double x,y,z,w;
 #endif
 public:
  void print(FILE *fptr);
  void scanf(FILE *fptr);
  void write(FILE *fptr);
  void read(FILE *fptr);
  inline rpoint& operator+=(rpoint p);
  inline rpoint& operator-=(rpoint p);
  rpoint& operator*=(double c);
  rpoint& operator/=(double c);
  void zero(); // sets coords to zero 
  rpoint& operator=(point p); 
  double distance(); 
  double distance_projection(double a, double b);
  double p_angle();

  double coord_x();
  double coord_y();
 #ifdef TWO_DIMENSIONS
  void assign(double xx,double yy);
 #endif
 #ifdef THREE_DIMENSIONS
  void assign(double xx,double yy,double zz);
  double coord_z();
 #endif
 #ifdef FOUR_DIMENSIONS
  void assign(double xx,double yy,double zz,double ww);
  double coord_z();
  double coord_w();
 #endif
}; 

rpoint operator*(double c,rpoint p);
rpoint operator*(rpoint p,double c);
rpoint operator/(rpoint p,double c);
rpoint operator+(rpoint p,rpoint q);
rpoint operator-(rpoint p,rpoint q);
rpoint operator*(rpoint p,rpoint q);
inline double inner_prod(rpoint* p,rpoint* q);

// excluded_class specifies the region the walk is excluded from.
// In the case of projecting the 3d walk onto 2d, it also specifies
// the x and y directions for the plane onto which we project.

class excluded_class {
 private:
 public:
  // region specifies the region the walk is excluded from
  // =0 : no excluded region, full plane 
  // =1 : half plane, walk must satisfy y>0
  // =2 : cut: walk not allowed to hit positive x-axis
  // =3 : quarter: walk not allowed to hit pos x-axis or pos y-axis
  // =4 : exterior of a disc: used by LERW
  int region;
  // xn and yn are unit vectors that define the x and y directions for the
  //   plane we project onto
  // We don't input xn and yn directly. They are defined by three angles 
  // theta, phi, xphi (one angle in two dims). 
  // After the angles are read in, xn and yn are computed.
  // This is done in parameters::scan()
  rpoint xn,yn;
  // for region=4, following specify the disk
  double radius;
  rpoint center;
}; 

class walk {
// stores a single walk
 private: 
  point* steps;
  long nsteps; // number of steps in the walk  
  // number of attempted pivots the walk has been through since it was a line 
  // in multiples of INNER_LOOP
  // So if INNER_LOOP=1,000,000, this is number in millions
  long niter; 
  long npivot; // number of accepted pivots not yet applied to walk
  long* ptime; // array of pivot times 
  long* igroup; // array of group elements  
  point* shift; // array of shifts
 public:
  int allocate(long n);
  void deallocate(); 
  void initialize(); 
  long get_nsteps();
  void set_nsteps(long n);
  long get_niter();
  void increment_niter();
  long get_npivot();
  void set_npivot(long n);
  point step_rval(long i);
  void assign_point(long i,point p);
  void print(FILE *fptr);
  void scan(FILE *fptr);
  void plot(FILE *fptr);
  void plot_projected(FILE *fptr,double rho);
  long pivot_strictly_saw(long pivot_loc,long isym,excluded_class excluded);
  long pivot_weakly_saw(long pivot_loc,long isym,double beta,double p,
    excluded_class excluded);
  walk& operator=(walk& w);
  void line_initialize(int direction);
  void add_pivot(long pivot_loc,long isym,point trans);
  void simplify(); //carry out the pivots, so npivot -> 0
  double turn_frac(); // fraction of sites with a turn in the walk
  void test(excluded_class* excluded); // checks walk is n.n., misses excluded
  int even(long pivot_loc); 

}; 

inline long find_segment(long itime,long npivot,long* ptime);

void setup_symmetry();


