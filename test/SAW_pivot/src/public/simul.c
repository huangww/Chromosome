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



// file: src/public/simul.c
//
// Runs the pivot algorithm to generate SAW's.
//

#include "local.h"

extern long group_inverse[NUM_SYM+1];

int main()
{
walk w;
excluded_class excluded;
long ntries,naccept,iseed;
long pivot_type,pivot_loc,iiter,inner;
char init_walk_fname[101],final_walk_fname[101],data_fname[101];
char wname[101]="walk.plt";
// sign of accept flag indicates walk was accepted (>0) or rejected (<0) 
// |accept_flag| is the number of tests done in the call to pivot routine
long accept_flag=0;
double accept_ratio=0.;
rpoint endpt;
double end_distance,total,inner_total;
FILE *fptr;
long nsteps,region,nsimplify,niter,nrecord,interaction,direction;
double beta,p;

setup_symmetry();

printf("Enter name of file containing initial walk, \n");
printf(" or 0 to have program generate a straight line for initial walk\n");
scanf("%s",init_walk_fname);
printf("Enter name of file in which to store final walk  \n");
scanf("%s",final_walk_fname);
printf("Enter name of file in which to store data  \n");
scanf("%s",data_fname);

printf("seed for random number generator \n");
scanf("%ld",&iseed);
printf("number of steps in the walk \n");
scanf("%ld",&nsteps);

printf("region \n");
scanf("%ld",&region);
excluded.region=region;
if (region < 0 || region>3)  
 {printf("bad value for region\n");
  exit(0);
 }

printf("1 for strictly or 2 for weakly or 3 for FWSAW \n");
scanf("%ld",&interaction);
switch (interaction) {
 case 1: 
 break;
 case 2: 
  printf("beta\n");
  scanf("%lf",&beta);
 break;
 case 3: 
  printf("beta\n");
  scanf("%lf",&beta);
  printf("p\n");
  scanf("%lf",&p);
 break;
 default:
  printf("interaction case not implemented in simul.c\n");
  exit(1);
 break;
} 

printf("nsimplify\n");
scanf("%ld",&nsimplify);
printf("niter\n");
scanf("%ld",&niter);
printf("nrecord\n");
scanf("%ld",&nrecord);

// nsimplify=0 => compute nsimplify
if (nsimplify==0) nsimplify=long(sqrt(double(nsteps/40)));

// seed random number generator
// if iseed is negative, we call a routine that uses system clock
#if RANDOM_NUMBER_GENERATOR==1
if (iseed<0) iseed=long(floor(mytime())); // hack
srand48(iseed);
#endif
#if RANDOM_NUMBER_GENERATOR==2
if (iseed<0) iseed=make_sprng_seed();  
init_sprng(iseed,SPRNG_DEFAULT);	
#endif

/////////////////////////////////////////////////////////////////
//                    allocate memory                          //
/////////////////////////////////////////////////////////////////

w.allocate(nsteps);

////////////////////////////////////////////////////////////
//  initialize walk to a line or  read walk from file     //
////////////////////////////////////////////////////////////

// if filename for initial walk starts with a 0 we generate a line
// for the initial walk. direction is its direction:
// 1=horizontal, 2=45 degs, 3=vertical, 4=60 degs 
direction=3;
if (init_walk_fname[0]=='0') w.line_initialize(direction);
else 
 {
  fptr=fopen(init_walk_fname,"r");
  if (fptr==NULL)
   {printf("FAILURE to open walk file %s, exiting \n",init_walk_fname);
    exit(0);
   }
  w.scan(fptr);
  fclose(fptr);
 }

// checks that initial walk is ok
w.test(&(excluded));
printf("Initial walk has turn fraction %f \n",w.turn_frac());

////////////////////////////////////////////////////////////
//      initialize running totals to zero                 //
////////////////////////////////////////////////////////////
ntries=0; 
naccept=0; 
total=0.;

for (iiter=1;iiter<=niter;iiter++) 
 {
 // inner_total will be the sum over the inner loop of the square of end 
 // to end distance 
 inner_total=0.;
 for (inner=1;inner<=INNER_LOOP;inner++)
  {
  ntries++;

  #if RANDOM_NUMBER_GENERATOR==1
  #define RNG_NAME drand48
  #endif
  #if RANDOM_NUMBER_GENERATOR==2
  #define RNG_NAME sprng
  #endif

  // a pivot location at nsteps will do nothing, but a pivot at 0 will
  pivot_loc=int(floor(nsteps*RNG_NAME()));

  ///////////////// MODEL DEPENDENT ////////////////////////////
  #if ALGORITHM_NUMBER==1
    // Manhattan lattice. Lattice operation is 
    //   reflect in +45 deg line (isym=1) if pivot point (x,y) has x+y even
    //   reflect in -45 deg line (isym=2) if pivot point (x,y) has x+y odd
    // so choice of lattice sym op is not random
  if (w.even(pivot_loc)==1) pivot_type=1;
  else pivot_type=2;
  #endif

  #if ALGORITHM_NUMBER==2
    // sq lattice with only reflections in 45 deg and 
    // -45 degree lines (this is supposedly ergodic)  
    pivot_type=int(floor(2*RNG_NAME()+1));
  #endif

  #if ALGORITHM_NUMBER==7
    // regular sq lattice with all 7 nontrivial lattice symmetries
    pivot_type=int(floor(7*RNG_NAME()+1));
  #endif

  #if ALGORITHM_NUMBER==5
    // hexagonal lattice with 5 nontrivial lattice symmetries
    pivot_type=int(floor(5*RNG_NAME()+1));
  #endif

  #if ALGORITHM_NUMBER==11
    // triangular lattice with all 11 nontrivial lattice symmetries
    pivot_type=int(floor(11*RNG_NAME()+1));
  #endif

  #if ALGORITHM_NUMBER==47
    // cubic lattice with all 47 nontrivial lattice symmetries
    pivot_type=int(floor(47*RNG_NAME()+1));
  #endif

  #if ALGORITHM_NUMBER==383
    // 4d hypercubic lattice with all 383 nontrivial lattice symmetries
    pivot_type=int(floor(383*RNG_NAME()+1));
  #endif

  switch (interaction) { 
  case 0:
    printf("interaction==0 not implemented, use NO_SAW \n");
    exit(1);
  break;
  case 1:
   accept_flag=w.pivot_strictly_saw(pivot_loc,pivot_type,excluded);
  break;
  case 2: // note we send p=0 to the routine here 
    accept_flag=w.pivot_weakly_saw(pivot_loc,pivot_type,beta,0.,excluded);
  break;
  case 3:
    accept_flag=w.pivot_weakly_saw(pivot_loc,pivot_type,beta,p,excluded);
  break;
  }
  if (w.get_npivot()==nsimplify)  w.simplify();
  if (accept_flag> 0) naccept++; 

  // compute our one random variable
  endpt=w.step_rval(nsteps);
  end_distance=endpt.distance();
  inner_total+=end_distance*end_distance;

  if (inner%PRINT_FREQ==0)
   {
    w.simplify();
    accept_ratio=100*double(naccept)/double(ntries);
    printf("%ldK iterations, turn=%f, accept ratio=%f\n", 
      (inner/1000+(iiter-1)*(INNER_LOOP/1000)),w.turn_frac(),accept_ratio);
    endpt.print(stdout);
    ntries=0;
    naccept=0;
   }

  } // end loop on inner

  total+=inner_total/double(INNER_LOOP);

  w.increment_niter();
  w.simplify();

  if (iiter%nrecord==0)
   {
    // divide by the number of times we did the inner loop 
    total/=double(nrecord);

    // record the "data"
    fptr=fopen(data_fname,"a");
    fprintf(fptr,"%14.10f\n",total);
    fclose(fptr);

    // record the walk itself
    fptr=fopen(final_walk_fname,"w");
    w.print(fptr);
    fclose(fptr);

    #ifdef TWO_DIMENSIONS
    // record the walk itself for plotting purposes
    fptr=fopen(wname,"w"); w.plot(fptr); fclose(fptr);
    #endif

    #ifdef THREE_DIMENSIONS
    // record the projected walk itself for plotting purposes
    // Note that at present we project using rho=0
    fptr=fopen(wname,"w"); w.plot_projected(fptr); fclose(fptr);
    #endif

    total=0.;
   }

} // end loop on iiter 

} // end main 

