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

// file name : src/public/lib1.c
//
// library for SAW routines - lib1.c contains routines needed for the 
// simulation itself. 

#include "local.h"

long group_product[NUM_SYM+1][NUM_SYM+1];
long group_inverse[NUM_SYM+1];

////////////////////////////////////////////////////////////////
//                      simple routines                       //
////////////////////////////////////////////////////////////////

inline double polar_angle(double x,double y)
// polar angle of (x,y) in degrees with range =[0,360]
{
double angle=atan2(y,x)*180./M_PI; if (angle<0.) angle+=360.;
return(angle);
} 

void setup_symmetry()
// this routine computes group_product[][] and group_inverse[] entries
// by trial and error 
{
long i,j,k;
point origin,p1,p2,q1,q2,temp;
#ifdef THREE_DIMENSIONS
point p3,q3;
#endif
#ifdef FOUR_DIMENSIONS
point p3,q3,p4,q4;
#endif

origin.zero();

for (i=0;i<=NUM_SYM;i++) for (j=0;j<=NUM_SYM;j++)  ///////////////////
 {group_product[i][j]=-1;
  #ifdef TWO_DIMENSIONS
  p1.assign(0,1);
  temp=p1; p1.euclidean_op(&temp,&origin,j);
  temp=p1; p1.euclidean_op(&temp,&origin,i);
  p2.assign(1,0);
  temp=p2; p2.euclidean_op(&temp,&origin,j);
  temp=p2; p2.euclidean_op(&temp,&origin,i);
  #endif
  #ifdef THREE_DIMENSIONS
  p1.assign(0,0,1);
  temp=p1; p1.euclidean_op(&temp,&origin,j);
  temp=p1; p1.euclidean_op(&temp,&origin,i);
  p2.assign(0,1,0);
  temp=p2; p2.euclidean_op(&temp,&origin,j);
  temp=p2; p2.euclidean_op(&temp,&origin,i);
  p3.assign(1,0,0);
  temp=p3; p3.euclidean_op(&temp,&origin,j);
  temp=p3; p3.euclidean_op(&temp,&origin,i);
  #endif
  #ifdef FOUR_DIMENSIONS
  p1.assign(0,0,0,1);
  temp=p1; p1.euclidean_op(&temp,&origin,j);
  temp=p1; p1.euclidean_op(&temp,&origin,i);
  p2.assign(0,0,1,0);
  temp=p2; p2.euclidean_op(&temp,&origin,j);
  temp=p2; p2.euclidean_op(&temp,&origin,i);
  p3.assign(0,1,0,0);
  temp=p3; p3.euclidean_op(&temp,&origin,j);
  temp=p3; p3.euclidean_op(&temp,&origin,i);
  p4.assign(1,0,0,0);
  temp=p4; p4.euclidean_op(&temp,&origin,j);
  temp=p4; p4.euclidean_op(&temp,&origin,i);
  #endif
  for (k=0;k<=NUM_SYM;k++)
  {
   #ifdef TWO_DIMENSIONS
   q1.assign(0,1);
   temp=q1; q1.euclidean_op(&temp,&origin,k);
   q2.assign(1,0);
   temp=q2; q2.euclidean_op(&temp,&origin,k);
   if (p1==q1 && p2==q2) group_product[i][j]=k;
   #endif
   #ifdef THREE_DIMENSIONS
   q1.assign(0,0,1);
   temp=q1; q1.euclidean_op(&temp,&origin,k);
   q2.assign(0,1,0);
   temp=q2; q2.euclidean_op(&temp,&origin,k);
   q3.assign(1,0,0);
   temp=q3; q3.euclidean_op(&temp,&origin,k);
   if (p1==q1 && p2==q2 && p3==q3) group_product[i][j]=k;
   #endif
   #ifdef FOUR_DIMENSIONS
   q1.assign(0,0,0,1);
   temp=q1; q1.euclidean_op(&temp,&origin,k);
   q2.assign(0,0,1,0);
   temp=q2; q2.euclidean_op(&temp,&origin,k);
   q3.assign(0,1,0,0);
   temp=q3; q3.euclidean_op(&temp,&origin,k);
   q4.assign(1,0,0,0);
   temp=q4; q4.euclidean_op(&temp,&origin,k);
   if (p1==q1 && p2==q2 && p3==q3 && p4==q4) group_product[i][j]=k;
   #endif
  }
  if (group_product[i][j]==-1)
   {printf("ERROR in setup_symmetry() %ld %ld\n",i,j);
    exit(1);
   } 
 }

for (i=0;i<=NUM_SYM;i++) 
 {
 group_inverse[i]=-1;
 for (k=0;k<=NUM_SYM;k++)
  {
   #ifdef TWO_DIMENSIONS
   p1.assign(0,1);
   temp=p1; p1.euclidean_op(&temp,&origin,i);
   temp=p1; p1.euclidean_op(&temp,&origin,k);
   q1.assign(0,1);
   p2.assign(1,0);
   temp=p2; p2.euclidean_op(&temp,&origin,i);
   temp=p2; p2.euclidean_op(&temp,&origin,k);
   q2.assign(1,0);
   if (p1==q1 && p2==q2) group_inverse[i]=k;
   #endif
   #ifdef THREE_DIMENSIONS
   p1.assign(0,0,1);
   temp=p1; p1.euclidean_op(&temp,&origin,i);
   temp=p1; p1.euclidean_op(&temp,&origin,k);
   q1.assign(0,0,1);
   p2.assign(0,1,0);
   temp=p2; p2.euclidean_op(&temp,&origin,i);
   temp=p2; p2.euclidean_op(&temp,&origin,k);
   q2.assign(0,1,0);
   p3.assign(1,0,0);
   temp=p3; p3.euclidean_op(&temp,&origin,i);
   temp=p3; p3.euclidean_op(&temp,&origin,k);
   q3.assign(1,0,0);
   if (p1==q1 && p2==q2 && p3==q3) group_inverse[i]=k;
   #endif
   #ifdef FOUR_DIMENSIONS
   p1.assign(0,0,0,1);
   temp=p1; p1.euclidean_op(&temp,&origin,i);
   temp=p1; p1.euclidean_op(&temp,&origin,k);
   q1.assign(0,0,0,1);
   p2.assign(0,0,1,0);
   temp=p2; p2.euclidean_op(&temp,&origin,i);
   temp=p2; p2.euclidean_op(&temp,&origin,k);
   q2.assign(0,0,1,0);
   p3.assign(0,1,0,0);
   temp=p3; p3.euclidean_op(&temp,&origin,i);
   temp=p3; p3.euclidean_op(&temp,&origin,k);
   q3.assign(0,1,0,0);
   p4.assign(0,1,0,0);
   temp=p4; p4.euclidean_op(&temp,&origin,i);
   temp=p4; p4.euclidean_op(&temp,&origin,k);
   q4.assign(0,1,0,0);
   if (p1==q1 && p2==q2 && p3==q3 && p4==q4) group_inverse[i]=k;
   #endif
  }
  if (group_inverse[i]==-1)
   {printf("ERROR in setup_symmetry()\n");
    exit(1);
   } 
 }
return;

} // end setup_symmetry()

double mytime() 
// returns time in secs using gettimeofday
// only used for timing, does not affect the computation
{
#include "sys/time.h"
double temp;
struct timeval tp;
struct timezone tzp;
gettimeofday(&tp,&tzp);
temp=tp.tv_sec+tp.tv_usec/1.e6;
return(temp);
} // end mytime()

////////////////////////////////////////////////////////////////
//             member functions for class rpoint               //
////////////////////////////////////////////////////////////////

void rpoint::print(FILE *fptr) 
{
#ifdef TWO_DIMENSIONS
fprintf(fptr,"%f %f \n",x,y);
#endif
#ifdef THREE_DIMENSIONS
fprintf(fptr,"%f %f %f\n",x,y,z);
#endif
#ifdef FOUR_DIMENSIONS
fprintf(fptr,"%f %f %f %f\n",x,y,z,w);
#endif
}

void rpoint::scanf(FILE *fptr) 
{
#ifdef TWO_DIMENSIONS
fscanf(fptr,"%lf %lf\n",&x,&y);
#endif
#ifdef THREE_DIMENSIONS
fscanf(fptr,"%lf %lf %lf\n",&x,&y,&z);
#endif
#ifdef FOUR_DIMENSIONS
fscanf(fptr,"%lf %lf %lf %lf\n",&x,&y,&z,&w);
#endif
}

void rpoint::write(FILE *fptr) 
{
fwrite(&x,sizeof(double),1,fptr);
fwrite(&y,sizeof(double),1,fptr);
#ifdef THREE_DIMENSIONS
fwrite(&z,sizeof(double),1,fptr);
#endif
#ifdef FOUR_DIMENSIONS
fwrite(&z,sizeof(double),1,fptr);
fwrite(&w,sizeof(double),1,fptr);
#endif
}

void rpoint::read(FILE *fptr)
{
fread(&x,sizeof(double),1,fptr);
fread(&y,sizeof(double),1,fptr);
#ifdef THREE_DIMENSIONS
fread(&z,sizeof(double),1,fptr);
#endif
#ifdef FOUR_DIMENSIONS
fread(&z,sizeof(double),1,fptr);
fread(&w,sizeof(double),1,fptr);
#endif
}

inline rpoint& rpoint::operator+=(rpoint p) 
{
x+=p.x;
y+=p.y;
#ifdef THREE_DIMENSIONS
z+=p.z;
#endif
#ifdef FOUR_DIMENSIONS
z+=p.z;
w+=p.w;
#endif
return *this;
} 

inline rpoint& rpoint::operator-=(rpoint p) 
{
x-=p.x;
y-=p.y;
#ifdef THREE_DIMENSIONS
z-=p.z;
#endif
#ifdef FOUR_DIMENSIONS
z-=p.z;
w-=p.w;
#endif
return *this;
} 

rpoint& rpoint::operator/=(double c) 
{
x/=c;
y/=c;
#ifdef THREE_DIMENSIONS
z/=c;
#endif
#ifdef FOUR_DIMENSIONS
z/=c;
w/=c;
#endif
return *this;
} 

rpoint& rpoint::operator*=(double c) 
{
x*=c;
y*=c;
#ifdef THREE_DIMENSIONS
z*=c;
#endif
#ifdef FOUR_DIMENSIONS
z*=c;
w*=c;
#endif
return *this;
} 

inline rpoint operator+(rpoint p,rpoint q) 
{
rpoint r=p; 
return r+=q;
} 

inline rpoint operator-(rpoint p,rpoint q) 
{
rpoint r=p; 
return r-=q;
} 


void rpoint::zero() 
{
x=0.; y=0.;
#ifdef THREE_DIMENSIONS
z=0.;
#endif
#ifdef FOUR_DIMENSIONS
z=0.;
w=0.;
#endif
} 


#ifdef TWO_DIMENSIONS
void rpoint::assign(double xx,double yy) 
{
x=xx; y=yy;
} 
#endif

#ifdef THREE_DIMENSIONS
void rpoint::assign(double xx,double yy,double zz) 
{
x=xx; y=yy; z=zz; 
} 
#endif

#ifdef FOUR_DIMENSIONS
void rpoint::assign(double xx,double yy,double zz,double ww) 
{
x=xx; y=yy; z=zz; w=ww;
} 
#endif

double rpoint::coord_x() 
{
return(x);
}

double rpoint::coord_y() 
{
return(y);
}

#ifdef THREE_DIMENSIONS
double rpoint::coord_z() 
{
return(z);
}
#endif

#ifdef FOUR_DIMENSIONS
double rpoint::coord_z() 
{
return(z);
}

double rpoint::coord_w() 
{
return(w);
}
#endif

rpoint operator*(double c,rpoint p) 
{
rpoint r=p; 
return r*=c;
} 

rpoint operator*(rpoint p ,double c)
{
rpoint r=p; 
return r*=c;
} 

inline double inner_prod(rpoint* p,rpoint* q) 
// computes inner product of the two rpoints
{
#ifdef TWO_DIMENSIONS
 return((*p).coord_x()*(*q).coord_x()+(*p).coord_y()*(*q).coord_y());
#endif
#ifdef THREE_DIMENSIONS
 return((*p).coord_x()*(*q).coord_x()+(*p).coord_y()*(*q).coord_y()
                             +(*p).coord_z()*(*q).coord_z());
#endif
#ifdef FOUR_DIMENSIONS
 return((*p).coord_x()*(*q).coord_x()+(*p).coord_y()*(*q).coord_y()
       +(*p).coord_z()*(*q).coord_z()+(*p).coord_w()*(*q).coord_w());
#endif
} 

rpoint operator/(rpoint p ,double c)
{
rpoint r=p; 
return r/=c;
} 

double rpoint::distance()  
// returns euclidean distance from point to origin
{
#ifdef TWO_DIMENSIONS
return(sqrt(x*x+y*y));
#endif
#ifdef THREE_DIMENSIONS
return(sqrt(x*x+y*y+z*z));
#endif
#ifdef FOUR_DIMENSIONS
return(sqrt(x*x+y*y+z*z+w*w));
#endif
} // end rpoint::distance()

double rpoint::distance_projection(double a, double b)  
// projects point onto plane spanned by y axis and (a,b), where (a,b) has 
// length 1, then returns euclidean distance from projection to origin
{
#ifdef THREE_DIMENSIONS
double c=x*a+z*b;
return(sqrt(c*c+y*y));
#endif
printf("error: distance_projection called, not in 3 dimensions\n"); 
exit(0);
return(0.);
} // end rpoint::distance_projection() 

double rpoint::p_angle() 
{
#ifdef TWO_DIMENSIONS
double angle=polar_angle(x,y);
return(angle);
#endif
printf("\nERROR: rpoint::p_angle() called in dimension other than 2\n");
exit(0);
return(0.); // to avoid compiler warning 
} // end rpoint::p_angle() 

rpoint& rpoint::operator=(point p)  
{
x=p.real_coord_x(); 
y=p.real_coord_y(); 
#ifdef THREE_DIMENSIONS
z=p.real_coord_z(); 
#endif
#ifdef FOUR_DIMENSIONS
z=p.real_coord_z(); 
w=p.real_coord_w(); 
#endif
return(*this);
} // end point::operator=()

////////////////////////////////////////////////////////////////
//             member functions for class point               //
////////////////////////////////////////////////////////////////

inline long point::coord_x() {return(x);} 
inline long point::coord_y() {return(y);} 
#ifdef THREE_DIMENSIONS
inline long point::coord_z() {return(z);} 
#endif
#ifdef FOUR_DIMENSIONS
inline long point::coord_z() {return(z);} 
inline long point::coord_w() {return(w);} 
#endif

// since the above routines are "inline", they can only be called from this 
// file. The following versions are for use in other files.

long point::coord_x_slow() {return(x);}
long point::coord_y_slow() {return(y);}
#ifdef THREE_DIMENSIONS
long point::coord_z_slow() {return(z);}
#endif
#ifdef FOUR_DIMENSIONS
long point::coord_z_slow() {return(z);}
long point::coord_w_slow() {return(w);}
#endif

#ifdef TWO_DIMENSIONS
void point::assign(long xx,long yy) {x=xx; y=yy;} 
#endif
#ifdef THREE_DIMENSIONS
void point::assign(long xx,long yy,long zz) {x=xx; y=yy; z=zz;} 
#endif
#ifdef FOUR_DIMENSIONS
void point::assign(long xx,long yy,long zz,long ww) {x=xx; y=yy; z=zz; w=ww;} 
#endif

inline point& point::operator+=(point p) 
{
x+=p.x;
y+=p.y;
#ifdef THREE_DIMENSIONS
z+=p.z;
#endif
#ifdef FOUR_DIMENSIONS
z+=p.z;
w+=p.w;
#endif
return *this;
} 

inline point& point::operator-=(point p) 
{
x-=p.x;
y-=p.y;
#ifdef THREE_DIMENSIONS
z-=p.z;
#endif
#ifdef FOUR_DIMENSIONS
z-=p.z;
w-=p.w;
#endif
return *this;
} 

inline point operator+(point p,point q)
{
point r=p; 
return r+=q;
} 

inline point operator-(point p,point q)
{
point r=p; 
return r-=q;
} 

void point::zero() 
{
x=0; y=0;
#ifdef THREE_DIMENSIONS
z=0;
#endif
#ifdef FOUR_DIMENSIONS
z=0;
w=0;
#endif
} 

double point::distance()  
// returns euclidean distance from point to origin. Used by walk::turn_frac() 
{
rpoint rp;
rp=(*this);
return(rp.distance());
} // end point::distance()

void point::print(FILE *fptr) 
{
#ifdef TWO_DIMENSIONS
fprintf(fptr,"%ld %ld \n",x,y);
#endif
#ifdef THREE_DIMENSIONS
fprintf(fptr,"%ld %ld %ld\n",x,y,z);
#endif
#ifdef FOUR_DIMENSIONS
fprintf(fptr,"%ld %ld %ld %ld\n",x,y,z,w);
#endif
}

double point::real_coord_x() 
// returns real x coordinate of point, rather than coords wrt lattice structure
{
#ifdef SQUARE_LATTICE
  return(x);
#endif

#ifdef HEXAGONAL_TRIANGULAR_LATTICE
  return(x-0.5*y);
#endif

#ifdef CUBIC_LATTICE
  return(x);
#endif

#ifdef HYPER_LATTICE
  return(x);
#endif

printf("error in real_coord_x - bad case \n");
exit(1);
} // end real_coord_x()

double point::real_coord_y()
// returns real y coordinate of point, rather than coords wrt lattice structure
{
#ifdef SQUARE_LATTICE
  return(y);
#endif

#ifdef HEXAGONAL_TRIANGULAR_LATTICE
  return(0.5*sqrt(3.)*y);
#endif

#ifdef CUBIC_LATTICE
  return(y);
#endif

#ifdef HYPER_LATTICE
  return(y);
#endif

printf("error in real_coord_y - bad case \n");
exit(1);
} // end real_coord_y()

#ifdef THREE_DIMENSIONS
double point::real_coord_z()  
// returns real z coordinate of point, rather than coords wrt lattice structure
{
#ifdef CUBIC_LATTICE
  return(z);
#endif

#ifdef HYPER_LATTICE
  return(z);
#endif

printf("error in real_coord_z - bad case \n");
exit(1);

} // end real_coord_z()
#endif

#ifdef FOUR_DIMENSIONS
double point::real_coord_z() 
// returns real z coordinate of point, rather than coords wrt lattice structure
{
#ifdef HYPER_LATTICE
  return(z);
#endif

printf("error in real_coord_z - bad case \n");
exit(1);

} // end real_coord_z()

double point::real_coord_w() 
// returns real w coordinate of point, rather than coords wrt lattice structure
{
#ifdef HYPER_LATTICE
  return(w);
#endif

printf("error in real_coord_w - bad case \n");
exit(1);

} // end real_coord_w()
#endif

point& point::operator*=(long c) 
{
x*=c;
y*=c;
#ifdef THREE_DIMENSIONS
z*=c;
#endif
#ifdef FOUR_DIMENSIONS
z*=c;
w*=c;
#endif
return *this;
} 

point operator*(long c,point p) 
{
point r=p; 
return r*=c;
} 

void point::scan(FILE *fptr) 
{
#ifdef TWO_DIMENSIONS
fscanf(fptr,"%ld %ld",&x,&y);
#endif
#ifdef THREE_DIMENSIONS
fscanf(fptr,"%ld %ld %ld",&x,&y,&z);
#endif
#ifdef FOUR_DIMENSIONS
fscanf(fptr,"%ld %ld %ld %ld",&x,&y,&z,&w);
#endif
}

int operator==(point p,point q)
{
#ifdef TWO_DIMENSIONS
if (p.coord_x()!=q.coord_x() || p.coord_y()!=q.coord_y()) return(0);
#endif
#ifdef THREE_DIMENSIONS
if (p.coord_x()!=q.coord_x() || p.coord_y()!=q.coord_y()
   || p.coord_z()!=q.coord_z()) return(0);
#endif
#ifdef FOUR_DIMENSIONS
if (p.coord_x()!=q.coord_x() || p.coord_y()!=q.coord_y()
 || p.coord_z()!=q.coord_z() || p.coord_w()!=q.coord_w()) return(0);
#endif
return(1);
} 

point& point::operator=(const point& p) 
{
x=p.x;
y=p.y;
#ifdef THREE_DIMENSIONS
z=p.z;
#endif
#ifdef FOUR_DIMENSIONS
z=p.z;
w=p.w;
#endif
return(*this);
} // end point::operator=()


inline long point::excluded_distance(excluded_class* excluded)
// Computes "distance" from point to excluded region.
// NB: this routine assumes the point is not inside the excluded region.
//  No error checking is done
// Ideally, distance is the minimum number of nearest neighbor steps needed 
// to get from the point to the excluded region, but we actually use 
// the l2 distance which is a lower bound on the number of nearest neighbor
// steps. We can round this distance up to the nearest integer. Hence
// the use of ceil() in the above. We use ceil(dist-1.e-10) to be sure
// a number doesn't get wrongly rounded up due to numerical error.
//
// cut plane: the excluded region is only the positive x-axis which may not
// include any points at all. So we widen the excluded region to be the 
// strip x>0, -0.5<=y<=0.5. Every bond that crosses the positive x-axis has
// an endpoint in this strip.
{
 double xx,yy,zz;
 zz=0.; // solely to avoid compiler warning
 rpoint rp;
 rp=*this;
 if ((*excluded).region>=0) 
  {xx=rp.coord_x();
   yy=rp.coord_y();
  }
 else 
  {xx=inner_prod(&rp,&((*excluded).xn));
   yy=inner_prod(&rp,&((*excluded).yn));
  }

 switch (abs((*excluded).region)) {
  case 0:  // no excluded region, full plane 
   return(LONG_MAX/2);
  break;
  case 1:  // half plane, walk must satisfy y>0
   if (yy>0.) return(long(ceil(yy-1.e-10))); 
   else return(0);
  break;
  case 2:  // cut: walk not allowed to hit positive x-axis
   if (xx<0.) return(long(ceil(sqrt(xx*xx+yy*yy)-1.e-10))); 
   else 
    {if (fabs(yy)<=0.5) return(0);
     else return(long(ceil(fabs(yy)-1.e-10)));
    }
  break;
  case 3:  // quarter: walk not allowed to hit pos x-axis or pos y-axis
   if (yy>0 && xx>0) 
    {if (yy>xx) return(long(ceil(xx-1.e-10))); 
     else return(long(ceil(yy-1.e-10)));
    }
   else return(0);
  break;
  case 4: // other quarter plane: walk restricted to z>=0, y>0
   #ifdef THREE_DIMENSIONS
   zz=rp.coord_z();
   if (yy>0 && zz>=0) 
    {if (yy>zz+1) return(long(ceil(zz+1.-1.e-10)));//+1 cause excluded is z<=-1
     else return(long(ceil(yy-1.e-10)));
    }
   else return(0);
   #endif
   printf("ERROR : excluded=4 and number of dimensions is 2\n");
   exit(1);
  break;
  default:
   printf("excluded case not implemented in excluded_distance \n");
   exit(10);
  break;
 }

printf("algorithm not implemented in excluded_distance \n");
exit(1);
} // end point::excluded_distance()

inline long point::lattice_distance() 
// computes minimum number of steps to get from point with integer
// coefs (x,y) to origin, or a lower bound on this number.
// NB : For hexagonal lattice we use the separation on the triangular 
//   lattice. This is valid since the triangular latttice separation will
//   be a lower bound on the hexagonal lattice separation, but it 
//   may not be the fastest.
{
long separation=0;

// in this routine we use if then statement in place of abs - its faster

#ifdef SQUARE_LATTICE
  if (x<0) separation-=x;
  else separation+=x;
  if (y<0) separation-=y;
  else separation+=y;
  return(separation);
#endif

// hexagonal just uses triangular lattice separation
#ifdef HEXAGONAL_TRIANGULAR_LATTICE
  // Point can be written in three ways:
  // x e1 + y e2
  // x (e2+e1) + (y-x) e2
  // (x-y) e1 + y (e1+e2)
  // dif1,dif2, dif3 are |  | of x,y,x-y. The min distance is min of 
  // dif1+dif2, dif1+dif3, dif2+dif3.

  long temp,dif1,dif2,dif3;
  if (x>0) dif1=x;
  else dif1=-x;
  if (y>0) dif2=y;
  else dif2=-y;
  dif3=x-y;
  if (dif3<0) dif3=-dif3;

  separation=dif1+dif2;
  temp=dif3+dif1;
  if (temp<separation) separation=temp;
  temp=dif3+dif2;
  if (temp<separation) separation=temp;
  return(separation);
#endif

#ifdef CUBIC_LATTICE
  if (x<0) separation-=x;
  else separation+=x;
  if (y<0) separation-=y;
  else separation+=y;
  if (z<0) separation-=z;
  else separation+=z;
  return(separation);
#endif

#ifdef HYPER_LATTICE
  if (x<0) separation-=x;
  else separation+=x;
  if (y<0) separation-=y;
  else separation+=y;
  if (z<0) separation-=z;
  else separation+=z;
  if (w<0) separation-=w;
  else separation+=w;
  return(separation);
#endif

printf("error in lattice_distance - bad case \n");
exit(1);

} // end point::lattice_distance()


inline void point::euclidean_op(point* q,point* p,long isym) 
// isym specifies the group element g. Routine compute gq+p
{
  
#ifdef SQUARE_LATTICE
  switch(isym) {
   case 0: x=(*q).x;  y=(*q).y;  break;  // identity
   case 1: x=(*q).y;  y=(*q).x;  break;  // reflect in 45 deg line
   case 2: x=-(*q).y; y=-(*q).x; break;  // reflect in -45 deg line
   case 3: x=(*q).y;  y=-(*q).x; break;  // rotate by 90 clockwise
   case 4: x=-(*q).x; y=-(*q).y; break;  // rotate by 180 clockwise
   case 5: x=-(*q).y; y=(*q).x;  break;  // rotate by 270 clockwise
   case 6: x= (*q).x; y=-(*q).y; break;  // reflect in x-axis 
   case 7: x=-(*q).x; y=(*q).y;  break;  // reflect in y-axis 
   default:
    printf("bad case square in euclidean_op, isym=%ld \n",isym);
    exit(0);
   break;
   } // end switch on isym
#endif

#ifdef HEXAGONAL_TRIANGULAR_LATTICE
  // rotations here are counterclockwise
  switch(isym) {
   case 0:   x=(*q).x;      y=(*q).y;          break;  // identity

   case 1:   x=-(*q).y;     y= (*q).x-(*q).y; break;   // R^2, rotation by 120
   case 2:   x=-(*q).x+(*q).y; y=-(*q).x;     break;   // R^4, rotation by 240

   case 3:   x= (*q).x-(*q).y; y=-(*q).y;     break;   // F_1, flip about e1 
   case 4:   x=-(*q).x;     y=-(*q).x+(*q).y; break;   // F_2, flip about e2
   case 5:   x= (*q).y;     y= (*q).x;     break;   // F_3, flip about e1+e2

   case 6:   x= (*q).x-(*q).y; y= (*q).x;     break;   // R, rotation by 60
   case 7:   x=-(*q).x;     y=-(*q).y;     break;   // R^3, rotation by 180
   case 8:   x= (*q).y;     y=-(*q).x+(*q).y; break;   // R^5, rotation by 300

   case 9:   x=-(*q).x+(*q).y; y= (*q).y;     break;   // -F_1
   case 10:  x= (*q).x;     y= (*q).x-(*q).y; break;   // -F_2
   case 11:  x=-(*q).y;     y=-(*q).x;     break;   // -F_3

   default:
    printf("bad case hex_tri in euclidean_op \n");
    exit(0);
   break;
   } // end switch on isym
#endif

#if ALGORITHM_NUMBER==47
  switch(isym%6) {
   case 0:   x=(*q).x;     y=(*q).y;     z=(*q).z;  break;   // identity
   case 1:   x=(*q).y;     y=(*q).x;     z=(*q).z;  break;   // (xy)
   case 2:   x=(*q).x;     y=(*q).z;     z=(*q).y;  break;   // (yz)
   case 3:   x=(*q).z;     y=(*q).y;     z=(*q).x;  break;   // (xz)
   case 4:   x=(*q).z;     y=(*q).x;     z=(*q).y;  break;   // (xyz)
   case 5:   x=(*q).y;     y=(*q).z;     z=(*q).x;  break;   // (xzy) 
   default:  
    printf("bad case cubic isym_perm in euclidean_op \n");
    exit(0);
   break;
   } // end switch on isym_perm
  switch(isym/6) {
   case 0:                        break;   // +++
   case 1:                 z=-z;  break;   // ++-
   case 2:          y=-y;         break;   // +-+
   case 3:          y=-y;  z=-z;  break;   // +--
   case 4:   x=-x;                break;   // -++
   case 5:   x=-x;         z=-z;  break;   // -+-
   case 6:   x=-x;  y=-y;         break;   // --+
   case 7:   x=-x;  y=-y;  z=-z;  break;   // ---
   default:
    printf("bad case cubic isym_sign in euclidean_op \n");
    exit(0);
   break;
   } // end switch on isym_sign
#endif

#if ALGORITHM_NUMBER==383
  switch(isym%24) {
   case 0 :   x=(*q).x;   y=(*q).y;   z=(*q).z;   w=(*q).w;  break; 
   case 1 :   x=(*q).x;   y=(*q).y;   z=(*q).w;   w=(*q).z;  break; 
   case 2 :   x=(*q).x;   y=(*q).z;   z=(*q).y;   w=(*q).w;  break; 
   case 3 :   x=(*q).x;   y=(*q).z;   z=(*q).w;   w=(*q).y;  break; 
   case 4 :   x=(*q).x;   y=(*q).w;   z=(*q).y;   w=(*q).z;  break; 
   case 5 :   x=(*q).x;   y=(*q).w;   z=(*q).z;   w=(*q).y;  break; 

   case 6 :   x=(*q).y;   y=(*q).x;   z=(*q).w;   w=(*q).z;  break; 
   case 7 :   x=(*q).y;   y=(*q).x;   z=(*q).z;   w=(*q).w;  break; 
   case 8 :   x=(*q).y;   y=(*q).z;   z=(*q).x;   w=(*q).w;  break; 
   case 9 :   x=(*q).y;   y=(*q).z;   z=(*q).w;   w=(*q).x;  break; 
   case 10:   x=(*q).y;   y=(*q).w;   z=(*q).x;   w=(*q).z;  break; 
   case 11:   x=(*q).y;   y=(*q).w;   z=(*q).z;   w=(*q).x;  break; 

   case 12:   x=(*q).z;   y=(*q).x;   z=(*q).y;   w=(*q).w;  break; 
   case 13:   x=(*q).z;   y=(*q).x;   z=(*q).w;   w=(*q).y;  break; 
   case 14:   x=(*q).z;   y=(*q).y;   z=(*q).x;   w=(*q).w;  break; 
   case 15:   x=(*q).z;   y=(*q).y;   z=(*q).w;   w=(*q).x;  break; 
   case 16:   x=(*q).z;   y=(*q).w;   z=(*q).x;   w=(*q).y;  break; 
   case 17:   x=(*q).z;   y=(*q).w;   z=(*q).y;   w=(*q).x;  break; 

   case 18:   x=(*q).w;   y=(*q).x;   z=(*q).y;   w=(*q).z;  break; 
   case 19:   x=(*q).w;   y=(*q).x;   z=(*q).z;   w=(*q).y;  break; 
   case 20:   x=(*q).w;   y=(*q).y;   z=(*q).x;   w=(*q).z;  break; 
   case 21:   x=(*q).w;   y=(*q).y;   z=(*q).z;   w=(*q).x;  break; 
   case 22:   x=(*q).w;   y=(*q).z;   z=(*q).x;   w=(*q).y;  break; 
   case 23:   x=(*q).w;   y=(*q).z;   z=(*q).y;   w=(*q).x;  break; 
   default:  
    printf("bad case cubic isym_perm in euclidean_op \n");
    exit(0);
   break;
   } // end switch on isym_perm
  switch(isym/24) {
   case 0:                              break; 
   case 1:                        w=-w; break; 
   case 2:                 z=-z;        break; 
   case 3:                 z=-z;  w=-w; break; 
   case 4:          y=-y;               break; 
   case 5:          y=-y;         w=-w; break; 
   case 6:          y=-y;  z=-z;        break; 
   case 7:          y=-y;  z=-z;  w=-w; break; 
   case 8:   x=-x;                      break; 
   case 9:   x=-x;                w=-w; break; 
   case 10:  x=-x;         z=-z;        break; 
   case 11:  x=-x;         z=-z;  w=-w; break; 
   case 12:  x=-x;  y=-y;               break; 
   case 13:  x=-x;  y=-y;         w=-w; break; 
   case 14:  x=-x;  y=-y;  z=-z;        break; 
   case 15:  x=-x;  y=-y;  z=-z;  w=-w; break; 
   default:
    printf("bad case hyper isym_sign in euclidean_op \n");
    exit(0);
   break;
   } // end switch on isym_sign
#endif

x+=(*p).x;
y+=(*p).y;
#ifdef THREE_DIMENSIONS
z+=(*p).z;
#endif
#ifdef FOUR_DIMENSIONS
z+=(*p).z;
w+=(*p).w;
#endif

} // end point::euclidean_op()


////////////////////////////////////////////////////////////////
//             member functions for class walk                //
////////////////////////////////////////////////////////////////

long walk::get_nsteps() {return(nsteps);} 
void walk::set_nsteps(long n) {nsteps=n;} 
long walk::get_niter() {return(niter);} 
void walk::increment_niter() {niter++;} 
long walk::get_npivot() {return(npivot);} 
void walk::set_npivot(long n) {npivot=n;} 
point walk::step_rval(long i) 
// This routine computes the ith point on the walk from the complicated 
// data structure. Often we duplicate this code for speed rather than call this
 {point temp;
  long iseg;
  iseg=find_segment(i,npivot,ptime);
  temp.euclidean_op(steps+i,shift+iseg,igroup[iseg]);
  return(temp);
 } // end step_rval


void walk::initialize() 
// NB this does not create a walk. It only initializes things used by the 
// implicit pivot data structure
{
npivot=0;
igroup[0]=0;
shift[0].zero();
// segment number iseg corresponds to times [ptime[iseg],ptime[iseg+1])
// So when npivot=0 we should have : 
ptime[0]=0;
ptime[1]=nsteps+1;
}

void walk::line_initialize(int direction)
{ // initializes walk to be straight line
  // 1 is horizontal 
  // 2 is 45 degs
  // 3 is vertical
  // 4 is 60 degs
long i,i1,i2,i3,i4;

initialize();
niter=0;
i1=0; i2=0; i3=0; i4=0;
for (i=0;i<=nsteps;i++) 
 {
  #ifdef TWO_DIMENSIONS
    steps[i].assign(i1,i2);
  #endif
  #ifdef THREE_DIMENSIONS
    steps[i].assign(i1,i2,i3);
  #endif
  #ifdef FOUR_DIMENSIONS
    steps[i].assign(i1,i2,i3,i4);
  #endif

  #ifdef SQUARE_LATTICE
    switch(direction) {
     case 1: 
       i1++; 
     break; 
     case 2: 
       switch (i%2) {
        case 0: i1++; break;
        case 1: i2++; break;
       }
     break; 
     case 3: 
       i2++; 
     break; 
     default: printf("bad case in line_initialize() a\n"); exit(0); break;
    }
  #endif

  #if ALGORITHM_NUMBER==5
    switch(direction) {
     case 1: 
       switch (i%4) {
        case 0: i1++; break;
        case 1: i1++; i2++; break;
        case 2: i1++; break;
        case 3: i2--; break;
       }
     break; 
     case 3: 
       switch (i%2) {
        case 0: i2++; break;
        case 1: i1++; i2++; break;
       }
     break; 
     default: printf("bad case in line_initialize() b\n"); exit(0); break;
    }
  #endif

  #if ALGORITHM_NUMBER==11
    switch(direction) {
     case 1: 
       i1++; 
     break; 
     case 3: 
       switch (i%2) {
        case 0: i2++; break;
        case 1: i1++; i2++; break;
       }
     break; 
     case 4: 
       i1++; i2++;
     break; 
     default: printf("bad case in line_initialize() c\n"); exit(0); break;
    }
  #endif

  #ifdef CUBIC_LATTICE
    switch(direction) {
     case 1: i1++; break; 
     case 2: if (i%2==0) i2++; else i1++; break; 
     case 3: i2++; break; 
     default: printf("bad case in line_initialize() d\n"); exit(0); break;
    }
  #endif

  #ifdef HYPER_LATTICE
    switch(direction) {
     case 1: i1++; break; 
     case 2: if (i%2==0) i2++; else i1++; break; 
     case 3: i2++; break; 
     default: printf("bad case in line_initialize() e\n"); exit(0); break;
    }
  #endif

 }
} // end line_initialize()

int walk::allocate(long n) 
{
 steps = new point[n+1];
 // this may be one larger than needed 
 ptime = new long[MAX_NPIVOT+2]; 
 igroup = new long[MAX_NPIVOT+1]; 
 shift = new point[MAX_NPIVOT+1]; 
 if (steps==0 || ptime==0 || igroup==0 || shift==0) //allocation failed
   {printf("allocate error in walk::allocate \n");
    return(0);
   }
 nsteps=n; 
 return(1);
} // end walk::allocate()

void walk::deallocate() 
{
 delete [] steps;
 delete [] ptime;
 delete [] igroup;
 delete [] shift;
}; 

void walk::print(FILE *fptr) 
// this prints the coordinates of the walk wrt the two lattice vectors
// not the actual x,y coords
// If npivot is > 0, it always prints the arrays ptime, igroup, shift.
// This is for debugging. It is incompatible with walk::scanf()
{
long i;
fprintf(fptr,"%ld %ld\n",nsteps,niter);
for (i=0;i<=nsteps;i++) steps[i].print(fptr);
if (npivot>0 || igroup[0]!=0) 
 {fprintf(fptr,"npivot=%ld\n",npivot);
  fprintf(fptr,"iseg     ptime   igroup    shift \n");
  for (i=0;i<=npivot;i++) 
    {
      fprintf(fptr,"%4ld  %8ld   %2ld    ",i,ptime[i],igroup[i]); 
      shift[i].print(fptr);
    } // end loop on i 
 } // end if (npivot>0)
} // end walk::print()

void walk::scan(FILE *fptr) 
{
long i;
fscanf(fptr,"%ld %ld",&nsteps,&niter);
for (i=0;i<=nsteps;i++) steps[i].scan(fptr);
initialize();
} // end walk::scan()

void walk::plot(FILE *fptr) 
// this prints out the actual x,y,z,... coordinates of the walk so it 
// can be plotted
{ 
long i;
if (npivot!=0) 
  {printf("ERROR: walk::print() called with npivot!=0 \n");
   exit(1L);
  } 

for (i=0;i<=nsteps;i++)  
#ifdef TWO_DIMENSIONS
  fprintf(fptr,"%f %f \n",steps[i].real_coord_x(),
	  steps[i].real_coord_y());
#endif
#ifdef THREE_DIMENSIONS
  fprintf(fptr,"%f %f %f\n",steps[i].real_coord_x(),
	  steps[i].real_coord_y(),steps[i].real_coord_z());
#endif
#ifdef FOUR_DIMENSIONS
  fprintf(fptr,"%f %f %f %f\n",steps[i].real_coord_x(),
    steps[i].real_coord_y(),steps[i].real_coord_z(),steps[i].real_coord_w());
#endif
} // end walk::plot()

void walk::plot_projected(FILE *fptr,double rho) 
// this prints out the x,y coordinates of the walk projected using rho
{ 
long i;
double x,y;
rpoint rpt;

if (npivot!=0) 
  {printf("ERROR: walk::print() called with npivot!=0 \n");
   exit(1L);
  } 

#ifdef THREE_DIMENSIONS
double a=cos(M_PI*rho/180);
double b=sin(M_PI*rho/180);
#endif

for (i=0;i<=nsteps;i++) 
 {
  // compute x,y
  rpt=steps[i]; 
  #ifdef TWO_DIMENSIONS
  x=rpt.coord_x();
  # endif
  #ifdef THREE_DIMENSIONS
  x=a*rpt.coord_x()+b*rpt.coord_z();
  # endif
  y=rpt.coord_y();
  fprintf(fptr,"%f %f\n",x,y);
 } // end loop on i
} // end walk::plot_projected()

walk& walk::operator=(walk& w)
{
long i;
if (w.get_npivot()!=0)
 {printf("ERROR in walk::=, npivot!=0 \n");
  exit(1);
 } 
(*this).initialize();
nsteps=w.get_nsteps();
niter=w.get_niter();
for (i=0;i<=nsteps;i++) steps[i]=w.steps[i];
return(*this);
} 

void walk::test(excluded_class* excluded) 
{
// this tests that the walk is a n.n. walk and that it does not hit the 
// excluded region
// If it fails the test we just crash (exit).
long i,iseg;
point stepsi;

// n.n. check
for (i=1;i<=nsteps;i++) 
  if (fabs((steps[i-1]-steps[i]).distance()-1.) 
   > 1.e-5) 
    {printf("BAD WALK  BAD WALK   BAD WALK   BAD WALK   BAD WALK \n");
     printf("BAD WALK  BAD WALK   BAD WALK   BAD WALK   BAD WALK \n");
     printf("not nearest neighbor\n");
     exit(1);
    } 

for (i=1;i<=nsteps;) 
 {
  // stepsi is w[i]   
  iseg=find_segment(i,npivot,ptime);
  stepsi.euclidean_op(steps+i,shift+iseg,igroup[iseg]); 
  long separation=stepsi.excluded_distance(excluded);
  if (separation==0) 
    {printf("BAD WALK  BAD WALK   BAD WALK   BAD WALK   BAD WALK \n");
     printf("BAD WALK  BAD WALK   BAD WALK   BAD WALK   BAD WALK \n");
     printf("hits excluded region\n");
     exit(1);
    } 
  i+=separation;
 } // end loop on i
} // end walk::test()

double walk::turn_frac()
{
// Finds fraction of times at which walk turns.
// For hexagonal lattice we look at next nearest neighbors
long i,count;
count=0;

// hexagonal lattice needs special treatment: n.n.n. are always same distance
#if ALGORITHM_NUMBER==5
for (i=1;i<nsteps-1;i++) 
  if (fabs((steps[i-1]-steps[i+2]).distance()-sqrt(7.)) > 1.e-5) count++;
#else
for (i=1;i<nsteps;i++) 
  if (fabs((steps[i-1]-steps[i+1]).distance()-2.) > 1.e-5) count++;
#endif

return(double(count)/double(nsteps-1));
} // end turn_frac

long walk::pivot_strictly_saw(long pivot_loc,long isym,excluded_class excluded)
// This is the version that changes i and j "simultaneously"
// If the pivot is accepted, this routine carries it out. Otherwise it is
// leaves the walk unchanged. 
// Routine returns count if the new walk  is self-avoiding, 
// returns -count if the pivot produces a self-intersection or 
// intersects the excluded region, 
// where count is the number of distance computations done.
// This return value is only used to study how long the routine takes. 
{
long i,j,ip,jp,iseg,jseg,igroup_jseg,igroup_iseg,imin,imax,jmin,jmax;
long separation,min_separation,sep_mod;
point origin,transi,transj,pp,stepsp,stepsi,stepsj,shift_jseg,shift_iseg;
long count,changei_flag;

// pivot is given by w[t] -> g (w[t]-w[pivot_loc])+w[pivot_loc]
// this is equivalent to w[t] -> g w[t] + trans 
// where trans= w[pivot_loc] - g w[pivot_loc]
// stepsp=w[pivot_loc]

origin.zero();
iseg=find_segment(pivot_loc,npivot,ptime);
stepsp.euclidean_op(steps+pivot_loc,shift+iseg,igroup[iseg]);
transi.euclidean_op(&stepsp,&origin,isym);
transi=stepsp-transi;
transj.euclidean_op(&stepsp,&origin,group_inverse[isym]);
transj=stepsp-transj;
count=0;

/////////////////////////////////////////////////////////////////////////////
// following checks if pivoted walk hits the excluded region
/////////////////////////////////////////////////////////////////////////////

for (i=pivot_loc+1;i<=nsteps;) 
 {
  // stepsi is w[i]   
  iseg=find_segment(i,npivot,ptime);
  stepsi.euclidean_op(steps+i,shift+iseg,igroup[iseg]); 
  // pp is w[i] after pivot
  pp.euclidean_op(&stepsi,&transi,isym);
  separation=pp.excluded_distance(&excluded);
  count++;
  if (separation==0) return(-count);
  i+=separation;
 } // end loop on i

/////////////////////////////////////////////////////////////////////////////
// 
// "Simultaneous" changing of i and j. 
// Given j<pivot_loc<i, we assume omega(jp)!=omega(ip) for j<jp<pivot_loc<ip<i
// We then either decrease j or increase i.
// 
/////////////////////////////////////////////////////////////////////////////

jmin=0; jmax=pivot_loc-1;
imin=pivot_loc+1; imax=nsteps;

j=pivot_loc-2; 
i=pivot_loc+1; 
// if pivot_loc==0, j=-2 which would allow jp=-1
if (j<jmin-1) j=jmin-1;

if (pivot_loc>nsteps) i=imax+1;

if (!NO_SAW) while (i<=imax || j>=jmin)
 {
  // changei_flag=1 means we will increase i, =0 means we will increase j
  // We change the index that is closer to pivot_loc (if allowed).
  if (i-pivot_loc>pivot_loc-j) changei_flag=0;
  else changei_flag=1;
  if (i>imax) changei_flag=0;
  if (j<jmin) changei_flag=1;
  
  if (changei_flag)
   {// increase i. Need lower bound on distance from pivoted omega[i] to 
    // {omega[jp]: j<jp<pivot_loc}
    // This lower bound will be min_separation
    // stepsi is w[i]   
    iseg=find_segment(i,npivot,ptime);
    stepsi.euclidean_op(steps+i,shift+iseg,igroup[iseg]); 
    // pp is w[i] after pivot
    pp.euclidean_op(&stepsi,&transi,isym);
    min_separation=nsteps;
    jseg=npivot;   // can change to jseg=iseg; ?
    shift_jseg=shift[jseg]-pp;
    igroup_jseg=igroup[jseg];
    for (jp=jmax;jp>j;) // note that j is decreased
     {
      if (ptime[jseg]>jp) 
       {
        #ifdef USE_FIND_SEGMENT
        jseg=find_segment(jp,npivot,ptime);  
        #else 
        while (ptime[jseg]>jp) jseg--;
        #endif
        shift_jseg=shift[jseg]-pp;
        igroup_jseg=igroup[jseg];
       }
      // stepsj is w[jp]   
      stepsj.euclidean_op(steps+jp,&shift_jseg,igroup_jseg);
      separation=stepsj.lattice_distance();
      count++;
      if (separation==0) return(-count);
      if (separation>=min_separation)
         {jp-=separation-min_separation+1;
         } 
      else
         {
          #ifdef USE_THIRD
          sep_mod=separation%3;
          min_separation=(2*separation+sep_mod)/3;
          jp-=1+(separation-sep_mod)/3;
          #else 
          sep_mod=separation%2;
          min_separation=(separation+sep_mod)/2;
          jp-=1+(separation-sep_mod)/2;
          #endif
         }
     } // end loop on jp
    i+=min_separation;
    if (i>imax+1) i=imax+1;
   } // end increase i 

  else
   {// decrease j. Need lower bound on distance from omega[j] to 
    // pivoted {omega[ip]: pivot_loc<ip<i}
    // Equivalently we can use a lower bound on distance from inverse
    // pivoted omega[j] and {omega[ip]: pivot_loc<ip<i}
    // This lower bound will be min_separation
    // stepsj is w[j]   
    jseg=find_segment(j,npivot,ptime);
    stepsj.euclidean_op(steps+j,shift+jseg,igroup[jseg]); 
    pp.euclidean_op(&stepsj,&transj,group_inverse[isym]);
    min_separation=nsteps;
    iseg=0;   
    shift_iseg=shift[iseg]-pp;
    igroup_iseg=igroup[iseg];
    for (ip=imin;ip<i;) // note that i is increased
     {
      if (ptime[iseg+1]<=ip) // check this
       {
        #ifdef USE_FIND_SEGMENT
        iseg=find_segment(ip,npivot,ptime);  
        #else 
        while (ptime[iseg+1]<=ip) iseg++; 
        #endif
        shift_iseg=shift[iseg]-pp;
        igroup_iseg=igroup[iseg];
       }
      stepsi.euclidean_op(steps+ip,&shift_iseg,igroup_iseg);
      separation=stepsi.lattice_distance();
      count++;
      if (separation==0) return(-count);
      if (separation>=min_separation)
         {ip+=separation-min_separation+1;
         } 
      else
         {
          #ifdef USE_THIRD
          sep_mod=separation%3;
          min_separation=(2*separation+sep_mod)/3;
          ip+=1+(separation-sep_mod)/3;
          #else 
          sep_mod=separation%2;
          min_separation=(separation+sep_mod)/2;
          ip+=1+(separation-sep_mod)/2;
          #endif
         }
     } // end loop on ip
    j-=min_separation;
    if (j<jmin-1) j=jmin-1;
   } // end decrease j 
 } // end while 

// If we reach this point the walk is self-avoiding. 
// The pivot operations were not done to the walk. So we must do them.
add_pivot(pivot_loc,isym,transi);
return(count);

} // end pivot_strictly_saw()

long walk::pivot_weakly_saw(long pivot_loc,long isym,double beta,double p,
  excluded_class excluded) 
// returns 1 if the pivot is accepted 
// returns 0 if the pivot is rejected and restores walk to its original state
// We accept the pivot if xrandom<exp(-beta*(new_energy-old_energy)) 
// where xrandom is a uniform random number in [0,1], 
// and reject the pivot otherwise
// The energy of an intersection between w[i] and w[i] is beta/|i-j|^p.
// Taking p=0 gives the usual WSAW. Setting p>0 gives the "forgetful" WSAW.
{
long i,j,ip,jp,iseg,jseg,igroup_jseg,igroup_iseg,imin,imax,jmin,jmax,hit_flag;
long separation,min_separation,sep_mod;
point origin,trans,transi,transj,pp,stepsp,stepsi,stepsj,shift_jseg,shift_iseg;
long count,changei_flag;
double xrandom,new_energy_cutoff,old_energy,new_energy;

// this code is very close to that in pivot_strictly_saw
// See that code for more comments
origin.zero();
iseg=find_segment(pivot_loc,npivot,ptime);
stepsp.euclidean_op(steps+pivot_loc,shift+iseg,igroup[iseg]);

trans.euclidean_op(&stepsp,&origin,isym);
trans=stepsp-trans;

transi.euclidean_op(&stepsp,&origin,isym);
transi=stepsp-transi;
transj.euclidean_op(&stepsp,&origin,group_inverse[isym]);
transj=stepsp-transj;

count=0;

/////////////////////////////////////////////////////////////////////////////
// following checks if pivoted walk hits the excluded region
/////////////////////////////////////////////////////////////////////////////

for (i=pivot_loc+1;i<=nsteps;) 
 {
  // stepsi is w[i]   
  iseg=find_segment(i,npivot,ptime);
  stepsi.euclidean_op(steps+i,shift+iseg,igroup[iseg]); 
  // pp is w[i] after pivot
  pp.euclidean_op(&stepsi,&transi,isym);
  separation=pp.excluded_distance(&excluded);
  count++;
  if (separation==0) return(-count);
  i+=separation;
 } // end loop on i

/////////////////////////////////////////////////////////////////////////////
// compute old energy
/////////////////////////////////////////////////////////////////////////////

old_energy=0.;
jmin=0; jmax=pivot_loc-1;
imin=pivot_loc+1; imax=nsteps;

j=pivot_loc-2; 
i=pivot_loc+1; 
// if pivot_loc==0, j=-2 which would allow jp=-1
if (j<jmin-1) j=jmin-1;

if (pivot_loc>nsteps) i=imax+1;

if (!NO_SAW) while (i<=imax || j>=jmin)
 {
  // changei_flag=1 means we will increase i, =0 means we will increase j
  // We change the index that is closer to pivot_loc (if allowed).
  if (i-pivot_loc>pivot_loc-j) changei_flag=0;
  else changei_flag=1;
  if (i>imax) changei_flag=0;
  if (j<jmin) changei_flag=1;
  
  if (changei_flag)
   {// increase i. Need lower bound on distance from pivoted omega[i] to 
    // {omega[jp]: j<jp<pivot_loc}
    // This lower bound will be min_separation
    // stepsi is w[i]   
    iseg=find_segment(i,npivot,ptime);
    stepsi.euclidean_op(steps+i,shift+iseg,igroup[iseg]); 
    min_separation=nsteps;
    jseg=npivot;   // can change to jseg=iseg; ?
    shift_jseg=shift[jseg]-stepsi;
    igroup_jseg=igroup[jseg];
    hit_flag=0; 
    for (jp=jmax;jp>j;) // note that j is decreased
     {
      if (ptime[jseg]>jp) 
       {   
        #ifdef USE_FIND_SEGMENT
        jseg=find_segment(jp,npivot,ptime);  
        #else 
        while (ptime[jseg]>jp) jseg--;
        #endif
        shift_jseg=shift[jseg]-stepsi;
        igroup_jseg=igroup[jseg];
       }
      // stepsj is w[jp]   
      stepsj.euclidean_op(steps+jp,&shift_jseg,igroup_jseg);
      separation=stepsj.lattice_distance();
      count++;
      if (separation==0) 
       {hit_flag=1;
        if (p>0.) old_energy+=pow(double(abs(jp-i)),-p); 
        else old_energy+=1.;
        min_separation=nsteps;
        jp--;
       }	
      else // separation>0
       {if (separation>=min_separation)
         {jp-=separation-min_separation+1;
         } 
        else
         {
          #ifdef USE_THIRD
          sep_mod=separation%3;
          min_separation=(2*separation+sep_mod)/3;
          jp-=1+(separation-sep_mod)/3;
          #else 
          sep_mod=separation%2;
          min_separation=(separation+sep_mod)/2;
          jp-=1+(separation-sep_mod)/2;
          #endif
         }
       } // end else (separation==0)
     } // end loop on jp
    if (hit_flag) i++;
    else i+=min_separation;
    if (i>imax+1) i=imax+1;
   } // end increase i 

  else
   {// decrease j. Need lower bound on distance from omega[j] to 
    // pivoted {omega[ip]: pivot_loc<ip<i}
    // Equivalently we can use a lower bound on distance from inverse
    // pivoted omega[j] and {omega[ip]: pivot_loc<ip<i}
    // This lower bound will be min_separation
    // stepsj is w[j]   
    jseg=find_segment(j,npivot,ptime);
    stepsj.euclidean_op(steps+j,shift+jseg,igroup[jseg]); 
    min_separation=nsteps;
    iseg=0;   
    shift_iseg=shift[iseg]-stepsj;
    igroup_iseg=igroup[iseg];
    hit_flag=0; 
    for (ip=imin;ip<i;) // note that i is increased
     {
      if (ptime[iseg+1]<=ip) // check this
       {
        #ifdef USE_FIND_SEGMENT
        iseg=find_segment(ip,npivot,ptime);  
        #else 
        while (ptime[iseg+1]<=ip) iseg++; 
        #endif
        shift_iseg=shift[iseg]-stepsj;
        igroup_iseg=igroup[iseg];
       }
      stepsi.euclidean_op(steps+ip,&shift_iseg,igroup_iseg);
      separation=stepsi.lattice_distance();
      count++;
      if (separation==0) 
       {hit_flag=1;
        if (p>0.) old_energy+=pow(double(abs(j-ip)),-p); 
        else old_energy+=1.;
        min_separation=nsteps;
        ip++;
       }
      else
       {
        if (separation>=min_separation)
         {ip+=separation-min_separation+1;
         } 
        else
         {
          #ifdef USE_THIRD
          sep_mod=separation%3;
          min_separation=(2*separation+sep_mod)/3;
          ip+=1+(separation-sep_mod)/3;
          #else 
          sep_mod=separation%2;
          min_separation=(separation+sep_mod)/2;
          ip+=1+(separation-sep_mod)/2;
          #endif
         }
       } // end else (separation==0)
     } // end loop on ip
    if (hit_flag) j--;
    else j-=min_separation;
    if (j<jmin-1) j=jmin-1;
   } // end decrease j 
 } // end while 

/////////////////////////////////////////////////////////////////////////////
// compute new energy
/////////////////////////////////////////////////////////////////////////////

// carry out the pivot and compute the new energy as we go. 
// Since the new energy can only increase, we can reject the pivot
// once xrandom>exp(-beta*(new_energy-old_energy))  i.e., 
// once old_energy-log(xrandom)/beta<new_energy
#if RANDOM_NUMBER_GENERATOR==1
xrandom=drand48();
#endif
#if RANDOM_NUMBER_GENERATOR==2
xrandom=sprng();
#endif
if (beta>1.e-10) new_energy_cutoff=old_energy-log(xrandom)/beta;
else new_energy_cutoff=1.e10;

new_energy=0.;
jmin=0; jmax=pivot_loc-1;
imin=pivot_loc+1; imax=nsteps;

j=pivot_loc-2; 
i=pivot_loc+1; 
// if pivot_loc==0, j=-2 which would allow jp=-1
if (j<jmin-1) j=jmin-1;

if (pivot_loc>nsteps) i=imax+1;

if (!NO_SAW) while (i<=imax || j>=jmin)
 {
  // changei_flag=1 means we will increase i, =0 means we will increase j
  // We change the index that is closer to pivot_loc (if allowed).
  if (i-pivot_loc>pivot_loc-j) changei_flag=0;
  else changei_flag=1;
  if (i>imax) changei_flag=0;
  if (j<jmin) changei_flag=1;
  
  if (changei_flag)
   {// increase i. Need lower bound on distance from pivoted omega[i] to 
    // {omega[jp]: j<jp<pivot_loc}
    // This lower bound will be min_separation
    // stepsi is w[i]   
    iseg=find_segment(i,npivot,ptime);
    stepsi.euclidean_op(steps+i,shift+iseg,igroup[iseg]); 
    // pp is w[i] after pivot
    pp.euclidean_op(&stepsi,&transi,isym);
    min_separation=nsteps;
    jseg=npivot;   // can change to jseg=iseg; ?
    shift_jseg=shift[jseg]-pp;
    igroup_jseg=igroup[jseg];
    hit_flag=0; 
    for (jp=jmax;jp>j;) // note that j is decreased
     {
      if (ptime[jseg]>jp) 
       {   
        #ifdef USE_FIND_SEGMENT
        jseg=find_segment(jp,npivot,ptime);  
        #else 
        while (ptime[jseg]>jp) jseg--;
        #endif
        shift_jseg=shift[jseg]-pp;
        igroup_jseg=igroup[jseg];
       }
      // stepsj is w[jp]   
      stepsj.euclidean_op(steps+jp,&shift_jseg,igroup_jseg);
      separation=stepsj.lattice_distance();
      count++;
      if (separation==0) 
       {hit_flag=1;
        if (p>0.) new_energy+=pow(double(abs(jp-i)),-p); 
        else new_energy+=1.;
        if (new_energy>new_energy_cutoff) return(-1);
        min_separation=nsteps;
        jp--;
       }	
      else // separation>0
       {if (separation>=min_separation)
         {jp-=separation-min_separation+1;
         } 
        else
         {
          #ifdef USE_THIRD
          sep_mod=separation%3;
          min_separation=(2*separation+sep_mod)/3;
          jp-=1+(separation-sep_mod)/3;
          #else 
          sep_mod=separation%2;
          min_separation=(separation+sep_mod)/2;
          jp-=1+(separation-sep_mod)/2;
          #endif
         }
       } // end else (separation==0)
     } // end loop on jp
    if (hit_flag) i++;
    else i+=min_separation;
    if (i>imax+1) i=imax+1;
   } // end increase i 

  else
   {// decrease j. Need lower bound on distance from omega[j] to 
    // pivoted {omega[ip]: pivot_loc<ip<i}
    // Equivalently we can use a lower bound on distance from inverse
    // pivoted omega[j] and {omega[ip]: pivot_loc<ip<i}
    // This lower bound will be min_separation
    // stepsj is w[j]   
    jseg=find_segment(j,npivot,ptime);
    stepsj.euclidean_op(steps+j,shift+jseg,igroup[jseg]); 
    pp.euclidean_op(&stepsj,&transj,group_inverse[isym]);
    min_separation=nsteps;
    iseg=0;   
    shift_iseg=shift[iseg]-pp;
    igroup_iseg=igroup[iseg];
    hit_flag=0; 
    for (ip=imin;ip<i;) // note that i is increased
     {
      if (ptime[iseg+1]<=ip) // check this
       {
        #ifdef USE_FIND_SEGMENT
        iseg=find_segment(ip,npivot,ptime);  
        #else 
        while (ptime[iseg+1]<=ip) iseg++; 
        #endif
        shift_iseg=shift[iseg]-pp;
        igroup_iseg=igroup[iseg];
       }
      stepsi.euclidean_op(steps+ip,&shift_iseg,igroup_iseg);
      separation=stepsi.lattice_distance();
      count++;
      if (separation==0) 
       {hit_flag=1;
        if (p>0.) new_energy+=pow(double(abs(j-ip)),-p); 
        else new_energy+=1.;
        if (new_energy>new_energy_cutoff) return(-1);
        min_separation=nsteps;
        ip++;
       }
      else
       {
        if (separation>=min_separation)
         {ip+=separation-min_separation+1;
         } 
        else
         {
          #ifdef USE_THIRD
          sep_mod=separation%3;
          min_separation=(2*separation+sep_mod)/3;
          ip+=1+(separation-sep_mod)/3;
          #else 
          sep_mod=separation%2;
          min_separation=(separation+sep_mod)/2;
          ip+=1+(separation-sep_mod)/2;
          #endif
         }
       } // end else (separation==0)
     } // end loop on ip
    if (hit_flag) j--;
    else j-=min_separation;
    if (j<jmin-1) j=jmin-1;
   } // end decrease j 
 } // end while 

// If we reach this point we accept the pivot. 
// The pivot operations were not done to the walk. So we must do them.
add_pivot(pivot_loc,isym,trans);

return(1);

} // end pivot_weakly_saw()

void walk::add_pivot(long pivot_loc,long isym,point trans) 
{
long iseg,ipivot;
point pp;

if (npivot>MAX_NPIVOT-1) 
 {printf("number of implicit pivots in walk exceeds MAX_NPIVOT \n");
  exit(1);
 } 

// first, we add the pivot time to the list of pivot times. 
iseg=find_segment(pivot_loc,npivot,ptime);
if (pivot_loc!=ptime[iseg])
 {ptime[npivot+2]=ptime[npivot+1]; 
  for (ipivot=npivot;ipivot>=iseg;ipivot--)
   {
    ptime[ipivot+1]=ptime[ipivot]; 
    igroup[ipivot+1]=igroup[ipivot]; 
    shift[ipivot+1]=shift[ipivot]; 
   } // end loop on ipivot
  ptime[iseg+1]=pivot_loc;
  npivot++;
  iseg++; // pivot will be applied to segments iseg to npivot
 }

// second, we update igroup and shift
for (ipivot=iseg;ipivot<=npivot;ipivot++)
 {
  pp=shift[ipivot];
  shift[ipivot].euclidean_op(&pp,&trans,isym); 
  igroup[ipivot]=group_product[isym][igroup[ipivot]]; 
 } // end loop on ipivot

} // add_pivot()


void walk::simplify() 
// carry out the pivot operations implicit in the walk, so npivot -> 0 
{
long ipivot,itime,igroup_ipivot;
point shift_ipivot,pp;

// even on the 0th segment there may be something to do 
for (ipivot=0;ipivot<npivot;ipivot++)
 {
  shift_ipivot=shift[ipivot];
  igroup_ipivot=igroup[ipivot];
  for (itime=ptime[ipivot];itime<ptime[ipivot+1];itime++)
   {pp=steps[itime];
    steps[itime].euclidean_op(&pp,&shift_ipivot,igroup_ipivot);
   }
 } // end loop on ipivot

// npivot segment is different
for (itime=ptime[npivot];itime<=nsteps;itime++)
  {pp=steps[itime];
   steps[itime].euclidean_op(&pp,shift+npivot,igroup[npivot]);
  }

initialize();
} // end walk::simplify()

////////////////////////////////////////////////////////////////
//             misc 
////////////////////////////////////////////////////////////////

inline long find_segment(long itime,long npivot,long* ptime) 
// Finds the segment number for a pivot time. 
// segment number iseg corresponds to times [ptime[iseg],ptime[iseg+1])
// Note that the left endpt, ptime[iseg], is included in segment iseg.
// Presently, routine is a slow search. This could be done better, but
// I am not sure it is worth it. 
{
long iseg,isegl,isegu;

if (itime>=ptime[npivot]) return(npivot);

isegl=0; isegu=npivot; iseg=0;

while (isegu>isegl+1)
 {
  iseg=(isegl+isegu)/2;
  if (itime<ptime[iseg]) isegu=iseg;
  else isegl=iseg;
 }
return(isegl);
} // find_segment()



