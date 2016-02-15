
!!!-----------------------------------------------------------
!! This program is for simulation of chromosome movements(BD). 
!! Author: wenwen@pks.mpg.de
!! update: Fri Mar  7 10:58:06 CET 2014
!! Scaling:
!! length l --- a   
!! time   t --- xi*a**2/(kb*tem) 
!! force  f --- kb*tem/a
!!    a : rod length
!!   xi : friction coefficent
!!   kb : Boltzmann constant
!!   tem: absolute temperature
!!!-----------------------------------------------------------

MODULE Vars
implicit none
integer, parameter :: ps=kind(1.0d0)
integer, parameter :: bead_num=10,dimen=3,max_step=2E+5
real(kind=ps), parameter :: dt=1.0d-4,PI=3.14159265359d0
real(kind=ps), parameter :: rod_len=1.0d-7
real(kind=ps) :: r(dimen,bead_num),rs(dimen,bead_num)
real(kind=ps) :: u(dimen,bead_num),tension(bead_num)
real(kind=ps) :: f_uc(dimen,bead_num),fc(dimen,bead_num)
real :: fb(dimen,bead_num)
END MODULE Vars

MODULE Random
! This is a module generates random numbers
CONTAINS

FUNCTION Ran(idum)
implicit none
integer, parameter :: k4b=selected_int_kind(9) 
integer(k4b), intent(inout) :: idum
real :: ran
integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836 
real, save :: am
integer(k4b), save :: ix=-1,iy=-1,k
if (idum <= 0 .or. iy < 0) then
  am=nearest(1.0,-1.0)/im 
  iy=ior(ieor(888889999,abs(idum)),1)
  ix=ieor(777755555,abs(idum)) 
  idum=abs(idum)+1
end if
ix=ieor(ix,ishft(ix,13)) 
ix=ieor(ix,ishft(ix,-17)) 
ix=ieor(ix,ishft(ix,5))
k=iy/iq
iy=ia*(iy-k*iq)-ir*k
if (iy < 0) iy=iy+im 
ran=am*ior(iand(im,ieor(ix,iy)),1) 
END FUNCTION Ran

FUNCTION RanNorm(idum)
implicit none
integer, parameter :: k4b=selected_int_kind(9) 
integer(k4b), intent(inout) :: idum
real :: rannorm
real :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472, &
r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
do
  u = ran(idum)
  v = ran(idum)
  v = 1.7156 * (v - 0.5)
  x = u - s
  y = abs(v) - t
  q = x**2 + y*(a*y - b*x)
  if (q < r1) exit
  if (q > r2) cycle
  if (v**2 < -4.0*log(u)*u**2) exit
end do
rannorm = v/u
RETURN
END FUNCTION RanNorm

SUBROUTINE RanMvnorm(idum, n, f, x, first)
implicit none
integer, parameter :: k4b=selected_int_kind(9) 
integer(k4b), intent(inout) :: idum
integer, intent(in) :: n
real, intent(inout) :: f(:) !f(n*(n+1)/2)
real, intent(out) :: x(:)
logical, intent(in) :: first
integer :: j, i, m
real :: d(n*(n+1)/2) 
real :: y, v 
integer, save :: n2


if (first) then                        
  d = 0.0
  do i = 1,n
    d(i*(i-1)/2+i)=1.0
  end do
  n2 = 2*n
  f(1) = sqrt(d(1))
  y = 1.0/f(1)
  do j = 2,n
    f(j) = d(1+j*(j-1)/2) * y
  end do
  do i = 2,n
    v = d(i*(i-1)/2+i)
    do m = 1,i-1
      v = v - f((m-1)*(n2-m)/2+i)**2
    end do
    v = sqrt(v)
    y = 1.0/v
    f((i-1)*(n2-i)/2+i) = v
    do j = i+1,n
      v = d(j*(j-1)/2+i)
      do m = 1,i-1
        v = v - f((m-1)*(n2-m)/2+i)*f((m-1)*(n2-m)/2 + j)
      end do 
      f((i-1)*(n2-i)/2 + j) = v*y
    end do 
  end do 
end if

x = 0.0
do j = 1,n
  y = RanNorm(idum)
  do i = j,n
    x(i) = x(i) + f((j-1)*(n2-j)/2 + i) * y
  end do 
end do 
RETURN
END SUBROUTINE RanMvnorm

END MODULE Random

MODULE Constraint
use Vars
CONTAINS

SUBROUTINE Eqs(a,b)
implicit none
real(kind=ps), intent(out) :: a(:,:),b(:) 
integer :: i,j,k,n

n = size(b)

! specify the metric matrix and the rhs of eqs.
do i=1,n
  j=i+1
  k=i-1
  if (j>n) j=j-n
  if (k<1) k=n
  a(i,j) = -dot_product(u(:,i),u(:,j))
  a(i,k) = -dot_product(u(:,i),u(:,k))
  a(i,i) = 2*dot_product(u(:,i),u(:,i))
  b(i) = dot_product(u(:,i),f_uc(:,j)-f_uc(:,i))
end do

RETURN
END SUBROUTINE Eqs

SUBROUTINE Lubksb(a,indx,b)
implicit none
integer, intent(in) :: indx(:) 
real(kind=ps), intent(in) :: a(:,:) 
real(kind=ps), intent(inout) :: b(:)
integer :: i,n,ii,ll
real(kind=ps) :: summ 
n=size(a,1)
ii=0
do i=1,n
  ll=indx(i) 
  summ=b(ll) 
  b(ll)=b(i)
  if (ii /= 0) then
    summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1)) 
  else if (summ /= 0.0) then
    ii=i 
  end if
  b(i)=summ 
end do
do i=n,1,-1
  b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
end do
END SUBROUTINE Lubksb

SUBROUTINE Ludcmp(a,indx,d)
implicit none
real(kind=ps), intent(inout) :: a(:,:)
integer, intent(out) :: indx(:)
real(kind=ps), intent(out) :: d
real(kind=ps) :: vv(size(a,1)),dum(size(a,1))
real(kind=ps), parameter :: tiny=1.0d-20 
integer :: j,n,imax
integer,dimension(1) :: imaxloc
n=size(a,1)
d=1.0  
vv=maxval(abs(a),dim=2) 
if (any(vv == 0.0)) write(*,*)'Singular matrix in ludcmp'
vv=1.0d0/vv 
do j=1,n
  imaxloc = maxloc(vv(j:n)*abs(a(j:n,j)))
  imax=(j-1)+imaxloc(1) 
  if (j /= imax) then
    dum=a(imax,:)
    a(imax,:)=a(j,:)
    a(j,:)=dum
    d=-d
    vv(imax)=vv(j)
  end if
  indx(j)=imax
  if (a(j,j) == 0.0) a(j,j)=TINY
  a(j+1:n,j)=a(j+1:n,j)/a(j,j) 
  a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
end do
END SUBROUTINE Ludcmp

FUNCTION Outerprod(a,b)
implicit none
real(kind=ps), dimension(:), intent(in) :: a,b
real(kind=ps), dimension(size(a),size(b)) :: outerprod 
outerprod = spread(a,dim=2,ncopies=size(b)) * &
spread(b,dim=1,ncopies=size(a))
END FUNCTION Outerprod

SUBROUTINE Newt(x,check)
implicit none
real(kind=ps),intent(inout) :: x(:)
logical, intent(out) :: check
real(kind=ps) :: a(size(x),size(x)),b(size(x))
integer :: indx(size(x)) 
real(kind=ps) :: d

a = 0.0d0
b = 0.0d0
call eqs(a,b)
call ludcmp(a,indx,d)
call lubksb(a,indx,b)
x = b
check = .true.

RETURN 
END SUBROUTINE Newt

END MODULE Constraint

MODULE Force
use Vars
use Constraint
CONTAINS

SUBROUTINE CalFc
implicit none
integer :: i,j
logical :: check

do i=1,bead_num
  j=i+1
  if (j>bead_num) j=j-bead_num
  u(:,i)=r(:,j)-r(:,i)
  u(:,i)=u(:,i)/sqrt(dot_product(u(:,i),u(:,i)))
end do
call Newt(tension,check)
fc(:,1) = tension(1)*u(:,1)-tension(bead_num)*u(:,bead_num)
do i=2,bead_num
  fc(:,i) = tension(i)*u(:,i)-tension(i-1)*u(:,i-1)
end do

RETURN
END SUBROuTINE CalFc 

SUBROUTINE LennardJones(f_lj)
implicit none
integer :: i,j
real(kind=ps),intent(inout) :: f_lj(dimen,bead_num)
real(kind=ps) :: eps,a,b,r_ij(dimen),r6

a = 0.6d0
eps = 6.0d0
b = a**6
do j=1,bead_num
  f_lj(:,j)=0.0d0
  do i=1,bead_num
    if(i/=j) then
    r_ij=r(:,i)-r(:,j)
    r6=dot_product(r_ij,r_ij)**3
    if(r6<=2*b) then
      f_lj(:,j)=f_lj(:,j)+r_ij*4*eps*(12*b*b/r6*r6-6*b/r6)/dot_product(r_ij,r_ij)
    end if
  end if
  end do
end do

RETURN
END SUBROUTINE LennardJones

END MODULE Force


PROGRAM Main
use Vars
use Random
use Force
implicit none
real(kind=ps) :: t
integer(kind=4) :: idum,task_id
integer :: i,j,k,step
real :: random_num(dimen*bead_num),temp(size(random_num)*(size(random_num)+1)/2),theta,phi
real(kind=ps) :: rend(dimen)
character(len=30) :: arg

!---------------------------------------------------
! Passing argument to get task_ID when calculating in cluster
!call get_command_argument(1,arg)
!read(arg,*)task_id
!idum = task_id
!open(10,file='beadrod'//trim(adjustl(arg))//'.dat')

open(10,file='beadrod.dat')
open(11,file='rend.dat')
idum=3

!---------------------------------------------------
! Initialization
call RanMvnorm(idum, dimen*bead_num, temp, random_num, .true.)

do i=1,bead_num
  r(1,i)=cos(PI*2*(i-1)/bead_num)/(2*sin(PI/bead_num))
  r(2,i)=sin(PI*2*(i-1)/bead_num)/(2*sin(PI/bead_num))
  r(3,i)=1.0d0*(j-1)
end do
fb=0.0d0
fc=0.0d0

! Output initial state
!do i=1,bead_num
!do j=1,dimen
!write(10,*)r(j,i)
!end do
!end do

!---------------------------------------------------
! Evolution of the dynamical system
t=0.0d0
do step=1,max_step
  rs = r
  ! Step I: calculate known forces, i.e. unconstraint force, at time t
  f_uc = 0.0d0 

  ! 1:Calculate random force which is using during the whole time step
  call RanMvnorm(idum, dimen*bead_num, temp, random_num, .false.)
  fb = sqrt(2/dt)*reshape(random_num,(/dimen,bead_num/))

  ! 2:Calculate Lennard-Jones collision force
  call LennardJones(f_uc)

  ! 3:External force would be added after some equilibrium time
  if (t>10.0d0) then
    f_uc(1,1)=f_uc(1,1)+5.0d+2*sin(10.0d0*t)
  end if
  ! Sum them to f_uc
  f_uc = f_uc + fb
  
  ! Step II: calculate constraint force at time t
  call calFc

  ! Step III: calculate mid-step r(t+1/2dt)
  r = r+(f_uc+fc)*dt/2

  ! Step IV: calculate constraint force at time t+1/2dt
  f_uc = 0.0d0
  call LennardJones(f_uc)
  if (t>10.0d0) then
    f_uc(1,1)=f_uc(1,1)+5.0d+2*sin(10.0d0*(t+dt/2))
  end if
  f_uc = f_uc + fb
  call CalFc

  ! Step V: calculate final r(t+dt)
  r = rs+(f_uc+fc)*dt

  ! Step VI: reset r(t+dt) if the rod length change too much
  do i=1,bead_num
    j = i+1
    if (j>bead_num) j=j-bead_num
    u(:,i)=r(:,j)-r(:,i)
    u(:,i)=u(:,i)/sqrt(dot_product(u(:,i),u(:,i)))
  end do
  rend=sum(u,dim=2)/bead_num
  do i=1,bead_num-1
    r(:,i+1)=r(:,i)+u(:,i)-rend
  end do

  ! Output bead position
  if (t>10.0d0) then
    if(mod(step,int(1E2))==0) then
      do i=1,bead_num
        do j=1,dimen
          write(10,*)r(j,i)
        end do
      end do
      !write(*,*)t,sum(tension)
    end if
  end if

  t=t+dt
end do

close(10)
STOP
END PROGRAM Main

