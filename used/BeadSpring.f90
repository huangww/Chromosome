!!!-----------------------------------------------------------
!! This program is for simulation of chromosome movements(BD). 
!! Author: wenwen@pks.mpg.de
!! update: Tue Aug  5 11:27:52 CEST 2014
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
integer,parameter :: chain_num=1,bead_num=50,spring_num=50
integer, parameter :: dimen=3,max_step=1E8
real(kind=ps), parameter :: dt=1.0d-4,PI=3.14159265359d0
integer :: link(2,spring_num)
real(kind=ps) :: r(dimen,bead_num)
real(kind=ps) :: q(dimen,spring_num),qs(dimen,spring_num)
real(kind=ps) :: fa(dimen,bead_num)
real(kind=ps) :: fb(dimen,bead_num)
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

SUBROUTINE Network 
implicit none
integer :: i,j,k,i0,chain,pair
integer :: i1,j1,i2,j2

open(99,file='network.dat')

!! One single ring
do i=1,spring_num
  i0=i+1
  if(i0>bead_num) i0=1
  link(1,i)=i
  link(2,i)=i0
end do

! Output the network topology
do i=1,spring_num
  write(99,*)link(1,i),link(2,i)
end do

close(99)
RETURN
END SUBROUTINE Network

END MODULE Constraint

MODULE Force
use Vars
use Constraint
CONTAINS

FUNCTION F_Hookean(q_spring)
implicit none
real(kind=ps) :: Hs
real(kind=ps) :: q_spring(dimen),F_Hookean(dimen)
Hs = 1.0d0
F_Hookean = Hs*q_spring
RETURN 
END FUNCTION F_Hookean

FUNCTION F_FENE(q_spring)
implicit none
real(kind=ps) :: Hs,q0
real(kind=ps) :: q_spring(dimen),F_FENE(dimen)
Hs = 1.0d0
q0 = 1.0d0
F_FENE = Hs*q_spring/(1.0d0-dot_product(q_spring,q_spring)/(q0*q0))
RETURN 
END FUNCTION F_FENE

SUBROUTINE Spring(f_sp)
implicit none
integer :: i,j,k
real(kind=ps),intent(inout) :: f_sp(dimen,bead_num)
real(kind=ps) :: eps

eps = 1.0d2
f_sp = 0.0d0
q = 0.0d0
qs = 0.0d0
do k=1,spring_num
  q(:,k) = r(:,link(2,k))-r(:,link(1,k))
  qs(:,k) = q(:,k)/sqrt(dot_product(q(:,k),q(:,k)))
  f_sp(:,link(1,k)) = f_sp(:,link(1,k))+eps*(q(:,k)-qs(:,k))
  f_sp(:,link(2,k)) = f_sp(:,link(2,k))-eps*(q(:,k)-qs(:,k))
end do

RETURN
END SUBROUTINE Spring

SUBROUTINE LennardJones(f_lj)
implicit none
integer :: i,j
real(kind=ps),intent(inout) :: f_lj(dimen,bead_num)
real(kind=ps) :: eps,a,b,r_ij(dimen),r6

a = 0.75d0
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

SUBROUTINE Boundary(f_bo)
implicit none
integer :: i,j
real(kind=ps),intent(inout) :: f_bo(dimen,bead_num)
real(kind=ps) :: eps,a,b,r_xyz(dimen),r6,r_cy,len_cy

r_cy = 10.0d0
len_cy = 100.0d0
eps = 1.0d3
do j=1,bead_num
  f_bo(:,j)=0.0d0
  r_xyz = (/r(1,j),r(2,j),0.0d0/)
  if(r(1,j)*r(1,j)+r(2,j)*r(2,j)>=r_cy*r_cy) then
    f_bo(:,j)=f_bo(:,j)-eps*r_xyz/sqrt(r(1,j)*r(1,j)+r(2,j)*r(2,j))
  end if
  r_xyz = (/0.0d0,0.0d0,r(3,j)/)
  if(abs(r(3,j))>=len_cy) then
    f_bo(:,j)=f_bo(:,j)-eps*r_xyz/sqrt(r(3,j)*r(3,j))
  end if
end do

RETURN
END SUBROUTINE Boundary

END MODULE Force


PROGRAM Main
use Vars
use Random
use Constraint
use Force
implicit none
integer :: idum,task_id,T_eff
integer :: i,j,k,step,i0,j0,chain
real :: temp(dimen*bead_num*(dimen*bead_num+1)/2)
real :: random_num(dimen*bead_num)
real :: start_time,end_time
real(kind=ps) :: rend(dimen),r0(dimen),v0,t,t0
character(len=60) :: arg,filename
logical :: have_input,lj_switch,pf_switch


!---------------------------------------------------
! Step I: initialization
call cpu_time(start_time)
! 1). Initial configuration
r = 0.0d0
inquire(file="initial.in", exist=have_input)
if(have_input) then
  open(8,file='initial.in')
  do i=1,dimen*bead_num
    read(8,*)random_num(i)
  end do
  r=reshape(random_num,(/dimen,bead_num/))
else
  r(1,2)=r(1,1)+cos(PI/6)
  r(2,2)=r(2,1)+sin(PI/6)
  r(1,bead_num)=r(1,1)+cos(-PI/6)
  r(2,bead_num)=r(2,1)+sin(-PI/6)
  do i=1,(bead_num-1)/2-1
    r(1,i+2)=r(1,2)+i
    r(2,i+2)=r(2,2)
    r(1,bead_num-i)=r(1,bead_num)+i
    r(2,bead_num-i)=r(2,bead_num)
  end do
  if (mod(bead_num,2)==0) then
    r(1,bead_num/2+1)=r(1,2)+(bead_num-1)/2-1+sqrt(3.0d0)/2
  end if
end if

! 2). Passing arguments if have any
call get_command_argument(1,arg)
read(arg,*)task_id

! 3). Initializing other variables or vectors
fa = 0.0d0
fb = 0.0d0
r0 = 0.0d0
t0 = 0.0d3
v0 = 0.02d0
!idum=3
idum = task_id
if(v0==0.0d0) then
  T_eff = 0
else
  T_eff = int(1.0d0/v0)
end if
lj_switch = .true.
pf_switch = .false.

! 4). Open a file for data output 
write(filename,'(4(A,I0),A)')'r_',dimen,'D_N',bead_num,'_T',T_eff,'_',task_id,'.dat'
open(10,file=filename)

! 5). Initializing random generator and rod network
call Network
!call RanMvnorm(idum, dimen*bead_num, temp, random_num, .true.)

! 6). Output initial settings
write(*,*)'Ring Configuration'
write(*,*)'bead_num = ',bead_num
write(*,*)'spring_num = ',spring_num
write(*,*)'dimen = ',dimen
write(*,*)'dt = ',dt
write(*,*)'max_step = ',max_step
write(*,*)'T_eff = ',T_eff
write(*,*)'Recording start time t0 = ',t0
write(*,*)'Lennard Jones potential ON:',lj_switch
write(*,*)'Pseudo Force ON:',pf_switch


! 7). Output initial state
do i=1,bead_num
	do j=1,dimen
		write(10,*)r(j,i)
	end do
end do
write(*,*)'Output initial condition done!'
!---------------------------------------------------
! Step II: Evolution of the dynamical system
t=0.0d0
do step=1,max_step
  fa = 0.0d0
! Generate random force
  do i=1,dimen*bead_num
    random_num(i) = RanNorm(idum)
  end do
	fb=sqrt(2.0d0/dt)*reshape(random_num,(/dimen,bead_num/))
  fa = fa + fb

! Calculate Lennard-Jones collision force
  !call LennardJones(fb)
  !fa = fa + fb

! Calculate spring potential force
  call Spring(fb)
  fa = fa + fb
  
! Add external force 
  fa(:,1) = fa(:,1)-2.0d3*(r(:,1)-r0)
  fa(1,2:bead_num) = fa(1,2:bead_num)+1.0d0*v0

! Calculate next step position vector
  r = r+fa*dt

! Output every second
	if (t>t0) then
    if(mod(step,int(1E4))==0) then
			do i=1,bead_num
				do j=1,dimen
					write(10,*)r(j,i)
				end do
			end do
      if(mod(step,max_step/100)==0) then
        write(*,*)step/(max_step/100),'% Done'
      end if
    end if
	end if

  t=t+dt
end do


close(10)
close(8)
call cpu_time(end_time)
write(*,*)'Timming:',(end_time-start_time)/3600,'hours'
STOP
END PROGRAM Main

