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
integer, parameter :: chain_num=1,bead_num=10,rod_num=10
integer, parameter :: dimen=3,max_step=1E5
real(kind=ps), parameter :: dt=1.0d-4,PI=3.14159265359d0
integer :: link(2,rod_num),g(rod_num,rod_num)
real(kind=ps) :: r(dimen,bead_num),rs(dimen,bead_num)
real(kind=ps) :: tension(rod_num),u(dimen,rod_num),b(dimen,rod_num)
real(kind=ps) :: fa(dimen,bead_num),fc(dimen,bead_num)
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

! specify the link matrix
!k=0
!do k = 1,rod_num/2-1
	!link(1,k) = k
	!link(2,k) = k+1
!end do
!link(1,rod_num/2) = rod_num/2
!link(2,rod_num/2) = 1
!link(1,rod_num/2+1) = 1
!link(2,rod_num/2+1) = rod_num/2+1
!do k = rod_num/2+2,bead_num
	!link(1,k) = k-1
	!link(2,k) = k
!end do
!link(1,rod_num) = bead_num
!link(2,rod_num) = 1

do i = 1, bead_num-1
  link(1,i) = i
  link(2,i) = i+1
end do
link(1,bead_num) = bead_num
link(2,bead_num) = 1

! Output the network topology
do i=1,rod_num
  write(99,*)link(1,i),link(2,i)
end do

!specify metric matrix corresponding to Adj
g = 0
do i=1,rod_num
  g(i,i) = -2
  do j=1,rod_num 
    if(j==i) cycle
    i1=link(1,i)
    i2=link(2,i)
    j1=link(1,j)
    j2=link(2,j)
    if ((i1-j2)*(i2-j1)==0) then
      g(i,j) = 1
    else if ((i1-j1)*(i2-j2)==0) then
      g(i,j) = -1
    end if
  end do
end do

close(99)
RETURN
END SUBROUTINE Network

SUBROUTINE Picard(x) 
implicit none
integer, parameter :: MAXITS=1E+3
real(kind=ps), intent(inout) :: x(:) 
integer :: i,j,k,step,n
integer :: indx(size(x))
real(kind=ps) :: d,xold(size(x)),p(size(x))   
real(kind=ps) :: a(size(x),size(x)),temp(dimen)
real(kind=ps) :: eps=1.0d-8,err(size(x))
integer :: info,ipiv(size(x))

n = size(x)

a = 0.0d0
do j=1,n
  do i=1,n
    a(i,j)=g(i,j)*dot_product(b(:,i),u(:,j))
  end do
end do

!call Ludcmp(a,indx,d)
call DGETRF(n, n, a, n, ipiv, info)


do step=1,MAXITS
  xold = x
  do i=1,n
    temp = 0.0d0
    do j=1,n
      temp = temp + g(i,j)*x(j)*u(:,j)
    end do
    p(i)=(1.0d0-dot_product(b(:,i),b(:,i)))/(2*dt) &
         -dt*dot_product(temp,temp)/2
  end do
  !call Lubksb(a,indx,p)
  call DGETRS('N',n, 1, a, n, ipiv, p, n, info)
  x = p
  if(maxval(abs(x-xold))<eps)   RETURN
end do
write(*,*)'MAXITS exceeded in Picard'
RETURN
END SUBROUTINE Picard


END MODULE Constraint

MODULE Force
use Vars
use Constraint
CONTAINS

SUBROUTINE CalFc
implicit none
integer :: i,j,k

! Calculate u and b
do k=1,rod_num
	i=link(1,k)
	j=link(2,k)
	b(:,k) = rs(:,j)-rs(:,i)
end do

tension=0.0d0
call Picard(tension)

fc = 0.0d0 
do i=1,bead_num 
  do j=1,rod_num
    if(link(1,j)==i) then
      fc(:,i)=fc(:,i)+tension(j)*u(:,j)
    end if
    if(link(2,j)==i) then
      fc(:,i)=fc(:,i)-tension(j)*u(:,j)
    end if
  end do
end do

RETURN
END SUBROuTINE CalFc 

SUBROUTINE LennardJones(f_lj)
implicit none
integer :: i,j
real(kind=ps),intent(inout) :: f_lj(dimen,bead_num)
real(kind=ps) :: eps,a,b,r_ij(dimen),r6

a = 0.75d0
eps = 1.0d0
b = a**6
f_lj=0.0d0
do j=1,bead_num
  do i=j+1,bead_num
    r_ij=r(:,i)-r(:,j)
    r6=dot_product(r_ij,r_ij)**3
    if(r6<=2*b) then
      f_lj(:,i)=f_lj(:,i)+48.0d0*eps*((b/r6)*(b/r6)-0.5*(b/r6))*r_ij/dot_product(r_ij,r_ij)
      f_lj(:,j)=f_lj(:,j)-48.0d0*eps*((b/r6)*(b/r6)-0.5*(b/r6))*r_ij/dot_product(r_ij,r_ij)
    end if
  end do
end do

RETURN
END SUBROUTINE LennardJones

SUBROUTINE PseudoForce(f_pseudo)
implicit none
integer :: i,j,i1,j1,i2,j2,k
real(kind=ps),intent(inout) :: f_pseudo(dimen,bead_num)
real(kind=ps) :: metric(rod_num,rod_num),uij,pgr(dimen)
real(kind=ps) :: metric_inv(rod_num,rod_num)
real(kind=ps) :: work(size(metric,1))
integer :: ipiv(size(metric,1)),info,n
logical :: change

do j = 1,rod_num
  do i = 1,rod_num
    metric(i,j) = -g(i,j)*dot_product(u(:,i),u(:,j))
  end do
end do

metric_inv = metric
n = rod_num
call DGETRF(n, n, metric_inv, n, ipiv, info)
call DGETRI(n, metric_inv, n, ipiv, work, n, info)

do k = 1,bead_num
	f_pseudo(:,k) = 0.0d0
	do i=1,rod_num
		do j=i+1,rod_num
			if (abs(g(j,i))==1) then 
				i1 = link(1,i)
				i2 = link(2,i)
				j1 = link(1,j)
				j2 = link(2,j)
				if ((k-i1)*(k-i2)*(k-j1)*(k-j2)==0) then
					uij = dot_product(u(:,i),u(:,j))
					pgr = (merge(1,0,k==i2)-merge(1,0,k==i1)) &
				     	 *(u(:,j)-uij*u(:,i)) &
					     +(merge(1,0,k==j2)-merge(1,0,k==j1)) &
					     *(u(:,i)-uij*u(:,j))
					f_pseudo(:,k) = f_pseudo(:,k) + g(j,i)*metric_inv(j,i)*pgr
				end if
			end if
		end do
	end do
end do

! For one single ring
!do k = 1,bead_num
  !i = k
  !j = k+1
  !if (j>rod_num) j = j-rod_num
  !uij = dot_product(u(:,i),u(:,j))
  !f_pseudo(:,k)=metric_inv(j,i)*(uij*u(:,i)-u(:,j))
  !i = k-1
  !if (i<1) i = rod_num+i
  !j = k
  !uij = dot_product(u(:,i),u(:,j))
  !f_pseudo(:,k)=f_pseudo(:,k)+metric_inv(j,i)*(1+uij)*(u(:,j)-u(:,i))
  !i = k-2
  !if (i<1) i = rod_num+i
  !j = k-1
  !if (j<1) j = rod_num+j
  !uij = dot_product(u(:,i),u(:,j))
  !f_pseudo(:,k)=f_pseudo(:,k)+metric_inv(j,i)*(u(:,i)-uij*u(:,j))
!end do

RETURN
END SUBROUTINE PseudoForce

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
 

!inquire(file="initial.in", exist=have_input)
!if(have_input) then
  !open(8,file='initial.in')
  !do i=1,dimen*bead_num
    !read(8,*)random_num(i)
  !end do
  !r=reshape(random_num,(/dimen,bead_num/))
! else
 ! r(1,2)=r(1,1)+cos(PI/6)
  !r(2,2)=r(2,1)+sin(PI/6)
  !r(1,rod_num/2)=r(1,1)+cos(-PI/6)
  !r(2,rod_num/2)=r(2,1)+sin(-PI/6)
  !do i=1,(rod_num/2-1)/2-1
    !r(1,i+2)=r(1,2)+i
    !r(2,i+2)=r(2,2)
    !r(1,rod_num/2-i)=r(1,rod_num/2)+i
    !r(2,rod_num/2-i)=r(2,rod_num/2)
  !end do
  !if (mod(rod_num/2,2)==0) then
    !r(1,rod_num/4+1)=r(1,2)+(rod_num/2-1)/2-1+sqrt(3.0d0)/2
  !end if
	!r(1,rod_num/2+1:bead_num) = -r(1,2:rod_num/2)
 ! r(2,rod_num/2+1:bead_num) = r(2,2:rod_num/2)
!end if

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

! 2). Passing arguments if have any
!call get_command_argument(1,arg)
!read(arg,*)task_id

! 3). Initializing other variables or vectors
rs = r
fa = 0.0d0
fb = 0.0d0
fc = 0.0d0
r0 = 0.0d0
t0 = 0.0d3
v0 = 1.0d0
idum=3
!idum = task_id
!if(v0==0.0d0) then
  !T_eff = 0
!else
  !T_eff = int(1.0d0/v0)
!end if
!lj_switch = .true.
!pf_switch = .true.

! 4). Open a file for data output 
write(filename,'(4(A,I0),A)')'r_',dimen,'D_N',bead_num,'_T',T_eff,'_',task_id,'.dat'
open(10,file=filename)

! 5). Initializing random generator and rod network
!call RanMvnorm(idum, dimen*bead_num, temp, random_num, .true.)
call Network

! 6). Output initial settings
!write(*,*)'Ring Configuration'
!write(*,*)'bead_num = ',bead_num
!write(*,*)'rod_num = ',rod_num
!write(*,*)'dimen = ',dimen
!write(*,*)'dt = ',dt
!write(*,*)'max_step = ',max_step
!write(*,*)'T_eff = ',T_eff
!write(*,*)'Recording start time t0 = ',t0
!write(*,*)'Lennard Jones potential ON:',lj_switch
!write(*,*)'Pseudo Force ON:',pf_switch


! 7). Output initial state
do i=1,bead_num
	do j=1,dimen
		write(10,*)r(j,i)
	end do
end do
!write(*,*)'Output initial condition done!'
open(11,file='temp.dat')
!---------------------------------------------------
! Step II: Evolution of the dynamical system
t=0.0d0
do step=1,max_step
! Calculate rod unit vector u 
	do k=1,rod_num
		i=link(1,k)
		j=link(2,k)
		u(:,k) = r(:,j)-r(:,i)
		u(:,k) = u(:,k)/sqrt(dot_product(u(:,k),u(:,k)))
	end do

	fa = 0.0d0
! Generate random force
	!call RanMvnorm(idum, dimen*bead_num, temp, random_num, .false.)
	do i = 1,dimen*bead_num
		random_num(i) = RanNorm(idum)
	end do
	fb=sqrt(2.0d0/dt)*reshape(random_num,(/dimen,bead_num/))
	fa = fa + fb

! Calculate Lennard-Jones collision force
	call LennardJones(fb)
	fa = fa + fb

! Calculate pseudo force
	call PseudoForce(fb)
  do i = 1, bead_num
    write(11,*)fb(1,i),fb(2,i),fb(3,i)
  end do
	fa = fa + fb
	
! Add external force 
	fa(:,1) = fa(:,1)-2.0d3*(r(:,1)-r0)
	fa(1,2:bead_num) = fa(1,2:bead_num)+1.0d0*v0

! Predictive step
	rs = r+fa*dt

! Calculate constraint force
	call CalFc

! Corrective step
	r  = rs+fc*dt

! Output every second
	if (t>t0) then
		if(mod(step,int(1E3))==0) then
			do i=1,bead_num
				do j=1,dimen
					write(10,*)r(j,i)
				end do
			end do
			!if(mod(step,max_step/100)==0) then
				!write(*,*)step/(max_step/100),'% Done'
			!end if
		end if
	end if

	t=t+dt
end do


close(10)
close(8)
call cpu_time(end_time)
write(*,*)'Timming:',(end_time-start_time),'seconds'
STOP
END PROGRAM Main

