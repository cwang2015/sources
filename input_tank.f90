subroutine input

!use param
use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: property
double precision soil_surface,water_surface
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    

! Set nozzle and tank geometry parameters

!call tank%set(xl=3.34d0,yl=1.86d0,m=334,n=186)
!call tank%set(xl=3.34d0,yl=1.86d0,m=167,n=93)
!call tank%set(xl=3.32d0,yl=1.86d0,m=166,n=93)
!call tank%set(xl=3.4d0,yl=1.9d0,m=34,n=19)
!call tank%set(xl=3.32d0,yl=1.86d0,m=664,n=372)
!call tank%set(xl=3.74d0,yl=2.06d0,m=374,n=206)
!call tank%set(xl=0.44d0,yl=0.22d0,m=352,n=176)

!1call tank%set(xl=2.04d0,yl=1.52d0,m=102,n=76)
!2call tank%set(xl=2.02d0,yl=1.51d0,m=202,n=151)
call tank%set(xl = 0.00100d0,yl = 0.00105d0,m=40,n=42)



!call tank%set(xl=3.26d0,yl=1.82d0,m=32,n=18)

npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y
!      just for test 1 

!tank%zone = 2
!do i = 1, tank%m*tank%n
!   if(tank%x(i)<0.1.or.tank%x(i)>3.1.or.tank%y(i)<0.1) tank%zone(i) = 1
!   if(tank%zone(i)==1.and.tank%x(i)<0.1.and.tank%y(i)>0.1) tank%zone(i)=3
!   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<3.1) tank%zone(i)=4
!   if(tank%zone(i)==1.and.tank%x(i)>3.1) tank%zone(i)=5
!   if(tank%zone(i)==2.and.tank%x(i)>1.3)tank%zone(i)= 7
!   if(tank%zone(i)==2.and.tank%y(i)>0.7)tank%zone(i)= 7

!enddo


! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
   if(tank%y(i)>0.0010251) tank%zone(i) = 1   
   if(tank%y(i)<0.000025) tank%zone(i) = 3      
!   if(tank%y(i)>0.0010) tank%zone(i) = 1   
!   if(tank%y(i)<0.00005) tank%zone(i) = 3      
enddo

call parts%take_real(tank,2)
!call parts%take_boundary_for_tank_wedge3(tank)!1是原来的，2是加密的，3是楔形体也等距的
!3call parts%take_boundary_for_tank_wedge3(tank)!1是原来的，2是加密的，3是楔形体也等距的
call parts%take_virtual_Couette(tank,1)
call parts%take_virtual_Couette(tank,3)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = 2*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.001d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = 0.0!property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%x(2,i)>0.41)parts%p%r(i)=0.0
   if(parts%x(2,i)>0.41)parts%p%r(i)=0.0
!   if(parts%zone(i)/=2.and.parts%x(1,i)>0.11)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
!do i = 1,parts%ntotal +parts%nvirt 
!    if(tank%zone(i)==1) parts%vol%r(i) =  tank%dx*tank%dy/2.0d0
!    if(tank%zone(i)==3) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!    if(tank%zone(i)==4) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!enddo

do i = 1,parts%ntotal + parts%nvirt 
    if(parts%zone(i)/=2)then 
        parts%vol%r(i) = parts%vol%r(i)/4.0d0
!1        if(parts%x(1,i)<0.02.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)>2.02.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)==0.85) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
!1        if(parts%x(1,i)>1.32.and.parts%x(1,i)<1.34) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
!1        if(parts%x(1,i)>1.08.and.parts%x(1,i)<1.1) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0

!2        if(parts%x(1,i)<0.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!2        if(parts%x(1,i)>2.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!2        if(parts%x(1,i)==0.865) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!2        if(parts%x(1,i)==1.365) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!2        if(parts%x(1,i)>1.11.and.parts%x(1,i)<1.12) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
    endif
enddo

parts%mass = parts%vol * parts%rho

do i = parts%ntotal +1 ,parts%ntotal + parts%nvirt
if(parts%zone(i) ==1)then
    parts%n(1,i) = 0!(parts%x(2,i) - parts%x(2,i+1))/sqrt((parts%x(1,i) - parts%x(1,i+1))**2 + (parts%x(2,i) - parts%x(2,i+1))**2)
    parts%n(2,i) = -1.!(parts%x(1,i+1) - parts%x(1,i))/sqrt((parts%x(1,i) - parts%x(1,i+1))**2 + (parts%x(2,i) - parts%x(2,i+1))**2)
else
    parts%n(1,i) = 0!(parts%x(2,i) - parts%x(2,i+1))/sqrt((parts%x(1,i) - parts%x(1,i+1))**2 + (parts%x(2,i) - parts%x(2,i+1))**2)
    parts%n(2,i) = 1.!(parts%x(1,i+1) - parts%x(1,i))/sqrt((parts%x(1,i) - parts%x(1,i+1))**2 + (parts%x(2,i) - parts%x(2,i+1))**2)
endif
      !              nx = (parts%x(2,j) - parts%x(2,j+1))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
!              ny = (parts%x(1,j+1) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
enddo

return
end subroutine


!DEC$IF(.FALSE.)
    subroutine input

use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: property
double precision soil_surface,water_surface
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    

! Set nozzle and tank geometry parameters

!call tank%set(xl=3.34d0,yl=1.86d0,m=334,n=186)
!call tank%set(xl=3.34d0,yl=1.86d0,m=167,n=93)
!call tank%set(xl=3.32d0,yl=1.86d0,m=166,n=93)
!call tank%set(xl=3.4d0,yl=1.9d0,m=34,n=19)
!call tank%set(xl=3.32d0,yl=1.86d0,m=664,n=372)
!call tank%set(xl=3.74d0,yl=2.06d0,m=374,n=206)
!call tank%set(xl=0.44d0,yl=0.22d0,m=352,n=176)

!1call tank%set(xl=3.26d0,yl=1.82d0,m=163,n=91)
call tank%set(xl=3.26d0,yl=1.82d0,m=326,n=182)



!call tank%set(xl=3.26d0,yl=1.82d0,m=32,n=18)

npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y
!      just for test 1 

!tank%zone = 2
!do i = 1, tank%m*tank%n
!   if(tank%x(i)<0.1.or.tank%x(i)>3.1.or.tank%y(i)<0.1) tank%zone(i) = 1
!   if(tank%zone(i)==1.and.tank%x(i)<0.1.and.tank%y(i)>0.1) tank%zone(i)=3
!   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<3.1) tank%zone(i)=4
!   if(tank%zone(i)==1.and.tank%x(i)>3.1) tank%zone(i)=5
!   if(tank%zone(i)==2.and.tank%x(i)>1.3)tank%zone(i)= 7
!   if(tank%zone(i)==2.and.tank%y(i)>0.7)tank%zone(i)= 7

!enddo


! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
!1   if(tank%x(i)<0.02.or.tank%x(i)>3.24.or.tank%y(i)<0.02) tank%zone(i) = 1
!1   if(tank%zone(i)==1.and.tank%x(i)<0.02.and.tank%y(i)>0.02) tank%zone(i)=3
!1   if(tank%zone(i)==1.and.tank%x(i)>0.02.and.tank%x(i)<3.24) tank%zone(i)=4
!1   if(tank%zone(i)==1.and.tank%x(i)>3.24) tank%zone(i)=5
!1   if(tank%x(i)<0.02.and.tank%y(i)<0.02)  tank%zone(i)=6
!1   if(tank%x(i)>3.24.and.tank%y(i)<0.02)  tank%zone(i)=6
   
   if(tank%x(i)<0.01.or.tank%x(i)>3.25.or.tank%y(i)<0.01) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)<0.01.and.tank%y(i)>0.01) tank%zone(i)=3
   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<3.25) tank%zone(i)=4
   if(tank%zone(i)==1.and.tank%x(i)>3.25) tank%zone(i)=5
   if(tank%x(i)<0.01.and.tank%y(i)<0.01)  tank%zone(i)=6
   if(tank%x(i)>3.25.and.tank%y(i)<0.01)  tank%zone(i)=6
   
!   if(tank%zone(i)==3.and.tank%y(i)<0.08) tank%zone(i)=6
!   if(tank%zone(i)==4.and.tank%x(i)<0.08) tank%zone(i)=6
!   if(tank%zone(i)==4.and.tank%x(i)>3.18) tank%zone(i)=6
!   if(tank%zone(i)==5.and.tank%y(i)<0.08) tank%zone(i)=6

   !   if(tank%zone(i)==1.and.tank%y(i)>0.66)tank%zone(i)=3
!   if(tank%x(i)<0.06) tank%zone(i) = 1
!   if(tank%y(i)<0.06) tank%zone(i) = 4
!   if(tank%x(i)>3.28) tank%zone(i) = 6
!   if(tank%x(i)<0.06.and.tank%y(i)<0.06) tank%zone(i) = 3
!   if(tank%x(i)>3.28.and.tank%y(i)<0.06) tank%zone(i) = 5
   if(tank%zone(i)==2.and.tank%x(i)>1.21)tank%zone(i)= 8
   if(tank%zone(i)==2.and.tank%y(i)>0.61)tank%zone(i)= 8

!   if(tank%zone(i)==2.and.tank%x(i)>1.26)tank%zone(i)=4
!   if(tank%zone(i)==2.and.tank%y(i)>0.66)tank%zone(i)=4
enddo
!      write(*,*) tank%zone


call parts%take_real(tank,2)
!call parts%take_virtual(tank,1)
!call parts%take_virtual(tank,3)
!call parts%take_virtual(tank,4)
!call parts%take_virtual(tank,5)
!call parts%take_virtual(tank,6)

call parts%take_boundary_for_tank(tank)
!call parts%take_virtual(tank,4)
!call parts%take_virtual(tank,5)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = 1.5*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.61d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==3.and.parts%x(2,i)>0.62)parts%p%r(i)=0.0
   if(parts%x(1,i)>1.21)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
!do i = 1,parts%ntotal +parts%nvirt 
!    if(tank%zone(i)==1) parts%vol%r(i) =  tank%dx*tank%dy/2.0d0
!    if(tank%zone(i)==3) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!    if(tank%zone(i)==4) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!enddo

do i = 1,parts%ntotal + parts%nvirt 
    if(parts%zone(i)/=2)then 
        parts%vol%r(i) = parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)<0.02.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)>3.24.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0

        if(parts%x(1,i)<0.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(parts%x(1,i)>3.25.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
    endif
enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
end subroutine

!DEC$ENDIF



!DEC$IF(.FALSE.)
    
tank_wedge3
    subroutine input

!use param
use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: property
double precision soil_surface,water_surface
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    

! Set nozzle and tank geometry parameters

!call tank%set(xl=3.34d0,yl=1.86d0,m=334,n=186)
!call tank%set(xl=3.34d0,yl=1.86d0,m=167,n=93)
!call tank%set(xl=3.32d0,yl=1.86d0,m=166,n=93)
!call tank%set(xl=3.4d0,yl=1.9d0,m=34,n=19)
!call tank%set(xl=3.32d0,yl=1.86d0,m=664,n=372)
!call tank%set(xl=3.74d0,yl=2.06d0,m=374,n=206)
!call tank%set(xl=0.44d0,yl=0.22d0,m=352,n=176)

!1call tank%set(xl=2.04d0,yl=1.52d0,m=102,n=76)
!2call tank%set(xl=2.02d0,yl=1.51d0,m=202,n=151)
call tank%set(xl = 2.44d0,yl = 1.52d0,m=122,n=76)



!call tank%set(xl=3.26d0,yl=1.82d0,m=32,n=18)

npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y
!      just for test 1 

!tank%zone = 2
!do i = 1, tank%m*tank%n
!   if(tank%x(i)<0.1.or.tank%x(i)>3.1.or.tank%y(i)<0.1) tank%zone(i) = 1
!   if(tank%zone(i)==1.and.tank%x(i)<0.1.and.tank%y(i)>0.1) tank%zone(i)=3
!   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<3.1) tank%zone(i)=4
!   if(tank%zone(i)==1.and.tank%x(i)>3.1) tank%zone(i)=5
!   if(tank%zone(i)==2.and.tank%x(i)>1.3)tank%zone(i)= 7
!   if(tank%zone(i)==2.and.tank%y(i)>0.7)tank%zone(i)= 7

!enddo


! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
!1   if(tank%x(i)<0.02.or.tank%x(i)>2.02.or.tank%y(i)<0.02) tank%zone(i) = 1
!1   if(tank%zone(i)==1.and.tank%x(i)<0.02.and.tank%y(i)>0.02) tank%zone(i)=3
!1   if(tank%zone(i)==1.and.tank%x(i)>0.02.and.tank%x(i)<2.02) tank%zone(i)=4
!1   if(tank%zone(i)==1.and.tank%x(i)>2.02) tank%zone(i)=5
!1   if(tank%x(i)<0.02.and.tank%y(i)<0.02)  tank%zone(i)=6
!1   if(tank%x(i)>2.02.and.tank%y(i)<0.02)  tank%zone(i)=6
!1   if(tank%zone(i)==2.and.tank%y(i)>0.42)tank%zone(i)= 8
!1   if(tank%zone(i)==2.and.tank%x(i)>0.84.and.tank%y(i)<tank%x(i) - 0.84.and.tank%x(i)<1.1)tank%zone(i)= 8
!1   if(tank%zone(i)==2.and.tank%x(i)>1.1.and.tank%y(i)<1.33 - tank%x(i).and.tank%x(i)<1.34)tank%zone(i)= 8
!1   if(tank%zone(i)==2.and.tank%x(i)>2.02) tank%zone(i)= 8

!2   if(tank%x(i)<0.01.or.tank%x(i)>2.01.or.tank%y(i)<0.01) tank%zone(i) = 1
!2   if(tank%zone(i)==1.and.tank%x(i)<0.01.and.tank%y(i)>0.01) tank%zone(i)=3
!2   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<2.01) tank%zone(i)=4
!2   if(tank%zone(i)==1.and.tank%x(i)>2.01) tank%zone(i)=5
!2   if(tank%x(i)<0.01.and.tank%y(i)<0.01)  tank%zone(i)=6
!2   if(tank%x(i)>2.01.and.tank%y(i)<0.01)  tank%zone(i)=6
!2   if(tank%zone(i)==2.and.tank%x(i)>0.51)tank%zone(i)= 8
!2   if(tank%zone(i)==2.and.tank%y(i)>1.01)tank%zone(i)= 8
   
   if(tank%x(i)<0.02.or.tank%x(i)>2.42.or.tank%y(i)<0.02) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)<0.02.and.tank%y(i)>0.02) tank%zone(i)=3
   if(tank%zone(i)==1.and.tank%x(i)>0.02.and.tank%x(i)<2.42) tank%zone(i)=4
   if(tank%zone(i)==1.and.tank%x(i)>2.42) tank%zone(i)=5
   if(tank%x(i)<0.02.and.tank%y(i)<0.02)  tank%zone(i)=6
   if(tank%x(i)>2.42.and.tank%y(i)<0.02)  tank%zone(i)=6
!   if(tank%zone(i)==2.and.tank%x(i)>1.105.and.tank%x(i)<1.115.and.tank%y(i)>0.235) tank%zone(i)=7
!   if(tank%zone(i)==2.and.tank%x(i)>1.115) tank%zone(i)=8
   if(tank%zone(i)==2.and.tank%y(i)>0.42)tank%zone(i)= 8
   if(tank%zone(i)==2.and.tank%x(i)>0.87.and.tank%y(i)<tank%x(i) - 0.85.and.tank%x(i)<1.1)tank%zone(i)= 8
   if(tank%zone(i)==2.and.tank%x(i)>1.1.and.tank%y(i)<1.34 - tank%x(i).and.tank%x(i)<1.32)tank%zone(i)= 8
   if(tank%zone(i)==2.and.tank%x(i)>2.22) tank%zone(i)= 8
   
   
enddo

call parts%take_real(tank,2)

call parts%take_boundary_for_tank_wedge3(tank)!1是原来的，2是加密的，3是楔形体也等距的
!3call parts%take_boundary_for_tank_wedge3(tank)!1是原来的，2是加密的，3是楔形体也等距的
call parts%take_virtual(tank,7)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = 2*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.4d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%x(2,i)>0.41)parts%p%r(i)=0.0
   if(parts%x(2,i)>0.41)parts%p%r(i)=0.0
!   if(parts%zone(i)/=2.and.parts%x(1,i)>0.11)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
!do i = 1,parts%ntotal +parts%nvirt 
!    if(tank%zone(i)==1) parts%vol%r(i) =  tank%dx*tank%dy/2.0d0
!    if(tank%zone(i)==3) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!    if(tank%zone(i)==4) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!enddo

do i = 1,parts%ntotal + parts%nvirt 
    if(parts%zone(i)/=2)then 
        parts%vol%r(i) = parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)<0.02.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)>2.02.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)==0.85) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
!1        if(parts%x(1,i)>1.32.and.parts%x(1,i)<1.34) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
!1        if(parts%x(1,i)>1.08.and.parts%x(1,i)<1.1) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0

!2        if(parts%x(1,i)<0.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!2        if(parts%x(1,i)>2.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!2        if(parts%x(1,i)==0.865) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!2        if(parts%x(1,i)==1.365) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!2        if(parts%x(1,i)>1.11.and.parts%x(1,i)<1.12) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0

        if(parts%x(1,i)<0.02.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(parts%x(1,i)>2.22.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(parts%x(1,i)==0.87) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
        if(parts%x(1,i)>1.31.and.parts%x(1,i)<1.33) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
        if(parts%x(1,i)>1.09.and.parts%x(1,i)<1.1) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
    endif
enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
    end subroutine

!DEC$ENDIF



    

!DEC$IF(.FALSE.)
    subroutine input

use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: property
double precision soil_surface,water_surface
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    

! Set nozzle and tank geometry parameters

!call tank%set(xl=3.34d0,yl=1.86d0,m=334,n=186)
!call tank%set(xl=3.34d0,yl=1.86d0,m=167,n=93)
!call tank%set(xl=3.32d0,yl=1.86d0,m=166,n=93)
!call tank%set(xl=3.4d0,yl=1.9d0,m=34,n=19)
!call tank%set(xl=3.32d0,yl=1.86d0,m=664,n=372)
!call tank%set(xl=3.74d0,yl=2.06d0,m=374,n=206)
!call tank%set(xl=0.44d0,yl=0.22d0,m=352,n=176)

!1call tank%set(xl=3.26d0,yl=1.82d0,m=163,n=91)
call tank%set(xl=3.26d0,yl=1.82d0,m=326,n=182)



!call tank%set(xl=3.26d0,yl=1.82d0,m=32,n=18)

npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y
!      just for test 1 

!tank%zone = 2
!do i = 1, tank%m*tank%n
!   if(tank%x(i)<0.1.or.tank%x(i)>3.1.or.tank%y(i)<0.1) tank%zone(i) = 1
!   if(tank%zone(i)==1.and.tank%x(i)<0.1.and.tank%y(i)>0.1) tank%zone(i)=3
!   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<3.1) tank%zone(i)=4
!   if(tank%zone(i)==1.and.tank%x(i)>3.1) tank%zone(i)=5
!   if(tank%zone(i)==2.and.tank%x(i)>1.3)tank%zone(i)= 7
!   if(tank%zone(i)==2.and.tank%y(i)>0.7)tank%zone(i)= 7

!enddo


! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
!1   if(tank%x(i)<0.02.or.tank%x(i)>3.24.or.tank%y(i)<0.02) tank%zone(i) = 1
!1   if(tank%zone(i)==1.and.tank%x(i)<0.02.and.tank%y(i)>0.02) tank%zone(i)=3
!1   if(tank%zone(i)==1.and.tank%x(i)>0.02.and.tank%x(i)<3.24) tank%zone(i)=4
!1   if(tank%zone(i)==1.and.tank%x(i)>3.24) tank%zone(i)=5
!1   if(tank%x(i)<0.02.and.tank%y(i)<0.02)  tank%zone(i)=6
!1   if(tank%x(i)>3.24.and.tank%y(i)<0.02)  tank%zone(i)=6
   
   if(tank%x(i)<0.01.or.tank%x(i)>3.25.or.tank%y(i)<0.01) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)<0.01.and.tank%y(i)>0.01) tank%zone(i)=3
   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<3.25) tank%zone(i)=4
   if(tank%zone(i)==1.and.tank%x(i)>3.25) tank%zone(i)=5
   if(tank%x(i)<0.01.and.tank%y(i)<0.01)  tank%zone(i)=6
   if(tank%x(i)>3.25.and.tank%y(i)<0.01)  tank%zone(i)=6
   
!   if(tank%zone(i)==3.and.tank%y(i)<0.08) tank%zone(i)=6
!   if(tank%zone(i)==4.and.tank%x(i)<0.08) tank%zone(i)=6
!   if(tank%zone(i)==4.and.tank%x(i)>3.18) tank%zone(i)=6
!   if(tank%zone(i)==5.and.tank%y(i)<0.08) tank%zone(i)=6

   !   if(tank%zone(i)==1.and.tank%y(i)>0.66)tank%zone(i)=3
!   if(tank%x(i)<0.06) tank%zone(i) = 1
!   if(tank%y(i)<0.06) tank%zone(i) = 4
!   if(tank%x(i)>3.28) tank%zone(i) = 6
!   if(tank%x(i)<0.06.and.tank%y(i)<0.06) tank%zone(i) = 3
!   if(tank%x(i)>3.28.and.tank%y(i)<0.06) tank%zone(i) = 5
   if(tank%zone(i)==2.and.tank%x(i)>1.21)tank%zone(i)= 8
   if(tank%zone(i)==2.and.tank%y(i)>0.61)tank%zone(i)= 8

!   if(tank%zone(i)==2.and.tank%x(i)>1.26)tank%zone(i)=4
!   if(tank%zone(i)==2.and.tank%y(i)>0.66)tank%zone(i)=4
enddo
!      write(*,*) tank%zone


call parts%take_real(tank,2)
!call parts%take_virtual(tank,1)
!call parts%take_virtual(tank,3)
!call parts%take_virtual(tank,4)
!call parts%take_virtual(tank,5)
!call parts%take_virtual(tank,6)

call parts%take_boundary_for_tank(tank)
!call parts%take_virtual(tank,4)
!call parts%take_virtual(tank,5)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = 1.5*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.61d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==3.and.parts%x(2,i)>0.62)parts%p%r(i)=0.0
   if(parts%x(1,i)>1.21)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
!do i = 1,parts%ntotal +parts%nvirt 
!    if(tank%zone(i)==1) parts%vol%r(i) =  tank%dx*tank%dy/2.0d0
!    if(tank%zone(i)==3) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!    if(tank%zone(i)==4) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!enddo

do i = 1,parts%ntotal + parts%nvirt 
    if(parts%zone(i)/=2)then 
        parts%vol%r(i) = parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)<0.02.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!1        if(parts%x(1,i)>3.24.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0

        if(parts%x(1,i)<0.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(parts%x(1,i)>3.25.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
    endif
enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
end subroutine

!DEC$ENDIF


