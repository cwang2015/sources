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


call tank%set(xl = 38.1d0,yl = 1.05d0,m=762,n=21)



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
   
   if(tank%x(i)<0.05.or.tank%x(i)>38.05.or.tank%y(i)<0.05) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)<0.05)tank%zone(i)= 3
   if(tank%zone(i)==1.and.tank%x(i)>0.05.and.tank%x(i)<38.05)tank%zone(i)= 4
   if(tank%zone(i)==1.and.tank%x(i)>38.05)tank%zone(i)= 5
   if(tank%zone(i)==2.and.tank%x(i)>15.51)tank%zone(i)= 8
   if(tank%zone(i)==2.and.tank%y(i)>0.76)tank%zone(i)= 8
   
enddo

call parts%take_real(tank,2)

!call parts%take_boundary_for_tank_wedge3(tank)

call parts%take_boundary_for_tank_wedge_freesurface(tank)

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
water_surface = 1.01d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==3.and.parts%x(2,i)>1.01)parts%p%r(i)=0.0
   if(parts%x(1,i)>0.511)parts%p%r(i)=0.0
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
!1        if(parts%x(1,i)==0.85) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!1        if(parts%x(1,i)==1.33) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!1        if(parts%x(1,i)>1.08.and.parts%x(1,i)<1.1) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0

!2        if(parts%x(1,i)<0.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!2        if(parts%x(1,i)>2.01.and.parts%x(2,i)<0.01) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!2        if(parts%x(1,i)==0.865) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!2        if(parts%x(1,i)==1.365) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
!2        if(parts%x(1,i)>1.11.and.parts%x(1,i)<1.12) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0

        if(parts%x(1,i)<0.015.and.parts%x(2,i)<0.015) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(parts%x(1,i)>2.22.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(parts%x(1,i)==0.86) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
        if(parts%x(1,i)>1.21.and.parts%x(1,i)<1.22) parts%vol%r(i)=parts%vol%r(i)*3.0d0/4.0d0
        if(parts%x(1,i)>1.03.and.parts%x(1,i)<1.04) parts%vol%r(i)=parts%vol%r(i)*3.0d0/2.0d0
    endif
enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
end subroutine




