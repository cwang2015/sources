subroutine input

use param
use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: wass,sand
double precision soil_surface,water_surface
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    

! Set nozzle and tank geometry parameters

!call tank%set(xl=0.5d0,yl=0.3d0,m=100,n=60)
call tank%set(xl=0.5d0,yl=0.3d0,m=40,n=24)
npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y

! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
!   if(tank%x(i)<0.01.or.tank%x(i)>0.49.or.tank%y(i)<0.01) tank%zone(i) = 1
   if(tank%x(i)<0.025.or.tank%x(i)>0.475.or.tank%y(i)<0.025) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%y(i)>0.20)tank%zone(i)=5
   if(tank%zone(i)==5.and.tank%y(i)>0.25)tank%zone(i)=6
   if(tank%zone(i)==2.and.tank%y(i)>0.20)tank%zone(i)=3
   if(tank%zone(i)==3.and.tank%y(i)>0.25)tank%zone(i)=4
enddo
!      write(*,*) tank%zone

call parts%take_real(tank,2)
call parts%take_real(tank,3)
call parts%take_virtual(tank,1)
call parts%take_virtual(tank,5)
call parts%take_virtual(tank,6)

call soil%take_real(tank,2)
call soil%take_virtual(tank,1)

! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = tank%dx*1.2

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx = 0.d0

! ...Stress. You must define the free surface first.
water_surface = 0.25
wass => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p(i) = wass%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==6)parts%p(i)=0.0
enddo


parts%c = wass%c

! ...Volume fraction

do i = 1, parts%ntotal + parts%nvirt
   parts%vof(i) = 1.0
   if(parts%zone(i)==1.or.parts%zone(i)==2)parts%vof(i) = 0.3
enddo
if(single_phase) parts%vof = 1.0
if(.not.volume_fraction)  parts%vof = 1.0

! Density and Mass

call initial_density(parts)
parts%rho = parts%rho * parts%vof
parts%mass = parts%vol * parts%rho     

! ----------------------------------- For Soil --------------------------------
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

soil%vol = tank%dx*tank%dy
soil%hsml = tank%dx*1.2

! itype is positive for real particles, negative for virtual particles.
 
call soil%setup_itype

! Set initial conditions of particles
! ...Velocity

soil%vx = 0.d0

! ...Stress. You must define the free surface first.
soil_surface = 0.20
sand => soil%material
do i = 1,soil%ntotal+soil%nvirt
   soil%p(i) = (sand%rho0-wass%rho0)*gravity*(soil%x(2,i)-soil_surface)
enddo
do i = 1, soil%ntotal+soil%nvirt   
   soil%sxy(i) = 0.d0
   soil%sxx(i) =  0.1905*soil%p(i)
   soil%syy(i) = -0.381*soil%p(i)
   soil%p(i)   =  0.619*soil%p(i)
enddo

soil%c = sand%c

! Volume fraction

soil%vof = 0.7
if(.not.volume_fraction)soil%vof = 1.0

! Density and Mass

soil%rho  = sand%rho0*soil%vof
soil%mass = soil%vol*soil%rho

return
end subroutine



