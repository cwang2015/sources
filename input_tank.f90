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

call tank%set(xl=3.43424d0,yl=1.18052d0,m=96,n=33)
npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y
!      just for test 1 

! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
   if(tank%x(i)<0.11.or.tank%x(i)>3.32424.or.tank%y(i)<0.11) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)>1.30732)tank%zone(i)=3
   if(tank%zone(i)==2.and.tank%x(i)>1.30732)tank%zone(i)=4
   if(tank%zone(i)==2.and.tank%y(i)>0.70732)tank%zone(i)=4
enddo
!      write(*,*) tank%zone

call parts%take_real(tank,2)
call parts%take_virtual(tank,1)
call parts%take_virtual(tank,3)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.3
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==3)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
parts%mass = parts%vol * parts%rho     

return
end subroutine




