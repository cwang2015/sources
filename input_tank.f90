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
call tank%set(xl=3.34d0,yl=1.86d0,m=167,n=93)
!call tank%set(xl=3.32d0,yl=1.86d0,m=166,n=93)
!call tank%set(xl=3.4d0,yl=1.9d0,m=34,n=19)
!call tank%set(xl=3.32d0,yl=1.86d0,m=664,n=372)
!call tank%set(xl=3.74d0,yl=2.06d0,m=374,n=206)
!call tank%set(xl=0.44d0,yl=0.22d0,m=352,n=176)
npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y
!      just for test 1 

! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
!   if(tank%x(i)<0.06.or.tank%x(i)>3.28.or.tank%y(i)<0.06) tank%zone(i) = 1
!   if(tank%zone(i)==1.and.tank%x(i)>1.28)tank%zone(i)=3
!   if(tank%zone(i)==1.and.tank%y(i)>0.66)tank%zone(i)=3
   if(tank%x(i)<0.06) tank%zone(i) = 1
   if(tank%y(i)<0.06) tank%zone(i) = 4
   if(tank%x(i)>3.28) tank%zone(i) = 6
   if(tank%x(i)<0.06.and.tank%y(i)<0.06) tank%zone(i) = 3
   if(tank%x(i)>3.28.and.tank%y(i)<0.06) tank%zone(i) = 5
   if(tank%zone(i)==2.and.tank%x(i)>1.26)tank%zone(i)= 7
   if(tank%zone(i)==2.and.tank%y(i)>0.66)tank%zone(i)= 7
!   if(tank%zone(i)==2.and.tank%x(i)>1.26)tank%zone(i)=4
!   if(tank%zone(i)==2.and.tank%y(i)>0.66)tank%zone(i)=4
enddo
!      write(*,*) tank%zone

call parts%take_real(tank,2)
call parts%take_virtual(tank,1)
call parts%take_virtual(tank,3)
call parts%take_virtual(tank,4)
call parts%take_virtual(tank,5)
call parts%take_virtual(tank,6)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = 1.2*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.66d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==1.and.parts%x(2,i)>0.66)parts%p%r(i)=0.0
   if(parts%x(1,i)>1.26)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
parts%mass = parts%vol * parts%rho     

return
end subroutine




