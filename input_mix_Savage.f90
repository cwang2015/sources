subroutine input

!use param
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

!call tank%set(xl=0.2d0,yl=0.1d0,m=80,n=40)
!call tank%set(xl=0.25d0,yl=0.125d0,m=100,n=50)
call tank%set(xl=0.26d0,yl=0.13d0,m=52,n=26)
!call tank%set(xl=0.26d0,yl=0.13d0,m=104,n=52)
!call tank%set(xl=0.26d0,yl=0.13d0,m=208,n=104)
!call tank%set(xl=0.26d0,yl=0.13d0,m=416,n=208)
npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!tank%x = tank%x-0.005; tank%y = tank%y - 0.005
tank%x = tank%x-0.01; tank%y = tank%y - 0.01
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y

! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
   if(tank%x(i)<0.d0.or.tank%x(i)>0.24.or.tank%y(i)<0.d0) tank%zone(i) = 1
!   if(tank%x(i)<0.025.or.tank%x(i)>0.475.or.tank%y(i)<0.025) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)>0.06)tank%zone(i)=3
   if(tank%zone(i)==1.and.tank%y(i)>0.048)tank%zone(i)=5
   if(tank%zone(i)==5.and.tank%y(i)>0.1)tank%zone(i)=6
   if(tank%zone(i)==3.and.tank%y(i)>0.1)tank%zone(i)=6
   if(tank%zone(i)==2.and.tank%x(i)>0.06)tank%zone(i)=4
   if(tank%zone(i)==2.and.tank%y(i)>0.048)tank%zone(i)=4
   if(tank%zone(i)==4.and.tank%y(i)>0.1)tank%zone(i)=7
enddo
!      write(*,*) tank%zone

call parts%take_real(tank,2)
call parts%take_real(tank,4)
call parts%take_virtual(tank,1)
call parts%take_virtual(tank,3)
call parts%take_virtual(tank,5)
call parts%take_virtual(tank,6)

call soil%take_real(tank,2)
call soil%take_virtual(tank,1)
call soil%take_virtual(tank,3)
call soil%take_virtual(tank,5)
call soil%take_virtual(tank,6)

call parts%setup_ndim1
call soil%setup_ndim1

! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = tank%dx*1.2
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Stress. You must define the free surface first.
water_surface = 0.1
wass => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = wass%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==6)parts%p%r(i)=0.0
enddo
!write(*,*) 'adf', parts%p%ndim1
!parts%p%ndim1 = ntotal

parts%c = wass%c

! ...Volume fraction

do i = 1, parts%ntotal + parts%nvirt
   parts%vof%r(i) = 1.0
   if(parts%zone(i)==1.or.parts%zone(i)==2)parts%vof%r(i) = 0.45
enddo
if(single_phase) parts%vof%r = 1.0
if(.not.volume_fraction)  parts%vof%r = 1.0
!parts%vof%ndim1 = ntotal

! Density and Mass
!write(*,*) ntotal,parts%vof%ndim1,parts%p%ndim1

call initial_density(parts)

!write(*,*) parts%rho%ndim1,parts%vof%ndim1

parts%rho = parts%rho * parts%vof
parts%mass = parts%vol * parts%rho     
! ----------------------------------- For Soil --------------------------------
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

soil%vol = tank%dx*tank%dy
soil%hsml = tank%dx*1.2
soil%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call soil%setup_itype

! Set initial conditions of particles
! ...Velocity

soil%vx%x = 0.d0
soil%vx%y = 0.d0

! ...Stress. You must define the free surface first.
soil_surface = 0.048
sand => soil%material
do i = 1,soil%ntotal+soil%nvirt
   soil%p%r(i) = (sand%rho0-wass%rho0)*gravity*(soil%x(2,i)-soil_surface)
   if(soil%zone(i)==3.or.soil%zone(i)==5.or.soil%zone(i)==6)soil%p%r(i)=0.0
enddo
do i = 1, soil%ntotal+soil%nvirt   
   soil%str%xy%r(i) = 0.d0
   soil%str%x%r(i) =  0.1905*soil%p%r(i)
   soil%str%y%r(i) = -0.381*soil%p%r(i)
   soil%p%r(i)   =  0.619*soil%p%r(i)
enddo
!soil%p%ndim1 = bntotal
!soil%sxx%ndim1 = bntotal; soil%sxy%ndim1 = bntotal; soil%syy%ndim1 = bntotal

soil%c = sand%c

! Volume fraction

soil%vof%r = 0.55
if(.not.volume_fraction)soil%vof%r = 1.0
!soil%vof%ndim1 = bntotal

! Density and Mass
write(*,*) "rrrrrr", soil%rho%ndim1, soil%vof%ndim1,bntotal
soil%rho  = sand%rho0*soil%vof
write(*,*) "tttttt"
soil%mass = soil%vol*soil%rho

return
end subroutine



