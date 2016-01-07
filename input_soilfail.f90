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

!call tank%set(xl=0.61d0,yl=0.21d0,m=488,n=168)
!!call tank%set(xl=1.01d0,yl=0.31d0,m=808,n=248)
!call tank%set(xl=1.01d0,yl=0.21d0,m=808,n=168)
call tank%set(xl=1.01d0,yl=0.11d0,m=808,n=88)
npoint = tank%m*tank%n
tank%np = npoint
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
tank%x = tank%x - 0.01d0; tank%y = tank%y - 0.01d0
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y

! Zoning
tank%zone = 2
do i = 1, tank%m*tank%n
   if(tank%x(i)<0.0.or.tank%x(i)>1.0.or.tank%y(i)<0.0) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)>0.2)tank%zone(i)=3    !0.1
   if(tank%zone(i)==2.and.tank%x(i)>0.2)tank%zone(i)=4    !0.1
enddo
!      write(*,*) tank%zone

call parts%take_real(tank,2)
call parts%take_virtual(tank,1)
call parts%take_virtual(tank,3)

!      write(*,*) parts%ntotal, parts%nvirt

!call parts%setup_ndim1

! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = tank%dx*1.2
parts%dspp = tank%dx
!parts%dt = dt

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0; parts%vx%y = 0.d0

! ...Stress. You must define the free surface first.
soil_surface = 0.1     !0.2
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*parts%gravity*(parts%x(2,i)-soil_surface)
   if(parts%zone(i)==3)parts%p%r(i)=0.0
enddo
do i = 1, parts%ntotal+parts%nvirt   
   parts%str%xy%r(i) = 0.d0
   parts%str%x%r(i) =  0.1905*parts%p%r(i)
   parts%str%y%r(i) = -0.381*parts%p%r(i)
   parts%p%r(i)   =  0.619*parts%p%r(i)
enddo

parts%c = property%c

! Density and Mass

parts%rho  = property%rho0
parts%mass = parts%vol * parts%rho     

return
end subroutine



