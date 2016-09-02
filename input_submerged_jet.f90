subroutine input

use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: wass,sand
double precision soil_surface,water_surface,nozzle_height
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    

! Set nozzle and tank geometry parameters

call tank%set(xl=0.6d0,yl=0.4d0,m=240,n=160)
npoint = tank%m*tank%n
tank%np= npoint
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
tank%x = tank%x - 0.3
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y

!call nozzle%set(xl=0.02d0,yl=0.1d0,m=16,n=80)
!npoint = nozzle%m*nozzle%n
!nozzle%np = npoint
!allocate(nozzle%x(npoint),nozzle%y(npoint),nozzle%zone(npoint))
!call nozzle%cell_center
!nozzle%x = nozzle%x - 0.01
!nozzle%y = nozzle%y + 0.22


water_surface = 0.3
soil_surface = 0.10
nozzle_height=0.2

write(*,*) 'Geometry generated!'
! Zoning
tank%zone = 0
do i = 1, tank%m*tank%n
   if(tank%x(i)<0.0)cycle
   tank%zone(i) = 2

   if(tank%x(i)<-0.29.or.tank%x(i)>0.29.or.tank%y(i)<0.01) tank%zone(i) = 1
!   if(tank%x(i)<0.025.or.tank%x(i)>0.475.or.tank%y(i)<0.025) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%y(i)>water_surface)tank%zone(i)=3
   if(tank%zone(i)==1.and.tank%y(i)>soil_surface)tank%zone(i)=7
!   if(tank%zone(i)==2.and.dabs(tank%x(i))<0.005.and.dabs(tank%x(i))>0.0025.and.tank%y(i)>0.12)tank%zone(i)=8
   if(tank%zone(i)==2.and.dabs(tank%x(i))<0.015.and.dabs(tank%x(i))>0.01.and.tank%y(i)>nozzle_height)tank%zone(i)=8
   !if(tank%zone(i)==2.and.dabs(tank%x(i))<0.005.and.tank%y(i)>0.12.and.tank%y(i)<0.2)tank%zone(i)=10
   !if(tank%zone(i)==2.and.dabs(tank%x(i))<0.005.and.tank%y(i)>0.2)tank%zone(i)=4
   if(tank%zone(i)==2.and.dabs(tank%x(i))<0.01.and.tank%y(i)>nozzle_height)tank%zone(i)=4
   if(tank%zone(i)==2.and.tank%y(i)>water_surface)tank%zone(i)=5
   if(tank%zone(i)==2.and.tank%y(i)>soil_surface)tank%zone(i)=6
   !if(tank%zone(i)==6.and.dabs(tank%x(i))<0.01.and.tank%y(i)>0.15.and.tank%y(i)<0.20)tank%zone(i)=9
enddo
!      write(*,*) tank%zone
!nozzle%zone = 4
!do i = 1, nozzle%m*nozzle%n    !nozzle not used!

!   if(nozzle%x(i)<0.0)cycle
!   nozzle%zone(i) = 4

!   if(nozzle%x(i)<-0.005.or.nozzle%x(i)>0.005) nozzle%zone(i) = 8
!   if(nozzle%zone(i)==4.and.nozzle%y(i)<0.3) nozzle%zone(i) = 10
!enddo

call parts%take_real(tank,2)
call parts%take_real(tank,6)
!call parts%take_real(nozzle,4)
!call parts%take_real(nozzle,10)
call parts%take_real(tank,4)
!call parts%take_real(tank,10)
call parts%take_virtual(tank,1)
call parts%take_virtual(tank,3)
call parts%take_virtual(tank,7)
!call parts%take_virtual(nozzle,8)
call parts%take_virtual(tank,8)

call soil%take_real(tank,2)
call soil%take_virtual(tank,1)
call soil%take_virtual(tank,3)
call soil%take_virtual(tank,7)
write(*,*) 'Particles generated!'
!call parts%setup_ndim1
!call soil%setup_ndim1

! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
parts%hsml = tank%dx*1.2
parts%dspp = tank%dx
!parts%dt   = dt

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype
!write(*,*) 'asdfa',parts%itype(10613), parts%zone(10613)

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0
do i = 1, parts%ntotal+parts%nvirt
   if(parts%zone(i)==4) parts%vx%y%r(i) = inlet_velocity 
enddo

! ...Stress. You must define the free surface first.
!water_surface = 0.3
wass => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = wass%rho0*parts%gravity*(parts%x(2,i)-water_surface)
   !if(parts%zone(i)==3.or.parts%zone(i)==4)parts%p%r(i)=0.0
   if(parts%p%r(i)<0.d0) parts%p%r(i) =0.d0
   if(parts%zone(i)==4) parts%p%r(i) = 0.d0
   !if(parts%zone(i)==8) parts%p%r(i) = 0.d0
enddo


parts%c = wass%c

! ...Volume fraction

do i = 1, parts%ntotal + parts%nvirt
   parts%vof%r(i) = 0.5
   if(parts%zone(i)==4)parts%vof%r(i) = 1.d0
   if(parts%zone(i)==6)parts%vof%r(i) = 1.d0
   if(parts%zone(i)==8)parts%vof%r(i) = 1.d0
   !if(parts%zone(i)==10)parts%vof%r(i) = 1.d0
enddo
if(parts%single_phase) parts%vof = 1.0
if(.not.parts%volume_fraction)  parts%vof = 1.0

! Density and Mass

call initial_density(parts)
parts%rho = parts%rho * parts%vof
parts%mass = parts%vol * parts%rho     

! ----------------------------------- For Soil --------------------------------
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

soil%vol = tank%dx*tank%dy
soil%hsml = tank%dx*1.2
soil%dspp = tank%dx
soil%dt = parts%dt

! itype is positive for real particles, negative for virtual particles.
 
call soil%setup_itype

! Set initial conditions of particles
! ...Velocity

soil%vx%x = 0.d0; soil%vx%y = 0.d0

! ...Stress. You must define the free surface first.
!soil_surface = 0.10
sand => soil%material
do i = 1,soil%ntotal+soil%nvirt
   soil%p%r(i) = (sand%rho0-wass%rho0)*soil%gravity*(soil%x(2,i)-soil_surface)
   !if(soil%zone(i)==3)soil%p%r(i)=0.0
   if(soil%x(2,i)>soil_surface)soil%p%r(i) = 0.d0
enddo
do i = 1, soil%ntotal+soil%nvirt   
   soil%str%xy%r(i) = 0.d0
   soil%str%x%r(i) =  0.1905*soil%p%r(i)
   soil%str%y%r(i) = -0.381*soil%p%r(i)
   soil%p%r(i)   =  0.619*soil%p%r(i)
enddo

soil%c = sand%c

! Volume fraction

soil%vof = 0.5d0
!do i = 1, soil%ntotal+soil%nvirt
!   if(soil%x(2,i)>0.1)soil%vof%r(i) = 0.d0
!enddo
if(.not.soil%volume_fraction)soil%vof = 1.0

! Density and Mass

soil%rho  = sand%rho0*soil%vof
soil%mass = soil%vol*soil%rho

return
end subroutine


!----------------------------------------
    subroutine inlet_boundary    ! For all particles present.
!----------------------------------------
use declarations_sph
implicit none
integer i, k

do i = 1, parts%ntotal+parts%nvirt
if(parts%zone(i)==4.and.parts%x(2,i)<0.22) parts%zone(i) = 100
if(parts%zone(i)==4.and.parts%x(2,i)>=0.22)then
   parts%vx%x%r(i) = 0.d0; parts%vx%y%r(i) = inlet_velocity
   parts%p%r(i) = 0.d0;    parts%vof%r(i) = 1.d0
   parts%rho%r(i) = h2o%rho0*parts%vof%r(i)
endif
enddo
return
end subroutine

!----------------------------------------
    subroutine inlet_boundary2  ! Particles flow in one by one.
!----------------------------------------
use declarations_sph
implicit none
integer i, k, ntotal
double precision nozzle_height,nozzle_length

nozzle_height = 0.2
nozzle_length = 0.2

ntotal = parts%ntotal + parts%nvirt
k = 0 
do i = 1, parts%ntotal+parts%nvirt
   !if(parts%zone(i)==4.and.parts%x(2,i)<0.20)then
   if(parts%zone(i)==4.and.parts%x(2,i)<nozzle_height)then
      parts%zone(i) = 100
      k = k + 1 
      parts%x(1,ntotal+k) = parts%x(1,i)
      parts%x(2,ntotal+k) = parts%x(2,i) + nozzle_length
      parts%vx%x%r(ntotal+k)= 0.d0
      parts%vx%y%r(ntotal+k)= inlet_velocity
      parts%p%r(ntotal+k)   = 0
      parts%vof%r(ntotal+k) = 1.0
      parts%vol%r(ntotal+k) = parts%vol%r(i)
      parts%hsml(ntotal+k) = parts%hsml(i)
      parts%rho%r(ntotal+k) = h2o%rho0*parts%vof%r(ntotal+k)
      parts%mass%r(ntotal+k)  = parts%mass%r(i)
      parts%itype(ntotal+k) = 2
      parts%zone(ntotal+k) = 4
      parts%c%r(ntotal+k) = h2o%c

      parts%rho_min%r(ntotal+k) = parts%rho%r(ntotal+k)
      parts%vof_min%r(ntotal+k) = parts%vof%r(ntotal+k)
      parts%v_min%x%r(ntotal+k) = parts%vx%x%r(ntotal+k)
      parts%v_min%y%r(ntotal+k) = parts%vx%y%r(ntotal+k)

   endif

if(parts%zone(i)==4.and.parts%x(2,i)>nozzle_height)then
   parts%vx%x%r(i) = 0.d0; parts%vx%y%r(i) = inlet_velocity
   parts%p%r(i) = 0.d0
   parts%vof%r(i) = 1.d0
   parts%rho%r(i) = h2o%rho0*parts%vof%r(i)
endif
enddo
parts%ntotal = parts%ntotal + k
call parts%setup_ndim1

return
end subroutine


!----------------------------------------
    subroutine symmetry_boundary(parts)  ! Verticle symmetry plane
!----------------------------------------
use m_particles
implicit none
type(particles) parts 
integer i, k,ntotal

ntotal = parts%ntotal + parts%nvirt
k = 0 
do i = 1, parts%ntotal + parts%nvirt
   if(parts%x(1,i)<0.05)then
      k = k + 1 
      parts%x(1,ntotal+k) = -parts%x(1,i)
      parts%x(2,ntotal+k) =  parts%x(2,i)
      parts%vx%x%r(ntotal+k)= -parts%vx%x%r(i)
      parts%vx%y%r(ntotal+k)=  parts%vx%y%r(i)
      parts%p%r(ntotal+k)   = parts%p%r(i)
      parts%vof%r(ntotal+k) = parts%vof%r(i)
      if(associated(parts%vof2))parts%vof2%r(ntotal+k) = parts%vof2%r(i)
      parts%vol%r(ntotal+k) = parts%vol%r(i)
      parts%hsml(ntotal+k) =  parts%hsml(i)
      parts%rho%r(ntotal+k) = parts%rho%r(i)
      parts%mass%r(ntotal+k)  = parts%mass%r(i)
      parts%itype(ntotal+k) = parts%itype(i)
      parts%zone(ntotal+k) = parts%zone(i)
      parts%c%r(ntotal+k) =  parts%c%r(i)
      if(associated(parts%psi)) parts%psi%r(ntotal+k) =  parts%psi%r(i)

      !parts%rho_min%r(ntotal+k) = parts%rho_min%r(i)
      !parts%p_min%r(ntotal+k)   = parts%p_min%r(i)
      !parts%vof_min%r(ntotal+k) = parts%vof_min%r(i)
      !parts%v_min%x%r(ntotal+k) =  parts%v_min%x%r(i)
      !parts%v_min%y%r(ntotal+k) = -parts%v_min%y%r(i)

      if(associated(parts%str))then
      parts%str%x%r(ntotal+k) =  parts%str%x%r(i)    !!! -?
      parts%str%y%r(ntotal+k) =  parts%str%y%r(i)
      parts%str%xy%r(ntotal+k)= -parts%str%xy%r(i)
      endif

      !parts%str_min%x%r(ntotal+k) = -parts%str_min%x%r(i)
      !parts%str_min%y%r(ntotal+k) =  parts%str_min%y%r(i)
      !parts%str_min%xy%r(ntotal+k)= -parts%str_min%xy%r(i)
      
      if(associated(parts%epsilon_p)) parts%epsilon_p%r(ntotal+k)= -parts%epsilon_p%r(i)
   endif
enddo

parts%nsymm = k

parts%nvirt = parts%nvirt + k
call parts%setup_ndim1

return
end subroutine
