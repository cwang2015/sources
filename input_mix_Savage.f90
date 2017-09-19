subroutine input

!use param
use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: wass,sand
double precision soil_surface,water_surface,initial_water_volume_fraction
logical :: dbg = .false.
double precision element_size, soil_submerged_depth, column_length    

! Set nozzle and tank geometry parameters

!call tank%set(xl=0.2d0,yl=0.1d0,m=80,n=40)
!call tank%set(xl=0.25d0,yl=0.125d0,m=100,n=50)
!call tank%set(xl=0.26d0,yl=0.13d0,m=52,n=26)
if(parts%particle_size==1)then
  call tank%set(xl=0.26d0,yl=0.13d0,m=104,n=52)     !0.0025
elseif(parts%particle_size==2)then
  call tank%set(xl=0.26d0,yl=0.13d0,m=208,n=104)     !0.00125
endif
!call tank%set(xl=0.26d0,yl=0.13d0,m=208,n=104)
!call tank%set(xl=0.26d0,yl=0.13d0,m=416,n=208)
!call tank%set(xl=0.26d0,yl=0.13d0,m=832,n=416)

npoint = tank%m*tank%n
tank%np= npoint
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
call tank%cell_center
!tank%x = tank%x-0.005; tank%y = tank%y - 0.005
tank%x = tank%x-0.01; tank%y = tank%y - 0.01
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y

! Zoning

! Loose or Dense packing:

soil_surface = 0.08 !0.042 !0.048
column_length = 0.06
initial_water_volume_fraction = 0.40 !0.45 0.40

tank%zone = 2
do i = 1, tank%m*tank%n
   if(tank%x(i)<0.d0.or.tank%x(i)>0.24.or.tank%y(i)<0.d0) tank%zone(i) = 1
!   if(tank%x(i)<0.025.or.tank%x(i)>0.475.or.tank%y(i)<0.025) tank%zone(i) = 1
   if(tank%zone(i)==1.and.tank%x(i)>column_length)tank%zone(i)=3
   if(tank%zone(i)==1.and.tank%y(i)>soil_surface)tank%zone(i)=5    !0.08 !0.048
   if(tank%zone(i)==5.and.tank%y(i)>0.1)tank%zone(i)=6      !0.1
   if(tank%zone(i)==3.and.tank%y(i)>0.1)tank%zone(i)=6      !0.1
   if(tank%zone(i)==2.and.tank%x(i)>column_length)tank%zone(i)=4
   if(tank%zone(i)==2.and.tank%y(i)>soil_surface)tank%zone(i)=4     !0.08  !0.048
   if(tank%zone(i)==4.and.tank%y(i)>0.1)tank%zone(i)=7      !0.1
enddo
!      write(*,*) tank%zone

call parts%take_real(tank,2)
call parts%take_real(tank,4)
call parts%take_virtual(tank,1)
call parts%take_virtual(tank,3)
call parts%take_virtual(tank,5)
call parts%take_virtual(tank,6)

!call parts%take_virtual(tank,7)   !!! Rigid lip assumption!

call soil%take_real(tank,2)

!! If Slip boundary used, comment these.

call soil%take_virtual(tank,1)
call soil%take_virtual(tank,3)
call soil%take_virtual(tank,5)
call soil%take_virtual(tank,6)

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

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Stress. You must define the free surface first.
water_surface = 0.1   
wass => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = wass%rho0*parts%gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==6)parts%p%r(i)=0.0
enddo
!write(*,*) 'adf', parts%p%ndim1
!parts%p%ndim1 = ntotal

parts%c = wass%c

! ...Volume fraction

do i = 1, parts%ntotal + parts%nvirt
   parts%vof%r(i) = 1.0
   if(parts%zone(i)==1.or.parts%zone(i)==2)parts%vof%r(i) =     &
   initial_water_volume_fraction
enddo
if(parts%single_phase) parts%vof%r = 1.0
if(.not.parts%volume_fraction)  parts%vof%r = 1.0

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
soil%dt = parts%dt

! itype is positive for real particles, negative for virtual particles.
 
call soil%setup_itype

! Set initial conditions of particles
! ...Velocity

soil%vx%x = 0.d0
soil%vx%y = 0.d0

! ...Stress. You must define the free surface first.
!soil_surface = 0.048
!soil_surface = 0.08
sand => soil%material
do i = 1,soil%ntotal+soil%nvirt
   soil%p%r(i) = (sand%rho0-wass%rho0)*soil%gravity*(soil%x(2,i)-soil_surface)
   if(soil%zone(i)==3.or.soil%zone(i)==5.or.soil%zone(i)==6)soil%p%r(i)=0.0
enddo
do i = 1, soil%ntotal+soil%nvirt   
   soil%str%xy%r(i) = 0.d0
   soil%str%x%r(i) =  0.1905*soil%p%r(i)
   soil%str%y%r(i) = -0.381*soil%p%r(i)
   soil%str%z%r(i) =  0.1905*soil%p%r(i)
   soil%p%r(i)   =  0.619*soil%p%r(i)
enddo
!soil%p%ndim1 = bntotal
!soil%sxx%ndim1 = bntotal; soil%sxy%ndim1 = bntotal; soil%syy%ndim1 = bntotal

soil%c = sand%c

! Volume fraction

soil%vof%r = 1.0-initial_water_volume_fraction
if(.not.soil%volume_fraction)soil%vof%r = 1.0

! Density and Mass
soil%rho  = sand%rho0*soil%vof
soil%mass = soil%vol*soil%rho

! Pore water pressure at soil particles
do i = 1, soil%ntotal + soil%nvirt
soil%spp0%r(i) = wass%rho0*parts%gravity*(soil%x(2,i)-water_surface)   
if(soil%zone(i)==3.or.soil%zone(i)==5.or.soil%zone(i)==6)soil%spp0%r(i)=0.0
enddo

call soil%cambridge_ini

return
end subroutine


!----------------------------------------
    subroutine inlet_boundary2  
!----------------------------------------
use declarations_sph
implicit none


return
end subroutine


!----------------------------------------
    subroutine symmetry_boundary(parts)  
!----------------------------------------
use m_particles
implicit none
type(particles) parts

return
end subroutine

!----------------------------------------
    subroutine slip_boundary(parts)  ! Verticle symmetry plane
!----------------------------------------
use m_particles
implicit none
type(particles) parts 
integer i, k,ntotal

ntotal = parts%ntotal + parts%nvirt
k = 0 
do i = 1, parts%ntotal + parts%nvirt
   if(parts%x(1,i)<0.01)then
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
      parts%str%xy%r(ntotal+k)=  -parts%str%xy%r(i)
      endif

      !parts%str_min%x%r(ntotal+k) = -parts%str_min%x%r(i)
      !parts%str_min%y%r(ntotal+k) =  parts%str_min%y%r(i)
      !parts%str_min%xy%r(ntotal+k)= -parts%str_min%xy%r(i)
      
      if(associated(parts%epsilon_p)) parts%epsilon_p%r(ntotal+k)= -parts%epsilon_p%r(i)  !!!+-?
   endif
enddo

parts%nsymm = k

parts%nvirt = parts%nvirt + k

! Bottom wall
ntotal = parts%ntotal + parts%nvirt
k = 0 
do i = 1, parts%ntotal + parts%nvirt
   if(parts%x(2,i)<0.01)then
      k = k + 1 
      parts%x(1,ntotal+k) = parts%x(1,i)
      parts%x(2,ntotal+k) = -parts%x(2,i)
      parts%vx%x%r(ntotal+k)= parts%vx%x%r(i)
      parts%vx%y%r(ntotal+k)=  -parts%vx%y%r(i)
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
      parts%str%xy%r(ntotal+k)=  -parts%str%xy%r(i)
      endif

      !parts%str_min%x%r(ntotal+k) = -parts%str_min%x%r(i)
      !parts%str_min%y%r(ntotal+k) =  parts%str_min%y%r(i)
      !parts%str_min%xy%r(ntotal+k)= -parts%str_min%xy%r(i)
      
      if(associated(parts%epsilon_p)) parts%epsilon_p%r(ntotal+k)= -parts%epsilon_p%r(i)  !!!+-?
   endif
enddo

parts%nsymm = parts%nsymm + k

parts%nvirt = parts%nvirt + k

!! Right wall

ntotal = parts%ntotal + parts%nvirt
k = 0 
do i = 1, parts%ntotal + parts%nvirt
   if(parts%x(1,i)>0.23)then
      k = k + 1 
      parts%x(1,ntotal+k) = parts%x(1,ntotal+k)+2*(0.24-parts%x(1,i))
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
      
      if(associated(parts%epsilon_p)) parts%epsilon_p%r(ntotal+k)= -parts%epsilon_p%r(i)  !!!+-?
   endif
enddo

parts%nsymm = parts%nsymm + k

parts%nvirt = parts%nvirt + k

!!-------------------------------------

call parts%setup_ndim1

return
end subroutine
