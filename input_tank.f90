      subroutine input
      
      use param
      use declarations_sph
      implicit none     

      if (waterjet) call water_jet

      write(*,*)'  **************************************************'
      write(*,*)'  Initial particle configuration generated   '   
      write(*,*)'  Total number of real particles   ', parts%ntotal    
      write(*,*)'  Total number of virtual particles  ', parts%nvirt    
      write(*,*)'  **************************************************'

      if(.not.single_phase)then
         write(*,*) 'Total number of real soil particles:', soil%ntotal
         write(*,*) 'Total number of virtual particles:  ', soil%nvirt
      endif
   
      return
      end subroutine           
       

      subroutine water_jet

      use param
      use declarations_sph
      implicit none     
      
      integer  ntotal, bntotal,npoint
      integer i, j, d, k, ntotal_nozzle
      type(block) nozzle, tank
      type(material), pointer :: property
      double precision xleft,xright,ybottom,soil_surface,water_surface
      integer :: fluid_zone(10) = 0, fluid_boundary_zone(10) = 0
      integer ::  soil_zone(10) = 0,  soil_boundary_zone(10) = 0
      logical :: dbg = .true.
      double precision element_size, soil_submerged_depth    

      if(dbg) write(*,*) 'In water_jet...'      

!     Set nozzle and tank geometry parameters

      call tank%set(xl=0.5d0,yl=0.3d0,m=4,n=3)
      npoint = tank%m*tank%n
      allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
      call tank%cell_center
!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y

! Zoning
      tank%zone = 2
      do i = 1, tank%m*tank%n
         if(tank%x(i)<0.125.or.tank%x(i)>0.375.or.tank%y(i)<0.1) tank%zone(i) = 1
      enddo
!      write(*,*) tank%zone

      call parts%take_real(tank,2)
      call parts%take_virtual(tank,1)
      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

      parts%vol = tank%dx*tank%dy
      parts%hsml = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
      call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

      parts%vx = 0.d0

! ...Pressure. You must define the water-surface first.

      property => parts%material
      do i = 1,parts%ntotal+parts%nvirt
         parts%p(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
      enddo

      ! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
      call initial_density(parts)
      
      parts%mass = parts%vol * parts%rho     

      return
      end subroutine

      subroutine int_force1(parts)
!---------------------------------------------------------------------
!   Subroutine to calculate the internal forces on the right hand side 
!   of the Navier-Stokes equations, i.e. the pressure gradient and the
!   gradient of the viscous stress tensor, used by the time integration. 
!   Moreover the entropy production due to viscous dissipation, tds/dt, 
!   and the change of internal energy per mass, de/dt, are calculated. 

      use param
      use m_particles
      use m_sph_fo
      implicit none

      type(particles) parts

!      integer ntotal,niac
!      double precision,pointer,dimension(:) :: pair_i,pair_j
!      double precision,pointer,dimension(:,:) :: dwdx
!      double precision,pointer,dimension(:) :: sxx, sxy, syy,     &
!                                               szz, sxz, syz
!      double precision h, he 
!      integer d

      parts%dvx(1,:) = parts%dvx(1,:) - df(parts%p,'x',parts) +       &
                  df(parts%sxx,'x',parts) + df(parts%sxy,'y',parts)

      parts%dvx(2,:) = parts%dvx(2,:) - df(parts%p,'x',parts) +       &
                  df(parts%sxx,'x',parts) + df(parts%sxy,'y',parts)

      where (parts%rho.gt.0.0) parts%dvx(1,:) = parts%dvx(1,:)/parts%rho
      where (parts%rho.gt.0.0) parts%dvx(2,:) = parts%dvx(2,:)/parts%rho

      return
      end subroutine

