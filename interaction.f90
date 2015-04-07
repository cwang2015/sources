      subroutine darcy_law(water, soil)

      use param, only: volume_fraction
      use m_particles
      implicit none

      type(particles) water, soil
      double precision dx(3), ks, ns, gw, cf, sp, rrw
      type(material), pointer :: h2o,sio2  
      type(numerical), pointer :: numeric
      double precision gravity   
      integer i, j, k, d, dim

      h2o => water%material
      sio2=>  soil%material
      numeric => water%numeric
      gravity = numeric%gravity
      dim = water%dim

      gw = h2o%rho0*(-gravity); ns = sio2%porosity
      ks = sio2%permeability
      cf = gw*ns/ks
      !cf = 6.e6
       
      do  k=1,water%niac
          i = water%pair_i(k)
          j = water%pair_j(k)  
          rrw = water%w(k)/(water%rho(i)*soil%rho(j))

! For staturated soil
!        if(volume_fraction) cf = water%vof(i)*water%rho(i)*(-gravity)/ks
        if(volume_fraction) cf = water%vof(i)*soil%vof(j)*water%rho(i)*(-gravity)/ks

          do d=1,dim
             sp = cf*(water%vx(d,i)-soil%vx(d,j))*rrw
             water%dvx(d,i) = water%dvx(d,i) - soil%mass(j)*sp
             soil%dvx(d,j)  = soil%dvx(d,j) + water%mass(i)*sp   
          enddo
      enddo

      return
      end subroutine 

!-----------------------------------------------------
      subroutine pore_water_pressure(water,soil)
!-----------------------------------------------------
      use param
      use m_particles
      implicit none
      
      type(particles) water, soil
      double precision mprr
      integer i, j, k, d
      
      do k = 1, water%niac
         i = water%pair_i(k)   ! water
         j = water%pair_j(k)   ! soil
         mprr = water%mass(i)*water%p(i)/(water%rho(i)*soil%rho(j))
!         mprr = water%mass(i)*(water%p(i)+soil%p(j))/       &      Bui2014
!                (water%rho(i)*soil%rho(j))
         do d = 1, dim
            soil%dvx(d,j) = soil%dvx(d,j) + mprr*water%dwdx(d,k)  &  !+
                            *soil%vof(j)  
         enddo

! saturated soil
         if(volume_fraction)then
         do d = 1, dim
            water%dvx(d,i) = water%dvx(d,i) -                     & ! Must be -
            soil%mass(j)*water%p(i)*soil%vof(j)*water%dwdx(d,k)/  &
            (water%rho(i)*soil%rho(j))
         enddo
         endif
      enddo

      return
      end subroutine       

!-------------------------------------------------------------
      subroutine volume_fraction_water(water, soil)
!-------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) water, soil
      integer i,j,k,d, ntotal
      type(material), pointer :: sio2

      sio2 => soil%material
      ntotal = water%ntotal+water%nvirt

      water%vof2 = 0.d0
      do k = 1, water%niac
         i = water%pair_i(k)
         j = water%pair_j(k)
         water%vof2(i) = water%vof2(i)+soil%mass(j)*water%w(k)
      enddo

      do k = 1, water%ntotal+water%nvirt
         water%vof2(k) = 1.d0 - water%vof2(k)/sio2%rho0
      enddo

      return
      end subroutine 

!-------------------------------------------------------------
      subroutine volume_fraction_water2(water, soil)
!-------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) water, soil
      integer i,j,k,d, ntotal
      type(material), pointer :: sio2
      double precision dvx(dim),tmp

      sio2 => soil%material
      ntotal = water%ntotal+water%nvirt

!      water%dvof = 0.d0

      do k = 1, water%niac
         i = water%pair_i(k)
         j = water%pair_j(k)
         do d = 1, dim
            dvx(d) = water%vx(d,i)-soil%vx(d,j)
         enddo 
         tmp = dvx(1)*water%dwdx(1,k)+dvx(2)*water%dwdx(2,k)
         water%dvof(i) = water%dvof(i)-soil%mass(j)*tmp/sio2%rho0
      enddo

      return
      end subroutine 

!--------------------------------------------------------------------
      subroutine volume_fraction_soil(parts)
!--------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts
      integer i,j,k, ntotal
      type(material), pointer :: sio2

      sio2 => parts%material
      ntotal = parts%ntotal + parts%nvirt

      do i = 1, ntotal
         parts%vof(i) = parts%rho(i)/sio2%rho0
      enddo

      return
      end subroutine

!--------------------------------------------------------------------------
      subroutine drag_force(water,soil)
!--------------------------------------------------------------------------
! For fixed porous media, drag force to fluid  

      !use param
      use m_particles
      implicit none

      type(particles) water, soil
      double precision dx(3), ks, ns, gw, cf, sp, rrw
      type(material), pointer :: h2o,sio2  
      type(numerical), pointer :: numeric
      double precision gravity   
      integer i, j, k, d, dim

      h2o => water%material
      sio2=>  soil%material
      numeric => water%numeric
      gravity = numeric%gravity
      dim = water%dim

      gw = h2o%rho0*(-gravity) 
      !ns = sio2%porosity; ks = sio2%permeability
      ns = 0.3; ks = 0.05      

      cf = gw*ns/ks
      !cf = 6.e6
      !write(*,*) ns,ks
       
      do i=1,water%ntotal+water%nvirt
         do d=1,dim
      water%dvx(d,i) = water%dvx(d,i) - water%vx(d,i)*cf/water%rho(i)
         enddo
      enddo

      return
      end subroutine 

