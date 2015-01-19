!DEC$IF(.FALSE.)
      subroutine p_gas(rho, u, p, c)
      
c----------------------------------------------------------------------
c   Gamma law EOS: subroutine to calculate the pressure and sound  
 
c     rho    : Density                                              [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
          
      implicit none
      double precision rho, u, p, c   
      double precision gamma 
          
c      For air (idea gas)

      gamma=1.4
      p = (gamma-1) * rho * u     
      c = sqrt((gamma-1) * u) 
     
      end         
      
      subroutine p_art_water(rho, p, c)
      
c----------------------------------------------------------------------
c   Artificial equation of state for the artificial compressibility 

c     rho    : Density                                              [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
c     Equation of state for artificial compressibility   

      implicit none
      double precision rho, u, p, c
      double precision gamma, rho0, b

c     Artificial EOS, Form 1 (Monaghan, 1994)
      gamma=7.
      rho0=1000.       
      b = 14290 ! 1.013e5
      p = b*((rho/rho0)**gamma-1)      
c      c = 1480.
      c = 10

c     Artificial EOS, Form 2 (Morris, 1997)
c      c = 0.01
c      p = c**2 * rho      
      
      end

      function rho_initial(p)
c---------------------------------------------------------------------
      implicit none
      double precision rho_initial, p
      double precision gamma, rho0, b

      gamma = 7.
      rho0=1000.  !*0.3
      b = 14290 !1.013e5

      rho_initial = rho0*(p/b+1)**(1/gamma)
 
      return
      end function

c------------------------For soil--------------------------------------

      subroutine p_art_soil(rho, p)
c----------------------------------------------------------------------
      implicit none
      double precision rho, u, p, c
      double precision k, rho0, b,pmax,g,hmax

      rho0=2560.
      g = -9.8
      hmax = -0.6
      pmax = rho0*g*hmax
      k = 10*pmax       

      p = k*((rho/rho0)-1)      
      
      return
      end subroutine

      function rho_soil_initial(p)      
c----------------------------------------------------------------------
      implicit none
      double precision rho, u, p, c
      double precision k, rho0, b, pmax, g, hmax, rho_soil_initial

      rho0=2560.   !*0.7
      g = -9.8
      hmax = -0.6

      pmax = rho0*g*hmax

      k = 10*pmax       
      rho_soil_initial = rho0*(p/k+1)      
      
      return
      end function

!DEC$ENDIF

      subroutine initial_density(parts)
!---------------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) parts
      type(material), pointer :: water, soil
      integer ntotal
      
      ntotal = parts%ntotal+parts%nvirt

      if(parts%imaterial=='water')then

         water => parts%material
         parts%rho(1:ntotal) = water%rho0*(parts%p(1:ntotal)/water%b+1) &
                            **(1/water%gamma)

      elseif(parts%imaterial=='soil')then

         soil => parts%material
         !parts%rho(1:ntotal) = soil%rho0*(parts%p(1:ntotal)/soil%k+1)

      endif
 
      return
      end subroutine

!----------------------------------------------------------------------
                      subroutine pressure(parts)
!----------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts
      type(material), pointer :: water, soil
      integer ntotal, i
      
      ntotal = parts%ntotal+parts%nvirt

      if(parts%imaterial=='water')then

         water => parts%material
         parts%p(1:ntotal) = water%b*((parts%rho(1:ntotal)/(water%rho0  &
                            *parts%vof(1:ntotal))) &   !!! False density
                            **water%gamma-1)

! Tension instability
                              if(water_tension_instability==1)then
         do i = 1, ntotal
            if(parts%p(i)<0)then
               parts%p(i) = 0.d0
               !parts%rho(i) = water%rho0*parts%vof(i)
            endif
         enddo
                              endif

         parts%c(1:ntotal) = water%c         

      elseif(parts%imaterial=='soil')then

         soil => parts%material
         if(soil_pressure==1)          &
         parts%p(1:ntotal) = soil%k*(parts%rho(1:ntotal)/soil%rho0-1)
!         parts%p(1:ntotal) = parts%p(1:ntotal)
!     &                      -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
         parts%c(1:ntotal) = soil%c

      endif

      return      
      end subroutine

!--------------------------------------------------------------
            subroutine pressure_water(parts)
!--------------------------------------------------------------
use m_particles
implicit none

type(particles) parts
type(material), pointer :: water, soil
integer ntotal, i
      
ntotal = parts%ntotal+parts%nvirt

water => parts%material
parts%p(1:ntotal) = water%b*((parts%rho(1:ntotal)/(water%rho0  &
                   *parts%vof(1:ntotal))) &   !!! False density
                  **water%gamma-1)

parts%c(1:ntotal) = water%c         

return
end subroutine

!--------------------------------------------------------------
         subroutine tension_instability_water(parts)
!--------------------------------------------------------------
use m_particles
implicit none

type(particles) parts
type(material), pointer :: water, soil
integer ntotal, i
      
ntotal = parts%ntotal+parts%nvirt
water => parts%material

do i = 1, ntotal
   if(parts%p(i)<0)then
      parts%p(i) = 0.d0
      !parts%rho(i) = water%rho0*parts%vof(i)
   endif
enddo

return
end subroutine


!------------------------------------------------------------
              subroutine pressure_soil(parts)
!------------------------------------------------------------
use m_particles
implicit none

type(particles) parts
type(material), pointer :: water, soil
integer ntotal, i
      
ntotal = parts%ntotal+parts%nvirt

soil => parts%material
parts%p(1:ntotal) = soil%k*(parts%rho(1:ntotal)/soil%rho0-1)
! parts%p(1:ntotal) = parts%p(1:ntotal)  &
!                    -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
parts%c(1:ntotal) = soil%c

return      
end subroutine
   
