      subroutine viscosity(parts)
c----------------------------------------------------------------------
c   Subroutine to define the fluid particle viscosity
 
c     ntotal  : Number of particles                                 [in]
c     itype    : Type of particle                                   [in]
c     x       : Coordinates of all particles                        [in]
c     rho     : Density                                             [in]
c     eta     : Dynamic viscosity                                  [out]

      use param
      use m_particles
      implicit none

      type(particles) parts

      integer ntotal,i
      integer, pointer, dimension(:) :: itype
      double precision, pointer, dimension(:) :: eta
      type(material), pointer :: liquid

      ntotal = parts%ntotal + parts%nvirt
      itype => parts%itype
      eta   => parts%eta
      liquid=> parts%material

      do i=1,ntotal
        if (abs(itype(i)).eq.1) then
          eta(i)=0.
        else if (abs(itype(i)).eq.2) then
          eta(i)=liquid%viscosity
        endif  
      enddo  
 
      return
      end subroutine

      subroutine shear_modulus(soil)
c-------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) soil
      type(material), pointer :: property
      double precision e,v
      integer i

      property => soil%material

      e = property%E
      v = property%niu

      do i = 1, soil%ntotal+soil%nvirt
         if(abs(soil%itype(i))==3)then
            soil%eta(i) = e/(2.0*(1+v)) 
         endif
      enddo

      return
      end subroutine

