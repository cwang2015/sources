      subroutine av_vel(parts)

c----------------------------------------------------------------------
c     Subroutine to calculate the average velocity to correct velocity
c     for preventing.penetration (monaghan, 1992)

c     ntotal : Number of particles                                  [in]
c     mass   : Particle masses                                      [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     w      : Kernel for all interaction pairs                     [in]
c     vx     : Velocity of each particle                            [in]
c     rho    : Density of each particle                             [in]
c     av     : Average velocityof each particle                    [out]
   
      use param 
      use m_particles
      implicit none

      type(particles) parts

      integer ntotal, niac
      integer, pointer, dimension(:) :: pair_i, pair_j
      double precision, pointer, dimension(:) :: mass, w, rho
      double precision, pointer, dimension(:,:) :: vx, av

      double precision   vcc, dvx(dim), epsilon
      integer i,j,k,d       

      ntotal = parts%ntotal
      niac   = parts%niac

      pair_i => parts%pair_i
      pair_j => parts%pair_j
      mass   => parts%mass
      w      => parts%w
      vx     => parts%vx
      rho    => parts%rho
      av     => parts%av
      
c     epsilon --- a small constants chosen by experence, may lead to instability.
c     for example, for the 1 dimensional shock tube problem, the E <= 0.3

      epsilon = 0.001  !0.5 !0.3
      
!      do i = 1, ntotal
!        do d = 1, dim
!          av(d,i) = 0.
!        enddo 
!      enddo

      av = 0.d0
!      !$omp parallel 
!      !$omp do private(k,i,j,d)
      do k=1,niac       
         i = pair_i(k)
         j = pair_j(k)       
         do d=1,dim
            dvx(d) = vx(d,i) - vx(d,j)            
            av(d, i) = av(d,i) - 2*mass(j)*dvx(d)/(rho(i)+rho(j))*w(k)
            av(d, j) = av(d,j) + 2*mass(i)*dvx(d)/(rho(i)+rho(j))*w(k)                      
         enddo                    
      enddo  
 !     !$omp end do
 !     !$omp do private(i,d)  
      do i = 1, ntotal
         do d = 1, dim
            av(d,i) = epsilon * av(d,i)
         enddo 
      enddo             
 !     !$omp end do
 !     !$omp end parallel
      return
      end subroutine
