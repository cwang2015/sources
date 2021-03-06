      subroutine sum_density(parts) 
C----------------------------------------------------------------------
C   Subroutine to calculate the density with SPH summation algorithm.

      use param 
      use m_particles
      implicit none

      type(particles) parts

      integer ntotal, niac
      integer, pointer, dimension(:) :: pair_i, pair_j, itype 
      double precision, pointer, dimension(:) :: hsml, mass, w, rho
 
      integer i, j, k, d      
      double precision selfdens, hv(dim), r, wi(maxn)     

      ntotal = parts%ntotal + parts%nvirt
      niac   = parts%niac
      pair_i =>parts%pair_i
      pair_j =>parts%pair_j
      itype  =>parts%itype
      hsml   =>parts%hsml
      mass   =>parts%mass
      w      =>parts%w
      rho    =>parts%rho

c     wi(maxn)---integration of the kernel itself
        
      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      
c     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call parts%kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*mass(i)/rho(i)
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        wi(i) = wi(i) + mass(j)/rho(j)*w(k)
        wi(j) = wi(j) + mass(i)/rho(i)*w(k)
      enddo

c     Secondly calculate the rho integration over the space

      do i=1,ntotal
        call parts%kernel(r,hv,hsml(i),selfdens,hv)
        rho(i) = selfdens*mass(i)
      enddo

c     Calculate SPH sum for rho:
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        rho(i) = rho(i) + mass(j)*w(k)
        rho(j) = rho(j) + mass(i)*w(k)
      enddo

c     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
!      if (nor_density) then 
        do i=1, ntotal
          rho(i)=rho(i)/wi(i)
        enddo
!      endif 
 
      end

      subroutine norm_density(parts) 

C----------------------------------------------------------------------
C   Subroutine to calculate the density with SPH summation algorithm.

C     ntotal : Number of particles                                  [in]
C     hsml   : Smoothing Length                                     [in]
C     mass   : Particle masses                                      [in]
C     niac   : Number of interaction pairs                          [in]
C     pair_i : List of first partner of interaction pair            [in]
C     pair_j : List of second partner of interaction pair           [in]
C     w      : Kernel for all interaction pairs                     [in]
c     itype   : type of particles                                   [in]
c     x       : Coordinates of all particles                        [in]
c     rho    : Density                                             [out]
   
      use param 
      use m_particles
      implicit none

      type(particles) parts

      integer ntotal, niac
      integer, pointer, dimension(:) :: pair_i, pair_j, itype 
      double precision, pointer, dimension(:) :: hsml, mass, w, rho
 
      integer i, j, k, d      
      double precision selfdens, hv(dim), r, wi(maxn)     

      ntotal = parts%ntotal + parts%nvirt
      niac   = parts%niac
      pair_i =>parts%pair_i
      pair_j =>parts%pair_j
      itype  =>parts%itype
      hsml   =>parts%hsml
      mass   =>parts%mass
      w      =>parts%w
      rho    =>parts%rho

c     wi(maxn)---integration of the kernel itself
        
      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      
c     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call parts%kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*mass(i)/rho(i)
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        wi(i) = wi(i) + mass(j)/rho(j)*w(k)
        wi(j) = wi(j) + mass(i)/rho(i)*w(k)
      enddo

c     Secondly calculate the rho integration over the space

!      do i=1,ntotal
!        call kernel(r,hv,hsml(i),selfdens,hv)
!        rho(i) = selfdens*mass(i)
!      enddo

c     Calculate SPH sum for rho:
!      do k=1,niac
!        i = pair_i(k)
!        j = pair_j(k)
!        rho(i) = rho(i) + mass(j)*w(k)
!        rho(j) = rho(j) + mass(i)*w(k)
!      enddo

c     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
!      if (nor_density) then 
        do i=1, ntotal
          rho(i)=rho(i)/wi(i)
        enddo
!      endif 
 
      end

      
      subroutine con_density(parts)

c----------------------------------------------------------------------
c     Subroutine to calculate the density with SPH continuiity approach.

c     ntotal : Number of particles                                  [in]
c     mass   : Particle masses                                      [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : derivation of Kernel for all interaction pairs       [in]
c     vx     : Velocities of all particles                          [in]
c     itype   : type of particles                                   [in]
c     x      : Coordinates of all particles                         [in]
c     rho    : Density                                              [in]
c     drhodt : Density change rate of each particle                [out]   

      use param
      use m_particles 
      implicit none
  
      type(particles) parts      

      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j, itype    
      double precision, pointer, dimension(:) :: mass, rho, drhodt
      double precision, pointer, dimension(:,:) :: x, dwdx, vx
      integer i,j,k,d    
      double precision vcc, dvx(dim) 
      
      ntotal = parts%ntotal + parts%nvirt
      niac   = parts%niac
      pair_i => parts%pair_i
      pair_j => parts%pair_j
      itype  => parts%itype
      x      => parts%x
      mass   => parts%mass
      dwdx   => parts%dwdx
      vx     => parts%vx
      rho    => parts%rho
      drhodt => parts%drho

      do i = 1, ntotal
        drhodt(i) = 0.
      enddo
     
      do k=1,niac      
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j) 
        enddo        
        vcc = dvx(1)*dwdx(1,k)        
        do d=2,dim
          vcc = vcc + dvx(d)*dwdx(d,k)
        enddo    
        drhodt(i) = drhodt(i) + mass(j)*vcc
        drhodt(j) = drhodt(j) + mass(i)*vcc       
      enddo    

      end
