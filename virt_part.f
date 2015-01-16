       subroutine virt_part

c----------------------------------------------------------------------
c   Subroutine to determine the information of virtual particles
c   Here only the Monaghan type virtual particles for the 2D shear
c   cavity driven problem are generated.
c     itimestep : Current time step                                 [in]
c     ntotal : Number of particles                                  [in]
c     nvirt  : Number of virtual particles                         [out]
c     hsml   : Smoothing Length                                 [in|out]
c     mass   : Particle masses                                  [in|out]
c     x      : Coordinates of all particles                     [in|out]
c     vx     : Velocities of all particles                      [in|out]
c     rho    : Density                                          [in|out]
c     u      : internal energy                                  [in|out]
c     itype   : type of particles                               [in|out]

      use param
      use declarations_sph
      implicit none

      double precision, pointer, dimension(:) :: mass, vol, rho, p, u
      double precision, pointer, dimension(:,:) :: x, vx

      integer i, j, d, im, mp, n1, n2
      double precision xl, dx, v_inf, pore_water_pressure
      double precision, external :: rho_initial, rho_soil_initial

      x    => parts%x
      vol  => parts%vol
      mass => parts%mass
      vx   => parts%vx
      rho  => parts%rho
      p    => parts%p
      u    => parts%u

      if (vp_input) then          
                        
c        open(xv_vp,file="../../data/xv_vp.dat")
c        open(state_vp,file="../../data/state_vp.dat")
c        open(other_vp,file="../../data/other_vp.dat")            
c        read(xv_vp,*) nvirt
c        do j = 1, nvirt   
c          i = ntotal + j      
c          read(xv_vp,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
c          read(state_vp,*)im, mass(i), rho(i), p(i), u(i)        
c          read(other_vp,*)im, itype(i), hsml(i)                            
c        enddo  
c        close(xv_vp)
c        close(state_vp) 
c        close(other_vp) 
      
      else 
       
!        call parts%append_virtual_particles 
    
!        n1 = parts%ntotal+1
!        n2 = parts%ntotal+parts%nvirt

!        vx(1:2,n1:n2) = 0.
!        rho(n1:n2) = 1000.  !*0.3
!        p(n1:n2) = rho(n1:n2)*gravity*(x(2,n1:n2)+wasserjet%ds) !0.
!        u(n1:n2) = 357.1

!        do i = n1, n2
!           rho(i) = rho_initial(p(i))
!        enddo

       ! call initial_density(parts)
!        mass(n1:n2) = vol(n1:n2)*rho(n1:n2)
        

      endif   

c      if (mod(itimestep,save_step).eq.0) then
cc        open(1,file="../../data/xv_vp.dat")
cc        open(2,file="../../data/state_vp.dat")
cc        open(3,file="../../data/other_vp.dat")            
cc        write(1,*) nvirt
c        do i = ntotal + 1, ntotal + nvirt         
c         write(xv_vp,1001) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)              
c         write(state_vp,1002) i, mass(i), rho(i), p(i), u(i)
c         write(other_vp,1003) i, itype(i), hsml(i)                               
c        enddo       
c1001    format(1x, I6, 6(2x, e15.8))
c1002    format(1x, I6, 7(2x, e15.8)) 
c1003    format(1x, I6, 2x, I4, 2x, e15.8)
cc        close(1)
cc        close(2) 
cc        close(3) 
c      endif 

      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
         print *,' >> Statistics: Virtual boundary particles:'
         print *,'          Number of virtual particles:',parts%NVIRT
        endif     
      endif

c
c For soil
c
        if (single_phase)return

!        call soil%append_virtual_particles
    
!        n1 = soil%ntotal+1
!        n2 = soil%ntotal+soil%nvirt

!        soil%vx(1:2,n1:n2) = 0.
!        soil%rho(n1:n2) = 2560.  !*0.7
     
!        do i = n1,n2
!          pore_water_pressure = 1000.*gravity*  !*0.3
!     &                    (soil%x(2,i)+wasserjet%ds)
!          !pore_water_pressure = 0.
!          soil%p(i) = -pore_water_pressure+soil%rho(i)*gravity*
!     &                    (soil%x(2,i)+wasserjet%ds) !0.
!        enddo

!        soil%u(n1:n2) = 357.1
!        !soil%st(n1:n2)%xx = 0.

!        do i = n1, n2
!           soil%rho(i) = rho_soil_initial(soil%p(i))
!        enddo
!        soil%mass(n1:n2) = soil%vol(n1:n2)*soil%rho(n1:n2)

      return
      end subroutine
