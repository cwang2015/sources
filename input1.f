      subroutine input
      
c----------------------------------------------------------------------
c     Subroutine for loading or generating initial particle information

c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c     hsml-- smoothing lengths of particles                        [out]
c     ntotal-- total particle number                               [out]
      
      use param
      use declarations_sph
      implicit none     

      integer i, d, im       

c     load initial particle information from external disk file

      if(config_input) then    
                        
c        open(1,file="../f_xv.dat")
c        open(2,file="../f_state.dat")
c        open(3,file="../f_other.dat")        
      
c        write(*,*)'  **************************************************'
c        write(*,*)'      Loading initial particle configuration...   '       
c        read (1,*) ntotal 
c        write(*,*)'      Total number of particles   ', ntotal    	
c        write(*,*)'  **************************************************'	
c        do i = 1, ntotal         
c          read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
c          read(2,*)im, mass(i), rho(i), p(i), u(i)        
c          read(3,*)im, itype(i), hsml(i)                                        
c        enddo
       
c        close(1)
c        close(2)
c        close(3)         
          
      else 
          
c     generate particle configuration
       
c        if (shocktube) call shock_tube(x, vx, mass, rho, p, u, 
c     &                      itype, hsml, ntotal)                

c        if (shearcavity) call shear_cavity(x, vx, mass, rho, p, u, 
c     &                        itype, hsml, ntotal)

        if (waterjet) call water_jet

        write(*,*)'  **************************************************'
        write(*,*)'   Initial particle configuration generated   '   
        write(*,*)'   Total number of real particles   ', parts%ntotal    
        write(*,*)'   Total number of virtual particles  ', parts%nvirt    
        write(*,*)'  **************************************************'

        if(.not.single_phase)then
         write(*,*) 'Total number of real soil particles:', soil%ntotal
         write(*,*) 'Total number of virtual particles:  ', soil%nvirt
        endif

      endif
   
      return
      end subroutine           
       
       
      subroutine shock_tube(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =1   ideal gas
c     hsml-- smoothing lengths of particles                        [out]
c     ntotal-- total particle number                               [out]

      use param
      implicit none     
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, d
      double precision space_x     

      ntotal=400
      space_x=0.6/80.      
      
      do i=1,ntotal
        mass(i)=0.75/400.
        hsml(i)=0.015
        itype(i)=1
        do d = 1, dim
          x(d,i) = 0. 
          vx(d,i) = 0.
        enddo        
      enddo                
                
      do i=1,320
        x(1,i)=-0.6+space_x/4.*(i-1)
      enddo
      
      do i=320+1,ntotal
        x(1,i)=0.+space_x*(i-320)
      enddo               
                           
      do i=1,ntotal
        if (x(1,i).le.1.e-8) then
          u(i)=2.5
          rho(i)=1.
          p(i)=1.
        endif 
        if (x(1,i).gt.1.e-8)  then
          u(i)=1.795
          rho(i)=0.25
          p(i)=0.1795
        endif        
      enddo        
	               
      end
      
      subroutine shear_cavity(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2 d shear driven cavity probem with Re = 1
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =2   water
c     h-- smoothing lengths of particles                           [out]
c     ntotal-- total particle number                               [out]

      use param
      implicit none     
C      include 'param.inc'
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy

c     Giving mass and smoothing length as well as other data.

      m = 41
      n = 41
      mp = m-1
      np = n-1
      ntotal = mp * np
      xl = 1.e-3
      yl = 1.e-3
      dx = xl/mp
      dy = yl/np

      do i = 1, mp
         do j = 1, np
             k = j + (i-1)*np
             x(1, k) = (i-1)*dx + dx/2.
             x(2, k) = (j-1)*dy + dy/2.
        enddo
      enddo

      do i = 1, mp*np
         vx(1, i) = 0.
         vx(2, i) = 0.      
        rho (i) = 1000.   
        mass(i) = dx*dy*rho(i)  
        p(i)= 0.   
        u(i)=357.1
        itype(i) = 2
        hsml(i) = dx
      enddo  

      end	 

      subroutine water_jet

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2 d water jet impacting probem with Re = ?
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =2   water
c     h-- smoothing lengths of particles                           [out]
c     ntotal-- total particle number                               [out]

      use param
      use declarations_sph
      implicit none     
      
      integer  ntotal, bntotal
      integer i, j, d, k, ntotal_nozzle
      type(block) nozzle, tank
      type(material), pointer :: property
      double precision xleft,xright,ybottom
      logical :: dbg = .false.
      real AAA

      if(dbg) write(*,*) 'In water_jet...'      

c     Set nozzle and tank geometry parameters

      if(trim(particle_size)=='small') 
     &   call wasserjet%set(dr=0.2d0,L=2.d0,H=0.6d0,mdr=10,dj=10.d0,
     &        mdj=500,ds=0.4d0, mL=100,mH=30)
      if(trim(particle_size)=='large')
     &   call wasserjet%set(dr=0.2d0,L=2.d0,H=0.6d0,mdr=5,dj=10.d0,
     &        mdj=250,ds=0.4d0, mL=50,mH=15)
       if(trim(particle_size)=='smallnozzle')
     &   call wasserjet%set(dr=0.1d0,L=2.d0,H=0.6d0,mdr=5,dj=10.d0,
     &        mdj=500,ds=0.4d0, mL=100,mH=30)
        if(trim(particle_size)=='bui')
!     &   call wasserjet%set(dr=0.02d0,L=0.5d0,H=0.2d0,mdr=4,dj=2.d0,
!     &        mdj=400,ds=0.1d0, mL=100,mH=40)
     &   call wasserjet%set(dr=0.02d0,L=0.5d0,H=0.2d0,mdr=2,dj=2.d0,
     &        mdj=200,ds=0.1d0, mL=50,mH=20)

!     &   call wasserjet%set(dr=0.02d0,L=0.44d0,H=0.22d0,mdr=4,dj=2.d0,    !Koshi
!     &        mdj=400,ds=0.d0, mL=88,mH=44)


        if(trim(particle_size)=='bigbox')
     &   call wasserjet%set(dr=0.02d0,L=1.d0,H=0.2d0,mdr=2,dj=2.d0,
     &        mdj=200,ds=0.1d0, mL=100,mH=20)
        if(trim(particle_size)=='dambreak')
     &   call wasserjet%set(L=0.8d0,H=0.2d0, ds=0.1d0, mL=160,mH=40)
        if(trim(particle_size)=='damsmall')
!     &   call wasserjet%set(L=0.8d0,H=0.06d0, ds=0.24d0, mL=320,mH=24)  !smallest  4:1
!     &   call wasserjet%set(L=0.8d0,H=0.12d0, ds=0.18d0, mL=320,mH=48)  !smallest  2:1
     &   call wasserjet%set(L=1.2d0,H=0.11d0, ds=0.d0, mL=240,mH=22)  !coarse alumibar  2:1
!     &   call wasserjet%set(L=1.2d0,H=0.11d0, ds=0.d0, mL=4800,mH=44)  !fine alumibar  2:1
        

!     &   call wasserjet%set(L=0.8d0,H=0.06d0, ds=0.24d0, mL=160,mH=12)  !coarse     4:1
!     &   call wasserjet%set(L=0.8d0,H=0.03d0, ds=0.27d0, mL=320,mH=12)  !smallest  8:1
!     &   call wasserjet%set(L=0.8d0,H=0.11d0, ds=0.19d0, mL=160,mH=22) !cube       1:1


C     Nozzle particle entrance

      call nozzle%set(xl = wasserjet%dr, yl = wasserjet%dj,
     &                 m = wasserjet%mdr, n = wasserjet%mdj)
      k = nozzle%m*nozzle%n
      allocate(nozzle%x(k),nozzle%y(k))
      call nozzle%cell_center
      nozzle%x = nozzle%x - wasserjet%dr/2.
      nozzle%y = nozzle%y - wasserjet%immerse

C     Paticles in tank
 
      call tank%set(xl = wasserjet%L, yl = wasserjet%H,  
     &               m = wasserjet%mL, n = wasserjet%mH)
      k = tank%m*tank%n
      allocate(tank%x(k),tank%y(k))
      call tank%cell_center       
      tank%x = tank%x  - wasserjet%L/2.
      tank%y = (tank%y + wasserjet%ds)*(-1.0)

c -----------------------Phase I------------------------------

      xleft = -0.9; xright = 0.9; ybottom = -0.9
      if(trim(particle_size)=='bui')then
         xleft = -0.225; xright = 0.225; ybottom = -0.275
      !   xleft = -0.2 ;  xright = 0.2; ybottom = -0.2     !Koshi
      elseif(trim(particle_size)=='bigbox')then
         xleft = -0.475; xright = 0.475; ybottom = -0.275
      elseif(trim(particle_size)=='dambreak')then
         xleft = -0.1; xright = 0.1; ybottom = -0.275
      endif

      if(trim(particle_size)=='damsmall')then
!         xleft = -0.05; xright = 0.05; ybottom = -0.29   !cube
!         xleft = -0.1; xright = 0.1; ybottom = -0.295      ! smallest
          xleft = -0.2; xright = 0.2; ybottom = -0.1
      endif

      k = 0
      ntotal_nozzle = 0
                                          if(with_nozzle)then 
      ntotal_nozzle = nozzle%m*nozzle%n
      do i =1, ntotal_nozzle
         k = k + 1
         parts%x(1, k) = nozzle%x(i)
         parts%x(2, k) = nozzle%y(i)
      enddo
                                          endif !nozzle

      do i = 1, tank%m*tank%n
      if(tank%x(i)<xleft.or.tank%x(i)>xright.or.tank%y(i)<ybottom)then ! origional
!      if(tank%x(i)<xleft.or.tank%x(i)>-0.1.or.tank%y(i)<ybottom)then    ! water dam
         else
            k = k + 1
            parts%x(1,k) = tank%x(i)
            parts%x(2,k) = tank%y(i)
         endif
      enddo
      parts%ntotal = k

      if(trim(particle_size)=='dambreak')then
         xleft = -9999.0; xright = 9999.0
      endif
      
      if(trim(particle_size)=='damsmall')then   !!! IFORT is shit!
         xleft = -9999.0d0
         xright = 9999.0d0
      endif

!      if(trim(particle_size)=='damsmall')then ! alumibar
!      do i = 1, tank%m*tank%n
!         if(tank%x(i)>-0.01.and.tank%x(i)<0.and.tank%y(i)>ybottom)then
!            k = k + 1
!            parts%x(1,k) = tank%x(i)
!            parts%x(2,k) = tank%y(i)
!         endif
!      enddo
!      endif     

      do i = 1, tank%m*tank%n
      if(tank%x(i)<xleft.or.tank%x(i)>xright.or.tank%y(i)<ybottom)then
!       if(tank%x(i)>-0.01.and.tank%y(i)<ybottom)then  ! alumibar
            k = k + 1
            parts%x(1,k) = tank%x(i)
            parts%x(2,k) = tank%y(i)
        endif
      enddo
      parts%nvirt = k - parts%ntotal

c      parts%nvirt = 0
     
c Basic settings for particles
      ntotal = parts%ntotal; bntotal = parts%nvirt

      parts%vol(1:ntotal+bntotal) = tank%dx*tank%dy
      parts%hsml(1:ntotal+bntotal)= tank%dx
      parts%zone(1:ntotal_nozzle) = 1

      if(trim(parts%imaterial)=='water')then
         parts%itype(1:ntotal) = 2
         parts%itype(ntotal+1:ntotal+bntotal) = -2
      elseif(trim(parts%imaterial)=='soil')then
         parts%itype(1:ntotal) = 3
         parts%itype(ntotal+1:ntotal+bntotal) = -3
      endif

c      do i = 1, parts%ntotal
c         if(parts%x(1,i)<-0.9.or.parts%x(1,i)>0.9.or.parts%x(2,i)<-0.9)
c     &      parts%itype(i) = -2
c      enddo

c     Set initial conditions of particles

      property => parts%material
      do i = 1,parts%ntotal+parts%nvirt
         parts%vx(1,i) = 0.
         parts%vx(2,i) = 0.
         parts%rho(i) = property%rho0
         parts%p(i) = parts%rho(i)*gravity*(parts%x(2,i)+wasserjet%ds)  !0.
         
         !if(parts%x(1,i)>-0.1) parts%p(i) = 0.d0                     ! water dam
         !if(parts%x(1,i)>-0.1.or.parts%itype(i)<0) parts%p(i) = 0.d0  ! water dam

         if(parts%x(2,i)>-wasserjet%ds)then
            parts%vx(2,i) = inlet_velocity
            parts%p(i) = 0
         endif
         if(trim(particle_size)=='damsmall')then
!           if(abs(parts%x(1,i))>0.1) parts%p(i)=0.d0  ! cube 0.05 !smallest 0.1
           if(abs(parts%x(1,i))>0.2) parts%p(i)=0.d0  ! cube 0.05 !smallest 0.1
         endif
      enddo

      if(trim(parts%imaterial)=='soil')then
         do i = 1, parts%ntotal+parts%nvirt   
            !parts%sxx(i) = -0.3*parts%p(i)
            !parts%syy(i) = 0.d0
            parts%sxy(i) = 0.d0
            parts%sxx(i) = 0. !0.1905*parts%p(i)
            parts%syy(i) = 0. !-0.381*parts%p(i)
            parts%p(i)   = 0.  !0.619*parts%p(i)
         enddo
      endif

      call initial_density(parts)

      parts%vof(1:ntotal_nozzle) = 1.0
      do i = ntotal_nozzle+1, ntotal+bntotal
         parts%vof(i) = 0.5
         !if(parts%x(2,i)>-0.16)parts%vof(i) = 1           !!! 0.16
         !parts%vof(ntotal_nozzle+1:ntotal+bntotal) =0.7  !!! origional
      enddo
      if(single_phase) parts%vof(ntotal_nozzle+1:ntotal+bntotal) = 1.0
      if(.not.volume_fraction) 
     &   parts%vof(ntotal_nozzle+1:ntotal+bntotal) = 1.0

      parts%rho(1:ntotal+bntotal) = parts%rho(1:ntotal+bntotal)*
     &                              parts%vof(1:ntotal+bntotal)

      parts%mass(1:ntotal+bntotal) = parts%vol(1:ntotal+bntotal)*
     &                               parts%rho(1:ntotal+bntotal)

!      parts%mass(1:ntotal_nozzle) = 
!     &parts%vol(1:ntotal_nozzle)*parts%rho(1:ntotal_nozzle)
     
!      parts%mass(ntotal_nozzle+1:ntotal+bntotal) = 
!     &parts%vol(ntotal_nozzle+1:ntotal+bntotal)*0.7*
!     &parts%rho(ntotal_nozzle+1:ntotal+bntotal)

c-------------------------------soil input---------------------------------

      if(single_phase)return

      k = 0; ntotal = 0; bntotal = 0
      do i = 1, parts%ntotal+parts%nvirt
         if(parts%x(2,i)>-wasserjet%ds)cycle
         !if(parts%x(2,i)>-0.16) cycle                    !!!0.16
         k = k + 1
         soil%x(:,k)   = parts%x(:,i)
         soil%vol(k)   = parts%vol(i)
         soil%hsml(k)  = parts%hsml(i) 
c
         soil%p(k)     = parts%p(i)
c
         if(parts%itype(i)==-2)then
            bntotal = bntotal + 1
            soil%itype(k) = -3
         elseif(parts%itype(i)==2)then
            ntotal = ntotal + 1
            soil%itype(k) =  3
         endif
      enddo
      soil%ntotal = ntotal; soil%nvirt = bntotal

c set initial conditions for soil particles

      property => soil%material
      do i = 1, soil%ntotal+soil%nvirt
         soil%vx(:,i) = 0.
         soil%rho(i)  = property%rho0
         if(dry_soil)then
         soil%p(i) = soil%rho(i)*gravity*(soil%x(2,i)+wasserjet%ds)
         else
         soil%p(i) = soil%rho(i)*gravity*(soil%x(2,i)+wasserjet%ds)  
     &             - soil%p(i)                                        ! origional
!         soil%p(i) = (soil%rho(i)-1000)*gravity*
!     &               (soil%x(2,i)+wasserjet%ds+0.06)                   !!!0.16
         endif      

         !if(soil%x(1,i)>-0.1) soil%p(i) = 0.d0                     ! water dam
         !if(soil%x(1,i)>-0.1.or.soil%itype(i)<0) soil%p(i) = 0.d0   ! water dam

      enddo


      do i = 1, soil%ntotal+soil%nvirt   
         soil%sxx(i) = -0.3*soil%p(i)
         soil%syy(i) = 0.d0
         soil%sxy(i) = 0.d0
      enddo

      call initial_density(soil)       

      soil%vof(1:ntotal+bntotal) = 0.5  !0.1
      if(.not.volume_fraction)soil%vof(1:ntotal+bntotal) = 1.0 

      soil%rho(1:ntotal+bntotal) = soil%rho(1:ntotal+bntotal)*
     &                             soil%vof(1:ntotal+bntotal)      

      soil%mass(1:ntotal+bntotal) = soil%vol(1:ntotal+bntotal)    !*0.3
     &                            * soil%rho(1:ntotal+bntotal)


      if(dry_soil)then
         parts%ntotal = ntotal_nozzle   ! no water in the tank
         parts%nvirt  = 0
      endif
c--------------------------------------------------------------------

      return
      end subroutine


