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
      
      integer  ntotal, bntotal
      integer i, j, d, k, ntotal_nozzle
      type(block) nozzle, tank
      type(material), pointer :: property
      double precision xleft,xright,ybottom,soil_surface
      integer :: fluid_zone(10) = 0, fluid_boundary_zone(10) = 0
      integer ::  soil_zone(10) = 0,  soil_boundary_zone(10) = 0
      logical :: dbg = .true.
      double precision element_size, soil_submerged_depth

      if(dbg) write(*,*) 'In water_jet...'      

!     Set nozzle and tank geometry parameters

select case(trim(particle_size))
   case('small') 
      call wasserjet%set(dr=0.2d0,L=2.d0,H=0.6d0,dj=10.d0,ds=0.4d0)
      element_size = 0.02
   case('large')
      call wasserjet%set(dr=0.2d0,L=2.d0,H=0.6d0,dj=10.d0,ds=0.4d0)
      element_size = 0.04
   case('smallnozzle')
      call wasserjet%set(dr=0.1d0,L=2.d0,H=0.6d0,dj=10.d0,ds=0.4d0)
      element_size = 0.02
   case('bui')
!     &   call wasserjet%set(dr=0.02d0,L=0.5d0,H=0.2d0,mdr=4,dj=2.d0,
!     &        mdj=400,ds=0.1d0, mL=100,mH=40)
      call wasserjet%set(dr=0.02d0,L=0.5d0,H=0.2d0,dj=2.d0,ds=0.1d0)
      element_size = 0.005  ! 0.01
   case('koshi')
      call wasserjet%set(dr=0.02d0,L=0.44d0,H=0.22d0,dj=2.d0,ds=0.d0)
      element_size = 0.01   ! 0.005
   case('bigbox')
      call wasserjet%set(dr=0.02d0,L=1.d0,H=0.2d0,dj=2.d0,ds=0.1d0)
      element_size = 0.01
   case('dambreak')
      call wasserjet%set(L=0.8d0,H=0.2d0, ds=0.1d0)
      element_size = 0.005
   case('damsmall')
!     &   call wasserjet%set(L=0.8d0,H=0.06d0, ds=0.24d0, mL=320,mH=24)  !smallest  4:1
!     &   call wasserjet%set(L=0.8d0,H=0.12d0, ds=0.18d0, mL=320,mH=48)  !smallest  2:1
      call wasserjet%set(L=1.2d0,H=0.11d0, ds=0.d0)  !coarse alumibar  2:1
      element_size = 0.005
!     &   call wasserjet%set(L=1.2d0,H=0.11d0, ds=0.d0, mL=4800,mH=44)  !fine alumibar  2:1
        
!     &   call wasserjet%set(L=0.8d0,H=0.06d0, ds=0.24d0, mL=160,mH=12)  !coarse     4:1
!     &   call wasserjet%set(L=0.8d0,H=0.03d0, ds=0.27d0, mL=320,mH=12)  !smallest  8:1
!     &   call wasserjet%set(L=0.8d0,H=0.11d0, ds=0.19d0, mL=160,mH=22) !cube       1:1
end select

!     Nozzle particle entrance

      call nozzle%set(xl = wasserjet%dr, yl = wasserjet%dj)
      nozzle%m = nozzle%xl/element_size; nozzle%n = nozzle%yl/element_size
      k = nozzle%m*nozzle%n
      allocate(nozzle%x(k),nozzle%y(k))
      call nozzle%cell_center
      nozzle%x = nozzle%x - wasserjet%dr/2.
      nozzle%y = nozzle%y - wasserjet%immerse

!     Paticles in tank
 
      call tank%set(xl = wasserjet%L, yl = wasserjet%H)
      tank%m = tank%xl/element_size; tank%n = tank%yl/element_size
      k = tank%m*tank%n
      allocate(tank%x(k),tank%y(k))
      allocate(tank%zone(k))
      call tank%cell_center       
      tank%x = tank%x  - wasserjet%L/2.
      tank%y = (tank%y + wasserjet%ds)*(-1.0)

! Zoning...

      xleft = -0.9; xright = 0.9; ybottom = -0.9
      select case (trim(particle_size))
         case ('bui')
         xleft = -0.225; xright = 0.225; ybottom = -0.275
      !   xleft = -0.2 ;  xright = 0.2; ybottom = -0.2     !Koshi
         case ('bigbox')
         xleft = -0.475; xright = 0.475; ybottom = -0.275
         case ('dambreak')
         xleft = -0.1; xright = 0.1; ybottom = -0.275
         case ('damsmall')
!         xleft = -0.05; xright = 0.05; ybottom = -0.29   !cube
!         xleft = -0.1; xright = 0.1; ybottom = -0.295      ! smallest
          xleft = -0.2; xright = 0.2; ybottom = -0.1
      end select

                         select case (trim(operating_mode))
                             case('bui')   ! U-tube
      tank%zone = 0
      do i = 1, tank%m*tank%n
         if(tank%x(i)<0.23.and.tank%x(i)>-0.23.and.tank%y(i)>-0.28) tank%zone(i) = 1
         if(tank%x(i)<0.19.and.tank%x(i)>-0.19.and.tank%y(i)>-0.24) tank%zone(i) = 2
         if(tank%x(i)<0.17.and.tank%x(i)>-0.17.and.tank%y(i)>-0.22) tank%zone(i) = 3
         if(0.19<tank%x(i).and.tank%x(i)<0.23.and.tank%y(i)>-0.20) tank%zone(i) = 4
      enddo
                             case('koshi')
      tank%zone = 0
      soil_submerged_depth = 0.05
      do i = 1, tank%m*tank%n
         if(tank%x(i)<-0.2.or.tank%x(i)>0.2.or.tank%y(i)<-0.2) tank%zone(i) = -1
         if(tank%zone(i)==-1.and.tank%y(i)<-0.05) tank%zone(i) = -2
         if(tank%zone(i)==-2.and.tank%x(i)>-0.1) tank%zone(i) = -3
         if(tank%x(i)>-0.2.and.tank%x(i)<0.2.and.tank%y(i)>-0.2) tank%zone(i) = 1
         if(tank%zone(i)==1.and.tank%x(i)<-0.1.and.tank%y(i)<-0.05) tank%zone(i) = 2
      enddo
      parts%nreal_zone = 2
      parts%nvirtual_zone = 3
      parts%real_zone(1:2)=(/1,2/)
      parts%virtual_zone(1:3)=(/-1,-2,-3/)

      soil%nreal_zone = 1
      soil%nvirtual_zone = 2
      soil%real_zone(1) = 2
      soil%virtual_zone(1:2)=(/-2,-3/)


                             case('static_submerged')
                                     write(*,*) 'case here'
      tank%zone = 0
      soil_submerged_depth = 0.05
      do i = 1, tank%m*tank%n
         if(tank%x(i)<-0.2.or.tank%x(i)>0.2.or.tank%y(i)<-0.2) tank%zone(i) = -1
         if(tank%zone(i)==-1.and.tank%y(i)<-0.05) tank%zone(i) = -2
         if(tank%x(i)>-0.2.and.tank%x(i)<0.2.and.tank%y(i)>-0.2) tank%zone(i) = 1
         if(tank%zone(i)==1.and.tank%y(i)<-0.05) tank%zone(i) = 2
      enddo
      parts%nreal_zone = 2
      parts%nvirtual_zone = 2
      parts%real_zone(1:2) = (/1,2/)
      parts%virtual_zone(1:2) = (/-1,-2/)
     
      soil%nreal_zone = 1
      soil%nvirtual_zone = 1
      soil%real_zone(1) = 2
      soil%virtual_zone(1) = -2

                             case('static_mixture')
                                     write(*,*) 'case: static_mixture'
      tank%zone = 0
      soil_submerged_depth = 0.d0
      do i = 1, tank%m*tank%n
         if(tank%x(i)<-0.2.or.tank%x(i)>0.2.or.tank%y(i)<-0.2) tank%zone(i) = -1
         if(tank%x(i)>-0.2.and.tank%x(i)<0.2.and.tank%y(i)>-0.2) tank%zone(i) = 1
      enddo
      parts%nreal_zone = 1
      parts%nvirtual_zone = 1
      parts%real_zone(1) = 1
      parts%virtual_zone(1) = -1
     
      soil%nreal_zone = 1
      soil%nvirtual_zone = 1
      soil%real_zone(1) = 1
      soil%virtual_zone(1) = -1

                         end select

! -----------------------Phase I------------------------------
      k = 0
      parts%ntotal = 0; parts%nvirt = 0

      ntotal_nozzle = 0
                                          if(with_nozzle)then 
      ntotal_nozzle = nozzle%m*nozzle%n
      do i =1, ntotal_nozzle
         k = k + 1
         parts%x(1, k) = nozzle%x(i)
         parts%x(2, k) = nozzle%y(i)
      enddo
                                          endif !nozzle

      call parts%take_real_points(tank)
      call parts%take_virtual_points(tank)

!      parts%nvirt = 0
     
! Basic settings for particles (vol,hsml,itype)

      parts%vol = tank%dx*tank%dy
      parts%hsml = tank%dx
      call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

      do i = 1,parts%ntotal+parts%nvirt
         parts%vx(1,i) = 0.
         parts%vx(2,i) = 0.
         if(parts%x(2,i)>-wasserjet%ds)then
            parts%vx(2,i) = inlet_velocity
         endif
      enddo

! ...Pressure

      property => parts%material
      do i = 1,parts%ntotal+parts%nvirt
         parts%p(i) = property%rho0*gravity*(parts%x(2,i)+wasserjet%ds)  !0.
         if(parts%x(2,i)>-wasserjet%ds)then
            parts%p(i) = 0
         endif
      enddo      

!      if(trim(parts%imaterial)=='soil')then
!         do i = 1, parts%ntotal+parts%nvirt   
!            !parts%sxx(i) = -0.3*parts%p(i)
!            !parts%syy(i) = 0.d0
!            parts%sxy(i) = 0.d0
!            parts%sxx(i) = 0. !0.1905*parts%p(i)
!            parts%syy(i) = 0. !-0.381*parts%p(i)
!            parts%p(i)   = 0. !0.619*parts%p(i)
!         enddo
!      endif

! ...Voulume fraction

if(trim(operating_mode)=='koshi'.or.trim(operating_mode)=='static_submerged')then
      do i = 1, parts%ntotal+parts%nvirt
         parts%vof(i) = 0.3
         if(parts%zone(i)==-1.or.parts%zone(i)==1)parts%vof(i) = 1.0
      enddo
elseif(trim(operating_mode)=='static_mixture')then
      parts%vof = 0.3
endif

      if(single_phase) parts%vof = 1.0
      if(.not.volume_fraction)  parts%vof = 1.0

      call initial_density(parts)
      parts%rho = parts%rho * parts%vof
      parts%mass = parts%vol * parts%rho     

!-------------------------------soil input---------------------------------

      if(single_phase)return
   
      call soil%take_real_points(tank)
      call soil%take_virtual_points(tank)
      soil%vol = tank%dx*tank%dy
      soil%hsml= tank%dx
      call soil%setup_itype

! set initial conditions for soil particles

      property => soil%material
      do i = 1, soil%ntotal+soil%nvirt
         soil%vx(:,i) = 0.
         soil%rho(i)  = property%rho0
         if(dry_soil)then
         soil%p(i) = soil%rho(i)*gravity*(soil%x(2,i)+wasserjet%ds)
         else
!         soil%p(i) = soil%rho(i)*gravity*(soil%x(2,i)+wasserjet%ds)   & 
!                   - soil%p(i)                                        ! origional
         soil%p(i) = (soil%rho(i)-1000)*gravity*     &
                    (soil%x(2,i)+wasserjet%ds+soil_submerged_depth)    !!!0.16
         endif      

         if(trim(operating_mode)=='koshi')then
            if(soil%x(1,i)>-0.1) soil%p(i) = 0.d0  
         endif

      enddo


      do i = 1, soil%ntotal+soil%nvirt   
         !soil%sxx(i) = -0.3*soil%p(i)
         !soil%syy(i) = 0.d0
         soil%sxy(i) = 0.d0
            soil%sxx(i) = 0.1905*soil%p(i)
            soil%syy(i) = -0.381*soil%p(i)
            soil%p(i)   =  0.619*soil%p(i)
      enddo

      soil%vof = 0.7
      if(.not.volume_fraction)soil%vof = 1.0
      soil%rho  = soil%rho*soil%vof
      soil%mass = soil%vol*soil%rho


      if(dry_soil)then
         parts%ntotal = ntotal_nozzle   ! no water in the tank
         parts%nvirt  = 0
      endif
!--------------------------------------------------------------------

      return
      end subroutine


