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
      double precision xleft,xright,ybottom
      logical :: dbg = .false.
      double precision element_size

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
      element_size = 0.01

!     &   call wasserjet%set(dr=0.02d0,L=0.44d0,H=0.22d0,mdr=4,dj=2.d0,    !Koshi
!     &        mdj=400,ds=0.d0, mL=88,mH=44)

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
      call tank%cell_center       
      tank%x = tank%x  - wasserjet%L/2.
      tank%y = (tank%y + wasserjet%ds)*(-1.0)

! -----------------------Phase I------------------------------
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
! Select real particles in tank

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

      do i = 1, tank%m*tank%n
      if(tank%x(i)<xleft.or.tank%x(i)>xright.or.tank%y(i)<ybottom)cycle ! origional
!      if(tank%x(i)<xleft.or.tank%x(i)>-0.1.or.tank%y(i)<ybottom)then    ! water dam
            k = k + 1
            parts%x(1,k) = tank%x(i)
            parts%x(2,k) = tank%y(i)
      enddo
      parts%ntotal = k

! Boundary particles

      select case (trim(particle_size))
         case ('dambreak')
            xleft = -9999.0; xright = 9999.0
         case ('damsmall')
            xleft = -9999.0d0; xright = 9999.0d0
      end select

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

!      parts%nvirt = 0
     
! Basic settings for particles
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

!      do i = 1, parts%ntotal
!         if(parts%x(1,i)<-0.9.or.parts%x(1,i)>0.9.or.parts%x(2,i)<-0.9)
!     &      parts%itype(i) = -2
!      enddo

!     Set initial conditions of particles

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
            parts%sxx(i) = 0.1905*parts%p(i)
            parts%syy(i) = -0.381*parts%p(i)
            parts%p(i)   = 0.619*parts%p(i)
         enddo
      endif

      call initial_density(parts)

      parts%vof(1:ntotal_nozzle) = 1.0
      do i = ntotal_nozzle+1, ntotal+bntotal
         parts%vof(i) = 1
         !if(parts%x(2,i)>-0.16)parts%vof(i) = 1           !!! 0.16
         !parts%vof(ntotal_nozzle+1:ntotal+bntotal) =0.7  !!! origional
         if(parts%x(2,i)<-0.16.and.parts%x(1,i)<-0.18)parts%vof(i) = 0.7
!         if(parts%x(2,i)<-0.16.and.parts%x(1,i)>-0.18      &
!           .and.parts%itype(i)<0)parts%vof(i) = 0.7
      enddo
      if(single_phase) parts%vof(ntotal_nozzle+1:ntotal+bntotal) = 1.0
      if(.not.volume_fraction)                            &
        parts%vof(ntotal_nozzle+1:ntotal+bntotal) = 1.0

      parts%rho(1:ntotal+bntotal) = parts%rho(1:ntotal+bntotal)*    &
                                   parts%vof(1:ntotal+bntotal)

      parts%mass(1:ntotal+bntotal) = parts%vol(1:ntotal+bntotal)*   &
                                    parts%rho(1:ntotal+bntotal)

!      parts%mass(1:ntotal_nozzle) = 
!     &parts%vol(1:ntotal_nozzle)*parts%rho(1:ntotal_nozzle)
     
!      parts%mass(ntotal_nozzle+1:ntotal+bntotal) = 
!     &parts%vol(ntotal_nozzle+1:ntotal+bntotal)*0.7*
!     &parts%rho(ntotal_nozzle+1:ntotal+bntotal)

!-------------------------------soil input---------------------------------

      if(single_phase)return

      k = 0; ntotal = 0; bntotal = 0
      do i = 1, parts%ntotal+parts%nvirt
         if(parts%x(2,i)>-wasserjet%ds)cycle
         if(parts%x(2,i)>-0.16) cycle                    !!!0.16
         if(parts%x(1,i)>-0.18.and.parts%itype(i)>0)cycle
         k = k + 1
         soil%x(:,k)   = parts%x(:,i)
         soil%vol(k)   = parts%vol(i)
         soil%hsml(k)  = parts%hsml(i) 
!
         soil%p(k)     = parts%p(i)
!
         if(parts%itype(i)==-2)then
            bntotal = bntotal + 1
            soil%itype(k) = -3
         elseif(parts%itype(i)==2)then
            ntotal = ntotal + 1
            soil%itype(k) =  3
         endif
      enddo
      soil%ntotal = ntotal; soil%nvirt = bntotal

! set initial conditions for soil particles

      property => soil%material
      do i = 1, soil%ntotal+soil%nvirt
         soil%vx(:,i) = 0.
         soil%rho(i)  = property%rho0
         if(dry_soil)then
         soil%p(i) = soil%rho(i)*gravity*(soil%x(2,i)+wasserjet%ds)
         else
!         soil%p(i) = soil%rho(i)*gravity*(soil%x(2,i)+wasserjet%ds)  
!     &             - soil%p(i)                                        ! origional
         soil%p(i) = (soil%rho(i)-1000)*gravity*     &
                    (soil%x(2,i)+wasserjet%ds+0.06)                   !!!0.16
         endif      

         !if(soil%x(1,i)>-0.1) soil%p(i) = 0.d0                     ! water dam
         !if(soil%x(1,i)>-0.1.or.soil%itype(i)<0) soil%p(i) = 0.d0   ! water dam

      enddo


      do i = 1, soil%ntotal+soil%nvirt   
         !soil%sxx(i) = -0.3*soil%p(i)
         !soil%syy(i) = 0.d0
         soil%sxy(i) = 0.d0
            soil%sxx(i) = 0.1905*soil%p(i)
            soil%syy(i) = -0.381*soil%p(i)
            soil%p(i)   =  0.619*soil%p(i)
      enddo

      call initial_density(soil)       

      soil%vof(1:ntotal+bntotal) = 0.3  !0.1
      if(.not.volume_fraction)soil%vof(1:ntotal+bntotal) = 1.0 

      soil%rho(1:ntotal+bntotal) = soil%rho(1:ntotal+bntotal)*    &
                                  soil%vof(1:ntotal+bntotal)      

      soil%mass(1:ntotal+bntotal) = soil%vol(1:ntotal+bntotal)    &  !*0.3
                                 * soil%rho(1:ntotal+bntotal)


      if(dry_soil)then
         parts%ntotal = ntotal_nozzle   ! no water in the tank
         parts%nvirt  = 0
      endif
!--------------------------------------------------------------------

      return
      end subroutine


