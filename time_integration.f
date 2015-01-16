      subroutine time_integration
c----------------------------------------------------------------------
      use param
      use declarations_sph
      implicit none     
c  
      double precision, pointer, dimension(:) :: rho, p, u, drho, du
      double precision, pointer, dimension(:,:) :: x,vx,av,dvx 
c
      integer :: i, j, k, d, ntotal, it
      double precision :: temp_rho, temp_u
      type(particles), pointer :: pl

      ntotal = parts%ntotal
      x     => parts%x
      vx    => parts%vx
      rho   => parts%rho
      p     => parts%p
      u     => parts%u
      drho  => parts%drho
      dvx   => parts%dvx
      du    => parts%du
      av    => parts%av
              
      do it = 1, maxtimestep   
  
        itimestep = itimestep+1

        parts%itimestep = itimestep
        soil%itimestep  = itimestep

        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  
  
!--------------For water--------------------------------------------------
      
c     If not first time step, then update thermal energy, density and 
c     velocity half a time step  

        if (itimestep .ne. 1) then
           pl => parts
           call first_half
        
           !if(nor_density) call norm_density(pl)

!!DEC$IF(.FALSE.)
!!!Added on 3.Aug.2014
           if(trim(pl%imaterial)=='soil')then
              call drucker_prager_failure_criterion(pl)
           endif

                             ! For soil
                             if(.not.single_phase)then
                                pl => soil
                                call first_half
                              call drucker_prager_failure_criterion(pl)
                                !call volume_fraction_soil(soil)
                                !call volume_fraction_water(parts, soil)
                             endif  

!!DEC$ENDIF
        endif

c---  Definition of variables out of the function vector:    

        call single_step

        do i = 1, ntotal
           if(parts%zone(i)==1.and.x(2,i)>-wasserjet%immerse)then
              drho(i) = 0.
              du(i)   = 0.
              dvx(:,i)  = 0.
              parts%dvof(i) = 0.
           endif
        enddo

        if (itimestep .eq. 1) then
           pl => parts
           call first_step
           !if(nor_density) call norm_density(pl)

!!DEC$IF(.FALSE.)
!!!Added on 3.Aug.2014
           if(trim(pl%imaterial)=='soil')then
              if(plasticity==2)then
                 call plastic_or_not(pl)
                 call return_mapping
              endif
              call drucker_prager_failure_criterion(pl)
           endif

                                if(.not.single_phase)then
                                    pl => soil
                                    call first_step

                                   if(plasticity==2)then
                                      call plastic_or_not(pl)
                                      call return_mapping
                                   endif

                               call drucker_prager_failure_criterion(pl)
                                 !call volume_fraction_soil(soil)
                                 !call volume_fraction_water(parts, soil)
                                endif
!!DEC$ENDIF
        else
   
           pl => parts
           call second_half          
           !if(nor_density) call norm_density(pl)

!!DEC$IF(.FALSE.)
!!!Added on 3.Aug.2014
           if(trim(pl%imaterial)=='soil')then
              if(plasticity==2)then
                 call plastic_or_not(pl)
                 call return_mapping
              endif
              call drucker_prager_failure_criterion(pl)
           endif

                                if(.not.single_phase)then
                                   pl => soil
                                   call second_half

                                   if(plasticity==2)then
                                      call plastic_or_not(pl)
                                      call return_mapping
                                   endif

                               call drucker_prager_failure_criterion(pl)
                                 !call volume_fraction_soil(soil)
                                 !call volume_fraction_water(parts, soil)
                                endif

!!DEC$ENDIF
        endif 


! For saturated flow

!        if(.not.single_phase)then
!           call volume_fraction_soil(soil)
!           call volume_fraction_water(parts, soil)
!        endif        

!--------------------Velocity Inlet-----------------------

!DEC$IF(.FALSE.)
        k = 0
        do i = 1,ntotal
           if(parts%zone(i)==1.and.x(2,i)<-wasserjet%immerse)then
              parts%zone(i) = 0
              k = k + 1
              x(1,ntotal+k) = x(1,i)
              x(2,ntotal+k) = x(2,i)+parts%hsml(i)
              parts%vol(ntotal+k) = parts%vol(i)
              parts%hsml(ntotal+k) = parts%hsml(i)
              parts%itype(ntotal+k) = 2
              vx(1,ntotal+k) = 0.
              vx(2,ntotal+k) = inlet_velocity
              parts%p(ntotal+k) = 0.
              rho(ntotal+k) = 1000.
              parts%mass(ntotal+k) = parts%vol(ntotal+k)*rho(ntotal+k)
              parts%zone(ntotal+k) = 1

              parts%dvx(:,ntotal+k) = 0.
              parts%du(ntotal+k)  = 0.
              parts%drho(ntotal+k)= 0.
           endif
        enddo
        parts%ntotal = parts%ntotal + k
        ntotal = parts%ntotal
!---------------------------------------------------------

!DEC$ENDIF
        time = time + dt

!	if (itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0) then
!          call output
!	endif 

        if (mod(itimestep,print_step).eq.0) then
          write(*,*)
          write(*,101)'x','velocity', 'dvx'    
          write(*,100)x(1,moni_particle), vx(1,moni_particle), 
     &                dvx(1,moni_particle)    
        endif
        
101     format(1x,3(2x,a12))	 
100     format(1x,3(2x,e13.6))

      enddo

      return

c -------------------
      contains
c -------------------

      subroutine first_half
c ---------------------------------------------------------------------
      implicit none

      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
         pl%u_min(i) = pl%u(i)
         temp_u=0.
         !if (dim.eq.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)
         pl%u(i) = pl%u(i) + (dt/2.)* (pl%du(i)+temp_u)
         if(pl%u(i).lt.0)  pl%u(i) = 0.                 
            
         if (.not.summation_density) then    
            pl%rho_min(i) = pl%rho(i)
            temp_rho=0.
            !if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
            pl%rho(i) = pl%rho(i) +(dt/2.)*( pl%drho(i)+ temp_rho)
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof_min(i) = pl%vof(i)
            pl%vof(i) = pl%vof(i)+(dt/2.)*pl%dvof(i) 
         endif
 
         if(trim(pl%imaterial)=='soil')then
                  if(stress_integration==1)then
         !pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)   ! Nothing to do!
         !pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
         !pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
                  elseif(stress_integration==2)then
         pl%sxx_min(i) = pl%sxx(i)
         pl%sxy_min(i) = pl%sxy(i)
         pl%syy_min(i) = pl%syy(i)
         pl%p_min(i)   = pl%p(i)
         pl%sxx(i) = pl%sxx(i) + (dt/2.)*pl%dsxx(i)
         pl%sxy(i) = pl%sxy(i) + (dt/2.)*pl%dsxy(i)
         pl%syy(i) = pl%syy(i) + (dt/2.)*pl%dsyy(i)
         pl%p(i) = pl%p(i) + (dt/2.)*pl%dp(i)      !!!simultaneous pressure
                  endif
         endif
          
         if(pl%itype(i)<0)cycle
         !if(pl%itype(i)<0.and.pl%x(2,i)<-0.2)cycle
         do d = 1, dim
         !if(pl%itype(i)<0.and.d==1)cycle
            pl%v_min(d, i) = pl%vx(d, i)
            pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
         enddo
      enddo
 
      return
      end subroutine

      subroutine first_step
c -------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         temp_u=0.
         !if (dim.eq.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)        
         pl%u(i) = pl%u(i) + (dt/2.)*(pl%du(i) + temp_u)
         if(pl%u(i).lt.0)  pl%u(i) = 0.             
         
         if (.not.summation_density ) then
            temp_rho=0.
            !if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
            pl%rho(i) = pl%rho(i) + (dt/2.)* (pl%drho(i)+temp_rho)
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof(i) = pl%vof(i)+(dt/2.)*pl%dvof(i)
         endif

         if(trim(pl%imaterial)=='soil')then
                  if(stress_integration==1)then
         pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)
         pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
         pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
         pl%p(i)   = pl%p(i)   + dt*pl%dp(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%sxx(i) = pl%sxx(i) + (dt/2.)*pl%dsxx(i)
         pl%sxy(i) = pl%sxy(i) + (dt/2.)*pl%dsxy(i)
         pl%syy(i) = pl%syy(i) + (dt/2.)*pl%dsyy(i)
         pl%p(i)   = pl%p(i) + (dt/2.)*pl%dp(i)    !!! simultaneous pressure
                  endif
         endif
        
         if(pl%itype(i)<0)cycle
         !if(pl%itype(i)<0.and.pl%x(2,i)<-0.2)cycle
         do d = 1, dim        
         !if(pl%itype(i)<0.and.d==1)cycle
            pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i) 
     &                  + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
         enddo           
         !if(trim(pl%imaterial)=='soil')write(*,*) pl%dvx(1:2,i) 
      enddo 
      
      return
      end subroutine

      subroutine second_half
c ------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
         temp_u=0.
         !if (dim.eq.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)                       
         pl%u(i) = pl%u_min(i) + dt*(pl%du(i)+temp_u)
         if(pl%u(i).lt.0)  pl%u(i) = 0.          
            
         if (.not.summation_density ) then 
            temp_rho=0.
            !if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)        	           
            pl%rho(i) = pl%rho_min(i) + dt*(pl%drho(i)+temp_rho)
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof(i) = pl%vof_min(i)+dt*pl%dvof(i) 
         endif
                 
         if(trim(pl%imaterial)=='soil')then 
                  if(stress_integration==1)then 
         pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)
         pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
         pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
         pl%p(i)   = pl%p(i)   + dt*pl%dp(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%sxx(i) = pl%sxx_min(i) + dt*pl%dsxx(i)
         pl%sxy(i) = pl%sxy_min(i) + dt*pl%dsxy(i)
         pl%syy(i) = pl%syy_min(i) + dt*pl%dsyy(i)
         pl%p(i)   = pl%p_min(i) + dt*pl%dp(i)     !!! simultaneous pressure
                  endif
         endif

         if(pl%itype(i)<0)cycle
         !if(pl%itype(i)<0.and.pl%x(2,i)<-0.2)cycle
         do d = 1, dim                   
         !if(pl%itype(i)<0.and.d==1)cycle
            pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i) 
     &                  + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)                  
         enddo
         
         !if(trim(pl%imaterial)=='soil')write(*,*) pl%vx(1:2,i) 
      enddo

      return
      end subroutine

      subroutine return_mapping
c ------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal+pl%nvirt  ! origionally pl%ntotal            
                 
         !if(trim(pl%imaterial)=='soil')then
         if(pl%fail(i)==1)then 
                  if(stress_integration==1)then 
         pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx2(i)
         pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy2(i)
         pl%syy(i) = pl%syy(i) + dt*pl%dsyy2(i)
         pl%p(i)   = pl%p(i)   + dt*pl%dp2(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%sxx(i) = pl%sxx_min(i) + dt*pl%dsxx(i)
         pl%sxy(i) = pl%sxy_min(i) + dt*pl%dsxy(i)
         pl%syy(i) = pl%syy_min(i) + dt*pl%dsyy(i)
         pl%p(i)   = pl%p_min(i) + dt*pl%dp(i)     !!! simultaneous pressure
                  endif
         endif
         !endif
         
      enddo

      return
      end subroutine

      end subroutine
