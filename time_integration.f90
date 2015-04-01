      subroutine time_integration
!----------------------------------------------------------------------
      use param
      use declarations_sph
      implicit none     
!
      integer :: i, j, k, d, ntotal, it
      type(particles), pointer :: pl
              
      do it = 1, maxtimestep   
  
        itimestep = itimestep+1

        parts%itimestep = itimestep
        soil%itimestep  = itimestep

        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  
  
!--------------For water--------------------------------------------------
      
!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

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

!---  Definition of variables out of the function vector:    

        call single_step

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


        time = time + dt

!	if (itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0) then
!          call output
!	endif 

!        if (mod(itimestep,print_step).eq.0) then
!          write(*,*)
!          write(*,101)'x','velocity', 'dvx'    
!          write(*,100)x(1,moni_particle), vx(1,moni_particle), 
!     &                dvx(1,moni_particle)    
!        endif
        
!101     format(1x,3(2x,a12))	 
!100     format(1x,3(2x,e13.6))

      enddo

      return

! -------------------
      contains
! -------------------

      subroutine first_half
! ---------------------------------------------------------------------
      implicit none

      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            pl%rho_min(i) = pl%rho(i)
            pl%rho(i) = pl%rho(i) +(dt/2.)* pl%drho(i)
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
! -------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if (.not.summation_density ) then
            pl%rho(i) = pl%rho(i) + (dt/2.)* pl%drho(i)
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
            pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   &
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
         enddo           
         !if(trim(pl%imaterial)=='soil')write(*,*) pl%dvx(1:2,i) 
      enddo 
      
      return
      end subroutine

      subroutine second_half
! ------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if (.not.summation_density ) then 
            pl%rho(i) = pl%rho_min(i) + dt*pl%drho(i)
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
            pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)   &
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)                  
         enddo
         
         !if(trim(pl%imaterial)=='soil')write(*,*) pl%vx(1:2,i) 
      enddo

      return
      end subroutine

      subroutine return_mapping
! ------------------------------------------------------------------
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

!----------------------------------------------------------------------      
      subroutine time_integration_for_water
!----------------------------------------------------------------------
use param
use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it
type(particles), pointer :: pl
              
do it = 1, maxtimestep    
   itimestep = itimestep+1
   parts%itimestep = itimestep

!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

   if(itimestep .ne. 1)then
      pl => parts
      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            pl%rho_min(i) = pl%rho(i)
            pl%rho(i) = pl%rho(i) +(dt/2.)* pl%drho(i)
         endif
          
         if(pl%itype(i)<0)cycle

         do d = 1, dim
            pl%v_min(d, i) = pl%vx(d, i)
            pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
         enddo
      enddo
           
   endif

   call single_step_for_water

   if(itimestep .eq. 1) then
      pl => parts
      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if(.not.summation_density )then
            pl%rho(i) = pl%rho(i) + (dt/2.)* pl%drho(i)
         endif

         if(pl%itype(i)<0)cycle

         do d = 1, dim        
            pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   & 
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
         enddo           

      enddo 
           
   else

       pl => parts
       do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if(.not.summation_density )then 
            pl%rho(i) = pl%rho_min(i) + dt*pl%drho(i)
         endif

         if(pl%itype(i)<0)cycle

         do d = 1, dim                   
            pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)    & 
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)                  
         enddo
         
      enddo

   endif 

   time = time + dt

   if(mod(itimestep,print_step).eq.0)then
      write(*,*)'______________________________________________'
      write(*,*)'  current number of time step =',                &
                itimestep,'     current time=', real(time+dt)
      write(*,*)'______________________________________________'
   endif  

   if(itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0)then
      call output
   endif 

enddo

return
end subroutine      

!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_verlet
!----------------------------------------------------------------------
use param
use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
type(particles), pointer :: pl
              
real(dp), allocatable, dimension(:,:) :: lastvx   
real(dp), allocatable, dimension(:,:) :: temp1    
real(dp), allocatable, dimension(:)   :: lastrho  
real(dp), allocatable, dimension(:)   :: temp2    

allocate(lastvx(2,parts%maxn))
allocate(temp1(2,parts%maxn))
allocate(lastrho(parts%maxn))
allocate(temp2(parts%maxn))
    pl => parts
    lastvx = pl%vx
    lastrho = pl%rho

do it = 1, maxtimestep 
    itimestep = itimestep+1
   call single_step_for_water

if(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal
       do d = 1, dim
          pl%x(d,i) = pl%x(d,i) + dt * pl%vx(d,i) +(dt**2./2.)*pl%dvx(d,i)
          lastvx(d,i) = lastvx(d,i) + 2.*dt*pl%dvx(d,i)
       
       enddo
    enddo
else
    do i = 1, pl%ntotal
       do d = 1, dim
          pl%x(d,i) = pl%x(d,i) + dt * pl%vx(d,i) +(dt**2./2.)*pl%dvx(d,i)
          lastvx(d,i) = pl%vx(d,i) + dt*pl%dvx(d,i)
      enddo
    enddo
endif
temp1 = pl%vx
pl%vx = lastvx
lastvx = temp1

!if(mod(itimestep,50) .ne. 0) then
!    do i = 1, pl%ntotal+pl%nvirt
!          lastrho(i) = pl%rho(i) + dt * 2. * pl%drho(i)
!    enddo
!else
   do i = 1, pl%ntotal+pl%nvirt
          lastrho(i) = pl%rho(i) + dt* pl%drho(i)
   enddo
!endif
temp2 = pl%rho
pl%rho = lastrho
lastrho = temp2

   time = time + dt

   if(mod(itimestep,print_step).eq.0)then
      write(*,*)'______________________________________________'
      write(*,*)'  current number of time step =',                &
                itimestep,'     current time=', real(time+dt)
      write(*,*)'______________________________________________'
   endif  

   if(itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0)then
      call output
   endif 

enddo

deallocate(lastvx)
deallocate(lastrho)
deallocate(temp1)
deallocate(temp2)


return
end subroutine      
         

!----------------------------------------------------------------------
             subroutine time_integration_for_soil
!----------------------------------------------------------------------
use param
use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it
type(particles), pointer :: pl
              
do it = 1, maxtimestep     
   itimestep = itimestep+1
   parts%itimestep = itimestep
  
!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

        if (itimestep .ne. 1) then
           pl => parts
           call first_half
        
              call drucker_prager_failure_criterion(pl)


        endif

!---  Definition of variables out of the function vector:    

        call single_step

        if (itimestep .eq. 1) then
           pl => parts
           call first_step

              if(plasticity==2)then
                 call plastic_or_not(pl)
                 call return_mapping
              endif
              call drucker_prager_failure_criterion(pl)

        else
   
           pl => parts
           call second_half          
           !if(nor_density) call norm_density(pl)

              if(plasticity==2)then
                 call plastic_or_not(pl)
                 call return_mapping
              endif
              call drucker_prager_failure_criterion(pl)


        endif 

        time = time + dt

        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  
        

!	if (itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0) then
!          call output
!	endif 

!        if (mod(itimestep,print_step).eq.0) then
!          write(*,*)
!          write(*,101)'x','velocity', 'dvx'    
!          write(*,100)x(1,moni_particle), vx(1,moni_particle), 
!     &                dvx(1,moni_particle)    
!        endif
        
!101     format(1x,3(2x,a12))	 
!100     format(1x,3(2x,e13.6))

      enddo

      return

! -------------------
      contains
! -------------------

      subroutine first_half
! ---------------------------------------------------------------------
      implicit none

      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            pl%rho_min(i) = pl%rho(i)
            pl%rho(i) = pl%rho(i) +(dt/2.)* pl%drho(i)
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
! -------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if (.not.summation_density ) then
            pl%rho(i) = pl%rho(i) + (dt/2.)* pl%drho(i)
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
            pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   &
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
         enddo           
         !if(trim(pl%imaterial)=='soil')write(*,*) pl%dvx(1:2,i) 
      enddo 
      
      return
      end subroutine

      subroutine second_half
! ------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if (.not.summation_density ) then 
            pl%rho(i) = pl%rho_min(i) + dt*pl%drho(i)
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
            pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)   &
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)                  
         enddo
         
         !if(trim(pl%imaterial)=='soil')write(*,*) pl%vx(1:2,i) 
      enddo

      return
      end subroutine

      subroutine return_mapping
! ------------------------------------------------------------------
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
