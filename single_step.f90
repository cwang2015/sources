!----------------------------------------------------------------------      
                   subroutine single_step
!----------------------------------------------------------------------
!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
use param 
use declarations_sph
!use m_sph_fo
implicit none

integer  nphase, iphase
type(particles), pointer :: pl
logical :: dbg = .false.
integer i, ntotal
                 
if(dbg) write(*,*) 'In single_step...'

                                      nphase = 2
                                      if(single_phase) nphase=1

                                      do iphase = 1, nphase
                                         if(iphase==1) pl => parts  !!!
                                         if(iphase==2) pl => soil   !!!

pl%dvx = 0.d0; pl%drho = 0.d0
if(trim(pl%imaterial)=='water')pl%dvof = 0.d0
if(trim(pl%imaterial)=='soil')then
   pl%dsxx = 0.d0; pl%dsxy = 0.d0; pl%dsyy = 0.d0; pl%dp = 0.d0
endif
 
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
if (pl%numeric%nnps.eq.1) then 
   call direct_find(pl)
else if (pl%numeric%nnps.eq.2) then
   call link_list(pl)     
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
!        call link_list(itimestep, parts%ntotal+parts%nvirt,
!     &       parts%hsml(1),parts%x,parts%niac,parts%pair_i,
!     &       parts%pair_j,parts%w,parts%dwdx,parts%countiac)
else if (pl%numeric%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
endif         

if(mod(itimestep,print_step).eq.0.and.int_stat) then
   call pl%interaction_statistics
endif   

!--- Added by Wang
!if(nor_density) call norm_density(pl)

!---  Density approximation or change rate
     
!if(summation_density)then      
!if(mod(itimestep,30)==0) call sum_density(pl)
!else             
if(pl%nthreads==1)then
    call con_density(pl)   
        write(*,*),"asdf"
    stop
else
    call con_density_omp(pl)

endif
!endif
      
if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
    if(pl%nthreads==1)then
      call delta_sph(pl,pl%rho,pl%drho)
    else
      call delta_sph_omp(pl,pl%rho,pl%drho)
    endif  
   !endif
endif

!---  Dynamic viscosity:

if(trim(pl%imaterial)=='water')then
   !if (visc) call viscosity(pl)
elseif(trim(pl%imaterial)=='soil')then
   !call shear_modulus(pl)
endif
       
!---  Internal forces:
if(pl%nthreads==1)then
 call shear_strain_rate(pl) 
else
 call shear_strain_rate_omp(pl) 
endif
if(trim(pl%imaterial)=='soil')call velocity_divergence(pl)

call pressure(pl)

if(trim(pl%imaterial)=='water')then
   call newtonian_fluid(pl)
elseif(trim(pl%imaterial)=='soil')then

   if(yield_criterion == 1)then
      if(pl%nthreads==1)then    
       call mohr_coulomb_failure_criterion(pl)
      else 
       call mohr_coulomb_failure_criterion_omp(pl)
      endif 
   elseif(yield_criterion == 2)then
      if(pl%nthreads==1)then    
       call drucker_prager_failure_criterion(pl)
      else
       call drucker_prager_failure_criterion_omp(pl)
      endif 
   endif

   if(mod(itimestep,print_step).eq.0)    &
            write(*,*) 'Failured points: ', pl%nfail
endif

!Calculate internal force for water phase !! -phi_f Grad(p)
if(pl%imaterial=='water')then
   if(pl%nthreads==1)then
   pl%dvx(1,:) = -pl%vof*pl%df(pl%p,'x') + pl%df(pl%vof*pl%sxx,'x') + pl%df(pl%vof*pl%sxy,'y')
   pl%dvx(2,:) = -pl%vof*pl%df(pl%p,'y') + pl%df(pl%vof*pl%sxy,'x') + pl%df(pl%vof*pl%syy,'y')
   else
   pl%dvx(1,:) = -pl%vof*pl%df_omp(pl%p,'x') + pl%df_omp(pl%vof*pl%sxx,'x') + pl%df_omp(pl%vof*pl%sxy,'y')
   pl%dvx(2,:) = -pl%vof*pl%df_omp(pl%p,'y') + pl%df_omp(pl%vof*pl%sxy,'x') + pl%df_omp(pl%vof*pl%syy,'y')           
   endif        

   where (pl%rho.gt.0.0) pl%dvx(1,:) = pl%dvx(1,:)/pl%rho
   where (pl%rho.gt.0.0) pl%dvx(2,:) = pl%dvx(2,:)/pl%rho

else      
   !call int_force(pl)
!   if(pl%nthreads==1)then
   pl%dvx(1,:) = -pl%df3(pl%vof*pl%p,'x') + pl%df3(pl%vof*pl%sxx,'x') + pl%df3(pl%vof*pl%sxy,'y')
   pl%dvx(2,:) = -pl%df3(pl%vof*pl%p,'y') + pl%df3(pl%vof*pl%sxy,'x') + pl%df3(pl%vof*pl%syy,'y')
 !  else
 !  pl%dvx(1,:) = -pl%df3_omp(pl%vof*pl%p,'x') + pl%df3_omp(pl%vof*pl%sxx,'x') + pl%df3_omp(pl%vof*pl%sxy,'y')
 !  pl%dvx(2,:) = -pl%df3_omp(pl%vof*pl%p,'y') + pl%df3_omp(pl%vof*pl%sxy,'x') + pl%df3_omp(pl%vof*pl%syy,'y')
!   endif
   where (pl%rho.gt.0.0) pl%dvx(1,:) = pl%dvx(1,:)/pl%rho
   where (pl%rho.gt.0.0) pl%dvx(2,:) = pl%dvx(2,:)/pl%rho 
endif      

!if(trim(pl%imaterial)=='water'.and.water_tension_instability==2) &
!   call tension_instability(pl) 

! --- Plasticity flow rule   ! This was done before Jaummann_rate, because we 
!     we need txx,tyy,tzz, which was destroyed in Jaumann_rate!

if(trim(pl%imaterial)=='soil')then   
   if(plasticity==1)then 
      call plastic_flow_rule(pl)
   elseif(plasticity==2)then     
      call plastic_flow_rule2(pl)
   elseif(plasticity==3)then
    if(pl%nthreads==1)then
      call plastic_or_not(pl)
      call plastic_flow_rule3(pl)
    else
      call plastic_or_not_omp(pl)
      call plastic_flow_rule3_omp(pl)
    endif  
   endif
endif

! --- Jaumann rate  !When???

if(trim(pl%imaterial)=='soil')then
   if(pl%nthreads==1)then
     call Jaumann_rate(pl)
   else
     call Jaumann_rate_omp(pl)
   endif 
endif
!---  Artificial viscosity:

if (visc_artificial) then
   if(pl%nthreads==1)then
     call pl%art_visc
   else
     call pl%art_visc_omp
   endif  
endif   
if(trim(pl%imaterial)=='soil'.and.soil_artificial_stress)then
        !call art_stress(pl)
   if(pl%nthreads==1)then
   call pl%delta_sph(pl%p,pl%dp)
   call pl%delta_sph(pl%sxx,pl%dsxx)
   call pl%delta_sph(pl%sxy,pl%dsxy)
   call pl%delta_sph(pl%syy,pl%dsyy)
   else
   call pl%delta_sph_omp(pl%p,pl%dp)
   call pl%delta_sph_omp(pl%sxx,pl%dsxx)
   call pl%delta_sph_omp(pl%sxy,pl%dsxy)
   call pl%delta_sph_omp(pl%syy,pl%dsyy)
   endif
endif        
if(trim(pl%imaterial)=='water'.and.water_artificial_volume)then  
        !call art_volume_fraction_water2(pl)
   if(pl%nthreads==1)then
        call pl%delta_sph(pl%vof,pl%dvof)
   else
        call pl%delta_sph_omp(pl%vof,pl%dvof)
    endif
endif
!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
   if(pl%nthreads==1)then
       call pl%repulsive_force
   else
       call pl%repulsive_force_omp
   endif  
!      endif

pl%dvx(2,:) = pl%dvx(2,:) + gravity

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (average_velocity) then
   if(pl%nthreads==1)then
       call av_vel(pl)
   else    
       call av_vel_omp(pl) 
   endif
endif
!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

                                            enddo ! iphase

!call drag_force(parts,soil)   !!! Porous media


                if(.not.single_phase)then

!-------------------Water/soil interaction-------------------------------

!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
      !if (numeric%nnps.eq.1) then 
         call direct_find_2(parts,soil)
        !call direct_find(parts)
      !else if (numeric%nnps.eq.2) then
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      !else if (numeric%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      !endif         

      if(mod(itimestep,print_step).eq.0.and.int_stat) then
         call parts%interaction_statistics
      endif    

      call darcy_law(parts,soil)          
      call pore_water_pressure(parts,soil) 

      if(volume_fraction)then
         call volume_fraction_soil(soil)
         call volume_fraction_water2(parts,soil)
         call volume_fraction_water(parts,soil)  ! phi_f = 1- phi_s
         if(volume_fraction_renorm)then
            if(mod(itimestep,40).eq.0) then
               ntotal = parts%ntotal+parts%nvirt
               parts%rho(1:ntotal) = parts%rho(1:ntotal)/parts%vof(1:ntotal) 
               parts%vof = parts%vof2
               parts%rho(1:ntotal) = parts%rho(1:ntotal)*parts%vof(1:ntotal) 
            endif
         endif
      endif

!      call direct_find(parts) 
!      call art_volume_fraction_water2(parts)
     
             endif ! .not.single_phase

      if(itimestep>=save_step_from.and.   &
         mod(itimestep,save_step).eq.0)then
         call output
      endif 

return
end subroutine


!-------------------------------------------------
      subroutine single_step_for_water
!-------------------------------------------------
use param
use declarations_sph
!use m_sph_fo
implicit none
type(particles), pointer :: pl
type(material),pointer :: property

pl => parts
property => pl%material

pl%dvx = 0.d0; pl%drho = 0.d0

!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length

call pl%find_pairs

if(mod(itimestep,print_step).eq.0.and.int_stat) then
   call pl%interaction_statistics
endif   

!---  Density approximation or change rate
     
!if(summation_density) call sum_density(pl)
call sum_density(pl)
pl%drho = -pl%rho*(pl%df2(pl%vx(1,:),'x')+pl%df2(pl%vx(2,:),'y'))
      
if(artificial_density)then
   !call renormalize_density_gradient(pl)
   !call art_density(pl)
   call delta_sph(pl,pl%rho,pl%drho)
endif
       
!---  Internal forces:

!Calculate pressure

where(pl%rho>0.0) pl%p = property%b*((pl%rho/property%rho0)**property%gamma-1)

!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

pl%txx = 2./3.*(2.0*pl%df(pl%vx(1,:),'x')-pl%df(pl%vx(2,:),'y'))
pl%txy = pl%df(pl%vx(1,:),'y')+pl%df(pl%vx(2,:),'x')
pl%tyy = 2./3.*(2.0*pl%df(pl%vx(2,:),'y')-pl%df(pl%vx(1,:),'x'))

!Newtonian fluid

pl%sxx = property%viscosity*pl%txx
pl%syy = property%viscosity*pl%tyy
pl%sxy = property%viscosity*pl%txy

!Calculate internal force

pl%dvx(1,:) = - pl%df(pl%p,'x') + pl%df(pl%sxx,'x') + pl%df(pl%sxy,'y')
pl%dvx(2,:) = - pl%df(pl%p,'y') + pl%df(pl%sxy,'x') + pl%df(pl%syy,'y')

where (pl%rho.gt.0.0) pl%dvx(1,:) = pl%dvx(1,:)/pl%rho
where (pl%rho.gt.0.0) pl%dvx(2,:) = pl%dvx(2,:)/pl%rho
       
!if(water_tension_instability==2) call tension_instability(pl) 

!---  Artificial viscosity:

if (visc_artificial) call pl%art_visc
    
!---  External forces:

pl%dvx(2,:) = pl%dvx(2,:) + gravity

!call repulsive_force(pl)

! Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)
     
! Calculating average velocity of each partile for avoiding penetration

if (average_velocity) call av_vel(pl) 

if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

return
end subroutine


!----------------------------------------------------------------------      
               subroutine single_step_for_soil
!----------------------------------------------------------------------
!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
use param 
use declarations_sph
!use m_sph_fo
implicit none

type(particles), pointer :: pl
type(material), pointer :: property

pl => parts
property => pl%material

pl%dvx = 0.d0; pl%drho = 0.d0
pl%dsxx = 0.d0; pl%dsxy = 0.d0; pl%dsyy = 0.d0
pl%dp = 0.d0
 
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
call pl%find_pairs

if(mod(itimestep,print_step).eq.0.and.int_stat) then
   call pl%interaction_statistics
endif   

!---  Density approximation or change rate
     
if (summation_density) then      
    call sum_density(pl)
else             
    call con_density(pl)         
endif
      
if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      call delta_sph(pl,pl%rho,pl%drho)
   !endif
endif

!call shear_modulus(pl)
       
!---  Internal forces:

call shear_strain_rate(pl)   
call velocity_divergence(pl)

call pressure(pl)

if(yield_criterion == 1)then
   call mohr_coulomb_failure_criterion(pl)
elseif(yield_criterion == 2)then
   call drucker_prager_failure_criterion(pl)
endif

!call int_force1(pl)

pl%dvx(1,:) = - pl%df(pl%p,'x') + pl%df(pl%sxx,'x') + pl%df(pl%sxy,'y')
pl%dvx(2,:) = - pl%df(pl%p,'y') + pl%df(pl%sxy,'x') + pl%df(pl%syy,'y')

where (pl%rho.gt.0.0) pl%dvx(1,:) = pl%dvx(1,:)/pl%rho
where (pl%rho.gt.0.0) pl%dvx(2,:) = pl%dvx(2,:)/pl%rho       

! --- Plasticity flow rule   ! This was done before Jaummann_rate, because we 
!     we need txx,tyy,tzz, which was destroyed in Jaumann_rate!

if(plasticity==1)then 
   call plastic_flow_rule(pl)
elseif(plasticity==2)then     
   call plastic_flow_rule2(pl)
elseif(plasticity==3)then
   call plastic_or_not(pl)
   call plastic_flow_rule3(pl)
endif

! --- Jaumann rate  !When???

call Jaumann_rate(pl)

!---  Artificial viscosity:

if (visc_artificial) call pl%art_visc

!if(soil_artificial_stress) call art_stress(pl)
if(soil_artificial_stress)then
   call pl%delta_sph(pl%p,pl%dp)
   call pl%delta_sph(pl%sxx,pl%dsxx)
   call pl%delta_sph(pl%sxy,pl%dsxy)
   call pl%delta_sph(pl%syy,pl%dsyy)
endif   

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
!      endif

pl%dvx(2,:) = pl%dvx(2,:) + gravity

! Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)
     
! Calculating average velocity of each partile for avoiding penetration

if (average_velocity) call av_vel(pl) 

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

return
end subroutine



