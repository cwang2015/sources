!-------------------------------
     subroutine allocate_sph
!-------------------------------
use param
use declarations_sph
implicit none

!---------------------Real particles (water)---------------------
parts%dim  = dim
parts%maxn = maxn
parts%max_interaction = max_interaction

call allocate_particles(parts)

allocate(parts%zone(maxn)); parts%zone = 0
!allocate(parts%real_zone(parts%max_zone)); parts%real_zone = 0
!allocate(parts%virtual_zone(parts%max_zone)); parts%virtual_zone = 0

allocate(parts%w(max_interaction));             parts%w   = 0.d0
allocate(parts%dwdx(dim,max_interaction));      parts%dwdx= 0.d0
allocate(parts%pair_i(max_interaction));        parts%pair_i = 0
allocate(parts%pair_j(max_interaction));        parts%pair_j = 0
allocate(parts%countiac(maxn));                 parts%countiac=0
allocate(parts%grid(parts%maxngx,parts%maxngy,parts%maxngz)); parts%grid = 0
allocate(parts%xgcell(dim,maxn)); parts%xgcell = 0
allocate(parts%celldata(maxn)); parts%celldata = 0

! Fields variables

call allocate_particles_fields(parts)

! Working arrays outside particles type

allocate(txx(maxn), tyy(maxn), txy(maxn))
txx = 0.d0; tyy = 0.d0; txy = 0.d0

parts%txx => txx; parts%tyy => tyy; parts%txy => txy

if(trim(parts%imaterial)=='water')then
   parts%sxx => txx; parts%syy => tyy; parts%sxy => txy
elseif(trim(parts%imaterial)=='soil')then
   parts%wxy => txy
endif   

! Material assignment

if(trim(parts%imaterial)=='water') parts%material => H2O
if(trim(parts%imaterial)=='soil')  parts%material => SiO2

! Numerical parameters
parts%numeric => numeric


!--------------- Soil particles --------------------------------

if(single_phase) return

soil%dim = dim
soil%maxn= maxn
soil%max_interaction = max_interaction

call allocate_particles(soil)

allocate(soil%zone(maxn)); soil%zone = 0
!allocate(soil%real_zone(soil%max_zone)); soil%real_zone = 0
!allocate(soil%virtual_zone(soil%max_zone)); soil%virtual_zone = 0

soil%w      => parts%w
soil%dwdx   => parts%dwdx
soil%pair_i => parts%pair_i
soil%pair_j => parts%pair_j
allocate(soil%countiac(maxn));    soil%countiac=0
soil%grid   => parts%grid; soil%xgcell => parts%xgcell
soil%celldata => parts%celldata

! Fields variables
call allocate_particles_fields(soil)

soil%txx => txx; soil%tyy => tyy; soil%txy => txy
soil%wxy => txy

soil%material => SiO2
soil%numeric => numeric

return
end subroutine

!-------------------------------------------
     subroutine allocate_particles(parts)
!-------------------------------------------
use m_particles
implicit none

type(particles) parts
integer dim,maxn,max_interaction

dim  = parts%dim
maxn = parts%maxn
max_interaction = parts%max_interaction 

if(.not.associated(parts%itype)) allocate(parts%itype(maxn))
if(.not.associated(parts%x))     allocate(parts%x(dim,maxn))
if(.not.associated(parts%vol))   allocate(parts%vol(maxn))
if(.not.associated(parts%mass))  allocate(parts%mass(maxn))
if(.not.associated(parts%hsml))  allocate(parts%hsml(maxn))

return
end subroutine

!---------------------------------------------------------
     subroutine allocate_particles_fields(parts)
!---------------------------------------------------------
! Fields variables
!---------------------------------------------------------
use m_particles
implicit none

type(particles) parts
integer dim,maxn

dim  = parts%dim
maxn = parts%maxn

allocate(parts%vx(dim,maxn));  parts%vx  = 0.d0
allocate(parts%rho(maxn));     parts%rho = 0.d0
allocate(parts%p(maxn));       parts%p   = 0.d0
allocate(parts%u(maxn));       parts%u   = 0.d0
allocate(parts%c(maxn));       parts%c   = 0.d0
!allocate(parts%eta(maxn));     parts%eta = 0.d0 

! Old value
allocate(parts%v_min(dim,maxn));   parts%v_min   = 0.d0
allocate(parts%rho_min(maxn));     parts%rho_min = 0.d0
allocate(parts%u_min(maxn));       parts%u_min   = 0.d0

! Accelerations
allocate(parts%drho(maxn));    parts%drho= 0.d0  
allocate(parts%dvx(dim,maxn)); parts%dvx = 0.d0  
allocate(parts%du(maxn));      parts%du  = 0.d0 
allocate(parts%av(dim,maxn));  parts%av  = 0.d0

! Volume Fraction
allocate(parts%vof(maxn))
parts%vof = 0.d0

if(trim(parts%imaterial)=='water')then
   allocate(parts%vof2(maxn))          ! phi_f = 1- phi_s
   parts%vof2 = 0.d0
   allocate(parts%dvof(maxn))
   parts%dvof = 0.d0
   allocate(parts%vof_min(maxn))
   parts%vof_min = 0.d0
endif

if(trim(parts%imaterial)=='soil')then
   allocate(parts%sxx(maxn), parts%syy(maxn), parts%sxy(maxn))
   parts%sxx = 0.d0; parts%sxy = 0.d0; parts%syy = 0.d0
   allocate(parts%sxx_min(maxn))
   allocate(parts%syy_min(maxn))
   allocate(parts%sxy_min(maxn))
   allocate(parts%p_min(maxn))
   allocate(parts%vcc(maxn))    
   allocate(parts%szz(maxn))    !!! Intel Fortran Compiler is shit!!!
   parts%sxx_min = 0.d0; parts%sxy_min = 0.d0; parts%syy_min = 0.d0
   parts%p_min = 0.d0
   parts%vcc = 0.d0
   allocate(parts%dsxx(maxn),parts%dsyy(maxn),parts%dsxy(maxn))
   parts%dsxx = 0.d0; parts%dsxy = 0.d0; parts%dsyy = 0.d0
   allocate(parts%dp(maxn))
   parts%dp(maxn) = 0.d0
! For return mapping algorithm
   allocate(parts%dsxx2(maxn),parts%dsyy2(maxn),parts%dsxy2(maxn))
   parts%dsxx2 = 0.d0; parts%dsxy2 = 0.d0; parts%dsyy2 = 0.d0
   allocate(parts%dp2(maxn))
   parts%dp2(maxn) = 0.d0

   allocate(parts%fail(maxn)); parts%fail = 0

   allocate(parts%epsilon_p(maxn))
   parts%epsilon_p = 0.d0
endif

allocate(parts%drhodx(dim,maxn))
parts%drhodx = 0.d0

return
end subroutine
