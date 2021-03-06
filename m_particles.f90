!---------------------------
    module m_particles
!---------------------------
implicit none
integer, parameter :: dp = kind(0.d0)

type block             
     real(dp) xl, yl
     integer  m,  n
     real(dp) dx, dy
     real(dp), pointer, dimension(:) :: x,y
     integer,  pointer, dimension(:) :: zone
     type(particles), pointer :: p1=>null(), p2=>null()
   contains
     procedure :: set => set_block_sub
     procedure :: cell_center => get_cell_center_sub
end type

type material

! For water

   real(dp) rho0,b,gamma,c,viscosity

! For soil

   real(dp) k,porosity,permeability,G,E,niu,cohesion,phi
   
end type

type numerical

!     Nearest neighbor particle searching (nnps) method
!     nnps = 1 : Simplest and direct searching
!            2 : Sorting grid linked list
!            3 : Tree algorithm
integer :: nnps = 1

! Gravitational acceleration
  real(dp) :: gravity = -9.8

! Artificial viscosity
! alpha: shear viscosity; beta: bulk viscosity; etq: parameter to avoid sigularities   
   real(dp) :: alpha=0.1d0, beta=0.d0, etq=0.1d0

! Leonard_Johns repulsive force
   real(dp) :: dd = 0.1d0, p1 = 12, p2 = 4

! Delta-SPH
   real(dp) :: delta = 0.1d0

! Velocity average
  real(dp) :: epsilon = 0.001
   
end type

type array
   integer :: ndim1
   real(dp), pointer, dimension(:) :: r  => null()
   type(array), pointer :: x => null(), y => null(), z => null()
   type(array), pointer :: xy => null(), xz => null(), yz => null()
   type(particles), pointer :: parts => null()
   contains
       final :: array_final
end type

! Particles in SPH method
type particles

!Physics parameter
!real(dp) :: gravity = -9.8

!Geometry parameters

!dim : Dimension of the problem (1, 2 or 3)
integer :: dim = 2        
        
   integer :: niac  = 0
   integer :: ntotal = 0, nvirt = 0
   real(dp) dspp  ! initial particle interval  

!Parameters for the computational geometry,  
!x_maxgeom : Upper limit of allowed x-regime 
!x_mingeom : Lower limit of allowed x-regime 
!y_maxgeom : Upper limit of allowed y-regime 
!y_mingeom : Lower limit of allowed y-regime 
!z_maxgeom : Upper limit of allowed z-regime 
!z_mingeom : Lower limit of allowed z-regime 
real(dp) :: x_maxgeom = 0.3e0, x_mingeom = -0.05e0,  &
            y_maxgeom = 0.15e0, y_mingeom = -0.05e0,  &
            z_maxgeom = 10.e0, z_mingeom = -10.e0

!Parameter used for sorting grid cells in the link list algorithm
!maxngx  : Maximum number of sorting grid cells in x-direction
!maxngy  : Maximum number of sorting grid cells in y-direction
!maxngz  : Maximum number of sorting grid cells in z-direction
!Determining maximum number of sorting grid cells:
!(For an homogeneous particle distribution:)
!1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
!2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
!3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)
integer :: maxngx = 100,maxngy = 100, maxngz = 1

integer, pointer, dimension(:,:,:) :: grid => null()
integer, pointer, dimension(:,:) :: xgcell => null()
integer, pointer, dimension(:) :: celldata => null() 
integer  ngridx(3),ghsmlx(3)
real(dp) mingridx(3),maxgridx(3),dgeomx(3)

!Memory allocation

!maxn: Maximum number of particles
!max_interation : Maximum number of interaction pairs
integer :: maxn = 12000, max_interaction = 10 * 12000
  
!SPH algorithm

!Particle approximation (pa_sph)
!pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
!         2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
integer :: pa_sph = 2 

!Nearest neighbor particle searching (nnps) method
!nnps = 1 : Simplest and direct searching
!       2 : Sorting grid linked list
!       3 : Tree algorithm
integer :: nnps = 1 

!Smoothing length evolution (sle) algorithm
!sle = 0 : Keep unchanged,
!      1 : h = fac * (m/rho)^(1/dim)
!      2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
!      3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) ) 
integer :: sle = 0 

!Smoothing kernel function 
!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!    = 2, Gauss kernel   (Gingold and Monaghan 1981) 
!    = 3, Quintic kernel (Morris 1997)
!    = 4, Wendland    
integer :: skf = 4 

!Switches for different senarios

!summation_density = .TRUE. : Use density summation model in the code, 
!                    .FALSE.: Use continuiity equation
!average_velocity = .TRUE. : Monaghan treatment on average velocity,
!                   .FALSE.: No average treatment.
!config_input = .TRUE. : Load initial configuration data,
!               .FALSE.: Generate initial configuration.
!virtual_part = .TRUE. : Use vritual particle,
!               .FALSE.: No use of vritual particle.
!vp_input = .TRUE. : Load virtual particle information,
!           .FALSE.: Generate virtual particle information.
!visc = .true. : Consider viscosity,
!       .false.: No viscosity.
!ex_force =.true. : Consider external force,
!          .false.: No external force.
!visc_artificial = .true. : Consider artificial viscosity,
!                  .false.: No considering of artificial viscosity.
!heat_artificial = .true. : Consider artificial heating,
!                  .false.: No considering of artificial heating.
!self_gravity = .true. : Considering self_gravity,
!               .false.: No considering of self_gravity
!nor_density =  .true. : Density normalization by using CSPM,
!               .false.: No normalization.

!integer :: integrate_scheme = 1  ! =1, LF; =2, Verlet
!logical :: summation_density  = .false.         
!logical :: average_velocity  = .true.         
!logical :: config_input  = .false. 
!logical :: virtual_part  = .true. 
!logical :: vp_input  = .false.  
!logical :: visc  = .true.  
!logical :: ex_force  = .true.
!logical :: visc_artificial  = .true. 
!logical :: heat_artificial  = .false. 
!logical :: self_gravity  = .true.      
!logical :: nor_density  = .false.              

!integer :: soil_pressure = 2  ! =1, eos; =2, mean trace
!integer :: stress_integration = 1
!integer :: yield_criterion = 2
!integer :: plasticity = 3  ! =0 non; =1 Bui, =2 return mapping =3 Lopez
!logical :: artificial_density = .true.                  
!logical :: soil_artificial_stress = .true.

!logical :: volume_fraction = .true.
!logical :: water_artificial_volume = .true.
!logical :: volume_fraction_renorm = .true.

! 0 ignor; 1 negative pressure to zero; 2 artficial stress
!integer :: water_tension_instability = 0

! Symmetry of the problem
! nsym = 0 : no symmetry,
!      = 1 : axis symmetry,
!      = 2 : center symmetry.     
integer :: nsym = 0

! Control parameters for output 
! int_stat = .true. : Print statistics about SPH particle interactions.
!                     including virtual particle information.
! print_step: Print Timestep (On Screen)
! save_step : Save Timestep    (To Disk File)
! moni_particle: The particle number for information monitoring.
logical :: int_stat = .true.
integer :: print_step, save_step, moni_particle = 264

!Recorde time interval
integer :: save_step_from = 0, save_step_to = 100

   character(len=32) :: imaterial
   class(*), pointer :: material => null()

! Numerical parameters
   class(numerical), pointer :: numeric => null()

! Particle interaction pair
   integer :: maxp = 0, minp = 0
   integer :: maxiac = 0, miniac = 0
   integer :: sumiac = 0, noiac  = 0
   integer, pointer, dimension(:)  :: pair_i => null()
   integer, pointer, dimension(:)  :: pair_j => null()
   integer, pointer, dimension(:)  :: countiac=>null()

! Kernel and its derivative
   real(dp), pointer, dimension(:)   :: w    => null()
   real(dp), pointer, dimension(:,:) :: dwdx => null()

! Particle fundamental data
   integer,  pointer, dimension(:)   :: itype=> null()
   real(dp), pointer, dimension(:,:) :: x    => null()
   real(dp), pointer, dimension(:)   :: vol  => null()
   real(dp), pointer, dimension(:)   :: mass => null()
   real(dp), pointer, dimension(:)   :: hsml => null()
   integer,  pointer, dimension(:)   :: zone => null()

! Volume of Fraction
   real(dp), pointer, dimension(:)   :: vof  => null(), vof2=>null() ! phi_f= 1-phi_s
   real(dp), pointer, dimension(:)   :: dvof  => null()
   real(dp), pointer, dimension(:)   :: vof_min  => null()

! Field variables
   real(dp), pointer, dimension(:)   :: rho  => null()
   real(dp), pointer, dimension(:,:) :: vx   => null()   
   real(dp), pointer, dimension(:)   :: p    => null()
   !real(dp), pointer, dimension(:)   :: eta  => null()   
   real(dp), pointer, dimension(:)   :: c    => null()   
   real(dp), pointer, dimension(:)   :: u    => null()

! Stress tensor 
   !type(tensor), pointer, dimension(:)   :: st   => null() 
   real(dp), pointer, dimension(:) :: sxx => null(), syy => null()
   real(dp), pointer, dimension(:) :: szz => null(), sxy => null()
   real(dp), pointer, dimension(:) :: sxz => null(), syz => null()
 
! Shear strain rate 
   real(dp), pointer, dimension(:) :: txx => null(), tyy => null()
   real(dp), pointer, dimension(:) :: tzz => null(), txy => null()
   real(dp), pointer, dimension(:) :: txz => null(), tyz => null()
! Bulk strain rate, i.e. divergency of velocity
   real(dp), pointer, dimension(:) :: vcc => null()

! Accumulative deviatoric strain
   real(dp), pointer, dimension(:) :: epsilon_p => null()

! Rotation tensor
   real(dp), pointer, dimension(:) :: wxy => null()

! Old values
   real(dp), pointer, dimension(:)   :: rho_min => null()
   real(dp), pointer, dimension(:,:) :: v_min   => null()
   real(dp), pointer, dimension(:)   :: u_min   => null()
   real(dp), pointer, dimension(:)   :: p_min   => null()
   real(dp), pointer, dimension(:) :: sxx_min => null(), syy_min => null()
   real(dp), pointer, dimension(:) :: szz_min => null(), sxy_min => null()
   real(dp), pointer, dimension(:) :: sxz_min => null(), syz_min => null()

! Acceleration
   real(dp), pointer, dimension(:)     :: drho => null()
   real(dp), pointer, dimension(:,:)   :: dvx  => null()
   real(dp), pointer, dimension(:)     :: du   => null()
   real(dp), pointer, dimension(:)     :: dp   => null()
   real(dp), pointer, dimension(:,:)   :: av   => null()
!   type(tensor), pointer, dimension(:) :: dst  => null()
   real(dp), pointer, dimension(:) :: dsxx => null(), dsyy => null()
   real(dp), pointer, dimension(:) :: dszz => null(), dsxy => null()
   real(dp), pointer, dimension(:) :: dsxz => null(), dsyz => null()

! For plasticity corrector step (return mapping algorithm)
   real(dp), pointer, dimension(:) :: dsxx2 => null(), dsyy2 => null()
   real(dp), pointer, dimension(:) :: dsxy2 => null(),   dp2 => null()

! Density gradient
   real(dp), pointer, dimension(:,:)     :: drhodx => null()

! For soil failure
   integer :: nfail = 0
   integer, pointer, dimension(:) :: fail => null()



! Boundry particles defined as type particles
   type(particles), pointer :: bor

! Multiphase flow
!   type(particles), pointer :: water => null(), soil => null()

   integer itimestep

   integer :: nthreads = 1
   integer, pointer, dimension(:) :: niac_start, niac_end

!   procedure(df),pointer :: diff => null() 

   contains

       procedure :: write_particles
       procedure :: write_component 
       procedure :: append_virtual_particles
       procedure :: depend_virtual_particles
       procedure :: interaction_statistics
       procedure :: minimum_time_step
       procedure :: take_real => take_real_points1
       procedure :: take_virtual => take_virtual_points1
       procedure :: setup_itype
       procedure :: get_scale_k
       procedure :: direct_find
       procedure :: init_grid
       procedure :: grid_geom
       procedure :: link_list
       procedure :: find_pairs
       procedure :: kernel
       procedure :: initial_density
       procedure :: pressure
       procedure :: pressure_water
       procedure :: pressure_soil
       procedure :: tension_instability_water 
       procedure :: sum_density
       procedure :: sum_density_omp
       procedure :: con_density
       procedure :: con_density_omp
       procedure :: newtonian_fluid
       procedure :: art_visc
       procedure :: art_visc_omp
       procedure :: delta_sph
       procedure :: delta_sph_omp
       procedure :: av_vel
       procedure :: av_vel_omp
       procedure :: repulsive_force
       procedure :: repulsive_force_omp
       procedure :: df
       procedure :: df_omp
       procedure :: df2
       procedure :: df2_omp
       procedure :: df3
       procedure :: df3_omp
       procedure :: velocity_divergence
       procedure :: velocity_divergence_omp
!       procedure :: time_integration_for_water
!       procedure :: find_particle_nearest2_point
!       procedure :: grad_scalar
!       procedure :: grad_tensor
!       generic :: grad => grad_scalar, grad_tensor
      procedure :: shear_strain_rate 
      procedure :: shear_strain_rate_omp
      procedure :: Jaumann_rate
      procedure :: Jaumann_rate_omp
      procedure :: mohr_coulomb_failure_criterion
      procedure :: mohr_coulomb_failure_criterion_omp
      procedure :: drucker_prager_failure_criterion
      procedure :: drucker_prager_failure_criterion_omp
      procedure :: plastic_or_not
      procedure :: plastic_or_not_omp
      procedure :: plastic_flow_rule
      procedure :: plastic_flow_rule2
      procedure :: plastic_flow_rule3
      procedure :: plastic_flow_rule3_omp

      procedure :: get_num_threads
      procedure :: get_niac_start_end

end type 

interface write_field
   module procedure :: write_scalar_field
   module procedure :: write_vector2d_field
end interface

interface operator(+)
   module procedure :: array_add_array 
end interface

interface operator(-)
   module procedure :: array_sub_array
   module procedure :: array_minus
end interface

interface operator(*)
   module procedure :: array_mul_array
   module procedure :: array_mul_real
   module procedure :: array_mul_double_real
end interface

!=======
contains
!=======

!DEC$IF(.FALSE.)
!---------------------------------------------------------------------------------
   subroutine set_geo_sub(this,dl,dr,dh,ds,dj,D,L,H,immerse,       &
                              mdl,mdr,mdh,mds,mdj,mD,mL,mH)
!---------------------------------------------------------------------------------
   implicit none
   class(geo) this
   real(dp), optional :: dl,dr,dh,ds,dj,D,L,H,immerse
   integer,  optional :: mdl,mdr,mdh,mds,mdj,mD,mL,mH

   if(present(dl)) this%dl = dl  
   if(present(dr)) this%dr = dr  
   if(present(dh)) this%dh = dh  
   if(present(ds)) this%ds = ds  
   if(present(dj)) this%dj = dj  
   if(present(D)) this%D = D  
   if(present(L)) this%L = L  
   if(present(H)) this%H = H  
   if(present(immerse)) this%immerse = immerse

   if(present(mdl)) this%mdl = mdl  
   if(present(mdr)) this%mdr = mdr  
   if(present(mdh)) this%mdh = mdh  
   if(present(mds)) this%mds = mds  
   if(present(mdj)) this%mdj = mdj  
   if(present(mD)) this%mD = mD  
   if(present(mL)) this%mL = mL  
   if(present(mH)) this%mH = mH  
  
   return
   end subroutine

!DEC$ENDIF
!--------------------------------------------
   subroutine set_block_sub(this,xl,yl,m,n)
!--------------------------------------------
implicit none
class(block) this
real(dp),optional :: xl,yl
integer, optional :: m, n

if(present(xl)) this%xl = xl
if(present(yl)) this%yl = yl
if(present(m))  this%m  = m
if(present(n))  this%n  = n

return
end subroutine

!------------------------------------------------
   subroutine get_cell_center_sub(this)
!------------------------------------------------
implicit none
class(block) this
integer i,j,k

this%dx = this%xl/this%m
this%dy = this%yl/this%n

k = 0
do i = 1, this%m
   do j = 1, this%n
      k = k + 1
      this%x(k) = (i-1)*this%dx + this%dx/2.
      this%y(k) = (j-1)*this%dy + this%dy/2.
   enddo
enddo

return
end subroutine

!-------------------------------------------
!   subroutine copy_particle_sub(this,that)
!-------------------------------------------
!implicit none
!class(particle) this
!class(particle)  that
!
!this%x = that%x
!
!return
!end subroutine

!----------------------------------------
   subroutine write_particles(this) 
!----------------------------------------
implicit none
class(particles) this
integer i

write(*,*) 'Real particles:'
write(*,*) 'Num.  itype  x   y   vol  hsml'
do i = 1, this%ntotal
   write(*,*) i,  this%itype(i), this%x(1,i), this%x(2,i), this%vol(i), this%hsml(i) 
end do

write(*,*) 'Virtual particles:'
write(*,*) 'Num.  itype  x   y   vol  hsml'
do i = this%ntotal+1, this%ntotal+this%nvirt
   write(*,*) i,  this%itype(i), this%x(1,i), this%x(2,i), this%vol(i), this%hsml(i) 
end do

return
end subroutine

!---------------------------------------------
   subroutine write_component(this,vars) 
!---------------------------------------------
implicit none
class(particles) this
character(len=*) vars
integer ntotal, nvirt, i

ntotal = this%ntotal
nvirt  = this%nvirt

select case (trim(vars))
   case ('x')
      write(*,*) 'Particle position:'
      do i = 1, ntotal+nvirt
         write(*,*) i,this%x(1,i),this%x(2,i)
      enddo
   case ('mass')
      write(*,*) 'Particle mass:'
      do i = 1, ntotal+nvirt
         write(*,*) i, this%mass(i)
      enddo
   case('w')
      write(*,*) 'Kernels:'
      write(*,*) (i,this%w(i), i = 1, this%niac)
   case('dwdx')
      write(*,*) 'Derivatives of Kernels:'
      write(*,*) (i,this%dwdx(1,i),this%dwdx(2,i),i = 1, this%niac)
   case ('pairs')
      write(*,*) 'Interaction pairs:'
      write(*,*) (i, this%pair_i(i), this%pair_j(i),'|',i=1,this%niac)
   case default
      write(*,*) 'Write this compoent is not implemented yet!'
end select

return
end subroutine

! Statistics for the interaction
!-------------------------------------------------
     subroutine interaction_statistics(parts)
!-------------------------------------------------
implicit none
class(particles) parts
integer ntotal, i,j,d
integer  sumiac, maxiac, miniac, noiac, maxp, minp
logical :: dbg = .false.

ntotal   =   parts%ntotal + parts%nvirt
sumiac = 0
maxiac = 0
miniac = 1000
noiac  = 0

do i=1,ntotal
   sumiac = sumiac + parts%countiac(i)
   if (parts%countiac(i).gt.maxiac) then
      maxiac = parts%countiac(i)
      maxp = i
   endif
   if (parts%countiac(i).lt.miniac) then 
      miniac = parts%countiac(i)
      minp = i
   endif
   if (parts%countiac(i).eq.0) noiac  = noiac + 1
enddo

parts%maxp = maxp;        parts%minp = minp  
parts%maxiac = maxiac;    parts%miniac = miniac
parts%sumiac = sumiac;    parts%noiac = noiac

print *,' >> Statistics: interactions per particle:'
print *,'**** Particle:',parts%maxp,                 &
        ' maximal interactions:', parts%maxiac
print *,'**** Particle:',parts%minp,                 &
        ' minimal interactions:', parts%miniac
print *,'**** Average :',real(parts%sumiac)/         &
                         real(parts%ntotal+parts%nvirt)
print *,'**** Total pairs : ',parts%niac
print *,'**** Particles with no interactions:',parts%noiac

return
end subroutine

!-------------------------------------------------
     subroutine minimum_time_step(parts)
!-------------------------------------------------
implicit none
class(particles) parts

type(material), pointer :: mat
integer ntotal
double precision c, miu 
double precision mint1, mint21, mint22, mint3

ntotal = parts%ntotal

mat => parts%material
c   = mat%c
miu = mat%viscosity/mat%rho0

! CFL condition
mint1 = 0.25*minval(parts%hsml(1:ntotal)/c)
! Acceleration 
mint21 = 0.25*minval(sqrt(parts%hsml(1:ntotal)/parts%dvx(1,1:ntotal)))
mint22 = 0.25*minval(sqrt(parts%hsml(1:ntotal)/parts%dvx(2,1:ntotal)))
! Viscous diffustion time
mint3 = 0.125*minval(parts%hsml(1:ntotal)**2/miu)

write(*,*) 'Minimum time step is: '
write(*,*) ' mint1 = ', mint1
write(*,*) ' mint21 = ', mint21, ' mint22 = ', mint22
write(*,*) ' mint3 = ', mint3

return
end subroutine

!-------------------------------------------------
     subroutine particle_monitor(this)
!-------------------------------------------------
implicit none
class(particles) this

!          write(*,*)
!          write(*,*) '**** Information for particle ****', 
!     &        		moni_particle         
!          write(*,101)'internal a ','artifical a=',
!     &         		'external a ','total a '   
!          write(*,100)indvxdt%rr(1,moni_particle),
!     &                ardvxdt%rr(1,moni_particle),
!     &                exdvxdt%rr(1,moni_particle),
!     &                parts%dvx(1,moni_particle)          
!      endif
!101   format(1x,4(2x,a12))      
!100   format(1x,4(2x,e13.6))      

return
end subroutine

!-------------------------------------------------
   subroutine write_scalar_field(field,parts)
!-------------------------------------------------
implicit none
double precision, dimension(:) :: field
type(particles) parts
integer i,ntotal

ntotal = parts%ntotal+parts%nvirt

write(*,*) (i, field(i),i=1,ntotal)

return
end subroutine

!-------------------------------------------------
   subroutine write_vector2d_field(field,parts)
!-------------------------------------------------
implicit none
double precision, dimension(:,:) :: field
type(particles) parts
integer i,ntotal

ntotal = parts%ntotal+parts%nvirt

write(*,*) (i, field(1,i),field(2,i),i=1,ntotal)

return
end subroutine

!---------------------------------------------------
    subroutine append_virtual_particles(this)
!---------------------------------------------------
implicit none
class(particles) this

type(particles), pointer :: bor
integer ntotal, i, k

bor    =>this%bor

! First, assembly all the real particles

!ntotal = this%ntotal
!k = 0 
!do i = 1, bor%ntotal
!   if(bor%itype(i)<0) cycle   ! real particles
!   k = k + 1
!   this%x(:,ntotal+k)  = bor%x(:,i)
!   this%hsml(ntotal+k) = bor%hsml(i)
!   this%vol(ntotal+k)  = bor%vol(i)
!   this%itype(ntotal+k)= bor%itype(i)
!   this%mass(ntotal+k) = bor%mass(i)
!enddo
!this%ntotal = this%ntotal + k

! Then assembly all the virtual particles

ntotal = this%ntotal
k = 0
do i = 1, bor%ntotal
   if(bor%itype(i)>0) cycle  ! virtual particles
   k = k + 1
   this%x(:,ntotal+k)  = bor%x(:,i)
   this%hsml(ntotal+k) = bor%hsml(i)
   this%itype(ntotal+k)= bor%itype(i)
   this%vol(ntotal+k)  = bor%vol(i)
   this%mass(ntotal+k) = bor%mass(i)
enddo
this%nvirt = k
!write(*,*) 'Appended virtual particle number: ', this%nvirt

return
end subroutine

!---------------------------------------------------
    subroutine depend_virtual_particles(this)
!---------------------------------------------------
implicit none
class(particles) this

!type(particles), pointer :: bor
!integer ntotal,nvirt, i,j,k

!ntotal = this%ntotal
!nvirt  = this%nvirt

!bor => this%bor

!k = 0
!do i = ntotal+1, ntotal+nvirt

!!   if(this%zone(i)==?) cycle      ! Symetery boundary virtural particles is skipped! 

!   j = i - ntotal
!   k = k + 1
!   bor%x(:,j)  = this%x(:,i)
!   bor%hsml(j) = this%hsml(i)
!   bor%vol(j)  = this%vol(i)
!   bor%zone(j) = this%zone(i)
!enddo
!this%nbor = k 

return
end subroutine

!-------------------------------------------------------
!    subroutine point_to_water(p, property)
!-------------------------------------------------------
!implicit none
!class(*), target :: property
!type(h2O), pointer :: p

!select type (property)
!   type is (h2O)
!      p => property
!end select
!
!return
!end subroutine


!DEC$IF(.FALSE.)
!-----------------------------------------------
      subroutine take_real_points(this,tank)
!-----------------------------------------------
implicit none

class(particles) this
type(block) tank

integer i,j,k

! Take real particles in tank

k = this%ntotal
do j = 1, this%nreal_zone
   do i = 1, tank%m*tank%n
      if(tank%zone(i)==this%real_zone(j))then             
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
enddo
this%ntotal = k

return
end subroutine

!DEC$ENDIF

!--------------------------------------------------
      subroutine take_real_points1(this,tank,zone)
!--------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take real particles in tank

k = this%ntotal

   do i = 1, tank%m*tank%n
      if(tank%zone(i)== zone)then             
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo

this%ntotal = k

return
end subroutine


!DEC$IF(.FALSE.)
!-----------------------------------------------------
      subroutine take_virtual_points(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank

integer i,j,k

! Take virtual particles in tank

k = this%ntotal
do j = 1, this%nvirtual_zone
   do i = 1, tank%m*tank%n
      if(tank%zone(i)==this%virtual_zone(j))then             
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
enddo
this%nvirt = k-this%ntotal

return
end subroutine

!DEC$ENDIF

!-----------------------------------------------------
      subroutine take_virtual_points1(this,tank,zone)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i)== zone)then             
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
this%nvirt = k-this%ntotal

return
end subroutine

!----------------------------------------
      subroutine setup_itype(parts)
!----------------------------------------
implicit none
class(particles) parts
integer ntotal, bntotal

ntotal = parts%ntotal; bntotal = parts%nvirt

if(trim(parts%imaterial)=='water')then
   parts%itype(1:ntotal) = 2
   parts%itype(ntotal+1:ntotal+bntotal) = -2
elseif(trim(parts%imaterial)=='soil')then
   parts%itype(1:ntotal) = 3
   parts%itype(ntotal+1:ntotal+bntotal) = -3
endif

return
end subroutine

!-------------------------------------------------
     function get_scale_k(this) result(scale_k)
!-------------------------------------------------
implicit none
class(particles) this
integer scale_k

select case (this%skf)
  case(1)
     scale_k = 2
  case(2)
     scale_k = 3
  case(3)
     scale_k = 3
  case(4)
     scale_k = 2
  case default
     stop 'get_scale_k: No such kernel function!'
end select 

end function

!----------------------------------------------------------------------
       subroutine direct_find(parts)
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by directly comparing the particle distance 
!   with the corresponding smoothing length.
!-----------------------------------------------------------------------
      implicit none
      
      class(particles) parts

      integer ntotal, niac,i,j,d,scale_k
      real(dp) dxiac(3), driac, r, mhsml, tdwdx(3)
      logical :: dbg = .false.

      if(dbg) write(*,*) 'In direct_find...'

      ntotal   =   parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
     
      do i=1,ntotal
        parts%countiac(i) = 0
      enddo

      niac = 0

      do i=1,ntotal-1     
        do j = i+1, ntotal
          dxiac(1) = parts%x(1,i) - parts%x(1,j)
          driac    = dxiac(1)*dxiac(1)
          do d=2,parts%dim
            dxiac(d) = parts%x(d,i) - parts%x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (parts%hsml(i)+parts%hsml(j))/2.
          if (sqrt(driac).lt.scale_k*mhsml) then
            if (niac.lt.parts%max_interaction) then    

!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 

              niac = niac + 1
              parts%pair_i(niac) = i
              parts%pair_j(niac) = j
              r = sqrt(driac)
              parts%countiac(i) = parts%countiac(i) + 1
              parts%countiac(j) = parts%countiac(j) + 1

!     Kernel and derivations of kernel
              call parts%kernel(r,dxiac,mhsml,parts%w(niac),tdwdx)
              do d=1,parts%dim
                parts%dwdx(d,niac) = tdwdx(d)
              enddo                                      
            else
              print *,  ' >>> ERROR <<< : Too many interactions' 
              stop
            endif
          endif
        enddo
      enddo  

      parts%niac = niac
 
      end subroutine

!----------------------------------------------------------------------
       subroutine direct_find_2(parts,part2)
!----------------------------------------------------------------------
      implicit none
      
      class(particles) parts, part2
      integer ntotal, niac,i,j,d, scale_k
      real(dp) dxiac(3), driac, r, mhsml, tdwdx(3)
      logical :: dbg = .false.

      if(dbg) write(*,*) 'In direct_find2...'

      ntotal   =   parts%ntotal + parts%nvirt
      scale_k  =   parts%get_scale_k()
     
      do i=1,ntotal
        parts%countiac(i) = 0
      enddo
      part2%countiac = 0

      niac = 0

      do i=1,ntotal     
        do j = 1, part2%ntotal+part2%nvirt
          dxiac(1) = parts%x(1,i) - part2%x(1,j)
          driac    = dxiac(1)*dxiac(1)
          do d=2,parts%dim
            dxiac(d) = parts%x(d,i) - part2%x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (parts%hsml(i)+part2%hsml(j))/2.
          if (sqrt(driac).lt.scale_k*mhsml) then
            if (niac.lt.parts%max_interaction) then    

!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 

              niac = niac + 1
              parts%pair_i(niac) = i
              parts%pair_j(niac) = j
              r = sqrt(driac)
              parts%countiac(i) = parts%countiac(i) + 1
              part2%countiac(j) = part2%countiac(j) + 1

!     Kernel and derivations of kernel
              call parts%kernel(r,dxiac,mhsml,parts%w(niac),tdwdx)
              do d=1,parts%dim
                parts%dwdx(d,niac) = tdwdx(d)
              enddo                                     
            else
              print *, ' >>> ERROR <<< : Too many interactions' 
              stop
            endif
          endif
        enddo
      enddo  

      parts%niac = niac

      end subroutine

!   Subroutine to established a pair linked list by sorting grid cell.
!   It is suitable for a homogeneous particle distribution with the 
!   same smoothing length in an instant. A fixed number of particles
!   lie in each cell.       
!----------------------------------------------------------------------      
      subroutine init_grid(parts)
!----------------------------------------------------------------------      
      implicit none
      class(particles) parts
      real(dp) hsml 
      integer i, j, k, d, dim, ntotal
!     Averaged number of particles per grid cell
      real(dp), parameter :: nppg = 3.e0

      hsml = parts%hsml(1)   !!! hsml is same for all particles
      dim = parts%dim; ntotal = parts%ntotal + parts%nvirt
      
!     Range of sorting grid

      parts%maxgridx(1) = parts%x_maxgeom;  parts%mingridx(1) = parts%x_mingeom
      parts%maxgridx(2) = parts%y_maxgeom;  parts%mingridx(2) = parts%y_mingeom
      parts%maxgridx(3) = parts%z_maxgeom;  parts%mingridx(3) = parts%z_mingeom
      parts%dgeomx = parts%maxgridx - parts%mingridx

!     Number of grid cells in x-, y- and z-direction:

      if (dim.eq.1) then
        parts%ngridx(1) = min(int(ntotal/nppg) + 1,parts%maxngx)
      else if (dim.eq.2) then
        parts%ngridx(1) = &
        min(int(sqrt(ntotal*parts%dgeomx(1)/(parts%dgeomx(2)*nppg))) + 1,parts%maxngx)
        parts%ngridx(2) = &
        min(int(parts%ngridx(1)*parts%dgeomx(2)/parts%dgeomx(1)) + 1,parts%maxngy)
      else if (dim.eq.3) then
        parts%ngridx(1) = min(int((ntotal*parts%dgeomx(1)*parts%dgeomx(1)/   &
           (parts%dgeomx(2)*parts%dgeomx(3)*nppg))**(1.e0/3.e0)) + 1,parts%maxngx)
        parts%ngridx(2) =  &
           min(int(parts%ngridx(1)*parts%dgeomx(2)/parts%dgeomx(1)) + 1,parts%maxngy)
        parts%ngridx(3) =  &
           min(int(parts%ngridx(1)*parts%dgeomx(3)/parts%dgeomx(1)) + 1,parts%maxngz)
      endif

!     Smoothing Length measured in grid cells:
      do d=1,dim
         parts%ghsmlx(d) = int(real(parts%ngridx(d))*hsml/parts%dgeomx(d)) + 1
      enddo

      parts%grid = 0   ! Initialize grid

      end subroutine

!   Subroutine to calculate the coordinates (xgcell) of the cell of 
!   the sorting  grid, in which the particle with coordinates (x) lies.      
!-----------------------------------------------------------------------      
      subroutine grid_geom(parts,i,x,xgcell)
!-----------------------------------------------------------------------
      implicit none
      class(particles) parts
      integer i, xgcell(3)
      real(dp) x(3)
      integer d

      do d=1,3
        xgcell(d) = 1
      enddo

      do d=1,parts%dim
        if ((x(d).gt.parts%maxgridx(d)).or.(x(d).lt.parts%mingridx(d))) then
          print *,' >>> ERROR <<< : Particle out of range'
          print *,'    Particle position: x(',i,d,') = ',x(d)
          print *,'    Range: [xmin,xmax](',D,') =                    &
                 [',parts%mingridx(d),',',parts%maxgridx(d),']'
          stop
        else
          xgcell(d) = int(real(parts%ngridx(d))/parts%dgeomx(d)*      &
                      (x(d)-parts%mingridx(d)) + 1.e0)
        endif
      enddo

      end subroutine

!----------------------------------------------------------------------      
      subroutine link_list(parts)
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by using a sorting grid linked list  
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer i, j, d, scale_k, ntotal, niac    
      integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)
      real(dp) hsml,dr,r,dx(3),tdwdx(3)

!      INTEGER nthreads,n_per_threads, niac_per_threads, it, n_start(8), n_end(8)
!      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

      ntotal  = parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
      hsml = parts%hsml(1)             !!! hsml is same for all particles!

      do i=1,ntotal
         parts%countiac(i) = 0
      enddo

      call parts%init_grid
      
!     Position particles on grid and create linked list:
      
      do i=1,ntotal
        call parts%grid_geom(i,parts%x(:,i),gcell)
        do d=1,parts%dim
          parts%xgcell(d,i) = gcell(d)
        enddo
        parts%celldata(i) = parts%grid(gcell(1),gcell(2),gcell(3))
        parts%grid(gcell(1),gcell(2),gcell(3)) = i
      enddo

!     Determine interaction parameters:
!!$omp parallel
!      nthreads = omp_get_num_threads()
!      write(*,*) 'ooo', nthreads
!!$omp end parallel       
!      n_per_threads = ntotal/nthreads
!      niac_per_threads = parts%max_interaction/threads
!      do it = 1, nthreads
!         n_start(it)=(it-1)*n_per*threads+1
!         n_end(it) = it*n_per_threads
!         niac_start(it)=(it-1)*niac_per_threads
!         niac_end(it)=niac_start(it)
!      enddo   

!!$omp parallel do private(niac,i,d)
!do it = 1, nthreads
!       niac = niac_start(it)
!       do i = n_start(it),n_end(it)

      niac = 0
      do i=1,ntotal-1

!     Determine range of grid to go through:
         
        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,parts%dim
          dnxgcell(d) = parts%xgcell(d,i) - parts%ghsmlx(d)
          dpxgcell(d) = parts%xgcell(d,i) + parts%ghsmlx(d)
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),parts%ngridx(d))
        enddo

!     Search grid:
      
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = parts%grid(xcell,ycell,zcell)
 1            if (j.gt.i) then
                dx(1) = parts%x(1,i) - parts%x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,parts%dim
                  dx(d) = parts%x(d,i) - parts%x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                  if (niac.lt.parts%max_interaction) then

!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 

                    niac = niac + 1
                    parts%pair_i(niac) = i
                    parts%pair_j(niac) = j
                    r = sqrt(dr)
                    parts%countiac(i) = parts%countiac(i) + 1
                    parts%countiac(j) = parts%countiac(j) + 1
                           
!--- Kernel and derivations of kernel

                    call parts%kernel(r,dx,hsml,parts%w(niac),tdwdx)  !!!!!!!!
	            do d = 1, parts%dim
	              parts%dwdx(d,niac)=tdwdx(d)
                    enddo                  
                  else
                    print *, ' >>> Error <<< : too many interactions'
                    stop
                  endif
                endif
                j = parts%celldata(j)
                goto 1
              endif
            enddo
          enddo
        enddo
      enddo
      parts%niac = niac

      end subroutine

      
!----------------------------------------------------------------------      
      subroutine link_list_omp(parts)
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by using a sorting grid linked list  
!----------------------------------------------------------------------
      use ifport
      implicit none

      class(particles) parts
      integer i, j, d, scale_k, ntotal, niac,m ,n
      integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)
      real(dp) hsml,dr,r,dx(3),tdwdx(3)

      INTEGER nthreads,n_per_threads, niac_per_threads, it, n_start(8), n_end(8)
      Integer niac_start(8),niac_end(8), last
      real(dp) t1,t2,t3

      ntotal  = parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
      hsml = parts%hsml(1)             !!! hsml is same for all particles!

      do i=1,ntotal
         parts%countiac(i) = 0
      enddo

      call parts%init_grid
      
!     Position particles on grid and create linked list:
      
      do i=1,ntotal
        call parts%grid_geom(i,parts%x(:,i),gcell)
        do d=1,parts%dim
          parts%xgcell(d,i) = gcell(d)
        enddo
        parts%celldata(i) = parts%grid(gcell(1),gcell(2),gcell(3))
        parts%grid(gcell(1),gcell(2),gcell(3)) = i
      enddo

!     Determine interaction parameters:
!!!$omp parallel
!!      nthreads = omp_get_num_threads()
!!!$omp end parallel
      nthreads = parts%nthreads
!      allocate(n_start(nthreads),n_end(nthreads))

      n_per_threads = ntotal/nthreads
      niac_per_threads = parts%max_interaction/nthreads
      do it = 1, nthreads
         n_start(it)=(it-1)*n_per_threads+1
         n_end(it) = it*n_per_threads
         parts%niac_start(it)=(it-1)*niac_per_threads
         parts%niac_end(it)=parts%niac_start(it)
      enddo   
      n_end(nthreads) = n_end(nthreads) + mod(ntotal,nthreads)

      parts%pair_i = 0; parts%pair_j=0
!      t1 = rtc()
!$omp parallel 
!$omp do private(niac,i,d,minxcell,maxxcell,dnxgcell,dpxgcell,xcell,ycell,zcell,j,dx,dr,r,tdwdx)
do it = 1, nthreads
       niac = parts%niac_start(it)
       do i = n_start(it),n_end(it)

!      niac = 0
!      do i=1,ntotal-1

!     Determine range of grid to go through:
         
        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,parts%dim
          dnxgcell(d) = parts%xgcell(d,i) - parts%ghsmlx(d)
          dpxgcell(d) = parts%xgcell(d,i) + parts%ghsmlx(d)
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),parts%ngridx(d))
        enddo

!     Search grid:
      
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = parts%grid(xcell,ycell,zcell)
 1            if (j.gt.i) then
                dx(1) = parts%x(1,i) - parts%x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,parts%dim
                  dx(d) = parts%x(d,i) - parts%x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                  if (niac.lt.parts%max_interaction) then
                  !if (niac.lt.niac_per_threads) then

!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 

                    niac = niac + 1
                    parts%pair_i(niac) = i
                    parts%pair_j(niac) = j
                    r = sqrt(dr)
                    !parts%countiac(i) = parts%countiac(i) + 1
                    !parts%countiac(j) = parts%countiac(j) + 1
                           
!--- Kernel and derivations of kernel

                    call parts%kernel(r,dx,hsml,parts%w(niac),tdwdx)  !!!!!!!!
                    do d = 1, parts%dim
                       parts%dwdx(d,niac)=tdwdx(d)
                    enddo                  
                  else
                    print *, ' >>> Error <<< : too many interactions', &
                    it,niac,niac_per_threads
                    stop
                  endif
                endif
                j = parts%celldata(j)
                goto 1
              endif
            enddo !xcell
          enddo !ycell
        enddo !zcell
      enddo !i

      parts%niac_start(it) = parts%niac_start(it)+1
      parts%niac_end(it) = niac
      enddo !it      
!$omp end do
!$omp end parallel
!      t2 = rtc()
!      write(*,*) t2-t1
      parts%niac = 0
      do it = 1, nthreads
         parts%niac = parts%niac + parts%niac_end(it)-parts%niac_start(it)+1
      enddo

      last = parts%max_interaction
      do while(parts%pair_i(last)==0)
         last = last - 1
      enddo  
      i = 1 
      do while(i<last)
         if(parts%pair_i(i)/=0)then
            i = i + 1
            cycle
         endif        
         parts%pair_i(i) = parts%pair_i(last)
         parts%pair_j(i) = parts%pair_j(last)
         parts%w(i)      = parts%w(last)
         parts%dwdx(1,i) = parts%dwdx(1,last)
         parts%dwdx(2,i) = parts%dwdx(2,last)
         last = last - 1
         do while(parts%pair_i(last)==0)
            last = last - 1
         enddo
         i = i + 1
      enddo   
      
call parts%get_niac_start_end

parts%niac_end(nthreads) = parts%niac_end(parts%nthreads)+n
!      t3 = rtc()
!      write(*,*) t1,t2,t3
!      write(*,*) t2-t1,t3-t2
!      pause
!      if(last/=parts%niac) stop 'adfasdf'
!      do i = 1, last
!         if(parts%pair_i(i)==0) stop 'adf'
!      enddo   

      end subroutine

!--------------------------------------------
     subroutine find_pairs(parts)
!--------------------------------------------
implicit none
class(particles) parts
integer nnps

nnps = parts%numeric%nnps
if(nnps == 1)then
   call direct_find(parts)
elseif(nnps == 2)then
   call link_list(parts)
!   call link_list(parts%itimestep, parts%ntotal+parts%nvirt,parts%hsml(1),  &
!   parts%x,parts%niac,parts%pair_i,parts%pair_j,parts%w,parts%dwdx,parts%countiac)
!   stop 'find_pairs: link_list method not implemented yet!'
elseif(nnps == 3)then
!    call tree_search(parts)
   stop 'find_pairs: tree_search method not implemented yet!'
endif

return
end subroutine

!----------------------------------------------------------------------
        subroutine kernel(this,r,dx,hsml,w,dwdx)   
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing kernel wij and its 
!   derivatives dwdxij.
!     if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!            = 2, Gauss kernel   (Gingold and Monaghan 1981) 
!            = 3, Quintic kernel (Morris 1997)
!            = 4, Wendland kernel
!----------------------------------------------------------------------
      implicit none
      class(particles) this
      real(dp) r, dx(3), hsml, w, dwdx(3)
      integer i, j, d, dim, skf      
      real(dp) q, dw, factor
      real(dp), parameter :: pi = 3.14159265358979323846 
      logical :: dbg = .false.

      if(dbg) write(*,*) 'In kernel...'

      dim = this%dim; skf = this%skf

      q = r/hsml 
      w = 0.e0
      do d=1,dim         
        dwdx(d) = 0.e0
      enddo   

      if (skf.eq.1) then     
        if (dim.eq.1) then
          factor = 1.e0/hsml
        elseif (dim.eq.2) then
          factor = 15.e0/(7.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 3.e0/(2.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif                                           
        if (q.ge.0.and.q.le.1.e0) then          
          w = factor * (2./3. - q*q + q**3 / 2.)
          do d = 1, dim
            dwdx(d) = factor * (-2.+3./2.*q)/hsml**2 * dx(d)       
          enddo   
        else if (q.gt.1.e0.and.q.le.2) then          
          w = factor * 1.e0/6.e0 * (2.-q)**3 
          do d = 1, dim
            dwdx(d) =-factor * 1.e0/6.e0 * 3.*(2.-q)**2/hsml * (dx(d)/r)        
          enddo              
	else
	  w=0.
          do d= 1, dim
            dwdx(d) = 0.
          enddo             
        endif     
                                    
      else if (skf.eq.2) then
      
        factor = 1.e0 / (hsml**dim * pi**(dim/2.))      
	if(q.ge.0.and.q.le.3) then
	  w = factor * exp(-q*q)
          do d = 1, dim
            dwdx(d) = w * ( -2.* dx(d)/hsml/hsml)
          enddo 
	else
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo 	   
	endif	       
	
      else if (skf.eq.3) then	
      
        if (dim.eq.1) then
          factor = 1.e0 / (120.e0*hsml)
        elseif (dim.eq.2) then
          factor = 7.e0 / (478.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif              
	if(q.ge.0.and.q.le.1) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * ( (-120 + 120*q - 50*q**2)      &
                              / hsml**2 * dx(d) )
          enddo 
	else if(q.gt.1.and.q.le.2) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4 + 30*(2-q)**4)     &
                             / hsml * (dx(d)/r) 
          enddo 
        else if(q.gt.2.and.q.le.3) then
          w = factor * (3-q)**5 
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4) / hsml * (dx(d)/r) 
          enddo 
        else   
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo  
        endif                      
     
      elseif(skf.eq.4)then    ! Wendland. refer to SPHysic manual

        if (dim.eq.1) then
          factor = 0.
        elseif (dim.eq.2) then
          factor = 7.e0 / (4.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 0.
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif              
	if(q.ge.0.and.q.le.2) then
          w = factor * ( (1-q/2)**4 *(1+2*q) )
          do d= 1, dim
            dwdx(d) = factor*(-5+15*q/2-15*q**2/4+5*q**3/8)/    &
                      hsml**2*dx(d)
          enddo 
	else   
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo  
        endif 
           
      endif 
		
      end subroutine

!--------------------------------------------------------------------------
      subroutine repulsive_force(parts)
!--------------------------------------------------------------------------
implicit none

class(particles) parts
double precision dx(3), rr, f, rr0, dd, p1, p2     
integer i, j, k, d, ii
           
rr0 = parts%dspp; dd = parts%numeric%dd
p1 = parts%numeric%p1; p2 = parts%numeric%p2
      
do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   if(parts%itype(i)*parts%itype(j)>0)cycle  
!   if(parts%itype(i).gt.0.and.parts%itype(j).lt.0)then  
      rr = 0.      
      do d=1,parts%dim
         dx(d) =  parts%x(d,i) -  parts%x(d,j)
         rr = rr + dx(d)*dx(d)
      enddo  
      rr = sqrt(rr)
      !if(rr.lt.rr0)then
      if(rr.gt.rr0)cycle
         f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
         
         ii = i
         if(parts%itype(i)<0)ii=j 

         do d = 1, parts%dim
            parts%dvx(d, ii) = parts%dvx(d, ii) + dd*dx(d)*f
         enddo
      !endif
   !endif        
enddo   
       
return
end subroutine
      
      !--------------------------------------------------------------------------
      subroutine repulsive_force_omp(parts)
!--------------------------------------------------------------------------
implicit none

class(particles) parts
double precision dx(3), rr, f, rr0, dd, p1, p2     
integer i, j, k, d, ii, it, ntotal
real(dp) local(2,parts%maxn,8)
ntotal = parts%ntotal + parts%nvirt
rr0 = parts%dspp; dd = parts%numeric%dd
p1 = parts%numeric%p1; p2 = parts%numeric%p2
  
!$omp parallel 
!$omp do private(i,it,j,d,k,rr,dx,f,ii)    
do it = 1,parts%nthreads
        do i =1,ntotal
          do d = 1,parts%dim
            local(d,i,it) = 0.
          enddo
        enddo
    do k=parts%niac_start(it),parts%niac_end(it)

   i = parts%pair_i(k)
   j = parts%pair_j(k)
   if(parts%itype(i)*parts%itype(j)>0)cycle  
!   if(parts%itype(i).gt.0.and.parts%itype(j).lt.0)then  
      rr = 0.      
      do d=1,parts%dim
         dx(d) =  parts%x(d,i) -  parts%x(d,j)
         rr = rr + dx(d)*dx(d)
      enddo  
      rr = sqrt(rr)
      !if(rr.lt.rr0)then
      if(rr.gt.rr0)cycle
         f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
         
         ii = i
         if(parts%itype(i)<0)ii=j 

         do d = 1, parts%dim
             local(d,ii,it) = local(d,ii,it) +dd*dx(d)*f
!            parts%dvx(d, ii) = parts%dvx(d, ii) + dd*dx(d)*f
         enddo
    enddo
      !endif
   !endif        
enddo   
!$omp enddo
!$omp barrier
!$omp do private(it)
     do i = 1,ntotal
         do it = 1,parts%nthreads
           do d = 1,parts%dim
             parts%dvx(d,i) = parts%dvx(d,i) +local(d,i,it)
           enddo
         enddo
     enddo
!$omp end do
!$omp end parallel

return
end subroutine


! Calculate partial derivatives of a field
!-------------------------------------------
          function df(parts,f,x)
!-------------------------------------------
implicit none

real(dp) f(:)
character(len=1) x
class(particles) parts
real(dp), allocatable :: df(:)
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k

allocate(df(size(f))); df = 0.

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f(i)+f(j))*dwdx(k)
   df(i) = df(i) + parts%mass(j)/parts%rho(j)*fwx
   df(j) = df(j) - parts%mass(i)/parts%rho(i)*fwx
enddo

end function

! Calculate partial derivatives of a field
!-------------------------------------------
          function df_omp(parts,f,x)
!-------------------------------------------
implicit none

real(dp) f(:)
character(len=1) x
class(particles) parts
real(dp), allocatable, dimension(:) :: df_omp
real(dp) df_local(parts%maxn,8)
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k, ntotal, it, nthreads

!write(*,*) 'In df_omp...'

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads
!write(*,*) 'sadf', nthreads

allocate(df_omp(ntotal))
!if(nthreads>1)then
!   allocate(df_local(ntotal,nthreads))
!   call parts%get_niac_start_end
!endif   

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

!$omp parallel
!$omp do private(i,j,k,fwx)
do it = 1, parts%nthreads
   do i = 1, ntotal
      df_local(i,it) = 0.d0
   enddo
   do k = parts%niac_start(it), parts%niac_end(it)


!do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f(i)+f(j))*dwdx(k)
   df_local(i,it) = df_local(i,it) + parts%mass(j)/parts%rho(j)*fwx
   df_local(j,it) = df_local(j,it) - parts%mass(i)/parts%rho(i)*fwx
   enddo !k
enddo !it
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1, ntotal
   df_omp(i) = 0.d0
   do it = 1, nthreads
      df_omp(i) = df_omp(i)+df_local(i,it)
   enddo
enddo   
!$omp end do
!$omp end parallel

end function

! Calculate partial derivatives of a field
!-------------------------------------------
          function df2(parts,f,x)
!-------------------------------------------
implicit none

real(dp) f(:)
character(len=1) x
class(particles) parts
real(dp), allocatable :: df2(:)
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k

allocate(df2(size(f))); df2 = 0.

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f(j)-f(i))*dwdx(k)
   df2(i) = df2(i) + parts%mass(j)*fwx
   df2(j) = df2(j) + parts%mass(i)*fwx
enddo

do i = 1, parts%ntotal + parts%nvirt 
   df2(i) = df2(i)/parts%rho(i)
enddo   

end function

!-------------------------------------------
          function df2_omp(parts,f,x)
!-------------------------------------------
implicit none

real(dp) f(:)
character(len=1) x
class(particles) parts
real(dp), allocatable :: df2_omp(:)
real(dp) df_local(parts%maxn,8)
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k,ntotal,nthreads,it

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads

allocate(df2_omp(size(f))); df2_omp = 0.
!if(nthreads>1)then
!   allocate(df_local(ntotal,nthreads))
!   call parts%get_niac_start_end
!endif  

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

!$omp parallel
!$omp do private(i,j,k,fwx)
do it = 1, parts%nthreads
   do i = 1, ntotal
      df_local(i,it) = 0.d0
   enddo
   
   do k = parts%niac_start(it), parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f(j)-f(i))*dwdx(k)
   df_local (i,it) = df_local(i,it) + parts%mass(j)*fwx
   df_local (j,it) = df_local(j,it) + parts%mass(i)*fwx
   enddo
enddo
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1, ntotal
   df2_omp(i) = 0.d0
   do it = 1, nthreads
      df2_omp(i) = df2_omp(i)+df_local(i,it)
   enddo
enddo   
!$omp end do
!$omp end parallel
!   df2(i) = df2(i) + parts%mass(j)*fwx
!   df2(j) = df2(j) + parts%mass(i)*fwx


do i = 1, parts%ntotal + parts%nvirt 
   df2_omp(i) = df2_omp(i)/parts%rho(i)
enddo   

end function
          
! Calculate partial derivatives of a field
!-------------------------------------------
          function df3(parts,f,x)
!------------------------------------------- 
implicit none

real(dp) f(:)
character(len=1) x
class(particles) parts
real(dp),allocatable :: df3(:)
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i,j,k

allocate(df3(size(f))); df3 = 0.

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = ((f(i)/parts%rho(i)**2)+(f(j)/parts%rho(j)**2))*dwdx(k)
   df3(i) = df3(i) + parts%mass(j)*fwx
   df3(j) = df3(j) - parts%mass(i)*fwx
enddo

do i = 1, parts%ntotal + parts%nvirt
   df3(i) = df3(i)*parts%rho(i)
enddo

end function

!-------------------------------------------
          function df3_omp(parts,f,x)
!------------------------------------------- 
implicit none

real(dp) f(:)
character(len=1) x
class(particles) parts
real(dp), allocatable, dimension(:) :: df3_omp
real(dp) df3_local(parts%maxn,8)
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i,j,k,it,ntotal,nthreads

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads

allocate(df3_omp(ntotal))

!if(nthreads>1)then
!   allocate(df3_local(ntotal,nthreads))
!   call parts%get_niac_start_end
!endif 

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

!$omp parallel
!$omp do private(i,j,k,fwx)
do it = 1, parts%nthreads
   do i = 1, ntotal
      df3_local(i,it) = 0.d0
   enddo
   
   do k = parts%niac_start(it), parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = ((f(i)/parts%rho(i)**2)+(f(j)/parts%rho(j)**2))*dwdx(k)
   df3_local (i,it)= df3_local (i,it) + parts%mass(j)*fwx
   df3_local (j,it)= df3_local (j,it) - parts%mass(i)*fwx
!   df3_omp(i) = df3_omp(i) + parts%mass(j)*fwx
!   df3_omp(j) = df3_omp(j) - parts%mass(i)*fwx
   enddo
enddo
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1, ntotal
   df3_omp(i) = 0.d0
   do it = 1, nthreads
      df3_omp(i) = df3_omp(i)+df3_local(i,it)
   enddo
enddo   
!$omp end do
!$omp end parallel

do i = 1, parts%ntotal + parts%nvirt
   df3_omp(i) = df3_omp(i)*parts%rho(i)
enddo


end function



!-----------------------------------------------------------------------
      subroutine velocity_divergence(parts)
!-----------------------------------------------------------------------
      implicit none

      class(particles) parts      
      integer, pointer, dimension(:) :: pair_i, pair_j
      double precision, pointer, dimension(:) :: mass, rho, vcc
      double precision, pointer, dimension(:,:) :: vx, dwdx
      integer i,j,k,d, dim, ntotal, niac
      double precision dvx(3), hvcc
      
      pair_i => parts%pair_i
      pair_j => parts%pair_j
      mass   => parts%mass
      rho    => parts%rho
      vcc    => parts%vcc
      vx     => parts%vx
      dwdx   => parts%dwdx

      ntotal = parts%ntotal + parts%nvirt
      niac   = parts%niac; dim = parts%dim 

      do i=1,ntotal
        vcc(i) = 0.e0         
      enddo
    
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,j) - vx(d,i) 
        enddo        
        hvcc = dvx(1)*dwdx(1,k)
        do d=2,dim
          hvcc = hvcc + dvx(d)*dwdx(d,k)
        enddo    
         vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
         vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
      enddo  
   
      return
      end subroutine

!-----------------------------------------------------------------------
      subroutine velocity_divergence_omp(parts)
!-----------------------------------------------------------------------
      implicit none

      class(particles) parts      
      integer, pointer, dimension(:) :: pair_i, pair_j
      double precision, pointer, dimension(:) :: mass, rho, vcc
      double precision, pointer, dimension(:,:) :: vx, dwdx
      real local(parts%maxn,8)
      integer i,j,k,d, dim, ntotal, niac,it
      double precision dvx(3), hvcc
      
      pair_i => parts%pair_i
      pair_j => parts%pair_j
      mass   => parts%mass
      rho    => parts%rho
      vcc    => parts%vcc
      vx     => parts%vx
      dwdx   => parts%dwdx

      ntotal = parts%ntotal + parts%nvirt
      niac   = parts%niac; dim = parts%dim 

      do i=1,ntotal
        vcc(i) = 0.e0         
      enddo
    
     !$omp parallel 
     !$omp do private(i,it,j,d,k,dvx,hvcc)    
     do it = 1,parts%nthreads
         do i =1,ntotal
             local(i,it) = 0.
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,j) - vx(d,i) 
        enddo        
        hvcc = dvx(1)*dwdx(1,k)
        do d=2,dim
          hvcc = hvcc + dvx(d)*dwdx(d,k)
        enddo    
        local(i,it) = local(i,it) + mass(j)*hvcc/rho(j)
        local(j,it) = local(j,it) + mass(i)*hvcc/rho(i)
!         vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
!         vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
      enddo  
     enddo
     !$omp end do
     !$omp barrier
     !$omp do private(it,i)
     do i = 1,ntotal
         do it = 1,parts%nthreads
             vcc(i) = vcc(i) + local(i,it)
         enddo
     enddo
     !$omp end do
     !$omp end parallel


   
      return
      end subroutine
!---------------------------------------------------------------------
      subroutine initial_density(parts)
!---------------------------------------------------------------------
      implicit none

      class(particles) parts

      type(material), pointer :: water, soil
      integer ntotal
      
      ntotal = parts%ntotal+parts%nvirt

      if(parts%imaterial=='water')then

         water => parts%material
         parts%rho(1:ntotal) = water%rho0*(parts%p(1:ntotal)/water%b+1) &
                            **(1/water%gamma)

      elseif(parts%imaterial=='soil')then

         soil => parts%material
         !parts%rho(1:ntotal) = soil%rho0*(parts%p(1:ntotal)/soil%k+1)

      endif
 
      return
      end subroutine

!----------------------------------------------------------------------
                      subroutine pressure(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water, soil
      integer ntotal, i
      
      ntotal = parts%ntotal+parts%nvirt

      if(parts%imaterial=='water')then

         water => parts%material
         parts%p(1:ntotal) = water%b*((parts%rho(1:ntotal)/(water%rho0  &
                            *parts%vof(1:ntotal))) &   !!! False density
                            **water%gamma-1)  

! Tension instability
!                              if(water_tension_instability==1)then
!         do i = 1, ntotal
!            if(parts%p(i)<0)then
!               parts%p(i) = 0.d0
!               parts%rho(i) = water%rho0*parts%vof(i)
!            endif
!         enddo
!                              endif

         !parts%c(1:ntotal) = water%c         
         parts%c(1:ntotal) = water%c*(parts%rho(1:ntotal)/(water%rho0*parts%vof(1:ntotal)))**3.0         

      elseif(parts%imaterial=='soil')then

         soil => parts%material
!         if(parts%soil_pressure==1)          &
!         parts%p(1:ntotal) = soil%k*(parts%rho(1:ntotal)/soil%rho0-1)
!!         parts%p(1:ntotal) = parts%p(1:ntotal)
!!     &                      -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
         parts%c(1:ntotal) = soil%c

      endif

      return      
      end subroutine

!--------------------------------------------------------------
            subroutine pressure_water(parts)
!--------------------------------------------------------------
implicit none

class(particles) parts
type(material), pointer :: water, soil
integer ntotal, i
      
ntotal = parts%ntotal+parts%nvirt

water => parts%material
parts%p(1:ntotal) = water%b*((parts%rho(1:ntotal)/(water%rho0  &
                   *parts%vof(1:ntotal))) &   !!! False density
                  **water%gamma-1)

parts%c(1:ntotal) = water%c*(parts%rho(1:ntotal)/(water%rho0*parts%vof(1:ntotal)))**3.0         

return
end subroutine

!--------------------------------------------------------------
         subroutine tension_instability_water(parts)
!--------------------------------------------------------------
implicit none

class(particles) parts
type(material), pointer :: water, soil
integer ntotal, i
      
ntotal = parts%ntotal+parts%nvirt
water => parts%material

do i = 1, ntotal
   if(parts%p(i)<0)then
      parts%p(i) = 0.d0
      !parts%rho(i) = water%rho0*parts%vof(i)
   endif
enddo

return
end subroutine


!------------------------------------------------------------
              subroutine pressure_soil(parts)
!------------------------------------------------------------
implicit none

class(particles) parts
type(material), pointer :: water, soil
integer ntotal, i
      
ntotal = parts%ntotal+parts%nvirt

soil => parts%material
parts%p(1:ntotal) = soil%k*(parts%rho(1:ntotal)/soil%rho0-1)
! parts%p(1:ntotal) = parts%p(1:ntotal)  &
!                    -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
parts%c(1:ntotal) = soil%c

return      
end subroutine

!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine sum_density(parts) 
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%maxn)     

      ntotal = parts%ntotal + parts%nvirt

!     wi(maxn)---integration of the kernel itself
        
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
!     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=selfdens*parts%mass(i)/parts%rho(i)
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        wi(i) = wi(i) + parts%mass(j)/parts%rho(j)*parts%w(k)
        wi(j) = wi(j) + parts%mass(i)/parts%rho(i)*parts%w(k)
      enddo

!     Secondly calculate the rho integration over the space

      do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        parts%rho(i) = selfdens*parts%mass(i)
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        parts%rho(i) = parts%rho(i) + parts%mass(j)*parts%w(k)
        parts%rho(j) = parts%rho(j) + parts%mass(i)*parts%w(k)
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
!      if (nor_density) then 
        do i=1, ntotal
          parts%rho(i)=parts%rho(i)/wi(i)
        enddo
!      endif 
 
      end subroutine

!----------------------------------------------------------------------
      subroutine sum_density_omp(parts) 
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d,it 
      real(dp) selfdens, hv(3), r, wi(parts%maxn)     
      real local(parts%maxn,8)
      ntotal = parts%ntotal + parts%nvirt

!     wi(maxn)---integration of the kernel itself
        
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
!     Firstly calculate the integration of the kernel over the space
     !$omp parallel
     !$omp do
      do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=selfdens*parts%mass(i)/parts%rho(i)
      enddo
      !$omp end do

     !$omp do private(i,it,j,d,k)    
     do it = 1,parts%nthreads
         do i =1,ntotal
             local(i,it) = 0.
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        local(i,it) = local(i,it) + parts%mass(j)/parts%rho(j)*parts%w(k)
        local(j,it) = local(j,it) + parts%mass(i)/parts%rho(i)*parts%w(k)
!        wi(i) = wi(i) + parts%mass(j)/parts%rho(j)*parts%w(k)
!        wi(j) = wi(j) + parts%mass(i)/parts%rho(i)*parts%w(k)
      enddo
     enddo
     !$omp enddo
     !$omp barrier
     !$omp do private(it)
     do i = 1,ntotal
         do it = 1,parts%nthreads
             wi(i) = wi(i) + local(i,it)
         enddo
     enddo
     !$omp end do


!     Secondly calculate the rho integration over the space
     !$omp do
      do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        parts%rho(i) = selfdens*parts%mass(i)
      enddo
     !$omp end do
!     Calculate SPH sum for rho:

     !$omp do private(i,it,j,d,k)   
     do it = 1,parts%nthreads
         do i =1,ntotal
             local(i,it) = 0.
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        local(i,it) = local(i,it) + parts%mass(j)*parts%w(k)
        local(j,it) = local(j,it) + parts%mass(i)*parts%w(k)
!        parts%rho(i) = parts%rho(i) + parts%mass(j)*parts%w(k)
!        parts%rho(j) = parts%rho(j) + parts%mass(i)*parts%w(k)
      enddo
     enddo
     !$omp enddo
     !$omp barrier
     !$omp do private(it)
     do i = 1,ntotal
         do it = 1,parts%nthreads
             parts%rho(i) = parts%rho(i) + local(i,it)
         enddo
     enddo
     !$omp end do

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
!      if (nor_density) then 
     !$omp do
        do i=1, ntotal
          parts%rho(i)=parts%rho(i)/wi(i)
        enddo
!      endif 
     !$omp end do
     !$omp end parallel
     
end subroutine


! Subroutine to calculate the density with SPH continuiity approach.
!----------------------------------------------------------------------      
      subroutine con_density(parts)
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      

      integer ntotal,i,j,k,d    
      double precision vcc, dvx(3) 
      
      ntotal = parts%ntotal + parts%nvirt

      do i = 1, ntotal
        parts%drho(i) = 0.
      enddo
     
      do k=1,parts%niac      
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        do d=1, parts%dim
          dvx(d) = parts%vx(d,i) - parts%vx(d,j) 
        enddo        
        vcc = dvx(1)* parts%dwdx(1,k)        
        do d=2, parts%dim
          vcc = vcc + dvx(d)*parts%dwdx(d,k)
        enddo    
        parts%drho(i) = parts%drho(i) + parts%mass(j)*vcc
        parts%drho(j) = parts%drho(j) + parts%mass(i)*vcc       
      enddo    

      end subroutine
! Subroutine to calculate the density with SPH continuiity approach.
!----------------------------------------------------------------------      
      subroutine con_density_omp(parts)
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      

      integer ntotal,i,j,k,d,it
      real  div_v_local(parts%maxn,8)
!      real, allocatable, dimension(:,:) :: div_v_local
      double precision vcc, dvx(3) 

      ntotal = parts%ntotal + parts%nvirt
 !      call parts%get_niac_start_end
      do i = 1, ntotal
        parts%drho(i) = 0.
      enddo
     !$omp parallel 
     !$omp do private(i,j,d,k,vcc,dvx)
     do it = 1,parts%nthreads
         do i =1,ntotal
            div_v_local(i,it) = 0.
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)      
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        do d=1, parts%dim
          dvx(d) = parts%vx(d,i) - parts%vx(d,j) 
        enddo        
        vcc = dvx(1)* parts%dwdx(1,k)        
        do d=2, parts%dim
          vcc = vcc + dvx(d)*parts%dwdx(d,k)
        enddo    
        div_v_local(i,it) = div_v_local(i,it) + parts%mass(j)*vcc
        div_v_local(j,it) = div_v_local(j,it) + parts%mass(i)*vcc       
      enddo    
     enddo
     !$omp end do
     !$omp barrier
     !$omp do private(it)
     do i = 1,ntotal
         do it = 1,parts%nthreads
             parts%drho(i) = parts%drho(i) + div_v_local(i,it)
         enddo
     enddo
     !$omp end do
     !$omp end parallel
 !    deallocate(div_v_local)
      end subroutine

!----------------------------------------------------------------------
      subroutine newtonian_fluid(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: liquid
      integer ntotal

      ntotal = parts%ntotal + parts%nvirt
      liquid => parts%material

      parts%sxx(1:ntotal) = liquid%viscosity*parts%txx(1:ntotal)
      parts%syy(1:ntotal) = liquid%viscosity*parts%tyy(1:ntotal)
      parts%sxy(1:ntotal) = liquid%viscosity*parts%txy(1:ntotal)

      return
      end subroutine

! Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
!----------------------------------------------------------------------      
      subroutine art_visc(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts

      type(numerical), pointer :: numeric
      real(dp) dx, dvx(3), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml
      integer i,j,k,d,dim,ntotal,niac

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim
            
      numeric  => parts%numeric      
      alpha = numeric%alpha; beta = numeric%beta; etq = numeric%etq
         
!     Calculate SPH sum for artificial viscosity
      
      do k=1,niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        mhsml= (parts%hsml(i)+parts%hsml(j))/2.
        vr = 0.e0
        rr = 0.e0
        do d=1,dim
          dvx(d) = parts%vx(d,i) - parts%vx(d,j)
          dx     = parts%x(d,i)  - parts%x(d,j)
          vr     = vr + dvx(d)*dx
          rr     = rr + dx*dx
        enddo

!     Artificial viscous force only if v_ij * r_ij < 0

        if (vr.lt.0.e0) then

!     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )
            
          muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)
          
!     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

          mc   = 0.5e0*(parts%c(i) + parts%c(j))
          mrho = 0.5e0*(parts%rho(i) + parts%rho(j))
          piv  = (beta*muv - alpha*mc)*muv/mrho              

!     Calculate SPH sum for artificial viscous force

          do d=1,dim
            h = -piv*parts%dwdx(d,k)
            parts%dvx(d,i) = parts%dvx(d,i) + parts%mass(j)*h
            parts%dvx(d,j) = parts%dvx(d,j) - parts%mass(i)*h
!            dedt(i) = dedt(i) - mass(j)*dvx(d)*h
!            dedt(j) = dedt(j) - mass(i)*dvx(d)*h
          enddo
        endif
      enddo

!     Change of specific internal energy:

      return
      end subroutine

      
! Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
!----------------------------------------------------------------------      
      subroutine art_visc_omp(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts

      type(numerical), pointer :: numeric
      real local(2,parts%maxn,8)
      real(dp) dx, dvx(3), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml
      integer i,j,k,d,dim,ntotal,niac,it

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim
            
      numeric  => parts%numeric      
      alpha = numeric%alpha; beta = numeric%beta; etq = numeric%etq
         
!     Calculate SPH sum for artificial viscosity
      
     !$omp parallel 
     !$omp do private(i,j,d,k,mhsml,dvx,dx,vr,rr,muv,mc,mrho,piv,h)    
     do it = 1,parts%nthreads
         do i =1,ntotal
             do d = 1,dim
             local(d,i,it) = 0.
             enddo
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)

        i = parts%pair_i(k)
        j = parts%pair_j(k)
        mhsml= (parts%hsml(i)+parts%hsml(j))/2.
        vr = 0.e0
        rr = 0.e0
        do d=1,dim
          dvx(d) = parts%vx(d,i) - parts%vx(d,j)
          dx     = parts%x(d,i)  - parts%x(d,j)
          vr     = vr + dvx(d)*dx
          rr     = rr + dx*dx
        enddo

!     Artificial viscous force only if v_ij * r_ij < 0

        if (vr.lt.0.e0) then

!     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )
            
          muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)
          
!     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

          mc   = 0.5e0*(parts%c(i) + parts%c(j))
          mrho = 0.5e0*(parts%rho(i) + parts%rho(j))
          piv  = (beta*muv - alpha*mc)*muv/mrho              

!     Calculate SPH sum for artificial viscous force

          do d=1,dim
            h = -piv*parts%dwdx(d,k)
            local(d,i,it) = local(d,i,it) + parts%mass(j)*h
            local(d,j,it) = local(d,j,it) - parts%mass(i)*h
!            parts%dvx(d,i) = parts%dvx(d,i) + parts%mass(j)*h
!            parts%dvx(d,j) = parts%dvx(d,j) - parts%mass(i)*h
!            dedt(i) = dedt(i) - mass(j)*dvx(d)*h
!            dedt(j) = dedt(j) - mass(i)*dvx(d)*h
          enddo
        endif
      enddo
      enddo
     !$omp enddo
     !$omp barrier
     !$omp do private(it,d)
     do i = 1,ntotal
         do d = 1,dim
             do it = 1,parts%nthreads
               parts%dvx(d,i) = parts%dvx(d,i) + local(d,i,it)
             enddo
         enddo
     enddo
     !$omp end do
     !$omp end parallel

!     Change of specific internal energy:

      return
      end subroutine


!-------------------------------------------------------------------
      subroutine delta_sph(parts,f,df)
!-------------------------------------------------------------------
      implicit none

      class(particles) parts
      real(dp), dimension(:) :: f,df

      real(dp) dx(3),delta, muv, rr, h
      integer i,j,k,d,ntotal,niac,dim

!      write(*,*) 'In art_density...'

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim            
      delta    = parts%numeric%delta
          
      do k=1,niac
         i = parts%pair_i(k)
         j = parts%pair_j(k)
         rr = 0.e0
         do d=1,dim
            dx(d)  =  parts%x(d,i) -  parts%x(d,j)
            rr     = rr + dx(d)*dx(d)
         enddo
            
         muv = 2.0*(f(i)-f(j))/rr

         h = 0.d0
         do d=1,dim
            !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
            h = h + dx(d)*muv*parts%dwdx(d,k)
         enddo
         df(i) = df(i) + delta*parts%hsml(i)*parts%c(i)*parts%mass(j)*h/parts%rho(j)
         df(j) = df(j) - delta*parts%hsml(j)*parts%c(j)*parts%mass(i)*h/parts%rho(i)
      enddo

      return
      end subroutine

!     Subroutine to calculate the average velocity to correct velocity
!     for preventing.penetration (monaghan, 1992)     

!-------------------------------------------------------------------
      subroutine delta_sph_omp(parts,f,df)
!-------------------------------------------------------------------
      implicit none

      class(particles) parts
      real(dp), dimension(:) :: f,df
      real local(parts%maxn,8)
      real(dp) dx(3),delta, muv, rr, h
      integer i,j,k,d,ntotal,niac,dim,it

!      write(*,*) 'In art_density...'

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim            
      delta    = parts%numeric%delta
          
     !$omp parallel 
     !$omp do private(i,it,j,d,k,rr,muv,dx,h)    
     do it = 1,parts%nthreads
         do i =1,ntotal
             local(i,it) = 0.
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)
         i = parts%pair_i(k)
         j = parts%pair_j(k)
         rr = 0.e0
         do d=1,dim
            dx(d)  =  parts%x(d,i) -  parts%x(d,j)
            rr     = rr + dx(d)*dx(d)
         enddo
            
         muv = 2.0*(f(i)-f(j))/rr

         h = 0.d0
         do d=1,dim
            !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
            h = h + dx(d)*muv*parts%dwdx(d,k)
         enddo
!          df(i) = df(i) + delta*parts%hsml(i)*parts%c(i)*parts%mass(j)*h/parts%rho(j)
!          df(j) = df(j) - delta*parts%hsml(j)*parts%c(j)*parts%mass(i)*h/parts%rho(i)
         local(i,it) = local(i,it) + delta*parts%hsml(i)*parts%c(i)*parts%mass(j)*h/parts%rho(j)
         local(j,it) = local(j,it) - delta*parts%hsml(j)*parts%c(j)*parts%mass(i)*h/parts%rho(i)
      enddo
     enddo
     !$omp enddo
     !$omp barrier
     !$omp do private(it,i)
     do i = 1,ntotal
         do it = 1,parts%nthreads
             df(i) = df(i) + local(i,it)
         enddo
     enddo
     !$omp end do
     !$omp end parallel
      return
      end subroutine
      
!     Subroutine to calculate the average velocity to correct velocity
!     for preventing.penetration (monaghan, 1992)      
!----------------------------------------------------------------------      
      subroutine av_vel(parts)
!----------------------------------------------------------------------   
      implicit none

      class(particles),target :: parts
      real(dp) vcc, dvx(3), epsilon, mrho
      integer i,j,k,d,ntotal,niac    

      ntotal = parts%ntotal
      niac   = parts%niac
      epsilon = parts%numeric%epsilon
      
      parts%av = 0.d0
     
      do k=1,niac       
         i = parts%pair_i(k)
         j = parts%pair_j(k)       
         mrho = (parts%rho(i)+parts%rho(j))/2.0
         do d=1,parts%dim
            dvx(d) = parts%vx(d,i) - parts%vx(d,j)            
            parts%av(d, i) = parts%av(d,i) - parts%mass(j)*dvx(d)/mrho*parts%w(k)
            parts%av(d, j) = parts%av(d,j) + parts%mass(i)*dvx(d)/mrho*parts%w(k)       
         enddo                    
      enddo  
        
      do i = 1, ntotal
         do d = 1, parts%dim
            parts%av(d,i) = epsilon * parts%av(d,i)
         enddo 
      enddo             

      return
      end subroutine

!----------------------------------------------------------------------      
      subroutine av_vel_omp(parts)
!----------------------------------------------------------------------   
      implicit none

      class(particles),target :: parts
      real(dp) vcc, dvx(3), epsilon, mrho
      integer i,j,k,d,ntotal,niac,it
      real local(2,parts%maxn,8)

      ntotal = parts%ntotal
      niac   = parts%niac
      epsilon = parts%numeric%epsilon
      
      parts%av = 0.d0
     
     !$omp parallel 
     !$omp do private(i,j,d,k,mrho,dvx)    
     do it = 1,parts%nthreads
         do d = 1,parts%dim
          do i =1,ntotal
             local(d,i,it) = 0.
          enddo
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)
    
         i = parts%pair_i(k)
         j = parts%pair_j(k)       
         mrho = (parts%rho(i)+parts%rho(j))/2.0
         do d=1,parts%dim
            dvx(d) = parts%vx(d,i) - parts%vx(d,j)            
            local(d,i,it) = local(d,i,it) - parts%mass(j)*dvx(d)/mrho*parts%w(k)
            local(d,j,it) = local(d,j,it) + parts%mass(i)*dvx(d)/mrho*parts%w(k)
!            parts%av(d, i) = parts%av(d,i) - parts%mass(j)*dvx(d)/mrho*parts%w(k)
!            parts%av(d, j) = parts%av(d,j) + parts%mass(i)*dvx(d)/mrho*parts%w(k)       
         enddo                    
      enddo  
     enddo
     !$omp enddo
     !$omp barrier
     !$omp do private(it,d)
     do i = 1,ntotal
         do d = 1,parts%dim
           do it = 1,parts%nthreads
             parts%av(d, i) = parts%av(d, i) + local(d,i,it)
           enddo
         enddo
     enddo
     !$omp end do
     !$omp end parallel
      do i = 1, ntotal
         do d = 1, parts%dim
            parts%av(d,i) = epsilon * parts%av(d,i)
         enddo 
      enddo             

      return
      end subroutine

!----------------------------------------------------------------------
      subroutine shear_strain_rate(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz
      integer i, j, k, d, dim, ntotal, niac

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac; dim = parts%dim

      parts%txx(1:ntotal) = 0.e0
      if(dim>=2)then
         parts%tyy(1:ntotal) = 0.e0
         parts%txy(1:ntotal) = 0.e0
      endif
      if(dim==3)then
         parts%tzz(1:ntotal) = 0.e0
         parts%txz(1:ntotal) = 0.e0
         parts%tyz(1:ntotal) = 0.e0
      endif
      
!     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

        do k=1,niac
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          do d=1,dim
             dvx(d) = parts%vx(d,j) - parts%vx(d,i)
          enddo
          if (dim.eq.1) then 
             hxx = 2.e0*dvx(1)*parts%dwdx(1,k)        
          else if (dim.eq.2) then           
             hxx = 2.e0*dvx(1)*parts%dwdx(1,k) - dvx(2)*parts%dwdx(2,k) 
             hxy = dvx(1)*parts%dwdx(2,k) + dvx(2)*parts%dwdx(1,k)
             hyy = 2.e0*dvx(2)*parts%dwdx(2,k) - dvx(1)*parts%dwdx(1,k)
          else if (dim.eq.3) then
             hxx = 2.e0*dvx(1)*parts%dwdx(1,k) - dvx(2)*parts%dwdx(2,k)   &
                   - dvx(3)*parts%dwdx(3,k) 
             hxy = dvx(1)*parts%dwdx(2,k) + dvx(2)*parts%dwdx(1,k)
             hxz = dvx(1)*parts%dwdx(3,k) + dvx(3)*parts%dwdx(1,k)          
             hyy = 2.e0*dvx(2)*parts%dwdx(2,k) - dvx(1)*parts%dwdx(1,k)   &    
                   - dvx(3)*parts%dwdx(3,k)
             hyz = dvx(2)*parts%dwdx(3,k) + dvx(3)*parts%dwdx(2,k)
             hzz = 2.e0*dvx(3)*parts%dwdx(3,k) - dvx(1)*parts%dwdx(1,k)   &
                   - dvx(2)*parts%dwdx(2,k)
          endif                              
          hxx = 2.e0/3.e0*hxx
          hyy = 2.e0/3.e0*hyy
          hzz = 2.e0/3.e0*hzz
          if (dim.eq.1) then 
             parts%txx(i) = parts%txx(i) + parts%mass(j)*hxx/parts%rho(j)
             parts%txx(j) = parts%txx(j) + parts%mass(i)*hxx/parts%rho(i)               
          else if (dim.eq.2) then           
             parts%txx(i) = parts%txx(i) + parts%mass(j)*hxx/parts%rho(j)
             parts%txx(j) = parts%txx(j) + parts%mass(i)*hxx/parts%rho(i)   
             parts%txy(i) = parts%txy(i) + parts%mass(j)*hxy/parts%rho(j)
             parts%txy(j) = parts%txy(j) + parts%mass(i)*hxy/parts%rho(i)            
             parts%tyy(i) = parts%tyy(i) + parts%mass(j)*hyy/parts%rho(j)
             parts%tyy(j) = parts%tyy(j) + parts%mass(i)*hyy/parts%rho(i)          
          else if (dim.eq.3) then
!             txx(i) = txx(i) + mass(j)*hxx/rho(j)
!             txx(j) = txx(j) + mass(i)*hxx/rho(i)   
!             txy(i) = txy(i) + mass(j)*hxy/rho(j)
!             txy(j) = txy(j) + mass(i)*hxy/rho(i) 
!             txz(i) = txz(i) + mass(j)*hxz/rho(j)
!             txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
!             tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
!             tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
!             tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
!             tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
!             tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
!             tzz(j) = tzz(j) + mass(i)*hzz/rho(i)                 
          endif                              
        enddo

      return
      end subroutine      

!----------------------------------------------------------------------
      subroutine shear_strain_rate_omp(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz
      integer i, j, k, d, dim, ntotal, niac,it
      real local_xx(parts%maxn,8)
      real local_xy(parts%maxn,8)
      real local_yy(parts%maxn,8)

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac; dim = parts%dim

      parts%txx(1:ntotal) = 0.e0
      if(dim>=2)then
         parts%tyy(1:ntotal) = 0.e0
         parts%txy(1:ntotal) = 0.e0
      endif
      if(dim==3)then
         parts%tzz(1:ntotal) = 0.e0
         parts%txz(1:ntotal) = 0.e0
         parts%tyz(1:ntotal) = 0.e0
      endif
      
!     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c
      
     !$omp parallel 
     !$omp do private(i,j,it,d,k,dvx,hxx,hxy,hyy)    
     do it = 1,parts%nthreads
         do i =1,ntotal 
             local_xx(i,it) = 0.
             local_xy(i,it) = 0.
             local_yy(i,it) = 0.
         enddo
      do k=parts%niac_start(it),parts%niac_end(it)

          i = parts%pair_i(k)
          j = parts%pair_j(k)
          do d=1,dim
             dvx(d) = parts%vx(d,j) - parts%vx(d,i)
          enddo
          if (dim.eq.1) then 
             hxx = 2.e0*dvx(1)*parts%dwdx(1,k)        
          else if (dim.eq.2) then           
             hxx = 2.e0*dvx(1)*parts%dwdx(1,k) - dvx(2)*parts%dwdx(2,k) 
             hxy = dvx(1)*parts%dwdx(2,k) + dvx(2)*parts%dwdx(1,k)
             hyy = 2.e0*dvx(2)*parts%dwdx(2,k) - dvx(1)*parts%dwdx(1,k)
          else if (dim.eq.3) then
             hxx = 2.e0*dvx(1)*parts%dwdx(1,k) - dvx(2)*parts%dwdx(2,k)   &
                   - dvx(3)*parts%dwdx(3,k) 
             hxy = dvx(1)*parts%dwdx(2,k) + dvx(2)*parts%dwdx(1,k)
             hxz = dvx(1)*parts%dwdx(3,k) + dvx(3)*parts%dwdx(1,k)          
             hyy = 2.e0*dvx(2)*parts%dwdx(2,k) - dvx(1)*parts%dwdx(1,k)   &    
                   - dvx(3)*parts%dwdx(3,k)
             hyz = dvx(2)*parts%dwdx(3,k) + dvx(3)*parts%dwdx(2,k)
             hzz = 2.e0*dvx(3)*parts%dwdx(3,k) - dvx(1)*parts%dwdx(1,k)   &
                   - dvx(2)*parts%dwdx(2,k)
          endif                              
          hxx = 2.e0/3.e0*hxx
          hyy = 2.e0/3.e0*hyy
          hzz = 2.e0/3.e0*hzz
          if (dim.eq.1) then 
             parts%txx(i) = parts%txx(i) + parts%mass(j)*hxx/parts%rho(j)
             parts%txx(j) = parts%txx(j) + parts%mass(i)*hxx/parts%rho(i)               
          else if (dim.eq.2) then           
             local_xx(i,it) = local_xx(i,it) + parts%mass(j)*hxx/parts%rho(j)
             local_xx(j,it) = local_xx(j,it) + parts%mass(i)*hxx/parts%rho(i)   
             local_xy(i,it) = local_xy(i,it) + parts%mass(j)*hxy/parts%rho(j)
             local_xy(j,it) = local_xy(j,it) + parts%mass(i)*hxy/parts%rho(i)            
             local_yy(i,it) = local_yy(i,it) + parts%mass(j)*hyy/parts%rho(j)
             local_yy(j,it) = local_yy(j,it) + parts%mass(i)*hyy/parts%rho(i) 
!             parts%txx(i) = parts%txx(i) + parts%mass(j)*hxx/parts%rho(j)
!             parts%txx(j) = parts%txx(j) + parts%mass(i)*hxx/parts%rho(i)   
!             parts%txy(i) = parts%txy(i) + parts%mass(j)*hxy/parts%rho(j)
!             parts%txy(j) = parts%txy(j) + parts%mass(i)*hxy/parts%rho(i)            
!             parts%tyy(i) = parts%tyy(i) + parts%mass(j)*hyy/parts%rho(j)
!             parts%tyy(j) = parts%tyy(j) + parts%mass(i)*hyy/parts%rho(i)          
          else if (dim.eq.3) then
!             txx(i) = txx(i) + mass(j)*hxx/rho(j)
!             txx(j) = txx(j) + mass(i)*hxx/rho(i)   
!             txy(i) = txy(i) + mass(j)*hxy/rho(j)
!             txy(j) = txy(j) + mass(i)*hxy/rho(i) 
!             txz(i) = txz(i) + mass(j)*hxz/rho(j)
!             txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
!             tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
!             tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
!             tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
!             tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
!             tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
!             tzz(j) = tzz(j) + mass(i)*hzz/rho(i)                 
          endif                              
      enddo
     enddo
     !$omp end do
     !$omp barrier
     !$omp do private(it)
     do i =1,ntotal
         do it = 1, parts%nthreads
             parts%txx(i) = parts%txx(i) + local_xx(i,it)
             parts%txy(i) = parts%txy(i) + local_xy(i,it)
             parts%tyy(i) = parts%tyy(i) + local_yy(i,it)
         enddo
     enddo
     !$omp end do
     !$omp end parallel

      return
      end subroutine      

      
!---------------------------------------------------------------------
      subroutine Jaumann_rate(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: soil
      real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz, G 
      integer i, j, k, d, dim, ntotal, niac

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac; dim = parts%dim
      soil => parts%material
      G = soil%e/(2.0*(1+soil%niu))

!     Hook's law

      do i = 1, ntotal
         parts%dsxx(i) = parts%dsxx(i)+G*parts%txx(i)   ! No accumulation origionaly
         parts%dsxy(i) = parts%dsxy(i)+G*parts%txy(i)
         parts%dsyy(i) = parts%dsyy(i)+G*parts%tyy(i)
      enddo

!         if(parts%soil_pressure==2)then
      do i = 1, ntotal
         parts%dp(i)   = parts%dp(i) - soil%k*parts%vcc(i)  !!! simultaneous pressure  
      enddo
!         endif         
     
!     spin tensor

      parts%wxy = 0.d0

      do k=1,niac
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          do d=1,dim
            dvx(d) = parts%vx(d,j) - parts%vx(d,i)
          enddo
          if (dim.eq.1) then 
            !hxx = 0.5e0*dvx(1)*dwdx(1,k)        
          else if (dim.eq.2) then           
            hxy = 0.5e0*(dvx(1)*parts%dwdx(2,k) - dvx(2)*parts%dwdx(1,k))
          else if (dim.eq.3) then
!            hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
!            hxz = dvx(1)*dwdx(3,k) + dvx(3)*dwdx(1,k)          
!            hyz = dvx(2)*dwdx(3,k) + dvx(3)*dwdx(2,k)
          endif                              
          if (dim.eq.1) then 
!            txx(i) = txx(i) + mass(j)*hxx/rho(j)
!            txx(j) = txx(j) + mass(i)*hxx/rho(i)                 
          else if (dim.eq.2) then           
            parts%wxy(i) = parts%wxy(i) + parts%mass(j)*hxy/parts%rho(j)
            parts%wxy(j) = parts%wxy(j) + parts%mass(i)*hxy/parts%rho(i)            
          else if (dim.eq.3) then
!            txy(i) = txy(i) + mass(j)*hxy/rho(j)
!            txy(j) = txy(j) + mass(i)*hxy/rho(i) 
!            txz(i) = txz(i) + mass(j)*hxz/rho(j)
!            txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
!            tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
!            tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
          endif                              
       enddo
   
!   Jaumann rate

      do i = 1, ntotal
         parts%dsxx(i) = parts%dsxx(i)+2.0*parts%sxy(i)*parts%wxy(i)
         parts%dsxy(i) = parts%dsxy(i)-(parts%sxx(i)-parts%syy(i))*parts%wxy(i)
         parts%dsyy(i) = parts%dsyy(i)-2.0*parts%sxy(i)*parts%wxy(i)
      enddo         

      return
      end subroutine
      

!---------------------------------------------------------------------
      subroutine Jaumann_rate_omp(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: soil
      real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz, G 
      real local(parts%maxn,8)
      integer i, j, k, d, dim, ntotal, niac,it

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac; dim = parts%dim
      soil => parts%material
      G = soil%e/(2.0*(1+soil%niu))

!     Hook's law
      !$omp parallel do
      do i = 1, ntotal
         parts%dsxx(i) = parts%dsxx(i)+G*parts%txx(i)   ! No accumulation origionaly
         parts%dsxy(i) = parts%dsxy(i)+G*parts%txy(i)
         parts%dsyy(i) = parts%dsyy(i)+G*parts%tyy(i)
      enddo
      !$omp end parallel do
!         if(parts%soil_pressure==2)then
      do i = 1, ntotal
         parts%dp(i)   = parts%dp(i) - soil%k*parts%vcc(i)  !!! simultaneous pressure  
      enddo
!         endif         
     
!     spin tensor

      parts%wxy = 0.d0
!   call parts%get_niac_start_end
   
     !$omp parallel 
     !$omp do private(i,it,j,d,k,dvx,hxy)    
     do it = 1,parts%nthreads
         do i =1,ntotal
             local(i,it) = 0.
         enddo
       do k=parts%niac_start(it),parts%niac_end(it)

          i = parts%pair_i(k)
          j = parts%pair_j(k)
          do d=1,dim
            dvx(d) = parts%vx(d,j) - parts%vx(d,i)
          enddo
          if (dim.eq.1) then 
            !hxx = 0.5e0*dvx(1)*dwdx(1,k)        
          else if (dim.eq.2) then           
            hxy = 0.5e0*(dvx(1)*parts%dwdx(2,k) - dvx(2)*parts%dwdx(1,k))
          else if (dim.eq.3) then
!            hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
!            hxz = dvx(1)*dwdx(3,k) + dvx(3)*dwdx(1,k)          
!            hyz = dvx(2)*dwdx(3,k) + dvx(3)*dwdx(2,k)
          endif                              
          if (dim.eq.1) then 
!            txx(i) = txx(i) + mass(j)*hxx/rho(j)
!            txx(j) = txx(j) + mass(i)*hxx/rho(i)                 
          else if (dim.eq.2) then     
              local(i,it) = local(i,it) + parts%mass(j)*hxy/parts%rho(j)
              local(j,it) = local(j,it) + parts%mass(i)*hxy/parts%rho(i)    
!            parts%wxy(i) = parts%wxy(i) + parts%mass(j)*hxy/parts%rho(j)
!            parts%wxy(j) = parts%wxy(j) + parts%mass(i)*hxy/parts%rho(i)            
          else if (dim.eq.3) then
!            txy(i) = txy(i) + mass(j)*hxy/rho(j)
!            txy(j) = txy(j) + mass(i)*hxy/rho(i) 
!            txz(i) = txz(i) + mass(j)*hxz/rho(j)
!            txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
!            tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
!            tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
          endif                              
       enddo
      enddo
       !$omp end do
       !$omp barrier
       !$omp do private(it,i)
       do i = 1,ntotal
          do it = 1,parts%nthreads
             parts%wxy(i) = parts%wxy(i) + local(i,it)
          enddo
       enddo
      !$omp end do
      !$omp end parallel
   
!   Jaumann rate
      !$omp parallel do
      do i = 1, ntotal
         parts%dsxx(i) = parts%dsxx(i)+2.0*parts%sxy(i)*parts%wxy(i)
         parts%dsxy(i) = parts%dsxy(i)-(parts%sxx(i)-parts%syy(i))*parts%wxy(i)
         parts%dsyy(i) = parts%dsyy(i)-2.0*parts%sxy(i)*parts%wxy(i)
      enddo         
      !$omp end parallel do
      return
      end subroutine

!------------------------------------------------------------
      subroutine mohr_coulomb_failure_criterion(soil)
!------------------------------------------------------------
      implicit none

      class(particles) soil
      type(material), pointer :: property
      double precision yield, phi, skale, cohesion, tmax
      integer i, k, ntotal

      ntotal = soil%ntotal+soil%nvirt
      property => soil%material
      cohesion = property%cohesion
      phi      = property%phi
   
      k = 0
      soil%fail = 0
      do i = 1, ntotal
         tmax  = sqrt(((soil%sxx(i)-soil%syy(i))/2)**2+soil%sxy(i)**2)
         if(tmax<1.e-6)cycle
         yield = cohesion*cos(phi)+soil%p(i)*sin(phi)

         if(yield<=0.)then    ! <
            yield=0.; soil%p(i)=-cohesion*tan(phi)**(-1.0)
            soil%rho(i) = property%rho0
         endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         if(tmax>yield)then
!         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = yield/tmax
            soil%sxx(i) = skale * soil%sxx(i)
            soil%sxy(i) = skale * soil%sxy(i)
            soil%syy(i) = skale * soil%syy(i)
         endif
         !sxx(i) = 0.; sxy(i) = 0.; syy(i) = 0.
      enddo
      soil%nfail = k

      return
      end subroutine      

!------------------------------------------------------------
      subroutine mohr_coulomb_failure_criterion_omp(soil)
!------------------------------------------------------------
      implicit none

      class(particles) soil
      type(material), pointer :: property
      double precision yield, phi, skale, cohesion, tmax
      integer i, k, ntotal

      ntotal = soil%ntotal+soil%nvirt
      property => soil%material
      cohesion = property%cohesion
      phi      = property%phi
   
      k = 0
      soil%fail = 0
     !$omp parallel do private(tmax,yield,skale) reduction(+:k)
      do i = 1, ntotal
         tmax  = sqrt(((soil%sxx(i)-soil%syy(i))/2)**2+soil%sxy(i)**2)
         if(tmax<1.e-6)cycle
         yield = cohesion*cos(phi)+soil%p(i)*sin(phi)

         if(yield<=0.)then    ! <
            yield=0.; soil%p(i)=-cohesion*tan(phi)**(-1.0)
            soil%rho(i) = property%rho0
         endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         if(tmax>yield)then
!         if(yield>0.and.tmax>yield)then

            k = k + 1

            soil%fail(i) = 1 
            skale = yield/tmax
            soil%sxx(i) = skale * soil%sxx(i)
            soil%sxy(i) = skale * soil%sxy(i)
            soil%syy(i) = skale * soil%syy(i)
         endif
         !sxx(i) = 0.; sxy(i) = 0.; syy(i) = 0.
      enddo
      !$omp end parallel do

      soil%nfail = k

      return
      end subroutine      
      
!------------------------------------------------------------
      subroutine drucker_prager_failure_criterion(soil)
!------------------------------------------------------------
      implicit none

      class(particles) soil
      type(material), pointer :: property
      double precision yield, phi, skale, cohesion, tmax, alpha1,I1,J2
      integer i, k, ntotal
        
      !if(soil%itimestep==82)write(*,*) 'in drucker...'
 
      ntotal = soil%ntotal+soil%nvirt
      property => soil%material
      cohesion = property%cohesion
      phi      = property%phi
      alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
   
      k = 0
      soil%fail = 0
      do i = 1, ntotal
         !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         !if(tmax<1.e-6)cycle
         !yield = cohesion*cos(phi)+p(i)*sin(phi)


         if(soil%p(i)<0.)then    ! <
            !k = k + 1
            !soil%fail(i) = 1
            soil%p(i)=0.d0
            !soil%rho(i) = property%rho0
         endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         I1 = 3.*soil%p(i) 

         J2 = soil%sxx(i)**2.+2.*soil%sxy(i)**2.+soil%syy(i)**2.+(soil%sxx(i)+soil%syy(i))**2.
         J2 = sqrt(J2/2.)+1.d-6
         !if(J2<1.e-6)cycle

         if(J2>alpha1*I1)then
!         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = alpha1*I1/J2
            soil%sxx(i) = skale * soil%sxx(i)
            soil%sxy(i) = skale * soil%sxy(i)
            soil%syy(i) = skale * soil%syy(i)
        endif

      enddo
      soil%nfail = k

      return
      end subroutine

!------------------------------------------------------------
      subroutine drucker_prager_failure_criterion_omp(soil)
!------------------------------------------------------------
      implicit none

      class(particles) soil
      type(material), pointer :: property
      double precision yield, phi, skale, cohesion, tmax, alpha1,I1,J2
      integer i, k, ntotal
        
      !if(soil%itimestep==82)write(*,*) 'in drucker...'
 
      ntotal = soil%ntotal+soil%nvirt
      property => soil%material
      cohesion = property%cohesion
      phi      = property%phi
      alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
   
      k = 0
      soil%fail = 0
      !$omp parallel do private(I1,J2,skale) reduction(+:k)
      do i = 1, ntotal
         !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         !if(tmax<1.e-6)cycle
         !yield = cohesion*cos(phi)+p(i)*sin(phi)


         if(soil%p(i)<0.)then    ! <
            !k = k + 1
            !soil%fail(i) = 1
            soil%p(i)=0.d0
            !soil%rho(i) = property%rho0
         endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         I1 = 3.*soil%p(i) 

         J2 = soil%sxx(i)**2.+2.*soil%sxy(i)**2.+soil%syy(i)**2.+(soil%sxx(i)+soil%syy(i))**2.
         J2 = sqrt(J2/2.)+1.d-6
         !if(J2<1.e-6)cycle

         if(J2>alpha1*I1)then
!         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = alpha1*I1/J2
            soil%sxx(i) = skale * soil%sxx(i)
            soil%sxy(i) = skale * soil%sxy(i)
            soil%syy(i) = skale * soil%syy(i)
        endif

      enddo
      !$omp end parallel do
      soil%nfail = k

      return
      end subroutine
      
!------------------------------------------------------------
      subroutine plastic_or_not(soil)
!------------------------------------------------------------
      implicit none

      class(particles) soil
      type(material), pointer :: property
      real(dp) yield, phi, skale, cohesion, tmax, alpha1,I1,J2
      integer i, k, ntotal
        
      !if(soil%itimestep==82)write(*,*) 'in drucker...'

      ntotal = soil%ntotal+soil%nvirt
      property => soil%material
      cohesion = property%cohesion
      phi      = property%phi
      alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
   
      k = 0
      soil%fail = 0
      do i = 1, ntotal
         !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         !if(tmax<1.e-6)cycle
         !yield = cohesion*cos(phi)+p(i)*sin(phi)


         !if(p(i)<0.)then    ! <
         !   !k = k + 1
         !   !soil%fail(i) = 1
         !   p(i)=0.d0
         !   !soil%rho(i) = property%rho0
         !endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         I1 = 3.*soil%p(i) 

         J2 = soil%sxx(i)**2.+2.*soil%sxy(i)**2.+soil%syy(i)**2.+(soil%sxx(i)+soil%syy(i))**2.
         J2 = sqrt(J2/2.)+1.d-6
         !if(J2<1.e-6)cycle

!         if(J2>alpha1*I1)then
         if(abs(J2-alpha1*I1)<1.0e-5)then
            k = k + 1
            soil%fail(i) = 1 
            !skale = alpha1*I1/J2
            !sxx(i) = skale * sxx(i)
            !sxy(i) = skale * sxy(i)
            !syy(i) = skale * syy(i)
        endif

      enddo
      soil%nfail = k

      return
      end subroutine

!------------------------------------------------------------
      subroutine plastic_or_not_omp(soil)
!------------------------------------------------------------
      implicit none

      class(particles) soil
      type(material), pointer :: property
      real(dp) yield, phi, skale, cohesion, tmax, alpha1,I1,J2
      integer i, k, ntotal
        
      !if(soil%itimestep==82)write(*,*) 'in drucker...'

      ntotal = soil%ntotal+soil%nvirt
      property => soil%material
      cohesion = property%cohesion
      phi      = property%phi
      alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
   
      k = 0
      soil%fail = 0
      !$omp parallel do private(I1,J2) reduction(+:k)
      do i = 1, ntotal
         !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         !if(tmax<1.e-6)cycle
         !yield = cohesion*cos(phi)+p(i)*sin(phi)


         !if(p(i)<0.)then    ! <
         !   !k = k + 1
         !   !soil%fail(i) = 1
         !   p(i)=0.d0
         !   !soil%rho(i) = property%rho0
         !endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         I1 = 3.*soil%p(i) 

         J2 = soil%sxx(i)**2.+2.*soil%sxy(i)**2.+soil%syy(i)**2.+(soil%sxx(i)+soil%syy(i))**2.
         J2 = sqrt(J2/2.)+1.d-6
         !if(J2<1.e-6)cycle

!         if(J2>alpha1*I1)then
         if(abs(J2-alpha1*I1)<1.0e-5)then
            k = k + 1
            soil%fail(i) = 1 
            !skale = alpha1*I1/J2
            !sxx(i) = skale * sxx(i)
            !sxy(i) = skale * sxy(i)
            !syy(i) = skale * syy(i)
        endif

      enddo
      !$omp end parallel do
      soil%nfail = k

      return
      end subroutine

!---------------------------------------------------------------------
      subroutine plastic_flow_rule(parts)
!---------------------------------------------------------------------
      implicit none

      class(particles) parts
      double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
      double precision, pointer, dimension(:) :: sxx, sxy, syy
      type(material), pointer :: property
      double precision alpha, phi, K, G,e,niu, J2, sde, dlambda
      double precision exx, exy, eyy                ! total strain rate
      double precision :: small_value = 1.d-10
      integer i, ntotal

      ntotal = parts%ntotal + parts%nvirt

      dsxx => parts%dsxx
      dsxy => parts%dsxy
      dsyy => parts%dsyy
      vcc  => parts%vcc
      sxx => parts%sxx
      sxy => parts%sxy
      syy => parts%syy

      property => parts%material
      phi = property%phi
      k   = property%k; e = property%E; niu = property%niu
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      G = e/(2.0*(1+niu))
      
      do i = 1, ntotal
                            !if(parts%fail(i)==1)then

      exx = parts%txx(i)/2.+parts%vcc(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%txy(i)/2.
      eyy = parts%tyy(i)/2.+parts%vcc(i)/3.

      sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
      J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
      J2 = J2 + small_value

      !G = parts%eta(i)
      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
      !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
      if(dlambda<0)cycle
      !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
      dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
      dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
      !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
      dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

      !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
      parts%dp(i) = parts%dp(i) + 3.*k*alpha*dlambda

! Accumulative deviatoric strain
 
      parts%epsilon_p(i) = parts%epsilon_p(i)    &
                         + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            !endif ! Fail
      enddo
      
      return
      end subroutine

!---------------------------------------------------------------------
      subroutine plastic_flow_rule2(parts)
!---------------------------------------------------------------------
      implicit none

      class(particles) parts
      double precision, pointer, dimension(:) :: dsxx2,dsxy2,dsyy2,vcc
      double precision, pointer, dimension(:) :: sxx, sxy, syy
      type(material), pointer :: property
      double precision alpha, phi, K, G,e,niu, J2, sde, dlambda
      double precision exx, exy, eyy                ! total strain rate
      double precision :: small_value = 1.d-10
      integer i, ntotal

      ntotal = parts%ntotal + parts%nvirt

      dsxx2 => parts%dsxx2
      dsxy2 => parts%dsxy2
      dsyy2 => parts%dsyy2
      vcc  => parts%vcc
      sxx => parts%sxx
      sxy => parts%sxy
      syy => parts%syy

      property => parts%material
      phi = property%phi
      k   = property%k; e = property%E; niu = property%niu
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      G = e/(2.0*(1+niu))
      
      do i = 1, ntotal
                            !if(parts%fail(i)==1)then

      exx = parts%txx(i)/2.+parts%vcc(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%txy(i)/2.
      eyy = parts%tyy(i)/2.+parts%vcc(i)/3.

      sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
      J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
      J2 = J2 + small_value

      !G = parts%eta(i)
      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)

      !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
      dsxx2(i) = -dlambda*(G/sqrt(J2)*sxx(i))
      dsxy2(i) = -dlambda*(G/sqrt(J2)*sxy(i))
      !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
      dsyy2(i) = -dlambda*(G/sqrt(J2)*syy(i))

      !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
      parts%dp2(i) =  3.*k*alpha*dlambda

! Accumulative deviatoric strain
 
!      parts%epsilon_p(i) = parts%epsilon_p(i) 
!     &                   + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            !endif ! Fail
      enddo
      
      return
      end subroutine

!---------------------------------------------------------------------
      subroutine plastic_flow_rule3(parts)
!---------------------------------------------------------------------
      implicit none

      class(particles) parts
      double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
      double precision, pointer, dimension(:) :: sxx, sxy, syy
      type(material), pointer :: property
      double precision alpha, phi, K, G,e,niu, J2, sde, dlambda
      double precision exx, exy, eyy                ! total strain rate
      double precision :: small_value = 1.d-10
      integer i, ntotal

      ntotal = parts%ntotal + parts%nvirt

      dsxx => parts%dsxx
      dsxy => parts%dsxy
      dsyy => parts%dsyy
      vcc  => parts%vcc
      sxx => parts%sxx
      sxy => parts%sxy
      syy => parts%syy

      property => parts%material
      phi = property%phi
      k   = property%k; e = property%E; niu = property%niu
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      G = e/(2.0*(1+niu))
      
      do i = 1, ntotal
                            if(parts%fail(i)==1)then

      exx = parts%txx(i)/2.+parts%vcc(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%txy(i)/2.
      eyy = parts%tyy(i)/2.+parts%vcc(i)/3.

      sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
      J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
      J2 = J2 + small_value

      !G = parts%eta(i)
      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
      !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
      if(dlambda<0)cycle
      !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
      dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
      dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
      !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
      dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

      !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
      parts%dp(i) = parts%dp(i) + 3.*k*alpha*dlambda

! Accumulative deviatoric strain
 
      parts%epsilon_p(i) = parts%epsilon_p(i)       & 
                         + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            endif ! Fail
      enddo
      
      return
      end subroutine

!---------------------------------------------------------------------
      subroutine plastic_flow_rule3_omp(parts)
!---------------------------------------------------------------------
      implicit none

      class(particles) parts
      double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
      double precision, pointer, dimension(:) :: sxx, sxy, syy
      type(material), pointer :: property
      double precision alpha, phi, K, G,e,niu, J2, sde, dlambda
      double precision exx, exy, eyy                ! total strain rate
      double precision :: small_value = 1.d-10
      integer i, ntotal

      ntotal = parts%ntotal + parts%nvirt

      dsxx => parts%dsxx
      dsxy => parts%dsxy
      dsyy => parts%dsyy
      vcc  => parts%vcc
      sxx => parts%sxx
      sxy => parts%sxy
      syy => parts%syy

      property => parts%material
      phi = property%phi
      k   = property%k; e = property%E; niu = property%niu
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      G = e/(2.0*(1+niu))
      !$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
      do i = 1, ntotal
                            if(parts%fail(i)==1)then

      exx = parts%txx(i)/2.+parts%vcc(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%txy(i)/2.
      eyy = parts%tyy(i)/2.+parts%vcc(i)/3.

      sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
      J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
      J2 = J2 + small_value

      !G = parts%eta(i)
      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
      !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
      if(dlambda<0)cycle
      !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
      dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
      dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
      !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
      dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

      !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
      parts%dp(i) = parts%dp(i) + 3.*k*alpha*dlambda

! Accumulative deviatoric strain
 
      parts%epsilon_p(i) = parts%epsilon_p(i)       & 
                         + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            endif ! Fail
      enddo
      !$omp end parallel do
      
      
      return
      end subroutine
      
!---------------------------------------------------------------
      subroutine darcy_law(water, soil)
!---------------------------------------------------------------
!      use param, only: volume_fraction
!      use m_particles
      implicit none

      type(particles) water, soil
      double precision dx(3), ks, ns, gw, cf, sp, rrw
      type(material), pointer :: h2o,sio2  
      type(numerical), pointer :: numeric
      double precision gravity   
      integer i, j, k, d, dim

      h2o => water%material
      sio2=>  soil%material
      numeric => water%numeric
      gravity = numeric%gravity
      dim = water%dim

      gw = h2o%rho0*(-gravity); ns = sio2%porosity
      ks = sio2%permeability
      cf = gw*ns/ks
      !cf = 6.e6
       
      do  k=1,water%niac
          i = water%pair_i(k)
          j = water%pair_j(k)  
          rrw = water%w(k)/(water%rho(i)*soil%rho(j))

! For staturated soil
!        if(volume_fraction) cf = water%vof(i)*water%rho(i)*(-gravity)/ks
!        if(water%volume_fraction) cf = water%vof(i)*soil%vof(j)*water%rho(i)*(-gravity)/ks
        cf = water%vof(i)*soil%vof(j)*water%rho(i)*(-gravity)/ks

          do d=1,dim
             sp = cf*(water%vx(d,i)-soil%vx(d,j))*rrw
             water%dvx(d,i) = water%dvx(d,i) - soil%mass(j)*sp
             soil%dvx(d,j)  = soil%dvx(d,j) + water%mass(i)*sp   
          enddo
      enddo

      return
      end subroutine 

!-----------------------------------------------------
      subroutine pore_water_pressure(water,soil)
!-----------------------------------------------------
!      use param
!      use m_particles
      implicit none
      
      type(particles) water, soil
      double precision mprr
      integer i, j, k, d
      
      do k = 1, water%niac
         i = water%pair_i(k)   ! water
         j = water%pair_j(k)   ! soil
         mprr = water%mass(i)*water%p(i)/(water%rho(i)*soil%rho(j))
!         mprr = water%mass(i)*(water%p(i)+soil%p(j))/       &      Bui2014
!                (water%rho(i)*soil%rho(j))
         do d = 1, water%dim
            soil%dvx(d,j) = soil%dvx(d,j) + mprr*water%dwdx(d,k)  &  !+
                            *soil%vof(j)  
         enddo

! saturated soil
         !if(water%volume_fraction)then
         !do d = 1, dim
            !water%dvx(d,i) = water%dvx(d,i) -                     & ! Must be -
            !soil%mass(j)*water%p(i)*soil%vof(j)*water%dwdx(d,k)/  &
            !(water%rho(i)*soil%rho(j))
         !enddo
         !endif
      enddo

      return
      end subroutine       

!-------------------------------------------------------------
      subroutine volume_fraction_water(water, soil)
!-------------------------------------------------------------
!      use param
!      use m_particles
      implicit none

      type(particles) water, soil
      integer i,j,k,d, ntotal
      type(material), pointer :: sio2

      sio2 => soil%material
      ntotal = water%ntotal+water%nvirt

      water%vof2 = 0.d0
      do k = 1, water%niac
         i = water%pair_i(k)
         j = water%pair_j(k)
         water%vof2(i) = water%vof2(i)+soil%mass(j)*water%w(k)
      enddo

      do k = 1, water%ntotal+water%nvirt
         water%vof2(k) = 1.d0 - water%vof2(k)/sio2%rho0
      enddo

      return
      end subroutine 

!-------------------------------------------------------------
      subroutine volume_fraction_water2(water, soil)
!-------------------------------------------------------------
!      use param
!      use m_particles
      implicit none

      type(particles) water, soil
      integer i,j,k,d, ntotal
      type(material), pointer :: sio2
      double precision dvx(3),tmp

      sio2 => soil%material
      ntotal = water%ntotal+water%nvirt

!      water%dvof = 0.d0

      do k = 1, water%niac
         i = water%pair_i(k)
         j = water%pair_j(k)
         do d = 1, water%dim
            dvx(d) = water%vx(d,i)-soil%vx(d,j)
         enddo 
         tmp = dvx(1)*water%dwdx(1,k)+dvx(2)*water%dwdx(2,k)
         water%dvof(i) = water%dvof(i)-soil%mass(j)*tmp/sio2%rho0
      enddo

      return
      end subroutine 

!--------------------------------------------------------------------
      subroutine volume_fraction_soil(parts)
!--------------------------------------------------------------------
!      use param
!      use m_particles
      implicit none

      type(particles) parts
      integer i,j,k, ntotal
      type(material), pointer :: sio2

      sio2 => parts%material
      ntotal = parts%ntotal + parts%nvirt

      do i = 1, ntotal
         parts%vof(i) = parts%rho(i)/sio2%rho0
      enddo

      return
      end subroutine

!--------------------------------------------------------------------------
      subroutine drag_force(water,soil)
!--------------------------------------------------------------------------
! For fixed porous media, drag force to fluid  

      !use param
!      use m_particles
      implicit none

      type(particles) water, soil
      double precision dx(3), ks, ns, gw, cf, sp, rrw
      type(material), pointer :: h2o,sio2  
      type(numerical), pointer :: numeric
      double precision gravity   
      integer i, j, k, d, dim

      h2o => water%material
      sio2=>  soil%material
      numeric => water%numeric
      gravity = numeric%gravity
      dim = water%dim

      gw = h2o%rho0*(-gravity) 
      !ns = sio2%porosity; ks = sio2%permeability
      ns = 0.3; ks = 0.05      

      cf = gw*ns/ks
      !cf = 6.e6
      !write(*,*) ns,ks
       
      do i=1,water%ntotal+water%nvirt
         do d=1,dim
      water%dvx(d,i) = water%dvx(d,i) - water%vx(d,i)*cf/water%rho(i)
         enddo
      enddo

      return
      end subroutine 


!---------------------------------------------------------------------
!   function find_particle_nearest2_point(parts,x,y ) result(index)
!---------------------------------------------------------------------
!implicit none
!class(particles) parts
!real(dp) x, y
!integer index
!integer ntotal, i
!real(dp) :: mindist = 10.d10, dist

!ntotal = parts%ntotal + parts%nvirt
!do i = 1, ntotal
!   dist = (parts%x(1,i)-x)**2.0 + (parts%x(2,i)-y)**2.0
!   if(dist<mindist)then
!      index = i
!      mindist = dist
!   endif        
!enddo

!end function

!------------------------------------------------------------------------
!   function find_particles_nearest2_points(parts,x,y ) result(indices)
!------------------------------------------------------------------------
!implicit none
!class(particles) parts
!real(dp):: x(:), y(:)
!integer index
!integer ntotal, i
!real(dp) :: mindist = 10.d10, dist

!ntotal = parts%ntotal + parts%nvirt
!do i = 1, ntotal
!   dist = (parts%x(1,i)-x)**2.0 + (parts%x(2,i)-y)**2.0
!   if(dist<mindist)then
!      index = i
!      mindist = dist
!   endif        
!enddo

!end function

!--------------------------------------------
    subroutine array_equal_array(a,b)
!--------------------------------------------        
implicit none
class(array),intent(INOUT) :: a
class(array),intent(IN)    :: b
integer ndim1

if(a%ndim1/=b%ndim1)stop 'Error: inequal length arrays!'

ndim1 = a%ndim1
a%r(1:ndim1) = b%r(1:ndim1)

return
end subroutine

!--------------------------------------------
       subroutine array_final(this)
!--------------------------------------------
implicit none
type(array) this

if(associated(this%r))deallocate(this%r)

return
end subroutine

!--------------------------------------------
    function array_add_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a, b
type(array), allocatable:: c
integer ndim1

if(a%ndim1/=b%ndim1)stop 'Cannot add arrays!'
if(.not.associated(a%parts,b%parts)) stop 'Cannot add arrays'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

c%r(1:ndim1) = a%r(1:ndim1) + b%r(1:ndim1)
c%parts => a%parts

end function

!--------------------------------------------
    function array_sub_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a, b
type(array), allocatable:: c
integer ndim1

if(a%ndim1/=b%ndim1)stop 'Cannot substract arrays!'
if(.not.associated(a%parts,b%parts)) stop 'Cannot substract arrays!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

c%r(1:ndim1) = a%r(1:ndim1) - b%r(1:ndim1)
c%parts => a%parts

end function

!--------------------------------------------
    function array_minus(a) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
type(array), allocatable:: c
integer ndim1

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

c%r(1:ndim1) = -a%r(1:ndim1)
c%parts => a%parts

end function

!--------------------------------------------
    function array_mul_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a,b
type(array), allocatable:: c
integer ndim1

if(a%ndim1/=b%ndim1)stop 'Cannot multiply arrays!'
if(.not.associated(a%parts,b%parts)) stop 'Cannot multiply arrays!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

c%r(1:ndim1) = a%r(1:ndim1)*b%r(1:ndim1)
c%parts => a%parts

end function

!--------------------------------------------
    function array_mul_real(a,r) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: r
type(array), allocatable:: c
integer ndim1

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

c%r(1:ndim1) = r*a%r(1:ndim1)
c%parts => a%parts

end function

!--------------------------------------------------
    function array_mul_double_real(a,r) result(c)
!--------------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: r
type(array), allocatable :: c
integer ndim1

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

c%r(1:ndim1) = r*a%r(1:ndim1)
c%parts => a%parts

end function

!--------------------------------------------------
    subroutine get_num_threads(this)
!--------------------------------------------------
implicit none
class(particles) this
integer i,j,k
integer, external :: OMP_GET_NUM_THREADS

!$omp parallel
this%nthreads = OMP_GET_NUM_THREADS()
!$omp end parallel

allocate(this%niac_start(this%nthreads),this%niac_end(this%nthreads))

end subroutine

!--------------------------------------------------
    subroutine get_niac_start_end(this)
!--------------------------------------------------
implicit none
class(particles) this
integer i,j,k,it

i = this%niac/this%nthreads
j = mod(this%niac,this%nthreads)

do it = 1, this%nthreads
   this%niac_start(it) = (it-1)*i + 1
   this%niac_end(it) = it*i
enddo

this%niac_end(this%nthreads) = this%niac_end(this%nthreads)+j

end subroutine

end module
