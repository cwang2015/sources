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

!---------------------------------------------------------
type parameters

!Math parameters
real(dp) :: pi = 3.14159265358979323846

!Physics parameter
real(dp) :: gravity = -9.8

!Geometry parameters

!dim : Dimension of the problem (1, 2 or 3)
integer :: dim = 2

!Parameters for the computational geometry,  
!x_maxgeom : Upper limit of allowed x-regime 
!x_mingeom : Lower limit of allowed x-regime 
!y_maxgeom : Upper limit of allowed y-regime 
!y_mingeom : Lower limit of allowed y-regime 
!z_maxgeom : Upper limit of allowed z-regime 
!z_mingeom : Lower limit of allowed z-regime 
real(dp) :: x_maxgeom = 10.e0, x_mingeom = -10.e0,  &
            y_maxgeom = 10.e0, y_mingeom = -10.e0,  &
            z_maxgeom = 10.e0, z_mingeom = -10.e0

!Memory allocation

!maxn: Maximum number of particles
!max_interation : Maximum number of interaction pairs
integer :: maxn = 12000, max_interaction = 100 * 12000
  
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

integer :: integrate_scheme = 1  ! =1, LF; =2, Verlet
logical :: summation_density  = .false.         
logical :: average_velocity  = .true.         
logical :: config_input  = .false. 
logical :: virtual_part  = .true. 
logical :: vp_input  = .false.  
logical :: visc  = .true.  
logical :: ex_force  = .true.
logical :: visc_artificial  = .true. 
logical :: heat_artificial  = .false. 
logical :: self_gravity  = .true.      
logical :: nor_density  = .false.              

integer :: soil_pressure = 2  ! =1, eos; =2, mean trace
integer :: stress_integration = 1
integer :: yield_criterion = 2
integer :: plasticity = 3  ! =0 non; =1 Bui, =2 return mapping =3 Lopez
logical :: artificial_density = .true.                  
logical :: soil_artificial_stress = .true.

logical :: volume_fraction = .true.
logical :: water_artificial_volume = .true.
logical :: volume_fraction_renorm = .true.

! 0 ignor; 1 negative pressure to zero; 2 artficial stress
integer :: water_tension_instability = 0

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
          
!Simulation cases

!shocktube = .true. : carry out shock tube simulation
!shearcavity = .true. : carry out shear cavity simulation
!logical shocktube, shearcavity, waterjet
logical :: shocktube  = .false. 
logical :: shearcavity  = .false. 
logical :: waterjet = .true.

!Recorde time interval
integer :: save_step_from = 0, save_step_to = 100

!Material parameters

!For equation of state (EOS) of water

!real(dp) rho0_f,b,gamma,cf,viscosity

!For soil

!real(dp) rho0_s,cs,k,porosity,permeability,G,E,niu,cohesion,phi
   
! Artificial viscosity
! alpha: shear viscosity; beta: bulk viscosity; etq: to avoid sigularities   
!real(dp) :: alpha=0.1d0, beta=0.d0, etq=0.1d0

! Leonard_Johns repulsive force
!real(dp) :: rr0 = 0.005d0, dd = 0.1d0, p1 = 12, p2 = 4

! Delta-SPH
!real(dp) :: delta = 0.1d0
        
end type parameters 
!-----------------------------------------------------------------

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
   real(dp) :: alpha=0.d0, beta=0.d0, etq=0.1d0

! Leonard_Johns repulsive force
   real(dp) :: rr0 = 0.005d0, dd = 0.1d0, p1 = 12, p2 = 4

! Delta-SPH
   real(dp) :: delta = 0.1d0
   
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

   integer :: dim   = 2
   integer :: maxn  = 0
!   integer :: nnps  = 1
   integer :: niac  = 0
   integer :: max_interaction = 0
   integer :: ntotal = 0, nvirt = 0
   integer :: max_zone = 10

   class(parameters), pointer :: params => null()

   character(len=32) :: imaterial
   class(*), pointer :: material => null()

! Background point
  type(block), pointer :: backpoint => null()
  integer nreal_zone, nvirtual_zone
  integer, pointer, dimension(:) :: real_zone => null(), virtual_zone => null()

! Numerical parameters
   class(numerical), pointer :: numeric => null()
   
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
   real(dp), pointer, dimension(:)   :: eta  => null()   
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

! Kernel and its derivative
   real(dp), pointer, dimension(:)   :: w    => null()
   real(dp), pointer, dimension(:,:) :: dwdx => null()

! Particle interaction pair
   integer :: maxp = 0, minp = 0
   integer :: maxiac = 0, miniac = 0
   integer :: sumiac = 0, noiac  = 0
   integer, pointer, dimension(:)  :: pair_i => null()
   integer, pointer, dimension(:)  :: pair_j => null()
   integer, pointer, dimension(:)  :: countiac=>null()

! Boundry particles defined as type particles
   type(particles), pointer :: bor

! Multiphase flow
!   type(particles), pointer :: water => null(), soil => null()

   integer itimestep

   contains

       procedure :: write_particles
       procedure :: write_component 
       procedure :: append_virtual_particles
       procedure :: depend_virtual_particles
       procedure :: interaction_statistics
       procedure :: minimum_time_step
       procedure :: take_real_points
       procedure :: take_real => take_real_points1
       procedure :: take_virtual_points
       procedure :: take_virtual => take_virtual_points1
       procedure :: setup_itype
       procedure :: find_pairs
       procedure :: repulsive_force

!       procedure :: find_particle_nearest2_point
!       procedure :: grad_scalar
!       procedure :: grad_tensor
!       generic :: grad => grad_scalar, grad_tensor
 
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

!-------------------------------------------------
     subroutine interaction_statistics(parts)
!-------------------------------------------------
implicit none
class(particles) parts

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
!   call link_list(parts)
   call link_list(parts%itimestep, parts%ntotal+parts%nvirt,parts%hsml(1),  &
   parts%x,parts%niac,parts%pair_i,parts%pair_j,parts%w,parts%dwdx,parts%countiac)
!   stop 'find_pairs: link_list method not implemented yet!'
elseif(nnps == 3)then
!    call tree_search(parts)
   stop 'find_pairs: tree_search method not implemented yet!'
endif

return
end subroutine

!--------------------------------------------------------------------------
      subroutine repulsive_force(parts)
!--------------------------------------------------------------------------
implicit none

class(particles) parts
double precision dx(3), rr, f, rr0, dd, p1, p2     
integer i, j, k, d
           
rr0 = parts%numeric%rr0; dd = parts%numeric%dd
p1 = parts%numeric%p1; p2 = parts%numeric%p2
      
do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)  
   if(parts%itype(i).gt.0.and.parts%itype(j).lt.0)then  
      rr = 0.      
      do d=1,parts%dim
         dx(d) =  parts%x(d,i) -  parts%x(d,j)
         rr = rr + dx(d)*dx(d)
      enddo  
      rr = sqrt(rr)
      if(rr.lt.rr0)then
         f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
         do d = 1, parts%dim
            parts%dvx(d, i) = parts%dvx(d, i) + dd*dx(d)*f
         enddo
      endif
   endif        
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

end module
