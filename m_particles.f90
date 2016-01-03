!---------------------------
    module m_particles
!---------------------------
use m_array
implicit none
!integer, parameter :: dp = kind(0.d0)

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
   real(dp) :: alpha=0.1d0, beta=0.d0, etq=0.1d0 !这里的alpha,beta的值值得商榷，在某一篇文献中，好像看到了类似的取值，但在LIU的书中，写的是1.0左右,Monaghan建议的值是0.1和0.0
!   real(dp) :: alpha=1.d0, beta=1d0, etq=0.1d0

! Leonard_Johns repulsive force
!   real(dp) :: dd = 0.1d0, p1 = 12, p2 = 4   !原来的数据,dd一般取最大速度的平方
   real(dp) :: dd = 1d0, p1 = 12, p2 = 4   !原来的数据,dd一般取最大速度的平方
!real(dp) :: dd = 0.01d0, p1 = 12, p2 = 4
!real(dp) :: dd = 0.001d0, p1 = 12, p2 = 4     !可以用的
!real(dp) :: dd = 0.0005d0, p1 = 12, p2 = 4
!   real(dp) :: dd = 1d0, p1 = 12, p2 = 4 
! real(dp) :: dd = 0.5d0, p1 = 12, p2 = 4 
 
! Delta-SPH
   real(dp) :: delta = 0.1d0,delta_alpha=0.1d0

! Velocity average epsilon是一个常数，通过施加临近粒子的影响使自身的运动速度与临近粒子的平均速度相近，在“改进的物理粘性SPH方法及其在溃坝问题中的应用”中的取值为0.3，书上的推荐值也是
!  real(dp) :: epsilon = 0.001
!  real(dp) :: epsilon = 0.1
!  real(dp) :: epsilon = 0.01
!   real(dp) :: epsilon = 0.001
!   real(dp) :: epsilon = 0.0003  !yuan
!real(dp) :: epsilon = 0.00005
!   real(dp) :: epsilon = 0.00001
!real(dp) :: epsilon = 0.000001
!   real(dp) :: epsilon = 0.0003   
!   real(dp) :: epsilon = 0.0003
   real(dp) :: epsilon = 0.0003
!   real(dp) :: epsilon = 0.0005
!  real(dp) :: epsilon = 0.000001
   
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

!real(dp) :: x_maxgeom = 38.5e0, x_mingeom = -0.5e0,  &
!            y_maxgeom = 2.10e0, y_mingeom = -0.5e0,  &
!            z_maxgeom = 10.e0, z_mingeom = -10.e0

real(dp) :: x_maxgeom = 3.85e0, x_mingeom = -0.05e0,  &
            y_maxgeom = 2.10e0, y_mingeom = -0.05e0,  &
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
integer :: maxn = 30000, max_interaction = 10 * 30000
  
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

   real(dp) res_F
   real(dp), pointer, dimension(:)   :: w    => null()
   real(dp), pointer, dimension(:)   :: ps    => null()
   real(dp), pointer, dimension(:)   :: rhos => null()
   real(dp), pointer, dimension(:)   :: gammaas => null()
   real(dp), pointer, dimension(:,:) :: dwdx => null()
   real(dp), pointer, dimension(:,:) :: n => null()
   real(dp), pointer, dimension(:,:) :: mons => null()
   real(dp), pointer, dimension(:,:) :: dgu => null()
! Particle fundamental data
   integer,  pointer, dimension(:)   :: itype=> null()
   real(dp), pointer, dimension(:,:) :: x    => null()
   type(array), pointer :: vol  => null()
   type(array), pointer :: mass => null()
   real(dp), pointer, dimension(:)   :: hsml => null()
   integer,  pointer, dimension(:)   :: zone => null()

! Volume of Fraction
   type(array), pointer :: vof  => null()
   type(array), pointer :: vof2=>null() ! phi_f= 1-phi_s
   type(array), pointer :: dvof  => null()
   type(array), pointer :: vof_min  => null()

! Field variables
   !real(dp), pointer, dimension(:)   :: rho  => null()
   type(array), pointer :: rho  => null()
!   type(array), pointer :: rhos  => null()
   type(array), pointer :: divvx  => null()
   !real(dp), pointer, dimension(:,:) :: vx   => null()   
   type(array), pointer :: vx   => null()   
   type(array), pointer :: mone   => null()   
!   type(array), pointer :: mons   => null()   
   !real(dp), pointer, dimension(:)   :: p    => null()
   type(array), pointer :: p => null()
   type(array), pointer :: gammaa => null()
!   type(array), pointer :: ps => null()
   !real(dp), pointer, dimension(:)   :: eta  => null()   
!   real(dp), pointer, dimension(:)   :: c    => null()   
   type(array), pointer :: c => null()
!   real(dp), pointer, dimension(:)   :: u    => null()

! Stress tensor 
   type(array), pointer :: str => null()
   !type(tensor), pointer, dimension(:)   :: st   => null() 
   !real(dp), pointer, dimension(:) :: sxx => null(), syy => null()
!   type(array), pointer :: sxx => null(), syy => null()
   !real(dp), pointer, dimension(:) :: szz => null(), sxy => null()
!   type(array), pointer :: szz => null(), sxy => null()
!   real(dp), pointer, dimension(:) :: sxz => null(), syz => null()
 
! Shear strain rate
   type(array), pointer :: tab => null()
!   type(array), pointer :: txx => null(), tyy => null()
!   type(array), pointer :: tzz => null(), txy => null()
!   type(array), pointer :: txz => null(), tyz => null()
!   real(dp), pointer, dimension(:) :: txx => null(), tyy => null()
!   real(dp), pointer, dimension(:) :: tzz => null(), txy => null()
!   real(dp), pointer, dimension(:) :: txz => null(), tyz => null()
! Bulk strain rate, i.e. divergency of velocity
!   real(dp), pointer, dimension(:) :: vcc => null()
   type(array), pointer :: vcc => null()
        
! Accumulative deviatoric strain
!   real(dp), pointer, dimension(:) :: epsilon_p => null()
   type(array), pointer :: epsilon_p => null()
! Rotation tensor
!   real(dp), pointer, dimension(:) :: wxy => null()
   type(array), pointer :: wxy => null()
! Old values
!   real(dp), pointer, dimension(:)   :: rho_min => null()
   type(array), pointer :: rho_min => null()
   type(array), pointer :: v_min => null()
!   real(dp), pointer, dimension(:,:) :: v_min   => null()
!   real(dp), pointer, dimension(:)   :: u_min   => null()
!   type(array), pointer :: u_min => null()
!   real(dp), pointer, dimension(:)   :: p_min   => null()
   type(array), pointer :: p_min   => null()
   type(array), pointer :: str_min => null()
!   real(dp), pointer, dimension(:) :: sxx_min => null(), syy_min => null()
!   real(dp), pointer, dimension(:) :: szz_min => null(), sxy_min => null()
!   real(dp), pointer, dimension(:) :: sxz_min => null(), syz_min => null()

! Acceleration
   type(array), pointer :: drho => null()
   !real(dp), pointer, dimension(:,:)   :: dvx  => null()
   type(array), pointer :: dvx  => null()
!   real(dp), pointer, dimension(:)     :: du   => null()
!   type(array), pointer :: du   => null()
!   real(dp), pointer, dimension(:)     :: dp   => null()
   type(array), pointer :: dp   => null()
!   real(dp), pointer, dimension(:,:)   :: av   => null()
   type(array), pointer :: av   => null()

!   type(tensor), pointer, dimension(:) :: dst  => null()
   type(array), pointer :: dstr => null()
!   real(dp), pointer, dimension(:) :: dsxx => null(), dsyy => null()
!   real(dp), pointer, dimension(:) :: dszz => null(), dsxy => null()
!   real(dp), pointer, dimension(:) :: dsxz => null(), dsyz => null()

! For plasticity corrector step (return mapping algorithm)
   type(array), pointer :: dstr2 => null()
!   real(dp), pointer, dimension(:) :: dsxx2 => null(), dsyy2 => null()
!   real(dp), pointer, dimension(:) :: dsxy2 => null(),  
!   real(dp), pointer, dimension(:) :: dp2 => null()
   type(array), pointer :: dp2 => null()

! Density gradient
!   type(array), pointer :: drhodx => null()

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
       procedure :: take_real_new_wedgecase
       procedure :: take_virtual => take_virtual_points1
       procedure :: take_virtual_Couette
       procedure :: take_boundary
       procedure :: take_boundary_for_tank
       procedure :: take_boundary_for_tank_jiami
       procedure :: take_boundary_for_tank_slip
       procedure :: take_boundary_for_tank_block
       procedure :: take_boundary_for_tank_block2
       procedure :: take_boundary_for_tank_wedge
       procedure :: take_boundary_for_tank_wedge2
       procedure :: take_boundary_for_tank_wedge3
       procedure :: take_boundary_for_tank_wedge_freesurface
       procedure :: take_boundary_for_tank_static
       procedure :: setup_itype
       procedure :: setup_ndim1
       procedure :: get_scale_k
       procedure :: direct_find
       procedure :: direct_find_Couette
       procedure :: init_grid
       procedure :: grid_geom
       procedure :: link_list
       procedure :: link_list2
       procedure :: find_pairs
       procedure :: kernel
       procedure :: initial_density
       procedure :: pressure
       procedure :: pressure_water
       procedure :: pressure_soil
       procedure :: pressure_nvirt
       procedure :: pressure_nvirt_unified
       procedure :: tension_instability_water 
       procedure :: sum_density
       procedure :: sum_density_MLS
       procedure :: get_gammaa
       procedure :: nvirt_density_unified
       procedure :: con_density
       procedure :: delta_gamma_unified2
       procedure :: delta_gamma_unified3
       procedure :: delta_gamma_unified4
       procedure :: delta_gamma_unified5
       procedure :: delta_gamma_unified6
       procedure :: real_density_unified
       procedure :: real_density_unified2
       procedure :: freesurface2 !用divr的方式
       procedure :: momentum_equation_unified       
       procedure :: momentum_equation_unified2    
       procedure :: momentum_equation_unified201    
       procedure :: momentum_equation_unified3
       procedure :: momentum_equation_unified4
       procedure :: momentum_equation_unified_Couette
       procedure :: momentum_equation_unified401 !delta p 求压强的方法
       procedure :: momentum_equation_unified402 !delta p和delta delta p 求压强的方法
       procedure :: momentum_equation_unified5
       procedure :: momentum_equation_unified6
       procedure :: freesurface
       procedure :: newtonian_fluid
       procedure :: art_visc
       procedure :: art_visc_omp
       procedure :: art_visc2
       procedure :: delta_sph
       procedure :: delta_sph_omp
       procedure :: delta_rho
       procedure :: delta_sph_vx       
       procedure :: av_vel
       procedure :: repulsive_force
       procedure :: repulsive_force_omp
       procedure :: df
       procedure :: df_omp
       procedure :: df_omp2
       procedure :: df2
       procedure :: df3
       procedure :: df3_omp
       procedure :: df4
       procedure :: df4_omp
       procedure :: psi
       procedure :: div
       procedure :: div_omp
       procedure :: div2
       procedure :: div4                     !师兄认为div4更适合两相流
       procedure :: div5                     !在连续性方程中使用，添加了XSPH 
       procedure :: velocity_divergence
!       procedure :: time_integration_for_water
!       procedure :: find_particle_nearest2_point
!       procedure :: grad_scalar
!       procedure :: grad_tensor
!       generic :: grad => grad_scalar, grad_tensor
      procedure :: shear_strain_rate 
      procedure :: Jaumann_rate
      procedure :: mohr_coulomb_failure_criterion
      procedure :: drucker_prager_failure_criterion
      procedure :: plastic_or_not
      procedure :: plastic_flow_rule
      procedure :: plastic_flow_rule2
      procedure :: plastic_flow_rule3

      procedure :: get_num_threads
      procedure :: get_niac_start_end

end type 

interface write_field
   module procedure :: write_scalar_field
   module procedure :: write_vector2d_field
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
   !write(*,*) i,  this%itype(i), this%x(1,i), this%x(2,i), this%vol(i), this%hsml(i) 
end do

write(*,*) 'Virtual particles:'
write(*,*) 'Num.  itype  x   y   vol  hsml'
do i = this%ntotal+1, this%ntotal+this%nvirt
   !write(*,*) i,  this%itype(i), this%x(1,i), this%x(2,i), this%vol(i), this%hsml(i) 
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
         write(*,*) i, this%mass%r(i)
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
!mint21 = 0.25*minval(sqrt(parts%hsml(1:ntotal)/parts%dvx(1,1:ntotal)))
!mint22 = 0.25*minval(sqrt(parts%hsml(1:ntotal)/parts%dvx(2,1:ntotal)))
! Viscous diffustion time
mint3 = 0.125*minval(parts%hsml(1:ntotal)**2/miu)

!write(*,*) 'Minimum time step is: '
!write(*,*) ' mint1 = ', mint1
!write(*,*) ' mint21 = ', mint21, ' mint22 = ', mint22
!write(*,*) ' mint3 = ', mint3

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
   this%x(1,ntotal+k)  = bor%x(1,i)
   this%x(2,ntotal+k)  = bor%x(2,i)
   this%hsml(ntotal+k) = bor%hsml(i)
   this%itype(ntotal+k)= bor%itype(i)
   this%vol%r(ntotal+k)  = bor%vol%r(i)
   this%mass%r(ntotal+k) = bor%mass%r(i)
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

!--------------------------------------------------
      subroutine take_real_new_wedgecase(this,tank,zone)
!--------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take real particles in tank

k = this%ntotal

   do i = 1, tank%m*tank%n
       if(tank%y(i) > (tank%x(i) - 0.865).and.tank%y(i)<(tank%x(i) - 0.855).and.tank%y(i) < 0.26.and.tank%y(i) > 0.08) then
         k = k + 1
         this%x(1,k) = tank%x(i) - tank%dx /6.
         this%x(2,k) = tank%y(i) - tank%dy * 5./6.
         this%zone(k) = tank%zone(i)
       elseif(tank%y(i) > (tank%x(i) - 0.865).and.tank%y(i)<(tank%x(i) - 0.855).and.tank%y(i) > 0.26.and.tank%y(i) < 0.28)then
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = tank%y(i) - tank%dy *8./9.
         this%zone(k) = tank%zone(i)
       elseif(tank%y(i)> (1.395 - tank%x(i)).and.tank%y(i)<(1.405 - tank%x(i)).and.tank%y(i) < 0.26.and.tank%y(i) >0.08) then
         k = k + 1
         this%x(1,k) = tank%x(i) + tank%dx /6.
         this%x(2,k) = tank%y(i) - tank%dy*5. /6.
         this%zone(k) = tank%zone(i)
       endif     
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


!-----------------------------------------------------
      subroutine take_virtual_Couette(this,tank,zone)
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
         this%x(1,k) = tank%x(i)/2 - 0.0005 
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i)== zone)then             
         k = k + 1
         this%x(1,k) = tank%x(i)/2 !+ 0.00025 
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i)== zone)then             
         k = k + 1
         this%x(1,k) = tank%x(i)/2 + 0.0005 
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i)== zone)then             
         k = k + 1
         this%x(1,k) = tank%x(i)/2 + 0.001 
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
this%nvirt = k-this%ntotal




return
end subroutine

!-----------------------------------------------------
      subroutine take_boundary(this,tank,zone)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == zone)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 0.01 + tank%dy * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

return
end subroutine
      
      
!-----------------------------------------------------
      subroutine take_boundary_for_tank(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
!1         this%x(2,k) = 1.81 - tank%dy * (k - this%ntotal - 1)
         this%x(2,k) = 1.815 - tank%dy * (k - this%ntotal - 1)

         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
!1         this%x(2,k) = 0.01 - tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = 0.005 - tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i) > 0.02.and.tank%x(i)<0.08)then
        k = k + 1
!1         this%x(1,k) = 0.03 + tank%dx * (k - this%ntotal - this%nvirt - 1)   
         this%x(1,k) = 0.015 + tank%dx * (k - this%ntotal - this%nvirt - 1)   
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then         
         k = k + 1
         this%x(2,k) = tank%y(i)
!1         this%x(1,k) = 0.03 + tank%dx * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = 0.015 + tank%dx * (k - this%ntotal - this%nvirt - 1)

         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

!k = this%ntotal+this%nvirt
!   do i = 1, tank%m*tank%n
!      if(tank%zone(i) == 6)then         
!         if(tank%x(i)>3.and.tank%x(i)<3.24)then
!         k = k + 1
!         this%x(2,k) = tank%y(i)    
!         this%x(1,k) = 3.19 + tank%dx * (k - this%ntotal - this%nvirt - 1)
!         this%zone(k) = tank%zone(i)
!         endif
!      endif
!   enddo
   
!this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>3.24)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
!1         this%x(2,k) = 0.01 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = 0.005 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
!1         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = 0.015 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = tank%x(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
      end subroutine


!-----------------------------------------------------
      subroutine take_boundary_for_tank_jiami(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 1.82 - tank%dy/2. * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
         this%x(2,k) = 0.01 - tank%dy/2. * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then         
         k = k + 1
         this%x(2,k) = tank%y(i)
         this%x(1,k) = 0.02 + tank%dx/2. * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>6.5)then
         k = k + 1
         this%x(1,k) = 3.26  
         this%x(2,k) = 0.01 + tank%dy * (k - this%ntotal - this%nvirt - 1)
!2         this%x(2,k) = 0.005 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
!         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = 0.02 + tank%dy/2. * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = 3.26
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
      end subroutine

      
!-----------------------------------------------------
      subroutine take_boundary_for_tank_slip(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3.and.tank%y(i)>0.2)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 0.99 - tank%dy * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3.and.tank%y(i)<0.2)then         
         k = k + 1
         this%x(1,k) = tank%x(i) + tank%dx * (k - this%ntotal - 1)
         this%x(2,k) = 0.2
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
         this%x(2,k) = 0.01 - tank%dy * (k - this%ntotal - this%nvirt - 1)  
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)<0.2)then          
           k = k + 1
           this%x(2,k) = tank%y(i) + tank%dx * (k - this%ntotal - 1)
           this%x(1,k) = 0.03         
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)>0.2)then          
           k = k + 1
           this%x(2,k) = tank%y(i) + tank%dx * (k - this%ntotal - 1)
           this%x(1,k) = 0.2 + tank%dx * (k - this%ntotal - 1) 
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>0.98)then
         k = k + 1
         this%x(2,k) = tank%y(i)    
         this%x(1,k) = 0.99 + tank%dx * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = tank%x(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
end subroutine
      
!-----------------------------------------------------
      subroutine take_boundary_for_tank_block(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3.and.tank%y(i)>0.2)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 0.99 - tank%dy * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3.and.tank%y(i)<0.2)then         
         k = k + 1
         this%x(1,k) = tank%x(i) + tank%dx * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = 0.19
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = 0.19
         this%x(2,k) = 0.19
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)<0.2)then          
           k = k + 1
           this%x(2,k) = 0.17 - tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%x(1,k) = 0.19         
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)>0.2)then          
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 0.21 + tank%dx * (k - this%ntotal - this%nvirt - 1) 
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>0.98)then
         k = k + 1
         this%x(2,k) = tank%y(i)    
         this%x(1,k) = 0.99 + tank%dx * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = tank%x(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
end subroutine

!-----------------------------------------------------
      subroutine take_boundary_for_tank_block2(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 0.99 - tank%dy * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = 0.01
         this%x(2,k) = 0.01
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)<0.4)then          
           k = k + 1
           this%x(2,k) = 0.01
           this%x(1,k) = 0.03 + tank%dx * (k - this%ntotal - this%nvirt - 1)        
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)>0.4.and.tank%x(i)<0.78)then          
           k = k + 1
           this%x(2,k) = 0.01 + tank%dy * (k - this%ntotal - this%nvirt - 1)      
           this%x(1,k) = 0.41         
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)>0.8.and.tank%x(i)<1.)then          
           k = k + 1
           this%x(2,k) = 0.39       
           this%x(1,k) = 0.41 + tank%dx * (k - this%ntotal - this%nvirt - 1)           
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)>1.02.and.tank%x(i)<1.4)then          
           k = k + 1
           this%x(2,k) = 0.37 - tank%dy * (k - this%ntotal - this%nvirt - 1)            
           this%x(1,k) = 0.59
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4.and.tank%x(i)>1.4)then          
           k = k + 1
           this%x(2,k) = 0.01      
           this%x(1,k) = 0.61 + tank%dx * (k - this%ntotal - this%nvirt - 1)          
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>0.98)then
         k = k + 1
         this%x(2,k) = 0.01 
         this%x(1,k) = 0.99 + tank%dx * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = 0.99
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
end subroutine
      
!-----------------------------------------------------
      subroutine take_boundary_for_tank_wedge(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 1.51 - tank%dy * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
         this%x(2,k) = 0.01 - tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)<0.84)then
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 0.03 + tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>0.84.and.tank%x(i)<1.1)then
           k = k + 1
           this%x(2,k) = 0.01 + tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%x(1,k) = 0.85 + tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>1.1.and.tank%x(i)<1.34)then
           k = k + 1
           this%x(2,k) = 0.23 - tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%x(1,k) = 1.11 + tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>1.34.and.tank%x(i)<2.02)then
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 1.35 + tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
         endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>2.02)then
         k = k + 1
         this%x(2,k) = tank%y(i)    
         this%x(1,k) = 2.03 + tank%dx * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = tank%x(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
      end subroutine
      

!-----------------------------------------------------
      subroutine take_boundary_for_tank_wedge2(this,tank)  !虚粒子的间距加密了
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 1.52 - tank%dy/2. * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
         this%x(2,k) = 0.01 - tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)<1.7)then
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 0.02 + tank%dx/2 * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>1.7.and.tank%x(i)<2.22)then
           k = k + 1
           this%x(2,k) = 0.01 + tank%dx/sqrt(2.)/2 * (k - this%ntotal - this%nvirt - 1)
           this%x(1,k) = 0.86 + tank%dx/sqrt(2.)/2 * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>2.22.and.tank%x(i)<2.72)then
           k = k + 1
           this%x(2,k) = 0.01 + tank%dx/sqrt(2.)/2 * 24 - tank%dx/sqrt(2.)/2 * (k - this%ntotal - this%nvirt - 1)
           this%x(1,k) = 0.86 + tank%dx/sqrt(2.)/2 * (k - this%ntotal - this%nvirt - 1 + 26)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>2.82.and.tank%x(i)<4.82)then
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 0.86 + tank%dx/sqrt(2.)/2 * 50 + tank%dx/2 * (k - this%ntotal - this%nvirt - 1 + 1)
           this%zone(k) = tank%zone(i)
         endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>2.42)then
         k = k + 1
         this%x(2,k) = tank%y(i)    
         this%x(1,k) = 0.86 + tank%dx/sqrt(2.)/2 * 50 + tank%dx/2 * 101
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.02 + tank%dy/2 * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = 0.86 + tank%dx/sqrt(2.)/2 * 50 + tank%dx/2 * 101
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
      end subroutine
!-----------------------------------------------------
      subroutine take_boundary_for_tank_wedge3(this,tank)   !这个设计是虚粒子与实粒子间距相等的设计
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 1.51 - tank%dy * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
         this%x(2,k) = 0.01 - tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)<0.86)then
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 0.03 + tank%dx * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>0.88.and.tank%x(i)<1.22)then
           k = k + 1
           this%x(2,k) = 0.01 + tank%dx/sqrt(2.) * (k - this%ntotal - this%nvirt - 1)
           this%x(1,k) = 0.87 + tank%dx/sqrt(2.) * (k - this%ntotal - this%nvirt - 1)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>1.22.and.tank%x(i)<1.54)then
           k = k + 1
           this%x(2,k) = 0.01 + tank%dx/sqrt(2.) * 15 - tank%dx/sqrt(2.) * (k - this%ntotal - this%nvirt - 1)
           this%x(1,k) = 0.87 + tank%dx/sqrt(2.) * (k - this%ntotal - this%nvirt - 1 + 17)
           this%zone(k) = tank%zone(i)
          endif 
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then   
         if(tank%x(i)>1.54)then
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 0.87 + tank%dx/sqrt(2.) * 32 + tank%dx * (k - this%ntotal - this%nvirt - 1 + 1)
           this%zone(k) = tank%zone(i)
         endif 
      endif
   enddo
  
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>2.42)then
         k = k + 1
         this%x(2,k) = tank%y(i)    
         this%x(1,k) = 0.87 + tank%dx/sqrt(2.) * 32 + tank%dx * 45
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = 0.87 + tank%dx/sqrt(2.) * 32 + tank%dx * 45
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
      end subroutine


!-----------------------------------------------------
      subroutine take_boundary_for_tank_wedge_freesurface(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 1.025 - tank%dy/2 * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 1.025 - 0.525 - 0.025 - tank%dy/2 * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then         
         if(tank%x(i)<25.55)then
         k = k + 1
         this%x(1,k) = 0.025 + tank%dx/2 * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then         
         if(tank%x(i)<25.55)then
         k = k + 1
         this%x(1,k) = 0.025 + 12.775 + 0.025 + tank%dx/2 * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = tank%y(i)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then         
         if(tank%x(i)>25.55.and.tank%x(i)<28.55)then
         k = k + 1
         this%x(1,k) = 0.025 + 12.775 + 0.025 + 9./sqrt(9.16)*tank%dx/2 * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = 0.025 + 1.2/sqrt(9.16)*tank%dy/2 * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then         
         if(tank%x(i)>25.55.and.tank%x(i)<28.55)then
         k = k + 1
         this%x(1,k) = 0.025 + 12.775 + 0.025 + 9./sqrt(9.16)*tank%dx/2 * 121
         this%x(2,k) = 0.025 + 1.2/sqrt(9.16)*tank%dy/2 * 121 
         this%zone(k) = tank%zone(i)
         endif
         exit
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then         
         if(tank%x(i)>28.55.and.tank%x(i)<31.55)then
         k = k + 1
         this%x(1,k) = 0.025 + 12.775 + 0.025 + 9./sqrt(9.16)*tank%dx/2 * 121 - 9./sqrt(9.16)*tank%dx/2 * (k - this%ntotal - this%nvirt - 1)
         this%x(2,k) = 0.025 + 1.2/sqrt(9.16)*tank%dy/2 * 121 - 1.2/sqrt(9.16)*tank%dy/2 * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
         exit
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.025 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = tank%x(i)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal


return
      end subroutine
      

      
!-----------------------------------------------------
      subroutine take_boundary_for_tank_static(this,tank)
!-----------------------------------------------------
implicit none

class(particles) this
type(block) tank
integer zone

integer i,j,k

! Take virtual particles in tank

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 3)then         
         k = k + 1
         this%x(1,k) = tank%x(i)
         this%x(2,k) = 0.99 - tank%dy * (k - this%ntotal - 1)
         this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)<0.02)then
         k = k + 1
         this%x(1,k) = tank%x(i)    
         this%x(2,k) = 0.01 - tank%dy * (k - this%ntotal - this%nvirt - 1)  
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal


k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 4)then          
           k = k + 1
           this%x(2,k) = tank%y(i)
           this%x(1,k) = 0.03 + tank%dx * (k - this%ntotal - this%nvirt - 1)           
           this%zone(k) = tank%zone(i)
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 6)then         
         if(tank%x(i)>0.98)then
         k = k + 1
         this%x(2,k) = tank%y(i)    
         this%x(1,k) = 0.99 + tank%dx * (k - this%ntotal - this%nvirt - 1)
         this%zone(k) = tank%zone(i)
         endif
      endif
   enddo
   
this%nvirt = k-this%ntotal

k = this%ntotal+this%nvirt
   do i = 1, tank%m*tank%n
      if(tank%zone(i) == 5)then         
         k = k + 1
         this%x(2,k) = 0.03 + tank%dy * (k - this%ntotal - this%nvirt - 1)
         this%x(1,k) = tank%x(i)
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

!----------------------------------------
      subroutine setup_ndim1(parts)
!----------------------------------------
implicit none
class(particles) parts
integer ntotal, bntotal

ntotal = parts%ntotal + parts%nvirt

if(associated(parts%rho))parts%rho%ndim1 = ntotal
!if(associated(parts%rhos))parts%rhos%ndim1 = ntotal
if(associated(parts%divvx))parts%divvx%ndim1 = ntotal
if(associated(parts%drho))parts%drho%ndim1 = ntotal
if(associated(parts%mass))parts%mass%ndim1 = ntotal
if(associated(parts%p))parts%p%ndim1 = ntotal
!if(associated(parts%ps))parts%p%ndim1 = ntotal
!if(associated(parts%wxy))parts%wxy%ndim1 = ntotal
if(associated(parts%vof))parts%vof%ndim1 = ntotal
if(associated(parts%vof2))parts%vof2%ndim1 = ntotal
if(associated(parts%dvof))parts%dvof%ndim1 = ntotal
if(associated(parts%vof_min))parts%vof_min%ndim1 = ntotal
if(associated(parts%str%x))parts%str%x%ndim1 = ntotal
if(associated(parts%str%xy))parts%str%xy%ndim1 = ntotal
if(associated(parts%str%y))parts%str%y%ndim1 = ntotal
if(associated(parts%dstr))then
   if(associated(parts%dstr%x))parts%dstr%x%ndim1 = ntotal
   if(associated(parts%dstr%xy))parts%dstr%xy%ndim1 = ntotal
   if(associated(parts%dstr%y))parts%dstr%y%ndim1 = ntotal
endif
if(associated(parts%dstr2))then
   if(associated(parts%dstr2%x))parts%dstr2%x%ndim1 = ntotal
   if(associated(parts%dstr2%xy))parts%dstr2%xy%ndim1 = ntotal
   if(associated(parts%dstr2%y))parts%dstr2%y%ndim1 = ntotal
endif
!if(associated(parts%tab%x))parts%tab%x%ndim1 = ntotal
!if(associated(parts%tab%xy))parts%tab%xy%ndim1 = ntotal
!if(associated(parts%tab%y))parts%tab%y%ndim1 = ntotal
if(associated(parts%vcc))parts%vcc%ndim1 = ntotal
if(associated(parts%v_min%x))parts%v_min%x%ndim1 = ntotal
if(associated(parts%v_min%y))parts%v_min%y%ndim1 = ntotal
if(associated(parts%epsilon_p))parts%epsilon_p%ndim1 = ntotal
if(associated(parts%vx%x))parts%vx%x%ndim1 = ntotal
if(associated(parts%vx%y))parts%vx%y%ndim1 = ntotal
if(associated(parts%mone%x))parts%mone%x%ndim1 = ntotal
if(associated(parts%mone%y))parts%mone%y%ndim1 = ntotal
!if(associated(parts%mons%x))parts%mons%x%ndim1 = ntotal
!if(associated(parts%mons%y))parts%mons%y%ndim1 = ntotal
if(associated(parts%dvx%x))parts%dvx%x%ndim1 = ntotal
if(associated(parts%dvx%y))parts%dvx%y%ndim1 = ntotal
if(associated(parts%vol))parts%vol%ndim1 = ntotal
if(associated(parts%rho_min))parts%rho_min%ndim1 = ntotal
!if(associated(parts%u_min))parts%u_min%ndim1 = ntotal
if(associated(parts%p_min))parts%p_min%ndim1 = ntotal
if(associated(parts%str_min))then
   if(associated(parts%str_min%x))parts%str_min%x%ndim1 = ntotal
   if(associated(parts%str_min%y))parts%str_min%y%ndim1 = ntotal
   if(associated(parts%str_min%xy))parts%str_min%xy%ndim1 = ntotal
endif
if(associated(parts%gammaa))parts%gammaa%ndim1 = ntotal
if(associated(parts%dp))parts%dp%ndim1 = ntotal
if(associated(parts%av%x))parts%av%x%ndim1 = ntotal
if(associated(parts%av%y))parts%av%y%ndim1 = ntotal
if(associated(parts%dp2))parts%dp2%ndim1 = ntotal
!if(associated(parts%drhodx))parts%drhodx%ndim1 = ntotal
if(associated(parts%c))parts%c%ndim1 = ntotal

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
  case(5)
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
       subroutine direct_find_Couette(parts)
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
        if(parts%x(1,i)>0.00078.and.parts%x(1,j)<0.00022.and.parts%itype(i)*parts%itype(j)>0)then
          dxiac(1) = parts%x(1,i) - (parts%x(1,j) + 0.001)
          driac    = dxiac(1)*dxiac(1)
          do d=2,parts%dim
            dxiac(d) = parts%x(d,i) - parts%x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (parts%hsml(i)+parts%hsml(j))/2.
!          write(*,*)'sasas'
          if (sqrt(driac).lt.scale_k*mhsml) then
           if (niac.lt.parts%max_interaction) then      
!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 
!          write(*,*)'sasas'
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
       elseif(parts%x(1,i)<0.00022.and.parts%x(1,j)>0.00078.and.parts%itype(i)*parts%itype(j)>0)then
          dxiac(1) = parts%x(1,i) - (parts%x(1,j) - 0.001)
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
        endif

        enddo
      enddo  

      parts%niac = niac
! write(*,*) 'niac=',niac
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
      use ifport
      implicit none

      class(particles) parts
      integer i, j, d, scale_k, ntotal, niac    
      integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)
      real(dp) hsml,dr,r,dx(3),tdwdx(3)

      INTEGER nthreads,n_per_threads, niac_per_threads, it
      integer,allocatable,dimension(:) :: n_start, n_end
      integer last
      real(dp) t1,t2,t3

      ntotal  = parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
      hsml = parts%hsml(1)             !!! hsml is same for all particles!

      do i=1,ntotal
         parts%countiac(i) = 0
      enddo

      call parts%init_grid
      
!     Position particles on grid and create linked list:
!      write(*,*) 'now is ', parts%itimestep,'in link_list'
!      if(parts%itimestep == 8694)then
!        write(*,*)'dfs'  
!      endif
      
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
      allocate(n_start(nthreads),n_end(nthreads))

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
      subroutine link_list2(parts,part2)
!--------------------------------------------
      use ifport
      implicit none
      
      class(particles) parts,part2
      integer i, j, d, scale_k, ntotal, niac    
      integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)
      real(dp) hsml,dr,r,dx(3),tdwdx(3)

      INTEGER nthreads,n_per_threads, niac_per_threads, it
      integer,allocatable,dimension(:) :: n_start, n_end
      integer last
      real(dp) t1,t2,t3
      
      ntotal  = parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
      hsml = parts%hsml(1)

      do i = 1,ntotal
         parts%countiac(i) = 0
      enddo

      call parts%init_grid
      parts%grid = 0; parts%xgcell = 0; parts%celldata = 0
      
!     Position particles on grid and create linked list:
      
      do i=1,ntotal
        call parts%grid_geom(i,parts%x(:,i),gcell)
        do d=1,parts%dim
          parts%xgcell(d,i) = gcell(d)
        enddo
        parts%celldata(i) = parts%grid(gcell(1),gcell(2),gcell(3))
        parts%grid(gcell(1),gcell(2),gcell(3)) = i
      enddo
     
      !call part2%init_grid
      part2%ngridx = parts%ngridx; part2%ghsmlx = parts%ghsmlx
      part2%grid = 0; part2%xgcell = 0; part2%celldata = 0

      
!     Position particles on grid and create linked list:
      
      do i=1,part2%ntotal+part2%nvirt
        call part2%grid_geom(i,part2%x(:,i),gcell)
        do d=1,part2%dim
          part2%xgcell(d,i) = gcell(d)
        enddo
        part2%celldata(i) = part2%grid(gcell(1),gcell(2),gcell(3))
        part2%grid(gcell(1),gcell(2),gcell(3)) = i
      enddo

      nthreads = parts%nthreads
      allocate(n_start(nthreads),n_end(nthreads))

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

!!$omp parallel 
!!$omp do private(niac,i,d,minxcell,maxxcell,dnxgcell,dpxgcell,xcell,ycell,zcell,j,dx,dr,r,tdwdx)
!do it = 1, nthreads
!       niac = parts%niac_start(it)
!       do i = n_start(it),n_end(it)

      niac = 0
      do i=1,ntotal   !ntotal-1

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
              j = part2%grid(xcell,ycell,zcell)
! 1            if (j.gt.i) then
 1            if (j/=0) then
                !dx(1) = parts%x(1,i) - part2%x(1,j)
                !dr    = dx(1)*dx(1)
                dr = 0.d0
                do d=1,parts%dim
                  dx(d) = parts%x(d,i) - part2%x(d,j)
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
                j = part2%celldata(j)
                goto 1
              endif
            enddo !xcell
          enddo !ycell
        enddo !zcell
      enddo !i

      parts%niac = niac

      parts%niac_start(it) = parts%niac_start(it)+1
      parts%niac_end(it) = niac
!      enddo !it      
!!$omp end do
!!$omp end parallel
return

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

! To check

      if(last/=parts%niac) stop 'adfasdf'
      do i = 1, last
         if(parts%pair_i(i)==0) stop 'adf'
      enddo   
end subroutine

!--------------------------------------------
      subroutine link_list2_omp(parts,part2)
!--------------------------------------------
      use ifport
      implicit none
      
      class(particles) parts,part2
      integer i, j, d, scale_k, ntotal, niac    
      integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)
      real(dp) hsml,dr,r,dx(3),tdwdx(3)

      INTEGER nthreads,n_per_threads, niac_per_threads, it
      integer,allocatable,dimension(:) :: n_start, n_end
      integer last
      real(dp) t1,t2,t3
      
      ntotal  = parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
      hsml = parts%hsml(1)

      do i = 1,ntotal
         parts%countiac(i) = 0
      enddo

      call parts%init_grid
      parts%grid = 0; parts%xgcell = 0; parts%celldata = 0
      
!     Position particles on grid and create linked list:
      
      do i=1,ntotal
        call parts%grid_geom(i,parts%x(:,i),gcell)
        do d=1,parts%dim
          parts%xgcell(d,i) = gcell(d)
        enddo
        parts%celldata(i) = parts%grid(gcell(1),gcell(2),gcell(3))
        parts%grid(gcell(1),gcell(2),gcell(3)) = i
      enddo
     
      !call part2%init_grid
      part2%ngridx = parts%ngridx; part2%ghsmlx = parts%ghsmlx
      part2%grid = 0; part2%xgcell = 0; part2%celldata = 0

      
!     Position particles on grid and create linked list:
      
      do i=1,part2%ntotal+part2%nvirt
        call part2%grid_geom(i,part2%x(:,i),gcell)
        do d=1,part2%dim
          part2%xgcell(d,i) = gcell(d)
        enddo
        part2%celldata(i) = part2%grid(gcell(1),gcell(2),gcell(3))
        part2%grid(gcell(1),gcell(2),gcell(3)) = i
      enddo

      nthreads = parts%nthreads
      allocate(n_start(nthreads),n_end(nthreads))

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
!      do i=1,ntotal   !ntotal-1

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
              j = part2%grid(xcell,ycell,zcell)
! 1            if (j.gt.i) then
 1            if (j/=0) then
                !dx(1) = parts%x(1,i) - part2%x(1,j)
                !dr    = dx(1)*dx(1)
                dr = 0.d0
                do d=1,parts%dim
                  dx(d) = parts%x(d,i) - part2%x(d,j)
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
                j = part2%celldata(j)
                goto 1
              endif
            enddo !xcell
          enddo !ycell
        enddo !zcell
      enddo !i

!      parts%niac = niac

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

! To check

      if(last/=parts%niac) stop 'adfasdf'
      do i = 1, last
         if(parts%pair_i(i)==0) stop 'adf'
      enddo   
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
           
    elseif(skf.eq.5)then    ! improved cubic spline kernel for a better pressure stability (M.B.Liu 2012)

        if (dim.eq.1) then
          factor = 1.e0/(7.e0*hsml)
        elseif (dim.eq.2) then
          factor = 1.e0 /(3.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 15.e0/(62.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif    
	if(q.eq.0) then
        w = factor * 6
        do d = 1, dim
            dwdx(d) = -6 * factor
        enddo
	else if(q.gt.0.and.q.le.1) then
          w = factor * ( q**3 - 6*q + 6 )
          do d= 1, dim
            dwdx(d) = factor * (3*q**2 - 6) /hsml*(dx(d)/r)
          enddo 
	else if(q.gt.1.and.q.le.2) then
          w = factor * ( 2 - q )**3
          do d= 1, dim
            dwdx(d) = factor * 3 * ( -q**2 + 4*q - 4) /hsml*(dx(d)/r)
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

         !do d = 1, parts%dim
         !   parts%dvx(d, ii) = parts%dvx(d, ii) + dd*dx(d)*f
         !enddo

         parts%dvx%x%r(ii) = parts%dvx%x%r(ii) + dd * dx(1)*f
         parts%dvx%y%r(ii) = parts%dvx%y%r(ii) + dd * dx(2)*f

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
real(dp),allocatable,dimension(:,:,:) :: local
integer i, j, k, d, ii,ntotal,nthreads,it,dim
           
rr0 = parts%dspp; dd = parts%numeric%dd
p1 = parts%numeric%p1; p2 = parts%numeric%p2
ntotal = parts%ntotal+parts%nvirt; dim = parts%dim
nthreads = parts%nthreads

allocate(local(dim,ntotal,nthreads))
call parts%get_niac_start_end

!$omp parallel
!$omp do private(i,j,rr,d,dx,f,ii,k) 
do it = 1,nthreads
   do i = 1,ntotal
      local(1:dim,i,it) = 0.d0
   enddo
do k=parts%niac_start(it),parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   if(parts%itype(i)*parts%itype(j)>0)cycle  
!   if(parts%itype(i).gt.0.and.parts%itype(j).lt.0)then  
      rr = 0.      
      do d=1,dim
         dx(d) =  parts%x(d,i) -  parts%x(d,j)
         rr = rr + dx(d)*dx(d)
      enddo  
      rr = sqrt(rr)
      !if(rr.lt.rr0)then
      if(rr.gt.rr0)cycle
         f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
         
         ii = i
         if(parts%itype(i)<0)ii=j 

         do d = 1,dim
            !parts%dvx(d, ii) = parts%dvx(d, ii) + dd*dx(d)*f
            local(d,ii,it) = local(d,ii,it) + dd * dx(d)*f
         enddo
      !endif
   !endif        
enddo  !k
enddo  !it
!$omp end do

!$omp do private(it)
do i = 1,ntotal
   do it = 1,nthreads
      parts%dvx%x%r(i) = parts%dvx%x%r(i) + local(1,i,it)
      parts%dvx%y%r(i) = parts%dvx%y%r(i) + local(2,i,it)
   enddo
enddo
!$omp end do
!$omp end parallel
return
end subroutine

      
! Calculate partial derivatives of a field
!--------------------------------------------
          function df(parts,f,x) result(res)
!--------------------------------------------
implicit none

type(array) :: f
character(len=1) x
class(particles) parts
type(array) :: res
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k, ntotal
ntotal = parts%ntotal +parts%nvirt

allocate(res%r(ntotal))
res%ndim1 = ntotal; res = 0.d0

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f%r(i)+f%r(j))*dwdx(k)
   res%r(i) = res%r(i) + parts%mass%r(j)/parts%rho%r(j)*fwx
   res%r(j) = res%r(j) - parts%mass%r(i)/parts%rho%r(i)*fwx
enddo

end function


! Calculate partial derivatives of a field
!-----------------------------------------------------
          function df_omp(parts,f,x)  result(fun)
!-----------------------------------------------------
implicit none

type(array) f
character(len=1) x
class(particles) parts
type(array) :: fun
real(dp), allocatable, dimension(:,:) :: local
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k, ntotal, it, nthreads

!write(*,*) 'In df_omp...'

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads

fun%ndim1 = ntotal
allocate(fun%r(ntotal))
if(nthreads>=1)then
   allocate(local(ntotal,nthreads))
   call parts%get_niac_start_end
endif   

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

!$omp parallel
!$omp do private(i,j,k,fwx) 
do it = 1, nthreads
   do i = 1, ntotal
      local(i,it) = 0.d0
   enddo 
   do k = parts%niac_start(it), parts%niac_end(it)
!do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f%r(i)+f%r(j))*dwdx(k)
   local(i,it) = local(i,it) + parts%mass%r(j)/parts%rho%r(j)*fwx
   local(j,it) = local(j,it) - parts%mass%r(i)/parts%rho%r(i)*fwx
   enddo !k
enddo !it
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1, ntotal
   fun%r(i) = 0.d0
  do it = 1, nthreads
      fun%r(i) = fun%r(i) + local(i,it)
  enddo
enddo   
!$omp end do
!$omp end parallel

end function

! Calculate partial derivatives of a field
!----------------------------------------------
      function df_omp2(parts,f,x) result(val)
!----------------------------------------------
implicit none

type(array) :: f
character(len=1) x
class(particles) parts
type(array) :: val
!real(dp), allocatable, dimension(:,:) :: df_local
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k, ntotal, it, nthreads, ii, jj

!write(*,*) 'In df_omp...'

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads
!write(*,*) 'sadf', nthreads

allocate(val%r(ntotal*nthreads))
!if(nthreads>=1)then
!   allocate(df_local(ntotal,nthreads))
   call parts%get_niac_start_end
!endif   

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

!$omp parallel
!$omp do private(i,j,k,fwx,ii,jj)
do it = 1, parts%nthreads
   do i = 1, ntotal
      ii = i + (it-1)*ntotal
      val%r(ii) = 0.d0
   enddo
   do k = parts%niac_start(it), parts%niac_end(it)


!do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f%r(i)+f%r(j))*dwdx(k)
   ii = i + (it-1)*ntotal; jj = j + (it-1)*ntotal
   val%r(ii) = val%r(ii) + parts%mass%r(j)/parts%rho%r(j)*fwx
   val%r(jj) = val%r(jj) - parts%mass%r(i)/parts%rho%r(i)*fwx   
!   df_local(i,it) = df_local(i,it) + parts%mass(j)/parts%rho(j)*fwx
!   df_local(j,it) = df_local(j,it) - parts%mass(i)/parts%rho(i)*fwx
   enddo !k
enddo !it
!$omp end do
!$omp barrier

!$omp do private(it,ii)
do i = 1, ntotal
   !val(i) = 0.d0
   do it = 2, nthreads
      ii = i + (it-1)*ntotal
      val%r(i) = val%r(i)+val%r(ii)
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

type(array) :: f
character(len=1) x
class(particles) parts
type(array) :: df2
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k, ntotal
ntotal = parts%ntotal +parts%nvirt

allocate(df2%r(ntotal))
df2%ndim1 = ntotal; df2 = 0.d0

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f%r(j)-f%r(i))*dwdx(k)
   df2%r(i) = df2%r(i) + parts%mass%r(j)*fwx
   df2%r(j) = df2%r(j) + parts%mass%r(i)*fwx
enddo

do i = 1, parts%ntotal + parts%nvirt 
   df2%r(i) = df2%r(i)/parts%rho%r(i)
enddo   

end function


! Calculate partial derivatives of a field
!-------------------------------------------
          function df3(parts,f,x)
!------------------------------------------- 
implicit none

type(array) :: f
character(len=1) x
class(particles) parts
type(array) :: df3
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i,j,k,ntotal
ntotal = parts%ntotal + parts%nvirt

allocate(df3%r(ntotal))
df3%ndim1 = ntotal; df3 = 0.d0

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = ((f%r(i)/parts%rho%r(i)**2)+(f%r(j)/parts%rho%r(j)**2))*dwdx(k)
   df3%r(i) = df3%r(i) + parts%mass%r(j)*fwx
   df3%r(j) = df3%r(j) - parts%mass%r(i)*fwx
enddo

do i = 1, parts%ntotal + parts%nvirt
   df3%r(i) = df3%r(i)*parts%rho%r(i)
enddo

end function

!-------------------------------------------
          function df3_omp(parts,f,x)
!------------------------------------------- 
implicit none

type(array) f
character(len=1) x
class(particles) parts
type(array) :: df3_omp
real(dp), allocatable, dimension(:,:) :: local
!real(dp) df3_local(parts%maxn,8)
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i,j,k,it,ntotal,nthreads

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads

df3_omp%ndim1 = ntotal
allocate(df3_omp%r(ntotal))

if(nthreads>=1)then
   allocate(local(ntotal,nthreads))
   call parts%get_niac_start_end
endif 

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

!$omp parallel
!$omp do private(i,j,k,fwx)
do it = 1, parts%nthreads
   do i = 1, ntotal
      local(i,it) = 0.d0
   enddo
   
   do k = parts%niac_start(it), parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = ((f%r(i)/parts%rho%r(i)**2)+(f%r(j)/parts%rho%r(j)**2))*dwdx(k)
   local(i,it)= local(i,it) + parts%mass%r(j)*fwx
   local(j,it)= local(j,it) - parts%mass%r(i)*fwx
!   df3_omp(i) = df3_omp(i) + parts%mass(j)*fwx
!   df3_omp(j) = df3_omp(j) - parts%mass(i)*fwx
   enddo
enddo
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1, ntotal
   df3_omp%r(i) = 0.d0
   do it = 1, nthreads
      df3_omp%r(i) = df3_omp%r(i)+local(i,it)
   enddo
enddo   
!$omp end do
!$omp end parallel

do i = 1, ntotal
   df3_omp%r(i) = df3_omp%r(i)*parts%rho%r(i)
enddo


end function

! Calculate partial derivatives of a field
!-------------------------------------------
          function df4(parts,f,x)
!-------------------------------------------
implicit none

type(array) :: f
character(len=1) x
class(particles) parts
type(array),allocatable :: df4
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k, ntotal
ntotal = parts%ntotal +parts%nvirt

allocate(df4); allocate(df4%r(ntotal))
df4%ndim1 = ntotal; ; df4 = 0.d0

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f%r(j)-f%r(i))*dwdx(k)
   df4%r(i) = df4%r(i) + parts%mass%r(j)/parts%rho%r(j)*fwx
   df4%r(j) = df4%r(j) + parts%mass%r(i)/parts%rho%r(i)*fwx
enddo

end function

!-------------------------------------------------
          function df4_omp(parts,f,x) result(fun)
!-------------------------------------------------
implicit none

type(array) f
character(len=1) x
class(particles) parts
type(array) :: fun
real(dp), allocatable, dimension(:,:) :: local
real(dp), pointer, dimension(:) :: dwdx
real(dp) fwx
integer i, j, k, ntotal, it, nthreads

!write(*,*) 'In df_omp...'

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads

fun%ndim1 = ntotal
allocate(fun%r(ntotal))
if(nthreads>=1)then
   allocate(local(ntotal,nthreads))
   call parts%get_niac_start_end
endif   

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

!$omp parallel
!$omp do private(i,j,k,fwx) 
do it = 1, parts%nthreads
   do i = 1, ntotal
      local(i,it) = 0.d0
   enddo  
   do k = parts%niac_start(it), parts%niac_end(it)
!do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   fwx = (f%r(j)-f%r(i))*dwdx(k)
   local(i,it) = local(i,it) + parts%mass%r(j)/parts%rho%r(j)*fwx
   local(j,it) = local(j,it) + parts%mass%r(i)/parts%rho%r(i)*fwx
   enddo !k
enddo !it
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1, ntotal
   fun%r(i) = 0.d0
   do it = 1, nthreads
      fun%r(i) = fun%r(i) + local(i,it)
   enddo
enddo   
!$omp end do
!$omp end parallel

end function

!-------------------------------------------
          function div(parts,f)
!-------------------------------------------
implicit none

type(array) :: f
class(particles) parts
type(array),allocatable :: div
real(dp) df(3)
real(dp) :: hdiv
real(dp),pointer,dimension(:,:) :: dwdx
integer, pointer, dimension(:) :: pair_i, pair_j
integer i,j,k,ntotal,nthreads,niac,dim,d
type(p2r) f_i(3),f_j(3)

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads
niac = parts%niac; dim = parts%dim
dwdx =>parts%dwdx

allocate(div);allocate(div%r(ntotal))
div%ndim1 = ntotal;div = 0.d0

do k = 1,niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   f_i = f%cmpt(i); f_j = f%cmpt(j)
   do d=1,dim
       df(d) = f_j(d)%p - f_i(d)%p
   enddo
   hdiv = df(1)*dwdx(1,k)
   do d = 2,dim
     hdiv = hdiv + df(d)*dwdx(d,k)
   enddo
   div%r(i) = div%r(i) + parts%mass%r(j)*hdiv/parts%rho%r(j)
   div%r(j) = div%r(j) + parts%mass%r(i)*hdiv/parts%rho%r(i)
enddo

end function

!--------------------------------------------------
          function div_omp(parts,f) result(fun)
!--------------------------------------------------
implicit none

type(array) :: f
class(particles) parts
type(array),allocatable :: fun
real(dp) df(3)
real(dp) :: hdiv
real(dp),pointer,dimension(:,:) :: dwdx
integer, pointer, dimension(:) :: pair_i, pair_j
real(dp), allocatable, dimension(:,:) :: local
integer i,j,k,ntotal,nthreads,niac,dim,d,it
type(p2r) f_i(3),f_j(3)

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads
niac = parts%niac; dim = parts%dim
dwdx =>parts%dwdx

allocate(fun);allocate(fun%r(ntotal))
fun%ndim1 = ntotal

if(nthreads>=1)then
   allocate(local(ntotal,nthreads))
   call parts%get_niac_start_end
endif   

!$omp parallel
!$omp do private(i,j,d,f_i,f_j,df,k,hdiv)
do it = 1, nthreads
   do i = 1, ntotal
      local(i,it) = 0.d0
   enddo   
   do k = parts%niac_start(it), parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   f_i = f%cmpt(i); f_j = f%cmpt(j)
   do d=1,dim
       df(d) = f_j(d)%p - f_i(d)%p
   enddo
   hdiv = df(1)*dwdx(1,k)
   do d = 2,dim
     hdiv = hdiv + df(d)*dwdx(d,k)
   enddo
   local(i,it) = local(i,it) + parts%mass%r(j)*hdiv/parts%rho%r(j)
   local(j,it) = local(j,it) + parts%mass%r(i)*hdiv/parts%rho%r(i)
   enddo
enddo
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1, ntotal
    fun%r(i) = 0.d0
    do it = 1,nthreads
      fun%r(i) =  fun%r(i) + local(i,it) 
    enddo  
enddo   
!$omp end do
!$omp end parallel

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
      type(p2r) vx_i(3), vx_j(3)
      
      pair_i => parts%pair_i
      pair_j => parts%pair_j
      mass   => parts%mass%r
      rho    => parts%rho%r
      vcc    => parts%vcc%r
      !vx     => parts%vx
      dwdx   => parts%dwdx

      ntotal = parts%ntotal + parts%nvirt
      niac   = parts%niac; dim = parts%dim 

      do i=1,ntotal
        vcc(i) = 0.e0         
      enddo
    
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
          !dvx(d) = vx(d,j) - vx(d,i) 
          dvx(d) = vx_j(d)%p - vx_i(d)%p 
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

!-------------------------------------------
          function div2(parts,f) result(res)
!-------------------------------------------
implicit none

type(array) :: f
class(particles) parts
type(array),allocatable :: res
real(dp) df(3)
real(dp) :: temp
real(dp),pointer,dimension(:,:) :: dwdx
integer, pointer, dimension(:) :: pair_i, pair_j
real(dp), allocatable, dimension(:,:) :: local
integer i,j,k,ntotal,nthreads,niac,dim,d, it
type(p2r) f_i(3),f_j(3)

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads
niac = parts%niac; dim = parts%dim
dwdx =>parts%dwdx

allocate(res);allocate(res%r(ntotal))
res%ndim1 = ntotal

allocate(local(ntotal,nthreads))
call parts%get_niac_start_end

!$omp parallel
!$omp do private(k,i,j,f_i,f_j,d,df,temp) 
do it = 1, nthreads
   do i = 1, ntotal
       local(i,it) = 0.d0
   enddo    
do k = parts%niac_start(it),parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   f_i = f%cmpt(i); f_j = f%cmpt(j)
   temp = 0.d0
   do d=1,dim
      df(d) = f_j(d)%p - f_i(d)%p
      temp = temp + df(d)*dwdx(d,k)
   enddo
   local(i,it) = local(i,it) + parts%mass%r(j)*temp
   local(j,it) = local(j,it) + parts%mass%r(i)*temp
enddo !k
enddo !it
!$omp end do

!$omp do private(it)
do i = 1, ntotal
    res%r(i) = 0
    do it = 1, nthreads
        res%r(i) = res%r(i) + local(i,it)
    enddo
    res%r(i) = res%r(i)/parts%rho%r(i)
enddo
!$omp end do
!$omp end parallel

end function
          
!船建师兄认为div4,适合两相流
!-------------------------------------------
          function div4(parts,f) result(res)
!-------------------------------------------
implicit none

type(array) :: f
class(particles) parts
type(array),allocatable :: res
real(dp) df(3)
real(dp) :: temp
real(dp),pointer,dimension(:,:) :: dwdx
integer, pointer, dimension(:) :: pair_i, pair_j
real(dp), allocatable, dimension(:,:) :: local
integer i,j,k,ntotal,nthreads,niac,dim,d, it
type(p2r) f_i(3),f_j(3)

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads
niac = parts%niac; dim = parts%dim
dwdx =>parts%dwdx

allocate(res);allocate(res%r(ntotal))
res%ndim1 = ntotal

allocate(local(ntotal,nthreads))
call parts%get_niac_start_end

!$omp parallel
!$omp do private(k,i,j,f_i,f_j,d,df,temp) 
do it = 1, nthreads
   do i = 1, ntotal
       local(i,it) = 0.d0
   enddo    
do k = parts%niac_start(it),parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   f_i = f%cmpt(i); f_j = f%cmpt(j)
   temp = 0.d0
   do d=1,dim
      df(d) = f_j(d)%p - f_i(d)%p
      temp = temp + df(d)*dwdx(d,k)
   enddo
   local(i,it) = local(i,it) + parts%mass%r(j)*temp/parts%rho%r(j)
   local(j,it) = local(j,it) + parts%mass%r(i)*temp/parts%rho%r(i)
enddo !k
enddo !it
!$omp end do

!$omp do private(it)
do i = 1, ntotal
    res%r(i) = 0
    do it = 1, nthreads
        res%r(i) = res%r(i) + local(i,it)
    enddo
!    res%r(i) = res%r(i)/parts%rho%r(i)
enddo
!$omp end do
!$omp end parallel

end function
          
!div2加上了XSPH的版本。          
!-------------------------------------------
          function div5(parts,f) result(res)
!-------------------------------------------
implicit none

type(array) :: f
class(particles) parts
type(array),allocatable :: res
real(dp) df(3)
real(dp) :: temp
real(dp),pointer,dimension(:,:) :: dwdx
integer, pointer, dimension(:) :: pair_i, pair_j
real(dp), allocatable, dimension(:,:) :: local
integer i,j,k,ntotal,nthreads,niac,dim,d, it
type(p2r) f_i(3),f_j(3),avi(3),avj(3)

ntotal = parts%ntotal + parts%nvirt
nthreads = parts%nthreads
niac = parts%niac; dim = parts%dim
dwdx =>parts%dwdx

allocate(res);allocate(res%r(ntotal))
res%ndim1 = ntotal

allocate(local(ntotal,nthreads))
call parts%get_niac_start_end

!$omp parallel
!$omp do private(k,i,j,f_i,f_j,d,df,temp) 
do it = 1, nthreads
   do i = 1, ntotal
       local(i,it) = 0.d0
   enddo    
do k = parts%niac_start(it),parts%niac_end(it)
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   f_i = f%cmpt(i); f_j = f%cmpt(j);avi = parts%av%cmpt(i);avj = parts%av%cmpt(j)
   temp = 0.d0
   do d=1,dim
      df(d) = (f_j(d)%p + avj(d)%p) - (f_i(d)%p + avi(d)%p)
      temp = temp + df(d)*dwdx(d,k)
   enddo
   local(i,it) = local(i,it) + parts%mass%r(j)*temp
   local(j,it) = local(j,it) + parts%mass%r(i)*temp
enddo !k
enddo !it
!$omp end do

!$omp do private(it)
do i = 1, ntotal
    res%r(i) = 0
    do it = 1, nthreads
        res%r(i) = res%r(i) + local(i,it)
    enddo
    res%r(i) = res%r(i)/parts%rho%r(i)
enddo
!$omp end do
!$omp end parallel

end function          
          
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
         parts%rho = water%rho0*(parts%p/water%b+1.d0)**(1/water%gamma)

      elseif(parts%imaterial=='soil')then

         soil => parts%material
         !parts%rho(1:ntotal) = soil%rho0*(parts%p(1:ntotal)/soil%k+1)

      endif
 
      return
      end subroutine

!---------------------------------------------------------------------
      subroutine get_norm(parts)
!---------------------------------------------------------------------
      implicit none

      class(particles) parts

      integer i
      
      do i = parts%ntotal +1 ,parts%ntotal + parts%nvirt-1
      parts%n(1,i) = (parts%x(2,i) - parts%x(2,i+1))/sqrt((parts%x(1,i) - parts%x(1,i+1))**2 + (parts%x(2,i) - parts%x(2,i+1))**2)
      parts%n(2,i) = (parts%x(1,i+1) - parts%x(1,i))/sqrt((parts%x(1,i) - parts%x(1,i+1))**2 + (parts%x(2,i) - parts%x(2,i+1))**2)
      
      !              nx = (parts%x(2,j) - parts%x(2,j+1))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
!              ny = (parts%x(1,j+1) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
      enddo

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
         parts%p = water%b*((parts%rho/(water%rho0*parts%vof))**water%gamma-1.d0) !!! False density  

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
         parts%c = water%c*(parts%rho/(water%rho0*parts%vof))**3.0         

      elseif(parts%imaterial=='soil')then

         soil => parts%material
!         if(parts%soil_pressure==1)          &
!         parts%p(1:ntotal) = soil%k*(parts%rho(1:ntotal)/soil%rho0-1)
!!         parts%p(1:ntotal) = parts%p(1:ntotal)
!!     &                      -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
         parts%c = soil%c

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
parts%p%r(1:ntotal) = water%b*((parts%rho%r(1:ntotal)/(water%rho0  &
                   *parts%vof%r(1:ntotal))) &   !!! False density
                  **water%gamma-1)

parts%c%r(1:ntotal) = water%c*(parts%rho%r(1:ntotal)/(water%rho0*parts%vof%r(1:ntotal)))**3.0         

return
end subroutine
!just for dam break
!----------------------------------------------------------------------
       subroutine pressure_nvirt(parts)
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by directly comparing the particle distance 
!   with the corresponding smoothing length.
!-----------------------------------------------------------------------
      implicit none
      
      class(particles) parts

      integer ntotal, niac,i,j,d,scale_k,k
      real(dp) dxiac(3), driac, r, mhsml, tdwdx(3),A,hv(3),selfdens
      integer,allocatable, dimension(:)  :: pair_i,pair_j
      real(dp),allocatable, dimension(:)  :: w,wi
      real(dp),allocatable, dimension(:,:)  :: dwdx,x
      real(dp),allocatable, dimension(:)  :: strxx,strxyy,strxyx,stryy
      logical :: dbg = .false.

      allocate(pair_i(parts%max_interaction))
      allocate(wi(parts%nvirt))
      allocate(strxx(parts%nvirt))
      allocate(strxyy(parts%nvirt))
      allocate(stryy(parts%nvirt))
      allocate(strxyx(parts%nvirt))
      allocate(pair_j(parts%max_interaction))
      allocate(w(parts%max_interaction))
      allocate(dwdx(parts%dim,parts%max_interaction))
      allocate(x(parts%dim,parts%max_interaction))
      if(dbg) write(*,*) 'In direct_find...'

      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
      ntotal   =   parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
      
      do i=1,ntotal
        parts%countiac(i) = 0
      enddo

      niac = 0


      
      do i= parts%ntotal+1,parts%ntotal+parts%nvirt
       if(parts%zone(i)==1.or.parts%zone(i)==4.or.parts%zone(i)==6) then 
 !       parts%p%r(i) = 0
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        parts%p%r(i) = selfdens*parts%p%r(i)*parts%mass%r(i)/parts%rho%r(i)
          wi(i-parts%ntotal) = 0
       endif
      enddo
      
      do i = parts%ntotal+1,parts%ntotal+parts%nvirt !for all nvirt particles
       if(parts%zone(i)==1.or.parts%zone(i)==4.or.parts%zone(i)==6) then 
          if(parts%zone(i)==1) then 
              x(1,i) = 0.12 - parts%x(1,i)
              x(2,i) = parts%x(2,i)
          elseif(parts%zone(i)==3) then
              write(*,*) 'parts%zone(i)==3 is nothing'
          elseif(parts%zone(i)==4) then
              x(1,i) = parts%x(1,i)
              x(2,i) = 0.12 - parts%x(2,i)
          elseif(parts%zone(i)==5) then
              write(*,*) 'parts%zone(i)==5 is nothing'
          elseif(parts%zone(i)==6) then
              x(1,i) = 3.28*2 - parts%x(1,i)
              x(2,i) = parts%x(2,i)
          endif
         do j = 1,parts%ntotal !for all real particles
          dxiac(1) = x(1,i) - parts%x(1,j)
          driac    = dxiac(1)*dxiac(1)
          do d=1,parts%dim
            dxiac(d) = x(d,i) - parts%x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (parts%hsml(i)+parts%hsml(j))/2.
          if (sqrt(driac).lt.scale_k*mhsml) then
            if (niac.lt.parts%max_interaction) then    
              niac = niac + 1
              pair_i(niac) = i
              pair_j(niac) = j
              r = sqrt(driac)
              parts%countiac(i) = parts%countiac(i) + 1
              parts%countiac(j) = parts%countiac(j) + 1
              call parts%kernel(r,dxiac,mhsml,w(niac),tdwdx)
              do d=1,parts%dim
                dwdx(d,niac) = tdwdx(d)
              enddo                                
                
                
!      do i=1,ntotal-1     
!        do j = i+1, ntotal
!          dxiac(1) = parts%x(1,i) - parts%x(1,j)
!          driac    = dxiac(1)*dxiac(1)
!          do d=2,parts%dim
!            dxiac(d) = parts%x(d,i) - parts%x(d,j)
!            driac    = driac + dxiac(d)*dxiac(d)
!          enddo
!          mhsml = (parts%hsml(i)+parts%hsml(j))/2.
!          if (sqrt(driac).lt.scale_k*mhsml) then
!            if (niac.lt.parts%max_interaction) then    

!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 

!              niac = niac + 1
!              parts%pair_i(niac) = i
!              parts%pair_j(niac) = j
!              r = sqrt(driac)
!              parts%countiac(i) = parts%countiac(i) + 1
!              parts%countiac(j) = parts%countiac(j) + 1

!     Kernel and derivations of kernel
!              call parts%kernel(r,dxiac,mhsml,parts%w(niac),tdwdx)
!              do d=1,parts%dim
!                parts%dwdx(d,niac) = tdwdx(d)
!              enddo                                      
            else
              print *,  ' >>> ERROR <<< : Too many interactions' 
              stop
            endif
          parts%p%r(i) = parts%p%r(i) + parts%p%r(j)*parts%mass%r(j)/parts%rho%r(j)*w(niac)
          strxx(i-parts%ntotal) = strxx(i-parts%ntotal) + (parts%str%x%r(i) + parts%str%x%r(j))*dwdx(1,niac)*parts%mass%r(j)/parts%rho%r(j)
          strxyx(i-parts%ntotal) = strxyx(i-parts%ntotal) + (parts%str%xy%r(i) + parts%str%xy%r(j))*dwdx(1,niac)*parts%mass%r(j)/parts%rho%r(j)
          strxyy(i-parts%ntotal) = strxyy(i-parts%ntotal) + (parts%str%xy%r(i) + parts%str%xy%r(j))*dwdx(2,niac)*parts%mass%r(j)/parts%rho%r(j)
          stryy(i-parts%ntotal) = stryy(i-parts%ntotal) + (parts%str%y%r(i) + parts%str%y%r(j))*dwdx(2,niac)*parts%mass%r(j)/parts%rho%r(j)
                   wi(i-parts%ntotal) = wi(i-parts%ntotal) + w(niac)*parts%mass%r(j)/parts%rho%r(j)
          endif
         enddo
       endif
      enddo  
      
      
!we do not konw weather to add itself
!      do i= parts%ntotal+1,parts%ntotal+parts%nvirt
!       if(parts%zone(i)==1.or.parts%zone(i)==4.or.parts%zone(i)==6) then 
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
!        parts%p%r(i) = parts%p%r(i)*selfdens*parts%mass%r(i)/parts%rho%r(i)
!       endif
!      enddo
      do i = 1,parts%nvirt
          if(wi(i)==0) wi(i)=1
      enddo
      
      do i= parts%ntotal+1,parts%ntotal+parts%nvirt
       if(parts%zone(i)==1) then 
           parts%p%r(i) = (parts%p%r(i) - (strxx(i-parts%ntotal) + strxyy(i-parts%ntotal))*(x(1,i) - 0.06)*2)/wi(i-parts%ntotal)
       endif
       if(parts%zone(i)==4) then 
           parts%p%r(i) = (parts%p%r(i) - 2*(x(2,i)-0.06)*(parts%rho%r(i)*parts%mass%r(i)*parts%numeric%gravity + strxyx(i-parts%ntotal) + stryy(i-parts%ntotal)))/wi(i-parts%ntotal)
       endif
       if(parts%zone(i)==6) then 
           parts%p%r(i) = (parts%p%r(i)+ (strxx(i-parts%ntotal) + strxyy(i-parts%ntotal))*(x(1,i) - 3.28)*2)/wi(i-parts%ntotal)
       endif
      enddo
      
      deallocate(pair_i)
      deallocate(wi)
      deallocate(strxx)
      deallocate(strxyy)
      deallocate(stryy)
      deallocate(strxyx)
      deallocate(pair_j)
      deallocate(w)
      deallocate(dwdx)
      deallocate(x)
      
!      parts%niac = niac
 
       end subroutine
            

!----------------------------------------------------------------------
       subroutine pressure_nvirt_unified(parts)
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by directly comparing the particle distance 
!   with the corresponding smoothing length.
!-----------------------------------------------------------------------
      implicit none
      
      class(particles) parts

      integer ntotal, niac,i,j,d,scale_k,k
      real(dp) dxiac(3), driac, r, mhsml, tdwdx(3),A,hv(3),selfdens
      real(dp)  q0,r1,r2,y1,y2,q00,zeta(2),nx,ny,hsml,cos1,cos2
      real(dp) , parameter :: pi = 3.1415926535898
      integer,allocatable, dimension(:)  :: pair_i,pair_j
      real(dp),allocatable, dimension(:)  :: w,wi,p,gammaas,gammaa,rho,rav1,rav2
      real(dp),allocatable, dimension(:,:)  :: dwdx,x
      real(dp),allocatable, dimension(:)  :: strxx,strxyy,strxyx,stryy
      type(material),pointer :: water
!      real(dp),allocatable, dimension(:)  :: strxx,strxyy,strxyx,stryy
      logical :: dbg = .false.

      allocate(pair_i(parts%max_interaction))
      allocate(wi(parts%nvirt))
      allocate(gammaas(parts%max_interaction))      
      allocate(gammaa(parts%nvirt))            
      allocate(p(parts%nvirt))
      allocate(rho(parts%nvirt))
      allocate(strxx(parts%nvirt))
      allocate(strxyy(parts%nvirt))
      allocate(stryy(parts%nvirt))
      allocate(strxyx(parts%nvirt))
      allocate(pair_j(parts%max_interaction))
      allocate(w(parts%max_interaction))
      allocate(dwdx(parts%dim,parts%max_interaction))
      allocate(x(parts%dim,parts%ntotal+parts%nvirt))
      if(dbg) write(*,*) 'In direct_find...'

      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      hsml = parts%hsml(1)
      ntotal   =   parts%ntotal + parts%nvirt
      scale_k = parts%get_scale_k()
    
      do i=1,ntotal
        parts%countiac(i) = 0
      enddo

      niac = 0
      gammaa = 0.
      rho = 0.
      strxx = 0.
      strxyy = 0.
      stryy = 0.
      strxyx = 0.
      
      do i= parts%ntotal+1,parts%ntotal+parts%nvirt
!       if(parts%zone(i)==3.or.parts%zone(i)==4.or.parts%zone(i)==5.or.parts%zone(i)==6) then 
        p(i-parts%ntotal) = selfdens*parts%p%r(i)*parts%mass%r(i)/parts%rho%r(i)
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        rho(i-parts%ntotal) = 0!selfdens*parts%mass%r(i)
!        parts%p%r(i) = selfdens*parts%p%r(i)*parts%mass%r(i)/parts%rho%r(i)
          wi(i-parts%ntotal) = selfdens*parts%mass%r(i)/parts%rho%r(i)!0
!       endif
      enddo
      
      do i = parts%ntotal+1,parts%ntotal+parts%nvirt !for all nvirt particles
!       if(parts%zone(i)==1.or.parts%zone(i)==4.or.parts%zone(i)==6) then 
          if(parts%zone(i)==3) then 
              x(1,i) = 0.04 - parts%x(1,i)
              x(2,i) = parts%x(2,i)  
!          elseif(parts%zone(i)==5) then
!              x(1,i) = 3.96 - parts%x(1,i)
!              x(2,i) = parts%x(2,i)                 
          elseif(i==1060) then
              x(1,i) = 0.04 - parts%x(1,i)
              x(2,i) = 0.04 - parts%x(2,i)
          elseif(1061<=i.and.i<=1079) then
              x(1,i) = parts%x(1,i)
              x(2,i) = 0.04 - parts%x(2,i)
          elseif(i==1080) then
              x(1,i) = 0.80 - parts%x(1,i)
              x(2,i) = 0.04 - parts%x(2,i)
          elseif(1081<=i.and.i<=1098) then
              x(1,i) = 0.80 - parts%x(1,i)
              x(2,i) = parts%x(2,i)
          elseif(i==1099) then
              x(1,i) = 0.80 - parts%x(1,i)
              x(2,i) = 0.80 - parts%x(2,i)
!              write(*,*)'i=',i
          elseif(1100<= i .and.i<=1107)then
              x(1,i) = parts%x(1,i)
              x(2,i) = 0.80 - parts%x(2,i)
!                            write(*,*)'i=',i
          elseif(i==1108) then
              x(1,i) = 1.20 - parts%x(1,i)
              x(2,i) = 0.80 - parts%x(2,i)              
          elseif(1109<=i.and.i<=1126)then
              x(1,i) = 1.20 - parts%x(1,i)
              x(2,i) = parts%x(2,i)
          elseif(i==1127) then
              x(1,i) = 1.20 - parts%x(1,i)
              x(2,i) = 0.04 - parts%x(2,i)     
          elseif(1128<=i.and.i<=1146) then
              x(1,i) = parts%x(1,i)
              x(2,i) = 0.04 - parts%x(2,i)   
          elseif(i==1147) then
              x(1,i) = 1.96 - parts%x(1,i)
              x(2,i) = 0.04 - parts%x(2,i)    
          elseif(1148<=i) then
              x(1,i) = 1.96 - parts%x(1,i)
              x(2,i) = parts%x(2,i)  
!                            write(*,*) 'i=',i

          endif
         do j = 1,parts%ntotal + parts%nvirt !for all real particles,但改成了所有粒子
          dxiac(1) = x(1,i) - parts%x(1,j)
          driac    = dxiac(1)*dxiac(1)
          do d=1,parts%dim
            dxiac(d) = x(d,i) - parts%x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (parts%hsml(i)+parts%hsml(j))/2.
          if (sqrt(driac).lt.scale_k*mhsml) then
            if (niac.lt.parts%max_interaction) then    
              niac = niac + 1
              pair_i(niac) = i
              pair_j(niac) = j
              r = sqrt(driac)
              parts%countiac(i) = parts%countiac(i) + 1
              parts%countiac(j) = parts%countiac(j) + 1
              call parts%kernel(r,dxiac,mhsml,w(niac),tdwdx)
              do d=1,parts%dim
                dwdx(d,niac) = tdwdx(d)
              enddo                                
                
                
!      do i=1,ntotal-1     
!        do j = i+1, ntotal
!          dxiac(1) = parts%x(1,i) - parts%x(1,j)
!          driac    = dxiac(1)*dxiac(1)
!          do d=2,parts%dim
!            dxiac(d) = parts%x(d,i) - parts%x(d,j)
!            driac    = driac + dxiac(d)*dxiac(d)
!          enddo
!          mhsml = (parts%hsml(i)+parts%hsml(j))/2.
!          if (sqrt(driac).lt.scale_k*mhsml) then
!            if (niac.lt.parts%max_interaction) then    

!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 

!              niac = niac + 1
!              parts%pair_i(niac) = i
!              parts%pair_j(niac) = j
!              r = sqrt(driac)
!              parts%countiac(i) = parts%countiac(i) + 1
!              parts%countiac(j) = parts%countiac(j) + 1

!     Kernel and derivations of kernel
!              call parts%kernel(r,dxiac,mhsml,parts%w(niac),tdwdx)
!              do d=1,parts%dim
!                parts%dwdx(d,niac) = tdwdx(d)
!              enddo                                      
            else
              print *,  ' >>> ERROR <<< : Too many interactions' 
              stop
            endif
            if(parts%itype(j)>0)then
            rav1 = parts%x(1,j) - parts%x(1,i)
            rav2 = parts%x(2,j) - parts%x(2,i)
          p(i-parts%ntotal) = p(i-parts%ntotal) + parts%p%r(j) *parts%mass%r(j)/parts%rho%r(j)*w(niac)
          strxx(i-parts%ntotal) = strxx(i-parts%ntotal) + (parts%str%x%r(i) + parts%str%x%r(j))*dwdx(1,niac)*parts%mass%r(j)/parts%rho%r(j)
          strxyx(i-parts%ntotal) = strxyx(i-parts%ntotal) + (parts%str%xy%r(i) + parts%str%xy%r(j))*dwdx(1,niac)*parts%mass%r(j)/parts%rho%r(j)
          strxyy(i-parts%ntotal) = strxyy(i-parts%ntotal) + (parts%str%xy%r(i) + parts%str%xy%r(j))*dwdx(2,niac)*parts%mass%r(j)/parts%rho%r(j)
          stryy(i-parts%ntotal) = stryy(i-parts%ntotal) + (parts%str%y%r(i) + parts%str%y%r(j))*dwdx(2,niac)*parts%mass%r(j)/parts%rho%r(j)
                   wi(i-parts%ntotal) = wi(i-parts%ntotal) + w(niac)*parts%mass%r(j)/parts%rho%r(j)
            endif
          if(parts%itype(j)<0)then
              nx = (parts%x(2,j) - parts%x(2,j+1))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
              ny = (parts%x(1,j+1) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
              q0 = sqrt(((x(1,i) - parts%x(1,j))*nx)**2 + ((x(2,i) - parts%x(2,j))*ny)**2)/hsml
              r1 = sqrt((x(1,i) - parts%x(1,j))**2 + (x(2,i) - parts%x(2,j))**2)
              r2 = sqrt((x(1,i) - parts%x(1,j+1))**2 + (x(2,i) - parts%x(2,j+1))**2)
              cos1 = -((x(1,i) - parts%x(1,j)) * (parts%x(1,j+1) - parts%x(1,j)) + (x(2,i) - parts%x(2,j)) * (parts%x(2,j+1) - parts%x(2,j)))/   &
                 (sqrt((x(1,i) - parts%x(1,j))**2 + (x(2,i) -parts%x(2,j))**2)*sqrt((parts%x(1,j+1) - parts%x(1,j))**2 + (parts%x(2,j+1) - parts%x(2,j))**2))
              cos2 = -((x(1,i) - parts%x(1,j+1)) * (parts%x(1,j+1) - parts%x(1,j)) + (x(2,i) - parts%x(2,j+1)) * (parts%x(2,j+1) - parts%x(2,j)))/  &
                 (sqrt((x(1,i) - parts%x(1,j+1))**2 + (x(2,i) -parts%x(2,j+1))**2)*sqrt((parts%x(1,j+1) - parts%x(1,j))**2 + (parts%x(2,j+1) - parts%x(2,j))**2))
              y1 = r1*cos1
              y2 = r2*cos2
              q00 = min(q0,2.)
              zeta(1) = min(sqrt(q00**2 + y1**2/hsml**2)/2,1.)
              zeta(2) = min(sqrt(q00**2 + y2**2/hsml**2)/2,1.)
              if(y1==0)then
              gammaas(niac) = -1./(4*pi)*(sign(1.,y2)*parts%psi(q0,zeta(2)))
              elseif(y2==0)then
              gammaas(niac) = -1./(4*pi)*( - sign(1.,y1)*parts%psi(q0,zeta(1)))
              else
              gammaas(niac) = -1./(4*pi)*(sign(1.,y2)*parts%psi(q0,zeta(2)) - sign(1.,y1)*parts%psi(q0,zeta(1)))    
              endif
              gammaa(i-parts%ntotal) = gammaa(i-parts%ntotal) + gammaas(niac)     !因为为负号的话，下面会>1、、、
          endif
             
          endif
         enddo
!       endif
      enddo  
      
      water => parts%material
      do i = parts%ntotal + 1,parts%ntotal + parts%nvirt
          
          if(wi(i-parts%ntotal)<0.1) wi(i-parts%ntotal)=1
          gammaa(i-parts%ntotal) = gammaa(i-parts%ntotal) + 1.
!          rho(i-parts%ntotal) = rho(i-parts%ntotal)/gammaa(i-parts%ntotal)
!          rho(i-parts%ntotal) = rho(i-parts%ntotal)/wi(i-parts%ntotal)
!          if(rho(i-parts%ntotal)>900)p(i-parts%ntotal) = water%b*((rho(i-parts%ntotal)/(water%rho0))**water%gamma-1.d0)
          p(i-parts%ntotal) = p(i-parts%ntotal)/ wi(i-parts%ntotal)
          if(parts%gammaa%r(i)>1.1d0)write(*,*) 'gammaa2 has something wrong'
      enddo
       
         
!we do not konw weather to add itself
!      do i= parts%ntotal+1,parts%ntotal+parts%nvirt
!       if(parts%zone(i)==1.or.parts%zone(i)==4.or.parts%zone(i)==6) then 
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
!        parts%p%r(i) = parts%p%r(i)*selfdens*parts%mass%r(i)/parts%rho%r(i)
!       endif
!      enddo
!      do i = 1,parts%nvirt
!          if(wi(i)<0.1) wi(i)=1
!      enddo
      
      do i= parts%ntotal+1,parts%ntotal+parts%nvirt
!       if(parts%zone(i)==1) then 
!           parts%p%r(i) = (parts%p%r(i) - (strxx(i-parts%ntotal) + strxyy(i-parts%ntotal))*(x(1,i) - 0.06)*2)/wi(i-parts%ntotal)
!           p(i-parts%ntotal) = p(i-parts%ntotal)/wi(i-parts%ntotal)
!       endif
!       write(*,*) 'gammaa(i-parts%ntotal)=',gammaa(i-parts%ntotal),'wi(i-parts%ntotal)=',wi(i-parts%ntotal),i
!       if(1060<=i<=1080) then 
!           p(i-parts%ntotal) = (p(i-parts%ntotal) - parts%hsml(1)/2.*(parts%rho%r(i)*parts%mass%r(i)*parts%numeric%gravity ))!/gammaa(i-parts%ntotal)
!       elseif(1099<=i<=1108)then
!           p(i-parts%ntotal) = (p(i-parts%ntotal) - parts%hsml(1)/2.*(parts%rho%r(i)*parts%mass%r(i)*parts%numeric%gravity ))!/gammaa(i-parts%ntotal)
!       elseif(1127<=i<=1147)then
!           p(i-parts%ntotal) = (p(i-parts%ntotal) - parts%hsml(1)/2.*(parts%rho%r(i)*parts%mass%r(i)*parts%numeric%gravity ))!/gammaa(i-parts%ntotal)
!       else
!           p(i-parts%ntotal) = p(i-parts%ntotal)!/gammaa(i-parts%ntotal)
!       endif
        if(1060<=i<=1080) then 
!             p(i-parts%ntotal) = p(i-parts%ntotal)!/gammaa(i-parts%ntotal)
            p(i-parts%ntotal) = (p(i-parts%ntotal) - parts%hsml(1)/2.*parts%rho%r(i)*(parts%numeric%gravity + strxyx(i-parts%ntotal) + stryy(i-parts%ntotal) ))!/gammaa(i-parts%ntotal)
        elseif(1099<=i<=1108)then
!                        p(i-parts%ntotal) = p(i-parts%ntotal)!/gammaa(i-parts%ntotal)
            p(i-parts%ntotal) = (p(i-parts%ntotal) - parts%hsml(1)/2.*parts%rho%r(i)*(parts%numeric%gravity + strxyx(i-parts%ntotal) + stryy(i-parts%ntotal) ))!/gammaa(i-parts%ntotal)
        elseif(1127<=i<=1147)then
!            p(i-parts%ntotal) = p(i-parts%ntotal)!/gammaa(i-parts%ntotal)
            p(i-parts%ntotal) = (p(i-parts%ntotal) - parts%hsml(1)/2.*parts%rho%r(i)*(parts%numeric%gravity + strxyx(i-parts%ntotal) + stryy(i-parts%ntotal)))!/gammaa(i-parts%ntotal)
        else
            p(i-parts%ntotal) = p(i-parts%ntotal) - parts%rho%r(i)*(strxx(i-parts%ntotal) + strxyy(i-parts%ntotal))*parts%hsml(1)/2.!/gammaa(i-parts%ntotal)
!                        p(i-parts%ntotal) = (p(i-parts%ntotal) - parts%hsml(1)/2.*(parts%rho%r(i)*parts%numeric%gravity ))!/gammaa(i-parts%ntotal)
        endif
      
      enddo
      
      do i= parts%ntotal+1,parts%ntotal+parts%nvirt
             parts%p%r(i) = p(i-parts%ntotal)
      enddo
                
      deallocate(pair_i)
      deallocate(wi)
      deallocate(strxx)
      deallocate(strxyy)
      deallocate(stryy)
      deallocate(strxyx)
      deallocate(pair_j)
      deallocate(w)
      deallocate(dwdx)
      deallocate(x)
      
!      parts%niac = niac
 
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
   if(parts%p%r(i)<0)then
      parts%p%r(i) = 0.d0
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
parts%p%r(1:ntotal) = soil%k*(parts%rho%r(1:ntotal)/soil%rho0-1)
! parts%p(1:ntotal) = parts%p(1:ntotal)  &
!                    -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
parts%c%r(1:ntotal) = soil%c

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
        wi(i)=selfdens*parts%mass%r(i)/parts%rho%r(i)
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        wi(i) = wi(i) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k)
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
      enddo

!     Secondly calculate the rho integration over the space

      do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        parts%rho%r(i) = selfdens*parts%mass%r(i)
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        parts%rho%r(i) = parts%rho%r(i) + parts%mass%r(j)*parts%w(k)
        parts%rho%r(j) = parts%rho%r(j) + parts%mass%r(i)*parts%w(k)
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
!      if (nor_density) then 
        do i=1, ntotal
          parts%rho%r(i)=parts%rho%r(i)/wi(i)
        enddo
!      endif 
 
      end subroutine
      

!Subroutine to calculate the gammaa .
!----------------------------------------------------------------------
      subroutine get_gammaa(parts) 
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%maxn),dvx(3) 
!      type(array),allocatable :: divvx
      real(dp), allocatable, dimension(:,:) :: dgua
      type(p2r) vx_i(3), vx_j(3) 

      r=0.
      hv=0.
      do i=1,parts%ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        parts%gammaa%r(i)=selfdens*parts%mass%r(i)/parts%rho%r(i)
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        parts%gammaa%r(i) = parts%gammaa%r(i) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k)
        parts%gammaa%r(j) = parts%gammaa%r(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
      enddo
      
      do i=1,parts%ntotal
        if(parts%gammaa%r(i)>0.99.and.parts%gammaa%r(i)<1.01) parts%gammaa%r(i)=1.d0
      enddo
     
      do i=parts%ntotal + 1,parts%ntotal + parts%nvirt
          parts%gammaa%r(i) = 0.d0
      enddo
      
      end subroutine


!-------------------------------------------------
          function psi(parts,q,zeta) result(val)
!-------------------------------------------------
implicit none

class(particles) parts
real(dp) zeta,q,val,a,b,c
integer i, j, k

!c = 4*zeta**2/q**4
!a = acosh(2*zeta/q) 
!b = atan(sqrt(4*zeta**2/q**2 -1))

val = q * sqrt((zeta**2 - q**2./4.)) * (-4./3.*zeta**5 + 7.*zeta**4 - (5./12.* q**2 +14.)*zeta**3 + 7./3.*(q**2+5.)*zeta**2   &
      - 1./4.*(5./8.*q**2 +21.)*q**2*zeta + 7./6.*q**4 +35./6.*q**2 - 7.) - (5./8.*q**2 + 21.)*q**5/16.*acosh(2*zeta/q)   &
      + 2*atan(sqrt(4*zeta**2/q**2 -1))
      
end function


!Subroutine to calculate the gammaa .
!----------------------------------------------------------------------
      subroutine get_gammaa2(parts) !j 和j+1   其实始终觉得这里的这种球get_gammaa2的方法，不太好，对于在流体内部的粒子，没有修正，要能是与get_gammaa的结果做一个调和平均就好了
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer  i, j, k, d      
      real(dp)  q0,r1,r2,y1,y2,q00,zeta(2),nx,ny,hsml,cos1,cos2
      real(dp) , parameter :: pi = 3.1415926535898
      
      hsml = parts%hsml(1)
      do i = 1, parts%ntotal
          parts%gammaa%r(i) = 0.d0
      enddo
      
      do k = 1,parts%niac
!          do d = 1,parts%dim
!              parts%dgu(d,k)=0.d0
!          enddo
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          if(parts%itype(i)*parts%itype(j)>0)cycle
              nx = (parts%x(2,j) - parts%x(2,j+1))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
              ny = (parts%x(1,j+1) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
              q0 = sqrt(((parts%x(1,i) - parts%x(1,j))*nx)**2 + ((parts%x(2,i) - parts%x(2,j))*ny)**2)/hsml
              r1 = sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) - parts%x(2,j))**2)
              r2 = sqrt((parts%x(1,i) - parts%x(1,j+1))**2 + (parts%x(2,i) - parts%x(2,j+1))**2)
              cos1 = -((parts%x(1,i) - parts%x(1,j)) * (parts%x(1,j+1) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j)) * (parts%x(2,j+1) - parts%x(2,j)))/   &
                 (sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) -parts%x(2,j))**2)*sqrt((parts%x(1,j+1) - parts%x(1,j))**2 + (parts%x(2,j+1) - parts%x(2,j))**2))
              cos2 = -((parts%x(1,i) - parts%x(1,j+1)) * (parts%x(1,j+1) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j+1)) * (parts%x(2,j+1) - parts%x(2,j)))/  &
                 (sqrt((parts%x(1,i) - parts%x(1,j+1))**2 + (parts%x(2,i) -parts%x(2,j+1))**2)*sqrt((parts%x(1,j+1) - parts%x(1,j))**2 + (parts%x(2,j+1) - parts%x(2,j))**2))
              y1 = r1*cos1
              y2 = r2*cos2
              q00 = min(q0,2.)
              zeta(1) = min(sqrt(q00**2 + y1**2/hsml**2)/2,1.)
              zeta(2) = min(sqrt(q00**2 + y2**2/hsml**2)/2,1.)
              parts%gammaas(k) = -1./(4*pi)*(sign(1.,y2)*parts%psi(q0,zeta(2)) - sign(1.,y1)*parts%psi(q0,zeta(1)))
              parts%gammaa%r(i) = parts%gammaa%r(i) + parts%gammaas(k)     !因为为负号的话，下面会>1、、、
      enddo           
 
      do i = 1, parts%ntotal
          parts%gammaa%r(i) = parts%gammaa%r(i) + 1.
          if(parts%gammaa%r(i)>1.1d0)write(*,*) 'gammaa2 has something wrong'
!          write(*,*)   parts%gammaa%r(60),parts%gammaa%r(120),parts%gammaa%r(180),parts%gammaa%r(240),parts%gammaa%r(300)
      enddo
!          write(*,*)   parts%gammaa%r(60),parts%gammaa%r(120),parts%gammaa%r(180),parts%gammaa%r(240),parts%gammaa%r(300)

      
      end subroutine

      
!----------------------------------------------------------------------
      subroutine get_gammaa4(parts)  !j和j-1
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer  i, j, k, d      
      real(dp)  q0,r1,r2,y1,y2,q00,zeta(2),nx,ny,hsml,cos1,cos2
      real(dp) , parameter :: pi = 3.1415926535898
      
      hsml = parts%hsml(1)
      do i = 1, parts%ntotal
          parts%gammaa%r(i) = 0.d0
      enddo
      
      do k = 1,parts%niac
!          do d = 1,parts%dim
!              parts%dgu(d,k)=0.d0
!          enddo
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          if(parts%itype(i)*parts%itype(j)>0)cycle
              nx = (parts%x(2,j) - parts%x(2,j-1))/sqrt((parts%x(1,j) - parts%x(1,j-1))**2 + (parts%x(2,j) - parts%x(2,j-1))**2)
              ny = (parts%x(1,j-1) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j-1))**2 + (parts%x(2,j) - parts%x(2,j-1))**2)
              q0 = sqrt(((parts%x(1,i) - parts%x(1,j))*nx)**2 + ((parts%x(2,i) - parts%x(2,j))*ny)**2)/hsml
              r1 = sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) - parts%x(2,j))**2)
              r2 = sqrt((parts%x(1,i) - parts%x(1,j-1))**2 + (parts%x(2,i) - parts%x(2,j-1))**2)
              cos1 = -((parts%x(1,i) - parts%x(1,j)) * (parts%x(1,j-1) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j)) * (parts%x(2,j-1) - parts%x(2,j)))/   &
                 (sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) -parts%x(2,j))**2)*sqrt((parts%x(1,j-1) - parts%x(1,j))**2 + (parts%x(2,j-1) - parts%x(2,j))**2))
              cos2 = -((parts%x(1,i) - parts%x(1,j-1)) * (parts%x(1,j-1) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j-1)) * (parts%x(2,j-1) - parts%x(2,j)))/  &
                 (sqrt((parts%x(1,i) - parts%x(1,j-1))**2 + (parts%x(2,i) -parts%x(2,j-1))**2)*sqrt((parts%x(1,j-1) - parts%x(1,j))**2 + (parts%x(2,j-1) - parts%x(2,j))**2))
              y1 = r1*cos1
              y2 = r2*cos2
              q00 = min(q0,2.)
              zeta(1) = min(sqrt(q00**2 + y1**2/hsml**2)/2,1.)
              zeta(2) = min(sqrt(q00**2 + y2**2/hsml**2)/2,1.)
              parts%gammaas(k) = -1./(4*pi)*(sign(1.,y2)*parts%psi(q0,zeta(2)) - sign(1.,y1)*parts%psi(q0,zeta(1)))
              parts%gammaa%r(i) = parts%gammaa%r(i) + parts%gammaas(k)     !因为为负号的话，下面会>1、、、
      enddo           
 
      do i = 1, parts%ntotal
          parts%gammaa%r(i) = parts%gammaa%r(i) + 1.
          if(parts%gammaa%r(i)>1.1d0)write(*,*) 'gammaa2 has something wrong'
      enddo
      
      
      
      end subroutine

      
!----------------------------------------------------------------------
      subroutine get_gammaa3(parts) 
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer  i, j, k, d ,m,k0,i0,j0     
      real(dp)  q0,r1,r2,y1,y2,q00,zeta(2),nx,ny,hsml,cos1,cos2,factor,dx
      integer , dimension(10) :: a  
      real(dp) , parameter :: pi = 3.1415926535898
      
      factor = 7.e0 / (4.e0*pi*hsml*hsml)
      dx = parts%hsml(1)/2.

!      do m = 1,10
!         a(m) = 30000
!      enddo
      
      hsml = parts%hsml(1)
      do i = 1, parts%ntotal
          parts%gammaa%r(i) = 0.d0
      enddo
      
      do k = 1,parts%niac
!          if(k == 37)then
!          do d = 1,parts%dim
!              parts%dgu(d,k)=0.d0
!          enddo
          i = parts%pair_i(k)
          j = parts%pair_j(k)
      do m = 1,10
         a(m) = 30000
      enddo
          m=0
          if(parts%itype(i)*parts%itype(j)>0)cycle
          if(parts%itype(i)>0)then
!              if(parts%zone(j) == 6)then
                 do k0 = 1, parts%niac
                     i0 = parts%pair_i(k0)
                     j0 = parts%pair_j(k0)
                     if(i0==i.and.j0>j)then
                         m = m+1
                         a(m) = j0
!                         write(*,*) 'a(m)',a(m),'m=',m
                     endif
                 enddo
!                 write(*,*) 'a(1)=',a(1),'a(2)=',a(2),'a(3)=',a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10)
                 j0 = min(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10))
                 if(j0/=30000)then
              nx = (parts%x(2,j) - parts%x(2,j0))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
              ny = (parts%x(1,j0) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
              q0 = sqrt(((parts%x(1,i) - parts%x(1,j))*nx)**2 + ((parts%x(2,i) - parts%x(2,j))*ny)**2)/hsml
              r1 = sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) - parts%x(2,j))**2)
              r2 = sqrt((parts%x(1,i) - parts%x(1,j0))**2 + (parts%x(2,i) - parts%x(2,j0))**2)
              cos1 = -((parts%x(1,i) - parts%x(1,j)) * (parts%x(1,j0) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j)) * (parts%x(2,j0) - parts%x(2,j)))/   &
                 (sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) -parts%x(2,j))**2)*sqrt((parts%x(1,j0) - parts%x(1,j))**2 + (parts%x(2,j0) - parts%x(2,j))**2))
              cos2 = -((parts%x(1,i) - parts%x(1,j0)) * (parts%x(1,j0) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j0)) * (parts%x(2,j0) - parts%x(2,j)))/  &
                 (sqrt((parts%x(1,i) - parts%x(1,j0))**2 + (parts%x(2,i) -parts%x(2,j0))**2)*sqrt((parts%x(1,j0) - parts%x(1,j))**2 + (parts%x(2,j0) - parts%x(2,j))**2))
              y1 = r1*cos1
              y2 = r2*cos2
              q00 = min(q0,2.)
              zeta(1) = min(sqrt(q00**2 + y1**2/hsml**2)/2,1.)
              zeta(2) = min(sqrt(q00**2 + y2**2/hsml**2)/2,1.)
              parts%gammaas(k) = -1./(4*pi)*(sign(1.,y2)*parts%psi(q0,zeta(2)) - sign(1.,y1)*parts%psi(q0,zeta(1)))
              parts%gammaa%r(i) = parts%gammaa%r(i) + parts%gammaas(k)     !因为为负号的话，下面会>1、、、
                 endif
          endif
!          endif
      enddo           

      do i = 1, parts%ntotal
          parts%gammaa%r(i) = parts%gammaa%r(i) + 1.
          if(parts%gammaa%r(i)>1.1d0)write(*,*) 'gammaa2 has something wrong'
      enddo
      
      
      
      end subroutine


      
      
      
      


!     wi(maxn)---integration of the kernel itself
      
!Subroutine to calculate the density with real_density_unified.
!----------------------------------------------------------------------
      subroutine real_density_unified(parts) 
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r,dvx(3) 
!      type(array),allocatable :: divvx
      real(dp), allocatable, dimension(:,:) :: dgua
      type(p2r) vx_i(3), vx_j(3) 

      ntotal = parts%ntotal + parts%nvirt
      allocate(dgua(parts%dim,parts%maxn))

!     Firstly calculate the integration of the kernel over the space

      do i=1,parts%ntotal
        parts%divvx%r(i)=0.d0 
      enddo

      do k=1,parts%niac
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
          do d=1,parts%dim
          dvx(d) = vx_i(d)%p - vx_j(d)%p
          enddo
           parts%divvx%r(i) = parts%divvx%r(i) - parts%mass%r(j)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))
           parts%divvx%r(j) = parts%divvx%r(j) - parts%mass%r(i)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))
          if(parts%itype(i)>0.and.parts%itype(j)<0)then
              if(parts%rhos(k)==0)cycle
              parts%divvx%r(i) = parts%divvx%r(i) + parts%rhos(k)*(parts%vx%x%r(i)*parts%dgu(1,k) + parts%vx%y%r(i)*parts%dgu(2,k))
          endif    
      enddo
      
      do i = 1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i) = 1
          parts%divvx%r(i) = parts%divvx%r(i)/(parts%rho%r(i) * (parts%gammaa%r(i)+10**-8))
      enddo

      end subroutine
      
!Subroutine to calculate the density with real_density_unified.
!----------------------------------------------------------------------
      subroutine real_density_unified2(parts)   !确保能量守恒的方式
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r,dvx(3) 
!      type(array),allocatable :: divvx
      real(dp), allocatable, dimension(:,:) :: dgua
      type(p2r) vx_i(3), vx_j(3) 

      ntotal = parts%ntotal + parts%nvirt
!      allocate(divvx);allocate(divvx%r(parts%maxn))
      allocate(dgua(parts%dim,parts%maxn))
!      divvx%ndim1 = ntotal

!     wi(maxn)---integration of the kernel itself
        
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
!     Firstly calculate the integration of the kernel over the space

      do i=1,parts%ntotal
        parts%divvx%r(i)=0.d0 
      enddo

      do k=1,parts%niac
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
          do d=1,parts%dim
          dvx(d) = vx_i(d)%p - vx_j(d)%p
          enddo
!          parts%drho%r(i) = parts%drho%r(i) + parts%mass%r(j)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))!书上的（4.34）
!          parts%drho%r(j) = parts%drho%r(j) + parts%mass%r(i)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))
           parts%divvx%r(i) = parts%divvx%r(i) + parts%mass%r(j)*parts%rho%r(i)*((-dvx(1))*parts%dwdx(1,k)+(-dvx(2))*parts%dwdx(2,k))
           parts%divvx%r(j) = parts%divvx%r(j) + parts%mass%r(i)*parts%rho%r(j)*(dvx(1)*(-parts%dwdx(1,k))+dvx(2)*(-parts%dwdx(2,k)))  !investigation of wall中的(38)
          if(parts%itype(i)>0.and.parts%itype(j)<0)then
              if(parts%rhos(k)==0)cycle
!              parts%drho%r(i) = parts%drho%r(i) - parts%rhos%r(j)*(parts%vx%x%r(i)*parts%dgu(1,k) + parts%vx%y%r(i)*parts%dgu(2,k))
              parts%divvx%r(i) = parts%divvx%r(i) - parts%rhos(k)*parts%rho%r(i)*((-parts%vx%x%r(i))*parts%dgu(1,k) + (-parts%vx%y%r(i))*parts%dgu(2,k))
          endif    
      enddo
      
      do i = 1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i) = 1
          parts%divvx%r(i) = parts%divvx%r(i)/(parts%rho%r(i)**2 * (parts%gammaa%r(i)))
      enddo

      end subroutine
      
!Subroutine to calculate the density with real_density_unified.
!----------------------------------------------------------------------
      subroutine freesurface2(parts)   !用Dr来求
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r ,dx(3),dxas(3)
      real(dp), allocatable, dimension(:,:) :: dgua
     type(array),allocatable :: divr

      ntotal = parts%ntotal + parts%nvirt
      allocate(dgua(parts%dim,parts%maxn))
      allocate(divr);allocate(divr%r(ntotal))
      divr%ndim1 = ntotal
        
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
!     Firstly calculate the integration of the kernel over the space

      do i=1,parts%ntotal
        divr%r(i)=0.d0 
      enddo

      do k=1,parts%niac
          i = parts%pair_i(k)
          j = parts%pair_j(k)
!          vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
          do d=1,parts%dim
!          dvx(d) = vx_i(d)%p - vx_j(d)%p
          dx(d) = parts%x(d,i) - parts%x(d,j)
          dxas(d) = parts%x(d,i) - (parts%x(d,j) + parts%x(d,j+1))/2
          enddo
!          parts%drho%r(i) = parts%drho%r(i) + parts%mass%r(j)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))!书上的（4.34）
!          parts%drho%r(j) = parts%drho%r(j) + parts%mass%r(i)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))
           divr%r(i) = divr%r(i) + parts%mass%r(j)*parts%rho%r(i)*((-dx(1))*parts%dwdx(1,k)+(-dx(2))*parts%dwdx(2,k))
           divr%r(j) = divr%r(j) + parts%mass%r(i)*parts%rho%r(j)*(dx(1)*(-parts%dwdx(1,k))+dx(2)*(-parts%dwdx(2,k)))  !investigation of wall中的(38)
          if(parts%itype(i)>0.and.parts%itype(j)<0)then
              if(parts%rhos(k)==0)cycle
!              parts%drho%r(i) = parts%drho%r(i) - parts%rhos%r(j)*(parts%vx%x%r(i)*parts%dgu(1,k) + parts%vx%y%r(i)*parts%dgu(2,k))
              divr%r(i) = divr%r(i) - parts%rhos(k)*parts%rho%r(i)*((-dxas(1))*parts%dgu(1,k) + (-dxas(2))*parts%dgu(2,k))
          endif    
      enddo
      
      do i = 1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i) = 1
          divr%r(i) = divr%r(i)/(parts%rho%r(i)**2 * (parts%gammaa%r(i)))
          if(divr%r(i)<1.5)parts%p%r(i)=0.!98
      enddo
      
      end subroutine
      

!Subroutine to calculate momentum equation by MLS.
!----------------------------------------------------------------------
      subroutine momentum_equation_unified(parts)  !dpre 用的Unified wall boundary 中的（28）
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt)
      type(p2r) vx_i(3), vx_j(3) 
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))
      
      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)*(parts%p%r(i)/parts%rho%r(i)**2+parts%p%r(j)/parts%rho%r(j)**2)*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)*(parts%p%r(j)/parts%rho%r(j)**2+parts%p%r(i)/parts%rho%r(i)**2)*parts%dwdx(d,k) !!这里是用-号吧？???
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - parts%rhos(k)*(parts%p%r(i)/parts%rho%r(i)**2+parts%ps(k)/parts%rhos(k)**2)*parts%dgu(d,k)
          enddo
        endif
      enddo
      
      do i =1,parts%ntotal
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)*parts%rho%r(i)/parts%gammaa%r(i)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
!      do i = parts%ntotal + 1, parts%ntotal+parts%nvirt
!          parts%p%r(i) = 0.
!      enddo


      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
           rr = 0.      
!       if(parts%x(1,i)>0.00012.and.parts%x(1,j)<0.00012.and.parts%itype(i)*parts%itype(j)>0)then
!           do d=1,parts%dim
!              dx(d) =  parts%x(d,i) - (parts%x(d,j) + 0.00025)
!              dvx(d) = vx_i(d)%p - vx_j(d)%p
!              rr = rr + dx(d)*dx(d)
!           enddo  
!       elseif(parts%x(1,i)<0.00012.and.parts%x(1,j)>0.00012.and.parts%itype(i)*parts%itype(j)>0)then
!           do d=1,parts%dim
!              dx(d) =  parts%x(d,i) - (parts%x(d,j) - 0.00025)
!              dvx(d) = vx_i(d)%p - vx_j(d)%p
!              rr = rr + dx(d)*dx(d)
!           enddo  
!       else
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) - parts%x(d,j)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
!              lap(d,i) = lap(d,i) + parts%mass%r(j)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/rr**2*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)    !这里的dwdx是对i吧？？？(33)
!              lap(d,j) = lap(d,j) - parts%mass%r(i)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/rr**2*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)    !改变一下优先级顺序
              lap(d,i) = lap(d,i) + parts%mass%r(j)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/rr**2*(dx(1)*parts%dwdx(1,k)+dx(2)*parts%dwdx(2,k))*dvx(d)    !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) - parts%mass%r(i)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/rr**2*(dx(1)*parts%dwdx(1,k)+dx(2)*parts%dwdx(2,k))*dvx(d)    !改变一下优先级顺序
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            do d = 1,parts%dim
                lap(d,i) = lap(d,i) - sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)*((water%viscosity*parts%divvx%r(i))*parts%n(d,j) + parts%mons(d,k))/parts%rho%r(i) 
            enddo
        endif
!        if(parts%itype(i)>0.and.parts%itype(j)<0)then
!            n(j-parts%ntotal) = n(j-parts%ntotal) + 1 
!        parts%p%r(j) = parts%p%r(j) + parts%p%r(i) - dpre(1,i) * dx(1) - dpre(2,i) *dx(2)
!        endif
!       endif
       
      enddo
    
!      do j = parts%ntotal + 1,parts%ntotal + parts%nvirt
!          if(n(j-parts%ntotal)==0) n(j-parts%ntotal)=1
!          parts%p%r(j)=parts%p%r(j)/n(j-parts%ntotal)
!          write(*,*) parts%p%r(j)
!      enddo
      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/parts%gammaa%r(i)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
!         parts%dvx%x%r(i) = lap(1,i)!这里的lap之前就除以rho了
!         parts%dvx%y%r(i) = 0. !这个地方的重力不知道是加号还是减号。。。。。。。。。。，         
      enddo
      
      end subroutine
      

      
!Subroutine to calculate momentum equation by MLS.
!----------------------------------------------------------------------
      subroutine momentum_equation_unified2(parts)  !dpre用的是Investigation of wall 中的能量守恒式，（37）
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt)
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))
      
      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo
        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/(parts%gammaa%r(i)+10**-8)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
     
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              dvx(1) =parts%vx%x%r(i) - parts%vx%x%r(j)
              dvx(2) =parts%vx%y%r(i) - parts%vx%y%r(j)
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2+0.01*parts%hsml(1)**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)    !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) - parts%mass%r(i)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2+0.01*parts%hsml(1)**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            do d = 1,parts%dim
                lap(d,i) = lap(d,i) - sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)*((water%viscosity*parts%divvx%r(i))*parts%n(d,j) + parts%mons(d,k))/parts%rho%r(i) 
            enddo
        endif
      enddo
      
     
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/(parts%gammaa%r(i)+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
      
!----------------------------------------------------------------------
      subroutine momentum_equation_unified201(parts)  !dpre用的是Investigation of wall 中的能量守恒式，（37）
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt)
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))
      
      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo
        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/(parts%gammaa%r(i)+10**-8)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
      do i = parts%ntotal + 1, parts%ntotal+parts%nvirt
          parts%p%r(i) = 0.
      enddo
      
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              dvx(1) =parts%vx%x%r(i) - parts%vx%x%r(j)
              dvx(2) =parts%vx%y%r(i) - parts%vx%y%r(j)
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2+0.01*parts%hsml(1)**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)    !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) - parts%mass%r(i)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2+0.01*parts%hsml(1)**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            do d = 1,parts%dim
                lap(d,i) = lap(d,i) - sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)*((water%viscosity*parts%divvx%r(i))*parts%n(d,j) + parts%mons(d,k))/parts%rho%r(i) 
            enddo
        endif
        
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            n(j-parts%ntotal) = n(j-parts%ntotal) + 1 
        parts%p%r(j) = parts%p%r(j) + parts%p%r(i) - dpre(1,i) * dx(1) - dpre(2,i) *dx(2)
        endif
      enddo
      
      do j = parts%ntotal + 1,parts%ntotal + parts%nvirt
          if(n(j-parts%ntotal)==0) n(j-parts%ntotal)=1
          parts%p%r(j)=parts%p%r(j)/n(j-parts%ntotal)
      enddo
      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/(parts%gammaa%r(i)+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
      

!Subroutine to calculate momentum equation by MLS.
!----------------------------------------------------------------------
      subroutine momentum_equation_unified3(parts)  !pre用的是Investigation of wall 中的能量守恒式，（37）,lap用的(18)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt)
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))

      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
          enddo
      enddo
      

      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo
        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/parts%gammaa%r(i)!+10**-8)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              rr = rr + dx(d)*dx(d)
           enddo  
              dvx(1) =parts%vx%x%r(i) - parts%vx%x%r(j)
              dvx(2) =parts%vx%y%r(i) - parts%vx%y%r(j)
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)    !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) - parts%mass%r(i)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            lapsx = 2.*sqrt(parts%dgu(1,k)**2 + parts%dgu(2,k)**2)*parts%mons(1,k)/parts%rho%r(i)
            lapsy = 2.*sqrt(parts%dgu(1,k)**2 + parts%dgu(2,k)**2)*parts%mons(2,k)/parts%rho%r(i)            
              lap(1,i) = lap(1,i) - lapsx
              lap(2,i) = lap(2,i) - lapsy
        endif
      enddo
      

      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/parts%gammaa%r(i)!+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
      
!----------------------------------------------------------------------
      subroutine momentum_equation_unified4(parts)  !pre用的是Investigation of wall 中的能量守恒式，（37）,lap用的DNS and LES中的(14)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt),rav1,rav2,zb
      type(p2r) vx_i(3), vx_j(3) 
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))

      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
!      write(*,*) 'this is momentum_equation_unified3'
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo

        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/parts%gammaa%r(i)!+10**-8)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)/parts%rho%r(j)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*parts%dwdx(d,k)) + (parts%dwdx(1,k) * dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * dvx(d))   !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) + parts%mass%r(i)/parts%rho%r(i)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*(-parts%dwdx(d,k))) + (parts%dwdx(1,k) *dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * (-dvx(d)))
           enddo
           
      if(parts%itype(i)>0.and.parts%itype(j)<0)then
           if(parts%rhos(k)==0) cycle
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2) - parts%hsml(1)/4
         lapsx = 2* sqrt(parts%dgu(1,k)**2 + parts%dgu(2,k)**2)*water%viscosity* (dvx(1) - parts%n(1,j)*(dvx(1)*parts%n(1,j) + dvx(2)*parts%n(2,j)))/(parts%rho%r(i)*zb) 
         lapsy = 2* sqrt(parts%dgu(1,k)**2 + parts%dgu(2,k)**2)*water%viscosity* (dvx(2) - parts%n(2,j)*(dvx(1)*parts%n(1,j) + dvx(2)*parts%n(2,j)))/(parts%rho%r(i)*zb)
             lap(1,i) = lap(1,i) - lapsx
             lap(2,i) = lap(2,i) - lapsy
       endif
           
           
           
!        if(parts%itype(i)>0.and.parts%itype(j)<0)then
!            if(parts%rhos(k)==0) cycle
!           lapsx = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(1,k)/parts%rhos(k)
!           lapsy = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(2,k)/parts%rhos(k)            
!             lap(1,i) = lap(1,i) - lapsx
!             lap(2,i) = lap(2,i) - lapsy
!        endif
      enddo

      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/parts%gammaa%r(i)!+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
       
!----------------------------------------------------------------------
      subroutine momentum_equation_unified_Couette(parts)  !pre用的是Investigation of wall 中的能量守恒式，（37）,lap用的DNS and LES中的(14)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt),y,zb,rav1,rav2
      type(p2r) vx_i(3), vx_j(3) 
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))

      do i = 1,parts%ntotal
          do d = 1,parts%dim
!              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
!      write(*,*) 'this is momentum_equation_unified3'


      
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
!              write(*,*)'dx='           
       if(parts%x(1,i)>0.00078.and.parts%x(1,j)<0.00022.and.parts%itype(i)*parts%itype(j)>0)then
!             write(*,*)'dx=',dx(1)
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) - (parts%x(d,j) + 0.001)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
       elseif(parts%x(1,i)<0.00022.and.parts%x(1,j)>0.00078.and.parts%itype(i)*parts%itype(j)>0)then
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) - (parts%x(d,j) - 0.001)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
       else
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) - parts%x(d,j)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
       endif
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)/parts%rho%r(j)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2&
             * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*parts%dwdx(d,k)) + (parts%dwdx(1,k) * dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * dvx(d))   !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) + parts%mass%r(i)/parts%rho%r(i)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*(-parts%dwdx(d,k))) + (parts%dwdx(1,k) *dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * (-dvx(d)))
           enddo
           
!      if(parts%itype(i)>0.and.parts%itype(j)<0)then
!           if(parts%rhos(k)==0) parts%rhos(k) =1000
!           y = abs(parts%x(2,i) - parts%x(2,j))-0.0000125
!         lapsx = 2*(water%viscosity* (vx_i(1)%p - vx_j(1)%p)/(parts%rho%r(i)*y))* sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)   
!         lap(1,i) = lap(1,i) - lapsx 

!w         lapsx = ((water%viscosity* (vx_i(1)%p - vx_j(1)%p)/(parts%rho%r(i)*y)) + parts%mons(1,k))* sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)
!w         lapsx =  2.*sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(1,k)/parts%rhos(k)
!w         lapsx =  sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * (parts%mons(1,k) + parts%rho%r(i)*water%viscosity*parts%vx%x%r(j))/parts%rhos(k)        
!w         lapsx =  sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * (parts%mons(1,k) + parts%rho%r(i)*water%viscosity*(parts%vx%x%r(i) - parts%vx%x%r(j)))/parts%rhos(k)        
!w         lap(1,i) = lap(1,i) - sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)*((water%viscosity*parts%divvx%r(i))*parts%n(1,j) + parts%mons(1,k))/parts%rho%r(i)       
!w         lap(1,i) =lap(1,i) + lapsx

!      endif

       if(parts%itype(i)>0.and.parts%itype(j)<0)then
           if(parts%rhos(k)==0) cycle
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2) !- parts%hsml(1)/4
         lapsx = 2* sqrt(parts%dgu(1,k)**2 + parts%dgu(2,k)**2)*water%viscosity* (dvx(1) - parts%n(1,j)*(dvx(1)*parts%n(1,j) + dvx(2)*parts%n(2,j)))/(parts%rho%r(i)*zb) 
!         lapsy = 2* sqrt(parts%dgu(1,k)**2 + parts%dgu(2,k)**2)*water%viscosity* (dvx(2) - parts%n(2,j)*(dvx(1)*parts%n(1,j) + dvx(2)*parts%n(2,j)))/(parts%rho%r(i)*zb)
             lap(1,i) = lap(1,i) - lapsx
!             lap(2,i) = lap(2,i) - lapsy
       endif
      
       
      enddo

      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/parts%gammaa%r(i)!+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) =  lap(1,i)!这里的lap之前就除以rho了
!         parts%dvx%y%r(i) = 0.!-dpre(2,i)/parts%rho%r(i) + lap(2,i) !+ parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
      
!----------------------------------------------------------------------
      subroutine momentum_equation_unified401(parts)  !pre用的是Investigation of wall 中的能量守恒式，（37）,lap用的DNS and LES中的(14)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt)
      type(p2r) vx_i(3), vx_j(3) 
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))
      n = 0.
      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
!      write(*,*) 'this is momentum_equation_unified3'
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo

        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/parts%gammaa%r(i)!+10**-8)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
      
      do i = parts%ntotal + 1, parts%ntotal+parts%nvirt !加上这个就是用DELTA P来求P了 
          parts%p%r(i) = 0.
      enddo
      
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)/parts%rho%r(j)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*parts%dwdx(d,k)) + (parts%dwdx(1,k) * dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * dvx(d))   !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) + parts%mass%r(i)/parts%rho%r(i)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*(-parts%dwdx(d,k))) + (parts%dwdx(1,k) *dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * (-dvx(d)))
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            lapsx = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(1,k)/parts%rhos(k)
            lapsy = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(2,k)/parts%rhos(k)            
              lap(1,i) = lap(1,i) - lapsx
              lap(2,i) = lap(2,i) - lapsy
        endif
        if(parts%itype(i)>0.and.parts%itype(j)<0)then   !用delta P的方式求
            n(j-parts%ntotal) = n(j-parts%ntotal) + 1 
        parts%p%r(j) = parts%p%r(j) + parts%p%r(i) - dpre(1,i) * dx(1) - dpre(2,i) *dx(2)
        endif
      enddo
    
      do j = parts%ntotal + 1,parts%ntotal + parts%nvirt!用delta P的方式求
          if(n(j-parts%ntotal)==0) n(j-parts%ntotal)=1
          parts%p%r(j)=parts%p%r(j)/n(j-parts%ntotal)
      enddo
      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/parts%gammaa%r(i)!+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
      

!----------------------------------------------------------------------
      subroutine momentum_equation_unified402(parts)  !pre用的是Investigation of wall 中的能量守恒式，（37）,lap用的DNS and LES中的(14),ddpre失败了
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
      real(dp), allocatable, dimension(:) :: lappre
      !      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt),lappreg
      type(p2r) vx_i(3), vx_j(3) 
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(lappre(parts%max_interaction))
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))
      n = 0.
      do i = 1,parts%ntotal
          lappre(i) = 0.d0
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
!      write(*,*) 'this is momentum_equation_unified3'
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo

        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/parts%gammaa%r(i)!+10**-8)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
           
!      do i = parts%ntotal + 1, parts%ntotal+parts%nvirt !加上这个就是用DELTA P来求P了 
!          parts%p%r(i) = 0.
!      enddo
      
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
              lappre(i) = lappre(i) + parts%mass%r(j)/parts%rho%r(j)* 2*(parts%p%r(i) - parts%p%r(j))/(rr**2)*(parts%dwdx(1,k)*dx(1) &
              + parts%dwdx(2,k)*dx(2)) 
              lappre(j) = lappre(j) + parts%mass%r(i)/parts%rho%r(i)* 2*(parts%p%r(j) - parts%p%r(i))/(rr**2)*(parts%dwdx(1,k)*dx(1) &
              + parts%dwdx(2,k)*dx(2))
              
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
              lappreg = 2. *(dpre(1,i) * parts%dgu(1,k) + dpre(2,i) * parts%dgu(2,k)) 
              lappre(i) = lappre(i) - lappreg 
        endif
      enddo
      
      do i = parts%ntotal + 1, parts%ntotal+parts%nvirt !加上这个就是用DELTA P来求P了 
          parts%p%r(i) = 0.
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)/parts%rho%r(j)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*parts%dwdx(d,k)) + (parts%dwdx(1,k) * dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * dvx(d))   !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) + parts%mass%r(i)/parts%rho%r(i)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*(-parts%dwdx(d,k))) + (parts%dwdx(1,k) *dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * (-dvx(d)))   
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            lapsx = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(1,k)/parts%rhos(k)
            lapsy = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(2,k)/parts%rhos(k)        
              lap(1,i) = lap(1,i) - lapsx
              lap(2,i) = lap(2,i) - lapsy
        endif
        if(parts%itype(i)>0.and.parts%itype(j)<0)then   !用delta P的方式求
            n(j-parts%ntotal) = n(j-parts%ntotal) + 1 
        parts%p%r(j) = parts%p%r(j) + parts%p%r(i) - dpre(1,i) * dx(1) - dpre(2,i) *dx(2) - lappre(i)/parts%gammaa%r(i)*rr**2
        endif
      enddo
    
      do j = parts%ntotal + 1,parts%ntotal + parts%nvirt!用delta P的方式求
          if(n(j-parts%ntotal)==0) n(j-parts%ntotal)=1
          parts%p%r(j)=parts%p%r(j)/n(j-parts%ntotal)
!    if(parts%p%r(j)/=0.)       write(*,*) parts%p%r(j)
      enddo
      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/parts%gammaa%r(i)!+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
      
!----------------------------------------------------------------------
      subroutine momentum_equation_unified5(parts)  !pre用的是Investigation of wall 中的能量守恒式，（37）,lap用的SPH in fluid dynamics中的(2.147)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,zb,rav1,rav2,das,n(parts%nvirt)
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
 !     allocate(divvx);allocate(divvx%r(ntotal))
 !     divvx%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))
      
      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
!      write(*,*) 'this is momentum_equation_unified3'
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo
        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/parts%gammaa%r(i)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
!      do i = parts%ntotal + 1, parts%ntotal+parts%nvirt !加上这个就是用DELTA P来求P了 
!          parts%p%r(i) = 0.
!      enddo
      
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              dvx(1) =parts%vx%x%r(i) - parts%vx%x%r(j)
              dvx(2) =parts%vx%y%r(i) - parts%vx%y%r(j)
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)    !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) - parts%mass%r(i)*water%viscosity*2/(parts%rho%r(i)*parts%rho%r(j))/(rr**2)*(dvx(1)*dx(1)+dvx(2)*dx(2))*parts%dwdx(d,k)
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2)
            das = max(parts%hsml(1)/2.,zb)
            lapsx = 2.*water%viscosity/parts%rho%r(i)*parts%vx%x%r(i)* sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)/das
            lapsy = 2.*water%viscosity/parts%rho%r(i)*parts%vx%y%r(i)* sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2)/das
              lap(1,i) = lap(1,i) - lapsx
              lap(2,i) = lap(2,i) - lapsy
        endif
!        if(parts%itype(i)>0.and.parts%itype(j)<0)then   !用delta P的方式求
!            n(j-parts%ntotal) = n(j-parts%ntotal) + 1 
!        parts%p%r(j) = parts%p%r(j) + parts%p%r(i) - dpre(1,i) * dx(1) - dpre(2,i) *dx(2)
!        endif
      enddo
      
!      do j = parts%ntotal + 1,parts%ntotal + parts%nvirt!用delta P的方式求
!          if(n(j-parts%ntotal)==0) n(j-parts%ntotal)=1
!          parts%p%r(j)=parts%p%r(j)/n(j-parts%ntotal)
!      enddo
      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= 0!lap(d,i)/(parts%gammaa%r(i)+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
      

!----------------------------------------------------------------------
      subroutine momentum_equation_unified6(parts)  !pre用的是Investigation of wall 中的能量守恒式，（37）,lap用的DNS and LES中的(14)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
      integer ntotal, i, j, k, d    
      real(dp), allocatable, dimension(:,:) :: dpre,lap
!      type(array),allocatable :: divvx
      type(array),allocatable :: dpre2x,dpre2y
      real(dp) selfdens, hv(3), r ,rr,dx(3),dvx(3),lapa,lapsx,lapsy,n(parts%nvirt)
      type(p2r) vx_i(3), vx_j(3) 
      
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
      parts%dvx%x%ndim1 = ntotal
      
 !     allocate(divvx);allocate(divvx%r(ntotal))
 
      allocate(dpre2x);allocate(dpre2x%r(ntotal))
      dpre2x%ndim1 = ntotal
      allocate(dpre2y);allocate(dpre2y%r(ntotal))
      dpre2y%ndim1 = ntotal
      allocate(dpre(parts%dim,parts%max_interaction))
      allocate(lap(parts%dim,parts%max_interaction))
!      allocate(laps(parts%dim,parts%max_interaction))
      n = 0.
      do i = 1,parts%ntotal
          do d = 1,parts%dim
              dpre(d,i) = 0.d0
              lap(d,i) = 0.d0
!              laps(d,i) = 0.d0
          enddo
      enddo
      
!      write(*,*) 'this is momentum_equation_unified3'
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) + parts%mass%r(j)/parts%rho%r(j)*(parts%rho%r(i)**2 * parts%p%r(j) + parts%rho%r(j)**2 * parts%p%r(i))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)   
           dpre(d,j) = dpre(d,j) - parts%mass%r(i)/parts%rho%r(i)*(parts%rho%r(j)**2 * parts%p%r(i) + parts%rho%r(i)**2 * parts%p%r(j))  &
                       /(parts%rho%r(i)*parts%rho%r(j))*parts%dwdx(d,k)
          enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            if(parts%rhos(k)==0)cycle
          do d = 1,parts%dim
           dpre(d,i) = dpre(d,i) - (parts%rho%r(i)**2 * parts%ps(k) + parts%rhos(k)**2 * parts%p%r(i))/(parts%rho%r(i)*parts%rhos(k))*parts%dgu(d,k)
          enddo

        endif
      enddo
      
      do i =1,parts%ntotal
!          if(parts%gammaa%r(i)==0) parts%gammaa%r(i) =1
          do d= 1,parts%dim
              dpre(d,i) = dpre(d,i)/parts%gammaa%r(i)!+10**-8)  !文章中的（28)式，还需要进步到（8）
          enddo
      enddo
      
      dpre2x = -parts%df(parts%p,'x')
      dpre2y = -parts%df(parts%p,'y')
      
      do i =1,parts%ntotal
         dpre(1,i) = dpre2x%r(i)
         dpre(2,i) = dpre2y%r(i)
      enddo
      
!      divvx  =  parts.div2(parts%vx)
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        if(parts%itype(i)>0.and.parts%itype(j)>0)then
           rr = 0.      
           do d=1,parts%dim
              dx(d) =  parts%x(d,i) -  parts%x(d,j)
              dvx(d) = vx_i(d)%p - vx_j(d)%p
              rr = rr + dx(d)*dx(d)
           enddo  
           rr = sqrt(rr)
           do d = 1,parts%dim
              lap(d,i) = lap(d,i) + parts%mass%r(j)/parts%rho%r(j)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*parts%dwdx(d,k)) + (parts%dwdx(1,k) * dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * dvx(d))   !这里的dwdx是对i吧？？？(33)
              lap(d,j) = lap(d,j) + parts%mass%r(i)/parts%rho%r(i)*(water%viscosity/parts%rho%r(i) + water%viscosity/parts%rho%r(j))/2  &
              * (((parts%dim + 2)/(rr**2)*(dvx(1)*dx(1) + dvx(2)*dx(2))*(-parts%dwdx(d,k))) + (parts%dwdx(1,k) *dx(1) &
              + parts%dwdx(2,k) * dx(2))/(rr**2) * (-dvx(d)))
           enddo
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            lapsx = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(1,k)/parts%rhos(k)
            lapsy = 2. * sqrt(parts%dgu(1,k)**2+parts%dgu(2,k)**2) * parts%mons(2,k)/parts%rhos(k)            
              lap(1,i) = lap(1,i) - lapsx
              lap(2,i) = lap(2,i) - lapsy
        endif
      enddo
      
      do i=1, parts%ntotal
!          if(parts%gammaa%r(i) == 0) parts%gammaa%r(i)=1
          do d=1, parts%dim
!             lap(d,i)= (lap(d,i)-laps(d,i)/parts%rho%r(i))/wi(i)
             lap(d,i)= lap(d,i)/parts%gammaa%r(i)!+10**-8)!这里用的（27）其实就已经除以了密度，是作业里面的（14）中的第二项
          enddo
      enddo
      
      do i = 1,parts%ntotal
         parts%dvx%x%r(i) = -dpre(1,i)/parts%rho%r(i) + lap(1,i)!这里的lap之前就除以rho了
         parts%dvx%y%r(i) = -dpre(2,i)/parts%rho%r(i) + lap(2,i) + parts%numeric%gravity  !这个地方的重力不知道是加号还是减号。。。。。。。。。。，
      enddo
      
      end subroutine
       
 
! Subroutine to calculate the delta gamma of unified condition.
!----------------------------------------------------------------------      
      subroutine delta_gamma_unified2(parts)    !算是自创的方法
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      
       integer k,i,j,d,m,k0,i0,j0
      real(dp) q0,q1,q2,q1c,q2c,Pq2,Pq1,factor,hsml,dx,rr,w,q,nx,ny,dr
      real(dp) , parameter :: pi = 3.1415926535898
      integer , dimension(40) :: a  
      hsml = parts%hsml(1)
      
      factor = 7.e0 / (4.e0*pi*hsml*hsml)
      dx = parts%hsml(1)/2.

      do m = 1,40
         a(m) = 30000
      enddo
            
      do k = 1, parts%niac
        parts%dgu(1,k) = 0.
        parts%dgu(2,k) = 0.
        parts%rhos(k) = 0.
        parts%mons(1,k) = 0.
        parts%mons(2,k) = 0.
        parts%ps(k) = 0.
          do d = 1,parts%dim
              parts%dgu(d,k)=0.d0
          enddo

          i = parts%pair_i(k)
          j = parts%pair_j(k)
          m=0
          if(parts%itype(i)*parts%itype(j)>0)cycle
          if(parts%itype(i)>0)then
!              if(parts%zone(j) == 6)then
                 do k0 = 1, parts%niac
                     i0 = parts%pair_i(k0)
                     j0 = parts%pair_j(k0)
                     if(i0==i.and.j0>j)then
                         m = m+1
                         a(m) = j0
!                         write(*,*) 'a(m)',a(m),'m=',m
                     endif
                 enddo
!                 write(*,*) 'a(1)=',a(1),'a(2)=',a(2),'a(3)=',a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10)
              j0 = min(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15),a(16),a(17),a(18),a(19),a(20))
              j0 = min(j0,a(21),a(22),a(23),a(24),a(25),a(26),a(27),a(28),a(29),a(30),a(31),a(32),a(33),a(34),a(35),a(36),a(37),a(38),a(39),a(40)) 
                 if(j0/=30000)then
!                    write(*,*) 'j0=',j0
                     rr = sqrt((parts%x(1,i) - parts%x(1,j0))**2 + (parts%x(2,i) - parts%x(2,j0))**2)
                     dr = sqrt((parts%x(1,j0) - parts%x(1,j))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                     q = rr/hsml
                     w = factor * ( (1-q/2)**4 *(1+2*q) )
                     nx = (parts%x(2,j) - parts%x(2,j0))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                     ny = (parts%x(1,j0) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                     parts%dgu(1,k) = (w + parts%w(k))/2 * dr * nx
                     parts%dgu(2,k) = (w + parts%w(k))/2 * dr * ny
                     parts%rhos(k) = (parts%rho%r(j) + parts%rho%r(j0))/2
                     parts%mons(1,k) = (parts%mone%x%r(j) + parts%mone%x%r(j0))/2
                     parts%mons(2,k) = (parts%mone%y%r(j) + parts%mone%y%r(j0))/2
                     parts%ps(k) = (parts%p%r(j)/parts%rho%r(j) + parts%p%r(j0)/parts%rho%r(j0))/2*parts%rhos(k)
!                 endif
                 a = 30000
              endif
          endif
      enddo
      
      end subroutine
      
! Subroutine to calculate the delta gamma of unified condition.
!----------------------------------------------------------------------      
      subroutine delta_gamma_unified3(parts)  !全部都是j与j+1
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      
      integer k,i,j,d,k0,i0,j0,m
      real(dp) q0,q1,q2,q1c,q2c,Pq2,Pq1,hsml,nx,ny,cos1,cos2
      real(dp) , parameter :: pi = 3.1415926535898
      integer , dimension(10) :: a  
      hsml = parts%hsml(1)
      
      do k = 1, parts%niac
          do d = 1,parts%dim
              parts%dgu(d,k)=0.d0
          enddo
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          if(parts%itype(i)*parts%itype(j)>0)cycle
              parts%dgu(1,k) = parts%n(1,j)*parts%w(k) * hsml/2
              parts%dgu(2,k) = parts%n(2,j)*parts%w(k) * hsml/2
              parts%rhos(k) = (parts%rho%r(j) + parts%rho%r(j+1))/2
              parts%mons(1,k) = (parts%mone%x%r(j) + parts%mone%x%r(j+1))/2
              parts%mons(2,k) = (parts%mone%y%r(j) + parts%mone%y%r(j+1))/2
              parts%ps(k) = (parts%p%r(j)/parts%rho%r(j) + parts%p%r(j+1)/parts%rho%r(j+1))/2*parts%rhos(k)

      enddo
      
      end subroutine
! Subroutine to calculate the delta gamma of unified condition.
!----------------------------------------------------------------------      
      subroutine delta_gamma_unified6(parts)  !全部都是j与j-1
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      
      integer k,i,j,d,k0,i0,j0,m
      real(dp) q0,q1,q2,q1c,q2c,Pq2,Pq1,hsml,nx,ny,cos1,cos2
      real(dp) , parameter :: pi = 3.1415926535898
      integer , dimension(10) :: a  
      hsml = parts%hsml(1)
      
      do k = 1, parts%niac
          do d = 1,parts%dim
              parts%dgu(d,k)=0.d0
          enddo
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          if(parts%itype(i)*parts%itype(j)>0)cycle
              q0 = sqrt(((parts%x(1,i) - parts%x(1,j-1))*parts%n(1,j-1))**2 + ((parts%x(2,i) - parts%x(2,j-1))*parts%n(2,j-1))**2)/hsml
              q1 = sqrt((parts%x(1,i) - parts%x(1,j-1))**2 + (parts%x(2,i) - parts%x(2,j-1))**2)/hsml
              q2 = sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) - parts%x(2,j))**2)/hsml
              cos1 = -((parts%x(1,i) - parts%x(1,j-1)) * (parts%x(1,j) - parts%x(1,j-1)) + (parts%x(2,i) - parts%x(2,j-1)) * (parts%x(2,j) - parts%x(2,j-1)))/   &
                 (sqrt((parts%x(1,i) - parts%x(1,j-1))**2 + (parts%x(2,i) -parts%x(2,j-1))**2)*sqrt((parts%x(1,j) - parts%x(1,j-1))**2 + (parts%x(2,j) - parts%x(2,j-1))**2))
              cos2 = -((parts%x(1,i) - parts%x(1,j)) * (parts%x(1,j) - parts%x(1,j-1)) + (parts%x(2,i) - parts%x(2,j)) * (parts%x(2,j) - parts%x(2,j-1)))/  &
                 (sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) -parts%x(2,j))**2)*sqrt((parts%x(1,j) - parts%x(1,j-1))**2 + (parts%x(2,j) - parts%x(2,j-1))**2))
              q1c = q1*cos1
              q2c = q2*cos2
              Pq2 = 7./192.*q2**5-21./64.*q2**4+35./32.*q2**3-35./24.*q2**2+7./4.+q0**2*(35./768.*q2**3-7./16.*q2**2+105./64.*q2-35./12.)+q0**4*(35./512.*q2-7./8.)
              Pq1 = 7./192.*q1**5-21./64.*q1**4+35./32.*q1**3-35./24.*q1**2+7./4.+q0**2*(35./768.*q1**3-7./16.*q1**2+105./64.*q1-35./12.)+q0**4*(35./512.*q1-7./8.)
              parts%dgu(1,k) = parts%n(1,j-1)*(q2c/pi*Pq2-q1c/pi*Pq1+q0**4/pi*(105./64.+35./512.*q0**2)*(sign(1.d0,q2c)*log((q2+abs(q2c))/abs(q0))-sign(1.d0,q1c)*log((q1+abs(q1c))/abs(q0))))/parts%hsml(j-1)
              parts%dgu(2,k) = parts%n(2,j-1)*(q2c/pi*Pq2-q1c/pi*Pq1+q0**4/pi*(105./64.+35./512.*q0**2)*(sign(1.d0,q2c)*log((q2+abs(q2c))/abs(q0))-sign(1.d0,q1c)*log((q1+abs(q1c))/abs(q0))))/parts%hsml(j-1)
              parts%rhos(k) = (parts%rho%r(j-1) + parts%rho%r(j))/2
              parts%mons(1,k) = (parts%mone%x%r(j-1) + parts%mone%x%r(j))/2
              parts%mons(2,k) = (parts%mone%y%r(j-1) + parts%mone%y%r(j))/2
              parts%ps(k) = (parts%p%r(j-1)/parts%rho%r(j-1) + parts%p%r(j)/parts%rho%r(j))/2*parts%rhos(k)
!          endif

      enddo
      
      end subroutine
      
! Subroutine to calculate the delta gamma of unified condition.
!----------------------------------------------------------------------      
      subroutine delta_gamma_unified4(parts)  !只用最标准的方案 
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      
      integer k,i,j,d,k0,i0,j0,m
      real(dp) q0,q1,q2,q1c,q2c,Pq2,Pq1,hsml,nx,ny,cos1,cos2
      real(dp) , parameter :: pi = 3.1415926535898
      integer , dimension(40) :: a  
      hsml = parts%hsml(1)
      
      do m = 1,40
         a(m) = 30000
      enddo
      do k = 1, parts%niac
          do d = 1,parts%dim
              parts%dgu(d,k)=0.d0
          enddo
          m=0
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          if(parts%itype(i)*parts%itype(j)>0)cycle
!          if(parts%zone(j)==6)then
              do k0 = 1, parts%niac
                 i0 = parts%pair_i(k0)
                 j0 = parts%pair_j(k0)
                    if(i0==i.and.j0>j)then
                       m = m+1
                       a(m) = j0
                    endif
              enddo
              j0 = min(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15),a(16),a(17),a(18),a(19),a(20))
              j0 = min(j0,a(21),a(22),a(23),a(24),a(25),a(26),a(27),a(28),a(29),a(30),a(31),a(32),a(33),a(34),a(35),a(36),a(37),a(38),a(39),a(40)) 
              if(j0/=30000)then
                nx = (parts%x(2,j) - parts%x(2,j0))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                ny = (parts%x(1,j0) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                q0 = sqrt(((parts%x(1,i) - parts%x(1,j))*nx)**2 + ((parts%x(2,i) - parts%x(2,j))*ny)**2)/hsml
                q1 = sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) - parts%x(2,j))**2)/hsml
                q2 = sqrt((parts%x(1,i) - parts%x(1,j0))**2 + (parts%x(2,i) - parts%x(2,j0))**2)/hsml
                cos1 = -((parts%x(1,i) - parts%x(1,j)) * (parts%x(1,j0) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j)) * (parts%x(2,j0) - parts%x(2,j)))/   &
                     (sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) -parts%x(2,j))**2)*sqrt((parts%x(1,j0) - parts%x(1,j))**2 + (parts%x(2,j0) - parts%x(2,j))**2))
                cos2 = -((parts%x(1,i) - parts%x(1,j0)) * (parts%x(1,j0) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j0)) * (parts%x(2,j0) - parts%x(2,j)))/  &
                     (sqrt((parts%x(1,i) - parts%x(1,j0))**2 + (parts%x(2,i) -parts%x(2,j0))**2)*sqrt((parts%x(1,j0) - parts%x(1,j))**2 + (parts%x(2,j0) - parts%x(2,j))**2))
                q1c = q1*cos1
                q2c = q2*cos2
                Pq2 = 7./192.*q2**5-21./64.*q2**4+35./32.*q2**3-35./24.*q2**2+7./4.+q0**2*(35./768.*q2**3-7./16.*q2**2+105./64.*q2-35./12.)+q0**4*(35./512.*q2-7./8.)
                Pq1 = 7./192.*q1**5-21./64.*q1**4+35./32.*q1**3-35./24.*q1**2+7./4.+q0**2*(35./768.*q1**3-7./16.*q1**2+105./64.*q1-35./12.)+q0**4*(35./512.*q1-7./8.)
                parts%dgu(1,k) = nx*(q2c/pi*Pq2-q1c/pi*Pq1+q0**4/pi*(105./64.+35./512.*q0**2)*(sign(1.d0,q2c)*log((q2+abs(q2c))/abs(q0))-sign(1.d0,q1c)*log((q1+abs(q1c))/abs(q0))))/parts%hsml(j)
                parts%dgu(2,k) = ny*(q2c/pi*Pq2-q1c/pi*Pq1+q0**4/pi*(105./64.+35./512.*q0**2)*(sign(1.d0,q2c)*log((q2+abs(q2c))/abs(q0))-sign(1.d0,q1c)*log((q1+abs(q1c))/abs(q0))))/parts%hsml(j)
                parts%rhos(k) = (parts%rho%r(j) + parts%rho%r(j0))/2
                parts%mons(1,k) =(parts%mone%x%r(j) + parts%mone%x%r(j0))/2
                parts%mons(2,k) =(parts%mone%y%r(j) + parts%mone%y%r(j0))/2
                parts%ps(k) = (parts%p%r(j)/parts%rho%r(j) + parts%p%r(j0)/parts%rho%r(j0))/2*parts%rhos(k)
             endif
                a = 30000
!          else 
!              nx = (parts%x(2,j) - parts%x(2,j+1))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
!              ny = (parts%x(1,j+1) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j+1))**2 + (parts%x(2,j) - parts%x(2,j+1))**2)
!              q0 = sqrt(((parts%x(1,i) - parts%x(1,j))*parts%n(1,j))**2 + ((parts%x(2,i) - parts%x(2,j))*parts%n(2,j))**2)/hsml
!              q1 = sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) - parts%x(2,j))**2)/hsml
!              q2 = sqrt((parts%x(1,i) - parts%x(1,j+1))**2 + (parts%x(2,i) - parts%x(2,j+1))**2)/hsml
!              cos1 = -((parts%x(1,i) - parts%x(1,j)) * (parts%x(1,j+1) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j)) * (parts%x(2,j+1) - parts%x(2,j)))/   &
!                 (sqrt((parts%x(1,i) - parts%x(1,j))**2 + (parts%x(2,i) -parts%x(2,j))**2)*sqrt((parts%x(1,j+1) - parts%x(1,j))**2 + (parts%x(2,j+1) - parts%x(2,j))**2))
!              cos2 = -((parts%x(1,i) - parts%x(1,j+1)) * (parts%x(1,j+1) - parts%x(1,j)) + (parts%x(2,i) - parts%x(2,j+1)) * (parts%x(2,j+1) - parts%x(2,j)))/  &
!                 (sqrt((parts%x(1,i) - parts%x(1,j+1))**2 + (parts%x(2,i) -parts%x(2,j+1))**2)*sqrt((parts%x(1,j+1) - parts%x(1,j))**2 + (parts%x(2,j+1) - parts%x(2,j))**2))
!              q1c = q1*cos1
!              q2c = q2*cos2
!              Pq2 = 7./192.*q2**5-21./64.*q2**4+35./32.*q2**3-35./24.*q2**2+7./4.+q0**2*(35./768.*q2**3-7./16.*q2**2+105./64.*q2-35./12.)+q0**4*(35./512.*q2-7./8.)
!              Pq1 = 7./192.*q1**5-21./64.*q1**4+35./32.*q1**3-35./24.*q1**2+7./4.+q0**2*(35./768.*q1**3-7./16.*q1**2+105./64.*q1-35./12.)+q0**4*(35./512.*q1-7./8.)
!              parts%dgu(1,k) = parts%n(1,j)*(q2c/pi*Pq2-q1c/pi*Pq1+q0**4/pi*(105./64.+35./512.*q0**2)*(sign(1.d0,q2c)*log((q2+abs(q2c))/abs(q0))-sign(1.d0,q1c)*log((q1+abs(q1c))/abs(q0))))/parts%hsml(j)
!              parts%dgu(2,k) = parts%n(2,j)*(q2c/pi*Pq2-q1c/pi*Pq1+q0**4/pi*(105./64.+35./512.*q0**2)*(sign(1.d0,q2c)*log((q2+abs(q2c))/abs(q0))-sign(1.d0,q1c)*log((q1+abs(q1c))/abs(q0))))/parts%hsml(j)
!              parts%rhos(k) = (parts%rho%r(j) + parts%rho%r(j+1))/2
!              parts%mons(1,k) = (parts%mone%x%r(j) + parts%mone%x%r(j+1))/2
!              parts%mons(2,k) = (parts%mone%y%r(j) + parts%mone%y%r(j+1))/2
!              parts%ps(k) = (parts%p%r(j)/parts%rho%r(j) + parts%p%r(j+1)/parts%rho%r(j+1))/2*parts%rhos(k)
!          endif

      enddo
      
      end subroutine

      
      
      
!----------------------------------------------------------------------      
      subroutine delta_gamma_unified5(parts)  
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      
       integer k,i,j,d,m,k0,i0,j0
      real(dp) q0,q1,q2,q1c,q2c,Pq2,Pq1,factor,hsml,dx,rr,w,q,nx,ny,dr
      real(dp) , parameter :: pi = 3.1415926535898
      integer , dimension(10) :: a  
      hsml = parts%hsml(1)
      
      factor = 7.e0 / (4.e0*pi*hsml*hsml)
      dx = hsml/2.

      do m = 1,10
         a(m) = 3000
      enddo
            
      do k = 1, parts%niac
          do d = 1,parts%dim
              parts%dgu(d,k)=0.d0
          enddo
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          m=0
          if(parts%itype(i)*parts%itype(j)>0)cycle
          if(parts%itype(i)>0)then
!              if(parts%zone(j) == 6)then
                 do k0 = 1, parts%niac
                     i0 = parts%pair_i(k0)
                     j0 = parts%pair_j(k0)
                     if(i0==i.and.j0>j)then
                         m = m+1
                         a(m) = j0
!                         write(*,*) 'a(m)',a(m),'m=',m
                     endif
                 enddo
!                 write(*,*) 'a(1)=',a(1),'a(2)=',a(2),'a(3)=',a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10)
                 j0 = min(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10))
                 if(j0/=3000)then
!                    write(*,*) 'j0=',j0
                     rr = sqrt((parts%x(1,i) - parts%x(1,j0))**2 + (parts%x(2,i) - parts%x(2,j0))**2)
                     dr = sqrt((parts%x(1,j0) - parts%x(1,j))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                     q = rr/hsml
                     w = factor * ( (1-q/2)**4 *(1+2*q) )
                     nx = (parts%x(2,j) - parts%x(2,j0))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                     ny = (parts%x(1,j0) - parts%x(1,j))/sqrt((parts%x(1,j) - parts%x(1,j0))**2 + (parts%x(2,j) - parts%x(2,j0))**2)
                     parts%dgu(1,k) = (w + parts%w(k))/2 * dr * nx
                     parts%dgu(2,k) = (w + parts%w(k))/2 * dr * ny
                     parts%rhos(k) = (parts%rho%r(j) + parts%rho%r(j0))/2
                     parts%mons(1,k) = (parts%mone%x%r(j) + parts%mone%x%r(j0))/2
                     parts%mons(2,k) = (parts%mone%y%r(j) + parts%mone%y%r(j0))/2
                     parts%ps(k) = (parts%p%r(j)/parts%rho%r(j) + parts%p%r(j0)/parts%rho%r(j0))/2*parts%rhos(k)
!                 endif
                 a = 3000
              endif
          endif
      enddo
      
      end subroutine
      
      
!Subroutine to calculate the density with SPH summation algorithm by MLS.
!----------------------------------------------------------------------
      subroutine sum_density_MLS(parts) 
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d, m, n,k0,p,s1,s2
      real(dp) selfdens, hv(3), r, W
      real(dp), allocatable, dimension(:) :: wi,phi0,phi1,phi2
      real(dp), allocatable, dimension(:,:,:) :: jz
      real(dp) jz2(3,3),l(3,3),u(3,3),lt(3,3),inverse_jz2(3,3),inverse_u(3,3),inverse_l(3,3),inverse_lt(3,3)
      integer i1,j1,k1,n1
      real(dp) sum1
     
      
      ntotal = parts%ntotal + parts%nvirt

      allocate(wi(parts%maxn))
      allocate(phi0(parts%maxn))
      allocate(phi1(parts%maxn))
      allocate(phi2(parts%maxn))
      allocate(jz(parts%maxn,3,3))
      
!     wi(maxn)---integration of the kernel itself
        
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
!     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        do m=1,3
            do n=1,3
                jz(i,m,n)=0
            enddo
        enddo
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        jz(i,1,1)=selfdens*parts%mass%r(i)/parts%rho%r(i)
      enddo
      
        do k0=1,parts%niac
          i = parts%pair_i(k0)
          j = parts%pair_j(k0)
          
          jz(i,1,1)=jz(i,1,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0) 
          jz(j,1,1)=jz(j,1,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)   
          jz(i,2,1)=jz(i,2,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))
          jz(j,2,1)=jz(j,2,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))
          jz(i,3,1)=jz(i,3,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(2,i)-parts%x(2,j))
          jz(j,3,1)=jz(j,3,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(2,j)-parts%x(2,i))
          jz(i,1,2)=jz(i,2,1)
          jz(j,1,2)=jz(j,2,1)
          jz(i,2,2)=jz(i,2,2) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))**2
          jz(j,2,2)=jz(j,2,2) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))**2
          jz(i,3,2)=jz(i,3,2) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))*(parts%x(2,i)-parts%x(2,j))
          jz(j,3,2)=jz(j,3,2) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))*(parts%x(2,j)-parts%x(2,i))
          jz(i,1,3)=jz(i,3,1)
          jz(j,1,3)=jz(j,3,1)
          jz(i,2,3)=jz(i,3,2)
          jz(j,2,3)=jz(j,3,2)
          jz(i,3,3)=jz(i,3,3) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(2,i)-parts%x(2,j))**2
          jz(j,3,3)=jz(j,3,3) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(2,j)-parts%x(2,i))**2    
        enddo    
        do i=1,ntotal
            p=0
          do m=1,3
            do n =1,3
                jz2(m,n )=jz(i,m,n)
!                if(jz2(m,n)==0)then
!                p=p+1
!                write(*,*) 'p=',p
!                endif
            enddo
          enddo
!         W = jz2(1,1)*jz2(2,2)*jz2(3,3)+jz2(1,2)*jz2(2,3)*jz2(3,1)+jz2(1,3)*jz2(2,3)*jz2(3,2)-jz2(3,1)*jz2(2,2)*jz2(1,3)-jz2(2,3)*jz2(3,2)*jz2(1,1)-jz2(1,2)*jz2(2,3)*jz2(3,3)
!         if(W/=0)then 
!            write(*,*)' A=0','i=',i
!             write(*,*)'1=',jz2(1,1),jz2(1,2),jz2(1,3)
!             write(*,*)'2=',jz2(2,1),jz2(2,2),jz2(2,3)
!             write(*,*)'3=',jz2(3,1),jz2(3,2),jz2(3,3)
             
!             write(*,*)'1=',jz(59,1,1),jz(59,1,2),jz(59,1,3)
!             write(*,*)'2=',jz(59,2,1),jz(59,2,2),jz(59,2,3)
!             write(*,*)'3=',jz(59,3,1),jz(59,3,2),jz(59,3,3)
             
 !            do k0 = 1, parts%niac
 !                s1 = parts%pair_i(k0)
 !                s2 = parts%pair_j(k0)
 !                if(s1==i.or.s2==i) then
 !                    p = p+1
 !                endif
 !            enddo
 !            write(*,*) 'p=',p
 !        endif
         call LU_break(jz2,3,l,u)
         call inverse_uptri_matrix(u,3,inverse_u)
         call transpose_matrix(l,3,3,lt)
         call inverse_uptri_matrix(lt,3,inverse_lt)
         call transpose_matrix(inverse_lt,3,3,inverse_l)
         call multiply_matrix(inverse_u,3,3,inverse_l,3,inverse_jz2)
         if((inverse_jz2(1,1)>-10**10.and.inverse_jz2(1,1)<10**10).and.(inverse_jz2(2,1) >-10**10.and.inverse_jz2(2,1) <10**10).and.(inverse_jz2(3,1)>-10**10.and.inverse_jz2(3,1) <10**10))then
         phi0(i) = inverse_jz2(1,1)  
         phi1(i) = inverse_jz2(2,1)
         phi2(i) = inverse_jz2(3,1) 
         else
         phi0(i) = 1
         phi1(i) = 0 
         phi2(i) = 0
         endif
         
       enddo 
        
        
      do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        parts%rho%r(i) = phi0(i)*selfdens*parts%mass%r(i)
      enddo
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        parts%rho%r(i) = parts%rho%r(i) + parts%mass%r(j)*(phi0(i)+phi1(i)*(parts%x(1,i)-parts%x(1,j))+phi2(i)*(parts%x(2,i)-parts%x(2,j)))*parts%w(k)
        parts%rho%r(j) = parts%rho%r(j) + parts%mass%r(i)*(phi0(j)+phi1(j)*(parts%x(1,j)-parts%x(1,i))+phi2(j)*(parts%x(2,j)-parts%x(2,i)))*parts%w(k)
      enddo

!      do i=1,ntotal
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
!        wi(i) = selfdens*parts%mass%r(i)/parts%rho%r(i)
!      enddo
!      do k=1,parts%niac
!        i = parts%pair_i(k)
!        j = parts%pair_j(k)
!        wi(i) = wi(i) + parts%mass%r(j)/parts%rho%r(j)*(phi0(i)+phi1(i)*(parts%x(1,i)-parts%x(1,j))+phi2(i)*(parts%x(2,i)-parts%x(2,j)))*parts%w(k)
!        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*(phi0(j)+phi1(j)*(parts%x(1,j)-parts%x(1,i))+phi2(j)*(parts%x(2,j)-parts%x(2,i)))*parts%w(k)
!      enddo
     
!        do i=1, ntotal
!          parts%rho%r(i)=parts%rho%r(i)/wi(i)
!        enddo
 
      end subroutine

      

!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine nvirt_density_unified(parts) !rho使用1阶修正,但是p用的是弱压缩性
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%maxn),vxi,vxj,dx,mu1,mu2,zb,rav1,rav2
      type(p2r) vx_i(3), vx_j(3) 
      type(material),pointer :: water
      type(array),allocatable :: rho
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
      mu1 = 0
      mu2 = 1
      
      allocate(rho);allocate(rho%r(ntotal))
      rho%ndim1 = ntotal
!     Firstly calculate the integration of the kernel over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)  !由于rho不能要自加，这里也要自加
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
        endif
      enddo
      
!     Secondly calculate the rho integration over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        parts%rho%r(i) = selfdens*parts%mass%r(i)   !rho不能为0
        rho%r(i) = 0       
        parts%mone%x%r(i) = 0
        parts%mone%y%r(i) = 0
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            rho%r(j) = rho%r(j) + parts%mass%r(i)*parts%w(k)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2)
                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%x%r(i)/zb*parts%w(k)
                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%y%r(i)/zb*parts%w(k)
        endif
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0)wi(i)=1
          rho%r(i) = rho%r(i)/wi(i)
          parts%mone%x%r(i) = parts%mone%x%r(i)/wi(i)
          parts%mone%y%r(i) = parts%mone%y%r(i)/wi(i)
!这个不注释的          if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
        enddo

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
            parts%p%r(i) =0.d0
        enddo
        
water => parts%material
do i = parts%ntotal + 1,parts%ntotal + parts%nvirt
parts%p%r(i) = water%b*((parts%rho%r(i)/(water%rho0))**water%gamma-1.d0)
enddo
        
        
      end subroutine
      
      

!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine nvirt_density_unified2(parts) !rho使用1阶修正，p用第一篇论文中的方法
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%maxn),vxi,vxj,dx,mu1,mu2,zb,rav1,rav2
      type(p2r) vx_i(3), vx_j(3) 
      type(material),pointer :: water
      type(array),allocatable :: rho
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
!      mu1 = 0
!      mu2 = 1
      

      hv = 0.
      r = 0.
      allocate(rho);allocate(rho%r(ntotal))
      rho%ndim1 = ntotal
!     Firstly calculate the integration of the kernel over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)  !由于rho不能要自加，这里也要自加
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
        endif
      enddo
      
!     Secondly calculate the rho integration over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        parts%rho%r(i) = selfdens*parts%mass%r(i)   !rho不能为0
        rho%r(i) = 0.!selfdens*parts%mass%r(i)       
        parts%mone%x%r(i) = 0
        parts%mone%y%r(i) = 0
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            rho%r(j) = rho%r(j) + parts%mass%r(i)*parts%w(k)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2)
!                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%x%r(i)/zb*parts%w(k)
!                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%y%r(i)/zb*parts%w(k)
                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%x%r(i)*parts%w(k)
                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%y%r(i)*parts%w(k)   !换一种方式
        endif
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0) wi(i)=1
          rho%r(i) = rho%r(i)/wi(i)
          parts%mone%x%r(i) = parts%mone%x%r(i)/wi(i)
          parts%mone%y%r(i) = parts%mone%y%r(i)/wi(i)
          if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
        enddo

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
            parts%p%r(i) =0.d0
        enddo
        
        do k=1,parts%niac
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            if(parts%itype(i)>0.and.parts%itype(j)<0)then
                vxi=0.d0
                vxj=0.d0
                vx_i = parts%vx%cmpt(i)
                vx_j = parts%vx%cmpt(j)
                do d=1, parts%dim
                  vxi = vxi + vx_i(d)%p*vx_i(d)%p  
                  vxj = vxj + vx_j(d)%p*vx_j(d)%p
                enddo     
                  dx = parts%x(2,i) - parts%x(2,j) 
                  parts%p%r(j) = parts%p%r(j) + (parts%p%r(i)/parts%rho%r(i) - parts%numeric%gravity*dx + (vxi - vxj)/2)*parts%w(k)*parts%mass%r(i)/parts%rho%r(i)
            endif!上面把g*dx的符号从负号改成了正号改了一下,感觉这个地方，还真的只能是负号
        enddo
        
        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0) wi(i)=1
          parts%p%r(i)=parts%p%r(i)/wi(i)
        enddo
        
        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          parts%p%r(i)=parts%p%r(i)*parts%rho%r(i)
        enddo
        

        
      end subroutine
      
!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine nvirt_density_unified2_Couette(parts) !rho使用1阶修正，p用第一篇论文中的方法
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%maxn),vxi,vxj,dx,mu1,mu2,zb,rav1,rav2
      type(p2r) vx_i(3), vx_j(3) 
      type(material),pointer :: water
      type(array),allocatable :: rho
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
!      mu1 = 0
!      mu2 = 1
      

      hv = 0.
      r = 0.
      allocate(rho);allocate(rho%r(ntotal))
      rho%ndim1 = ntotal
!     Firstly calculate the integration of the kernel over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)  !由于rho不能要自加，这里也要自加
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
        endif
      enddo
      
!     Secondly calculate the rho integration over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        parts%rho%r(i) = selfdens*parts%mass%r(i)   !rho不能为0
!        rho%r(i) = 0.!selfdens*parts%mass%r(i)       
        parts%mone%x%r(i) = 0
        parts%mone%y%r(i) = 0
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
!            rho%r(j) = rho%r(j) + parts%mass%r(i)*parts%w(k)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2) - 0.0000125
!               parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%x%r(i)/zb*parts%w(k)
!               parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%y%r(i)/zb*parts%w(k)   ! 这种方式很不好，这是1的方式
                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)*water%viscosity*(parts%vx%x%r(i) - parts%vx%x%r(j))*parts%w(k)
                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)*water%viscosity*(parts%vx%y%r(i) - parts%vx%y%r(j))*parts%w(k)   !换一种方式
!                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%x%r(i)*parts%w(k)
!                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%y%r(i)*parts%w(k)   !换一种方式
!                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)*water%viscosity*(-parts%vx%x%r(j))*parts%w(k)
!                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)*water%viscosity*(-parts%vx%y%r(j))*parts%w(k)   !换一种方式

        endif
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0) wi(i)=1
!          rho%r(i) = rho%r(i)/wi(i)
          parts%mone%x%r(i) = parts%mone%x%r(i)/wi(i)
          parts%mone%y%r(i) = parts%mone%y%r(i)/wi(i)
!          if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
        enddo



        
      end subroutine
      
!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine nvirt_density_unified201(parts) !rho使用1阶修正，p用第一篇论文中的方法
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%maxn),vxi,vxj,dx,mu1,mu2,zb,rav1,rav2
      type(p2r) vx_i(3), vx_j(3) 
      type(material),pointer :: water
      type(array),allocatable :: rho
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
!      mu1 = 0
!      mu2 = 1
      
      allocate(rho);allocate(rho%r(ntotal))
      rho%ndim1 = ntotal
!     Firstly calculate the integration of the kernel over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)  !由于rho不能要自加，这里也要自加
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
        endif
      enddo
      
!     Secondly calculate the rho integration over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        parts%rho%r(i) = selfdens*parts%mass%r(i)   !rho不能为0
        rho%r(i) = 0       
        parts%mone%x%r(i) = 0
        parts%mone%y%r(i) = 0
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            rho%r(j) = rho%r(j) + parts%mass%r(i)*parts%w(k)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2)
                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%x%r(i)/zb*parts%w(k)
                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%y%r(i)/zb*parts%w(k)
        endif
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0) wi(i)=1
          rho%r(i) = rho%r(i)/wi(i)
          parts%mone%x%r(i) = parts%mone%x%r(i)/wi(i)
          parts%mone%y%r(i) = parts%mone%y%r(i)/wi(i)
          if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
        enddo

!        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!            parts%p%r(i) =0.d0
!        enddo
        
!        do k=1,parts%niac
!            i = parts%pair_i(k)
!            j = parts%pair_j(k)
!            if(parts%itype(i)>0.and.parts%itype(j)<0)then
!                vxi=0.d0
!                vxj=0.d0
!                vx_i = parts%vx%cmpt(i)
!                vx_j = parts%vx%cmpt(j)
!                do d=1, parts%dim
!                  vxi = vxi + vx_i(d)%p*vx_i(d)%p  
!                  vxj = vxj + vx_j(d)%p*vx_j(d)%p
!                enddo     
!                  dx = parts%x(2,i) - parts%x(2,j) 
!                  parts%p%r(j) = parts%p%r(j) + (parts%p%r(i)/parts%rho%r(i) - parts%numeric%gravity*dx + (vxi - vxj)/2)*parts%w(k)*parts%mass%r(i)/parts%rho%r(i)
!            endif!上面把g*dx的符号从负号改成了正号改了一下,感觉这个地方，还真的只能是负号

            
            
 !       enddo
        
!        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!          if(wi(i)==0) wi(i)=1
!          parts%p%r(i)=parts%p%r(i)/wi(i)
!        enddo
        
!        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!          parts%p%r(i)=parts%p%r(i)*parts%rho%r(i)
!        enddo
        

        
      end subroutine 
      
!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine nvirt_density_unified202(parts) !rho使用1阶修正，p用第一篇论文中的方法
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d, m, n,k0,p,s1,s2,i1,j1,k1,n1
      real(dp) selfdens, hv(3), r,vxi,vxj,dx,mu1,mu2,zb,rav1,rav2, W, sum1
      type(p2r) vx_i(3), vx_j(3) 
      real(dp), allocatable, dimension(:) :: wi,phi0,phi1,phi2
      real(dp), allocatable, dimension(:,:,:) :: jz
      real(dp) jz2(3,3),l(3,3),u(3,3),lt(3,3),inverse_jz2(3,3),inverse_u(3,3),inverse_l(3,3),inverse_lt(3,3)
      type(material),pointer :: water
      type(array),allocatable :: rho
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
!      mu1 = 0
!      mu2 = 1
      

      allocate(wi(parts%maxn))
      allocate(phi0(parts%maxn))
      allocate(phi1(parts%maxn))
      allocate(phi2(parts%maxn))
      allocate(jz(parts%maxn,3,3))
      
!     wi(maxn)---integration of the kernel itself
        
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
      allocate(rho);allocate(rho%r(ntotal))
      rho%ndim1 = ntotal
      
!     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        do m=1,3
            do n=1,3
                jz(i,m,n)=0
            enddo
        enddo
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        jz(i,1,1)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)
      enddo
      
        do k0=1,parts%niac
          i = parts%pair_i(k0)
          j = parts%pair_j(k0)
          
          jz(i,1,1)=jz(i,1,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0) 
          jz(j,1,1)=jz(j,1,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)   
          jz(i,2,1)=jz(i,2,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))
          jz(j,2,1)=jz(j,2,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))
          jz(i,3,1)=jz(i,3,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(2,i)-parts%x(2,j))
          jz(j,3,1)=jz(j,3,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(2,j)-parts%x(2,i))
          jz(i,1,2)=jz(i,2,1)
          jz(j,1,2)=jz(j,2,1)
          jz(i,2,2)=jz(i,2,2) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))**2
          jz(j,2,2)=jz(j,2,2) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))**2
          jz(i,3,2)=jz(i,3,2) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))*(parts%x(2,i)-parts%x(2,j))
          jz(j,3,2)=jz(j,3,2) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))*(parts%x(2,j)-parts%x(2,i))
          jz(i,1,3)=jz(i,3,1)
          jz(j,1,3)=jz(j,3,1)
          jz(i,2,3)=jz(i,3,2)
          jz(j,2,3)=jz(j,3,2)
          jz(i,3,3)=jz(i,3,3) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(2,i)-parts%x(2,j))**2
          jz(j,3,3)=jz(j,3,3) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(2,j)-parts%x(2,i))**2    
        enddo    
        do i=1,ntotal
            p=0
          do m=1,3
            do n =1,3
                jz2(m,n)=jz(i,m,n)
            enddo
          enddo
         call LU_break(jz2,3,l,u)
         call inverse_uptri_matrix(u,3,inverse_u)
         call transpose_matrix(l,3,3,lt)
         call inverse_uptri_matrix(lt,3,inverse_lt)
         call transpose_matrix(inverse_lt,3,3,inverse_l)
         call multiply_matrix(inverse_u,3,3,inverse_l,3,inverse_jz2)
         if((inverse_jz2(1,1)>-10**10.and.inverse_jz2(1,1)<10**10).and.(inverse_jz2(2,1) >-10**10.and.inverse_jz2(2,1) <10**10).and.(inverse_jz2(3,1)>-10**10.and.inverse_jz2(3,1) <10**10))then
         phi0(i) = inverse_jz2(1,1)  
         phi1(i) = inverse_jz2(2,1)
         phi2(i) = inverse_jz2(3,1) 
         else
         phi0(i) = 1
         phi1(i) = 0 
         phi2(i) = 0
         endif
       enddo 
        
        
      do i=parts%ntotal + 1, parts%ntotal + parts%nvirt
         wi(i) = 0.
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
          wi(j)= wi(j) + parts%mass%r(i)/parts%rho%r(i)*(phi0(j)+phi1(j)*(parts%x(1,j)-parts%x(1,i))+phi2(j)*(parts%x(2,j)-parts%x(2,i)))*parts%w(k)        
        endif
      enddo

      do i=1,ntotal
        do m=1,3
            do n=1,3
                jz(i,m,n)=0
            enddo
        enddo
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        jz(i,1,1)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)
      enddo
      
        do k0=1,parts%niac
          i = parts%pair_i(k0)
          j = parts%pair_j(k0)
          
          jz(i,1,1)=jz(i,1,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0) 
          jz(j,1,1)=jz(j,1,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)   
          jz(i,2,1)=jz(i,2,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))
          jz(j,2,1)=jz(j,2,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))
          jz(i,3,1)=jz(i,3,1) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(2,i)-parts%x(2,j))
          jz(j,3,1)=jz(j,3,1) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(2,j)-parts%x(2,i))
          jz(i,1,2)=jz(i,2,1)
          jz(j,1,2)=jz(j,2,1)
          jz(i,2,2)=jz(i,2,2) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))**2
          jz(j,2,2)=jz(j,2,2) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))**2
          jz(i,3,2)=jz(i,3,2) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(1,i)-parts%x(1,j))*(parts%x(2,i)-parts%x(2,j))
          jz(j,3,2)=jz(j,3,2) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(1,j)-parts%x(1,i))*(parts%x(2,j)-parts%x(2,i))
          jz(i,1,3)=jz(i,3,1)
          jz(j,1,3)=jz(j,3,1)
          jz(i,2,3)=jz(i,3,2)
          jz(j,2,3)=jz(j,3,2)
          jz(i,3,3)=jz(i,3,3) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k0)*(parts%x(2,i)-parts%x(2,j))**2
          jz(j,3,3)=jz(j,3,3) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k0)*(parts%x(2,j)-parts%x(2,i))**2    
        enddo    
        do i=1,ntotal
            p=0
          do m=1,3
            do n =1,3
                jz2(m,n )=jz(i,m,n)
            enddo
          enddo

         call LU_break(jz2,3,l,u)
         call inverse_uptri_matrix(u,3,inverse_u)
         call transpose_matrix(l,3,3,lt)
         call inverse_uptri_matrix(lt,3,inverse_lt)
         call transpose_matrix(inverse_lt,3,3,inverse_l)
         call multiply_matrix(inverse_u,3,3,inverse_l,3,inverse_jz2)
         if((inverse_jz2(1,1)>-10**10.and.inverse_jz2(1,1)<10**10).and.(inverse_jz2(2,1) >-10**10.and.inverse_jz2(2,1) <10**10).and.(inverse_jz2(3,1)>-10**10.and.inverse_jz2(3,1) <10**10))then
         phi0(i) = inverse_jz2(1,1)  
         phi1(i) = inverse_jz2(2,1)
         phi2(i) = inverse_jz2(3,1) 
         else
         phi0(i) = 1
         phi1(i) = 0 
         phi2(i) = 0
         endif
       enddo 
        
        
      do i=parts%ntotal + 1, parts%ntotal + parts%nvirt
         rho%r(i) = selfdens*parts%mass%r(i)/parts%rho%r(i)
      enddo
      
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
          rho%r(j)= rho%r(j) + parts%mass%r(i)*(phi0(j)+phi1(j)*(parts%x(1,j)-parts%x(1,i))+phi2(j)*(parts%x(2,j)-parts%x(2,i)))*parts%w(k)        
        endif
      enddo
      
!      wi(2144)=0.
!      wi(2187)=0.
!      wi(2219)=0.
!      wi(2264)=0.
!       rho%r(2144) =0.
!       rho%r(2187) = 0.
!       rho%r(2219) = 0.
!       rho%r(2264) = 0.
       
!      do k=1,parts%niac
!        i = parts%pair_i(k)
!        j = parts%pair_j(k)
!        if(parts%itype(i)<0.and.parts%itype(j)<0)cycle
!        if(j==2144) then
!        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
!        write(*,*)  wi(j)
!        endif
!      enddo
      
!      do k=1,parts%niac
!        i = parts%pair_i(k)
!        j = parts%pair_j(k)
!        if(parts%itype(i)<0.and.parts%itype(j)<0)cycle
!        if(j==2144)then
!            rho%r(j) = rho%r(j) + parts%mass%r(i)*parts%w(k)
!            write(*,*)  rho%r(j),parts%w(k)
!        endif
!      enddo
      
      
      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!         if(i==2144.or.i==2187.or.i=2219.or.i=2264)cycle
         if(wi(i) == 0)wi(i)=1.
         rho%r(i)= rho%r(i)/wi(i)
         if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
      enddo
      
!     Firstly calculate the integration of the kernel over the space

!     Secondly calculate the rho integration over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        parts%rho%r(i) = selfdens*parts%mass%r(i)   !rho不能为0
!        rho%r(i) = 0.!selfdens*parts%mass%r(i)       
        parts%mone%x%r(i) = 0
        parts%mone%y%r(i) = 0
        if(wi(i) == 0)wi(i)=1.
!         parts%rho%r(i) =  parts%rho%r(i)/wi(i)
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
!            rho%r(j) = rho%r(j) + parts%mass%r(i)*parts%w(k)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2)
!                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%x%r(i)/zb*parts%w(k)
!                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%y%r(i)/zb*parts%w(k)
                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%x%r(i)*parts%w(k)
                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%y%r(i)*parts%w(k)   !换一种方式
        endif
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0) wi(i)=1
!          rho%r(i) = rho%r(i)/wi(i)
          parts%mone%x%r(i) = parts%mone%x%r(i)/wi(i)
          parts%mone%y%r(i) = parts%mone%y%r(i)/wi(i)
!          if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
        enddo

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
            parts%p%r(i) =0.d0
        enddo

        do k=1,parts%niac
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            if(parts%itype(i)>0.and.parts%itype(j)<0)then
                vxi=0.d0
                vxj=0.d0
                vx_i = parts%vx%cmpt(i)
                vx_j = parts%vx%cmpt(j)
                do d=1, parts%dim
                  vxi = vxi + vx_i(d)%p*vx_i(d)%p  
                  vxj = vxj + vx_j(d)%p*vx_j(d)%p
                enddo     
                  dx = parts%x(2,i) - parts%x(2,j) 
                  parts%p%r(j) = parts%p%r(j) + (parts%p%r(i)/parts%rho%r(i) - parts%numeric%gravity*dx + (vxi - vxj)/2)*parts%w(k)*parts%mass%r(i)/parts%rho%r(i)
            endif!上面把g*dx的符号从负号改成了正号改了一下,感觉这个地方，还真的只能是负号
        enddo

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0) wi(i)=1
          parts%p%r(i)=parts%p%r(i)/wi(i)
        enddo

        
        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          parts%p%r(i)=parts%p%r(i)*parts%rho%r(i)
        enddo
        
!      wi(2144)=0.
!      wi(2187)=0.
!      wi(2219)=0.
!      wi(2264)=0.
        parts%rho%r(2144) = (parts%rho%r(2143) + parts%rho%r(2145))/2
        parts%rho%r(2187) = (parts%rho%r(2186) + parts%rho%r(2188))/2
        parts%rho%r(2219) = (parts%rho%r(2218) + parts%rho%r(2220))/2
        parts%rho%r(2264) = (parts%rho%r(2263) + parts%rho%r(2265))/2
        parts%p%r(2144) = (parts%p%r(2143) + parts%p%r(2145))/2
        parts%p%r(2187) = (parts%p%r(2186) + parts%p%r(2188))/2
        parts%p%r(2219) = (parts%p%r(2218) + parts%p%r(2220))/2
        parts%p%r(2264) = (parts%p%r(2263) + parts%p%r(2265))/2

        
      end subroutine
      
!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine nvirt_density_unified3(parts) !rho和p使用的都是Investigation of wall bounded flows中的一阶的方法
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%ntotal+parts%nvirt),vxi,vxj,dx,mu1,mu2,rav1,rav2,zb
      type(p2r) vx_i(3), vx_j(3) 
      type(material),pointer :: water
      type(array),allocatable :: rho
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt
      mu1 = 0
      mu2 = 1
      
      allocate(rho);allocate(rho%r(ntotal))
      
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
      rho%ndim1 = ntotal
!     Firstly calculate the integration of the kernel over the space

      do i=parts%ntotal + 1,parts%ntotal + parts%nvirt
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)  !由于rho不能要自加，这里也要自加
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
!        wi(i) = wi(i) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k)
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
        endif
      enddo

!     Secondly calculate the rho integration over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
        rho%r(i) = 0.     
        parts%mone%x%r(i) = 0.
        parts%mone%y%r(i) = 0.
        parts%p%r(i) = 0.! selfdens*parts%p%r(i)*parts%mass%r(i)/parts%rho%r(i)
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            rho%r(j) = rho%r(j) + parts%mass%r(i)*parts%w(k)
            parts%p%r(j) = parts%p%r(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k) * (parts%p%r(i) - parts%rho%r(i)*(rav2 - (rav1 * parts%n(1,j) + rav2 * parts%n(2,j))*parts%n(2,j))*(-9.81))            
!            parts%p%r(j) = parts%p%r(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k) * (parts%p%r(i) - parts%rho%r(i)*(rav1 * parts%n(1,j) + rav2 * parts%n(2,j))*parts%n(2,j)*(-9.81))            
!            parts%p%r(j) = parts%p%r(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k) * parts%p%r(i)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2)
!                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%x%r(i)/zb*parts%w(k)
!                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%y%r(i)/zb*parts%w(k)
                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%x%r(i)*parts%w(k)
                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)*water%viscosity*parts%vx%y%r(i)*parts%w(k)   !换一种方式
        endif
      enddo


!      if(j>2187)      write(*,*)parts%n(1,j),parts%n(2,j)
!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          if(wi(i)==0)wi(i)=1
          rho%r(i) = rho%r(i)/wi(i) 
          parts%mone%x%r(i) = parts%mone%x%r(i)/wi(i)
          parts%mone%y%r(i) = parts%mone%y%r(i)/wi(i)
          if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
        enddo
        
        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!          if(wi(i)==0)wi(i)=1
          parts%p%r(i)=parts%p%r(i)/wi(i)! + selfdens*parts%mass%r(i)/parts%rho%r(i) )
!          if(2187<i.and.i<2219)write(*,*)wi(i), parts%p%r(i),i
        enddo
        
      end subroutine
      

!Subroutine to calculate the density with SPH summation algorithm.
!----------------------------------------------------------------------
      subroutine nvirt_density_unified4(parts) !仅p使用2阶修正,然而没有成功过
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r, wi(parts%ntotal+parts%nvirt),vxi,vxj,dx,rav,rav1,rav2,tempp,zb
      real(dp) r1vw(parts%ntotal+parts%nvirt),r2vw(parts%ntotal+parts%nvirt),r3vw(parts%ntotal+parts%nvirt),r4vw(parts%ntotal+parts%nvirt)
      type(p2r) vx_i(3), vx_j(3) 
      type(material),pointer :: water
      type(array),allocatable :: rho
      water => parts%material
      ntotal = parts%ntotal + parts%nvirt

      
      allocate(rho);allocate(rho%r(ntotal))
      rho%ndim1 = ntotal
!     Firstly calculate the integration of the kernel over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=0.!selfdens*parts%mass%r(i)/parts%rho%r(i)  !由于rho不能要自加，这里也要自加
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
        endif
      enddo
!     Secondly calculate the rho integration over the space

      do i=parts%ntotal +1,parts%ntotal + parts%nvirt
        rho%r(i) = 0.       
        parts%mone%x%r(i) = 0.
        parts%mone%y%r(i) = 0.
      enddo

!     Calculate SPH sum for rho:
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            rav1 = parts%x(1,i) - parts%x(1,j)
            rav2 = parts%x(2,i) - parts%x(2,j)
            rho%r(j) = rho%r(j) + parts%mass%r(i) * parts%w(k)
            zb = abs(parts%n(1,j)*rav1 + parts%n(2,j)*rav2)
                parts%mone%x%r(j) = parts%mone%x%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%x%r(i)/zb*parts%w(k)
                parts%mone%y%r(j) = parts%mone%y%r(j) + parts%mass%r(i)/parts%rho%r(i)*water%viscosity*parts%vx%y%r(i)/zb*parts%w(k)
        endif
      enddo

!     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
          parts%p%r(i) = 0.
        enddo
         
      r3vw = 0.e0
      r1vw = 0.e0
      r4vw = 0.e0
!      vw = 0.e0
      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
        rav  = (parts%x(1,j)  - parts%x(1,i))*parts%n(1,j) + (parts%x(2,j)  - parts%x(2,i))*parts%n(2,j) 
        r2vw(j-parts%ntotal) = r2vw(j-parts%ntotal) + rav**2 * parts%mass%r(i)/parts%rho%r(i) * parts%w(k)
        r4vw(j-parts%ntotal) = r4vw(j-parts%ntotal) + rav**4 * parts%mass%r(i)/parts%rho%r(i) * parts%w(k)
        endif
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)>0.and.parts%itype(j)<0)then
        rav1 = parts%x(1,i) - parts%x(1,j)
        rav2 = parts%x(2,i) - parts%x(2,j)
        rav  = (parts%x(1,i)  - parts%x(1,j))*parts%n(1,j) + (parts%x(2,i)  - parts%x(2,j))*parts%n(2,j) 
!        parts%p%r(j) = parts%p%r(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k) * (parts%p%r(i) - parts%rho%r(i)*(rav2 - (rav1 * parts%n(1,j) +  &
!                       rav2 * parts%n(2,j))*parts%n(2,j))*(-9.81)) * ( r4vw(j-parts%ntotal) - r2vw(j-parts%ntotal) * rav**2 )
        parts%p%r(j) = parts%p%r(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k) * (parts%p%r(i)  - parts%rho%r(i)* (rav1 * parts%n(1,j) +  &
                       rav2 * parts%n(2,j))*parts%n(2,j)*(-9.81)) * ( r4vw(j-parts%ntotal) - r2vw(j-parts%ntotal) * rav**2 )         
                       
        endif
     enddo

!     Secondly calculate the rho integration over the space

!      tempp = 0.e0
        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
            tempp = r4vw(i-parts%ntotal)*wi(i) - r2vw(i-parts%ntotal)*r2vw(i-parts%ntotal) 
          if(wi(i)==0)wi(i)=1
            rho%r(i) = rho%r(i)/wi(i)
            parts%mone%x%r(i) = parts%mone%x%r(i)/wi(i)
            parts%mone%y%r(i) = parts%mone%y%r(i)/wi(i)
          if(rho%r(i)/=0)parts%rho%r(i) = rho%r(i)
          if(tempp==0)tempp=1
            parts%p%r(i) = parts%p%r(i)/tempp
        enddo
      
      end subroutine
      
      
      
      
! Subroutine to calculate the density with SPH continuiity approach.
!----------------------------------------------------------------------      
      subroutine con_density(parts)
!----------------------------------------------------------------------
      implicit none
  
      class(particles) parts      

      integer ntotal,i,j,k,d    
      double precision vcc, dvx(3), dwdx(3) 
      type(p2r) vx_i(3), vx_j(3) 
      
      ntotal = parts%ntotal + parts%nvirt

      do i = 1, ntotal
        parts%drho%r(i) = 0.
      enddo
     
      do k=1,parts%niac      
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1, parts%dim
          dwdx(d) = parts%dwdx(d,k) 
          dvx(d) = vx_i(d)%p - vx_j(d)%p 
        enddo        
        !vcc = dvx(1)* parts%dwdx(1,k)        
        vcc = dvx(1)*dwdx(1)        
        do d=2, parts%dim
          !vcc = vcc + dvx(d)*parts%dwdx(d,k)
          vcc = vcc + dvx(d)*dwdx(d)
        enddo    
        parts%drho%r(i) = parts%drho%r(i) + parts%mass%r(j)*vcc
        parts%drho%r(j) = parts%drho%r(j) + parts%mass%r(i)*vcc       
      enddo    

      end subroutine
      
      
!见光滑粒子动力学SPH方法应力不稳定性的一种改进方案
!----------------------------------------------------------------------
      subroutine freesurface(parts) 
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      integer ntotal, i, j, k, d      
      real(dp) selfdens, hv(3), r!, wi(parts%maxn)     
      real(dp), allocatable, dimension(:) :: wi
      type(material),pointer :: water
      water => parts%material
      allocate(wi(parts%maxn))
      ntotal = parts%ntotal + parts%nvirt

      
      hv = 0.d0

!     Self density of each particle: Wii (Kernel for distance 0)
!     and take contribution of particle itself:

      r=0.d0
      
!     wi(maxn)---integration of the kernel itself

      do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=selfdens*parts%mass%r(i)/parts%rho%r(i)
      enddo

      do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        wi(i) = wi(i) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k)
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
      enddo
      
      do i = 1,ntotal 
          if(wi(i).lt.0.9e0) then
          parts%rho%r(i) = water%rho0
          parts%p%r(i) = 0.d0
          endif
      enddo      
      
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
      parts%tab%x%ndim1 = ntotal; parts%tab%y%ndim1 = ntotal; parts%tab%xy%ndim1 = ntotal

      parts%str%x = liquid%viscosity*parts%tab%x
      parts%str%y = liquid%viscosity*parts%tab%y
      parts%str%xy = liquid%viscosity*parts%tab%xy

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
      type(p2r) vx_i(3), vx_j(3)
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
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
          !dvx(d) = parts%vx(d,i) - parts%vx(d,j)
          dvx(d) = vx_i(d)%p - vx_j(d)%p
          dx     = parts%x(d,i)  - parts%x(d,j)
          vr     = vr + dvx(d)*dx
          rr     = rr + dx*dx
        enddo

!     Artificial viscous force only if v_ij * r_ij < 0

        if (vr.lt.0.e0) then

!     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )
            
          muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)
          
!     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

          mc   = 0.5e0*(parts%c%r(i) + parts%c%r(j))
          mrho = 0.5e0*(parts%rho%r(i) + parts%rho%r(j))
          piv  = (beta*muv - alpha*mc)*muv/mrho              

!     Calculate SPH sum for artificial viscous force

          !do d=1,dim
          !  h = -piv*parts%dwdx(d,k)
          !  parts%dvx(d,i) = parts%dvx(d,i) + parts%mass(j)*h
          !  parts%dvx(d,j) = parts%dvx(d,j) - parts%mass(i)*h
          !enddo

            h = -piv*parts%dwdx(1,k)
            parts%dvx%x%r(i) = parts%dvx%x%r(i) + parts%mass%r(j)*h
            parts%dvx%x%r(j) = parts%dvx%x%r(j) - parts%mass%r(i)*h
            if(dim.ge.2)then
            h = -piv*parts%dwdx(2,k)
            parts%dvx%y%r(i) = parts%dvx%y%r(i) + parts%mass%r(j)*h
            parts%dvx%y%r(j) = parts%dvx%y%r(j) - parts%mass%r(i)*h
            endif
            

        endif
      enddo

      return
      end subroutine

! Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
!----------------------------------------------------------------------      
      subroutine art_visc_omp(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts

      type(numerical), pointer :: numeric
      real(dp) dx, dvx(3), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml
      real(dp), allocatable, dimension(:,:,:) :: local
      type(p2r) vx_i(3), vx_j(3)
      integer i,j,k,d,dim,ntotal,niac,nthreads,it

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim
      nthreads =  parts%nthreads
            
      allocate(local(dim,ntotal,nthreads))
      call parts%get_niac_start_end
    
      numeric  => parts%numeric      
      alpha = numeric%alpha; beta = numeric%beta; etq = numeric%etq
         
!     Calculate SPH sum for artificial viscosity

!$omp parallel 
!$omp do private(i,j,k,mhsml,vr,rr,vx_i,vx_j,d,dvx,dx,muv,mc,mrho,piv,h) 
do it = 1, nthreads 
   do i = 1, ntotal
      local(1:dim,i,it) = 0.d0
   enddo
   do k=parts%niac_start(it),parts%niac_end(it)
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        mhsml= (parts%hsml(i)+parts%hsml(j))/2.
        vr = 0.e0
        rr = 0.e0
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
          !dvx(d) = parts%vx(d,i) - parts%vx(d,j)
          dvx(d) = vx_i(d)%p - vx_j(d)%p
          dx     = parts%x(d,i)  - parts%x(d,j)
          vr     = vr + dvx(d)*dx
          rr     = rr + dx*dx
        enddo

!     Artificial viscous force only if v_ij * r_ij < 0

        if (vr.lt.0.e0) then

!     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )
            
          muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)
          
!     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

          mc   = 0.5e0*(parts%c%r(i) + parts%c%r(j))
          mrho = 0.5e0*(parts%rho%r(i) + parts%rho%r(j))
          piv  = (beta*muv - alpha*mc)*muv/mrho              

!     Calculate SPH sum for artificial viscous force

          do d=1,dim
             h = -piv*parts%dwdx(d,k)
          !  parts%dvx(d,i) = parts%dvx(d,i) + parts%mass(j)*h
          !  parts%dvx(d,j) = parts%dvx(d,j) - parts%mass(i)*h
             local(d,i,it) = local(d,i,it) + parts%mass%r(j)*h
             local(d,j,it) = local(d,j,it) - parts%mass%r(i)*h
          enddo
        endif
  enddo !k
enddo !it
!$omp end do

!$omp do private(it)
do i = 1, ntotal
    do it = 1 ,nthreads
      parts%dvx%x%r(i) = parts%dvx%x%r(i) + local(1,i,it)
      if (dim.ge.2) parts%dvx%y%r(i) = parts%dvx%y%r(i) +local(2,i,it)
    enddo  
enddo
!$omp end do
!$omp end parallel
         
return
end subroutine


! Subroutine to calculate the artificial viscosity (Colagrossi, 2003)
      
!----------------------------------------------------------------------      
      subroutine art_visc2(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts

      type(numerical), pointer :: numeric
      real(dp) dx, dvx(3), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, hsml
      real(dp) dudxi,dudxj,dudyi,dvdxi,dudyj,dvdxj,dvdyi,dvdyj,EiEj,eplison,uij
      real(dp),allocatable,dimension(:) :: king
      type(p2r) vx_i(3), vx_j(3)
      type(array),allocatable :: dudx,dudy,dvdx,dvdy,divu
      integer i,j,k,d,dim,ntotal,niac

      allocate(king(parts%maxn))
      allocate(dudx);allocate(dudx%r(parts%maxn))
      allocate(dudy);allocate(dudy%r(parts%maxn))
      allocate(dvdx);allocate(dvdx%r(parts%maxn))
      allocate(dvdy);allocate(dvdy%r(parts%maxn))
      allocate(divu);allocate(divu%r(parts%maxn))
      dudx%ndim1 = parts%ntotal + parts%nvirt
      dudy%ndim1 = parts%ntotal + parts%nvirt      
      dvdx%ndim1 = parts%ntotal + parts%nvirt
      dvdy%ndim1 = parts%ntotal + parts%nvirt 
      divu%ndim1 = parts%ntotal + parts%nvirt 
      
      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim      
      numeric  => parts%numeric      
      alpha = 0.03                !0.005~0.03
      beta = numeric%beta; etq = numeric%etq
!      eplison = 0.1
      hsml= parts%hsml(1)
!     Calculate SPH sum for artificial viscosity
        dudx = parts%df(parts%vx%x,'x')
        dudy = parts%df(parts%vx%x,'y')
        dvdx = parts%df(parts%vx%y,'x')
        dvdy = parts%df(parts%vx%y,'y')
        
        divu=parts%div(parts%vx)

      do i = 1,parts%ntotal+parts%nvirt
          king(i) = 0 
      enddo
        
      do k=1,niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)

!        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        dudxi = parts%df(vx_i(1),x)
!        dudxj = parts%df(vx_j(1),x)
!        dudyi = parts%df(vx_i(1),y)
!        dvdxi = parts%df(vx_i(2),x)
!        dudyj = parts%df(vx_j(1),y)
!        dvdxj = parts%df(vx_j(2),x)
!        dvdyi = parts%df(vx_i(2),y)
!        dvdyj = parts%df(vx_j(2),y)xxczxc
        EiEj = sqrt(abs(dudx%r(i)*dudx%r(j) + (dudy%r(i) + dvdx%r(i))*(dudy%r(j) + dvdx%r(j))/2 +dvdy%r(i)*dvdy%r(j)))
        king(i)=king(i) + abs(divu%r(i))/(abs(divu%r(i)) + EiEj + parts%c%r(i)*1.0E-4/hsml)
        king(j)=king(j) + abs(divu%r(i))/(abs(divu%r(i)) + EiEj + parts%c%r(j)*1.0E-4/hsml)
        vr = 0.e0
        rr = 0.e0
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
          dvx(1) = parts%vx%x%r(i) - parts%vx%x%r(j)
          dvx(2) = parts%vx%y%r(i) - parts%vx%y%r(j)
          !dvx(d) = vx_i(d)%p - vx_j(d)%p
          dx     = parts%x(d,i)  - parts%x(d,j)
          vr     = vr + dvx(d)*dx
          rr     = rr + dx*dx
        enddo
        uij = hsml*(king(i) + king(j))*vr/(2*(rr+etq**2*hsml**2))
!     Artificial viscous force only if v_ij * r_ij < 0
        if (vr.lt.0.e0) then
          mc   = 0.5e0*(parts%c%r(i) + parts%c%r(j))
          mrho = 0.5e0*(parts%rho%r(i) + parts%rho%r(j))
          piv  = -alpha*uij*mc/mrho

         
!     Calculate SPH sum for artificial viscous force
         h = -piv*parts%dwdx(1,k)
         parts%dvx%x%r(i) = parts%dvx%x%r(i) + parts%mass%r(j)*h
         parts%dvx%x%r(j) = parts%dvx%x%r(j) - parts%mass%r(i)*h
         if(dim.ge.2)then
         h = -piv*parts%dwdx(2,k)
         parts%dvx%y%r(i) = parts%dvx%y%r(i) + parts%mass%r(j)*h
         parts%dvx%y%r(j) = parts%dvx%y%r(j) - parts%mass%r(i)*h
         endif
         
        endif
      enddo

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
         df(i) = df(i) + delta*parts%hsml(i)*parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
         df(j) = df(j) - delta*parts%hsml(j)*parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)
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
      type(array) f
      type(array) df
      real(dp), allocatable, dimension(:,:) :: local
      real(dp) dx(3),delta, muv, rr, h
      integer i,j,k,d,ntotal,niac,dim,it,nthreads

!      write(*,*) 'In art_density...'

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim; nthreads = parts%nthreads            
      delta    = parts%numeric%delta

      if(nthreads>=1)then
         allocate(local(ntotal,nthreads))
         call parts%get_niac_start_end
      endif   

     !$omp parallel 
     !$omp do private(i,j,d,k,rr,muv,dx,h)
     do it = 1, nthreads 
        do i =1,ntotal
           local(i,it) = 0.d0
        enddo     
        do k = parts%niac_start(it),parts%niac_end(it)
           i = parts%pair_i(k)
           j = parts%pair_j(k)
           rr = 0.e0
           do d=1,dim
              dx(d)  =  parts%x(d,i) -  parts%x(d,j)
              rr     = rr + dx(d)*dx(d)
           enddo
           muv = 2.0*(f%r(i)-f%r(j))/rr

         h = 0.d0
         do d=1,dim
            !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
            h = h + dx(d)*muv*parts%dwdx(d,k)
         enddo
!          df(i) = df(i) + delta*parts%hsml(i)*parts%c(i)*parts%mass(j)*h/parts%rho(j)
!          df(j) = df(j) - delta*parts%hsml(j)*parts%c(j)*parts%mass(i)*h/parts%rho(i)
         local(i,it) = local(i,it) + delta*parts%hsml(i)*parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
         local(j,it) = local(j,it) - delta*parts%hsml(j)*parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)
      enddo
     enddo
     !$omp end do
     !$omp barrier

     !$omp do private(it)
     do i = 1,ntotal
        do it = 1,nthreads
             df%r(i) =  df%r(i) + local(i,it)
        enddo
     enddo
     !$omp end do
     !$omp end parallel
return
      end subroutine
      
!-------------------------------------------------------------------
      subroutine delta_rho(parts,f,df)
!-------------------------------------------------------------------
      implicit none
      
      class(particles) parts
      type(array) f
      type(array) df
      real(dp), allocatable, dimension(:,:,:) :: w,rew
      real(dp), allocatable, dimension(:,:) :: d_rho
      real(dp), allocatable, dimension(:,:) :: local
      real(dp) dx(3),delta, muv, rr, h, m, n,deter
      integer i,j,k,d,ntotal,niac,dim,it,nthreads
    
      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac; dim = parts%dim            
      delta    = parts%numeric%delta
          
      allocate(w(ntotal,dim,dim))
      allocate(rew(ntotal,dim,dim))
      allocate(d_rho(dim,ntotal))

!      if(nthreads>=1)then
!         allocate(local(ntotal,nthreads))
!         call parts%get_niac_start_end
!      endif   
      
      do m =1,parts%dim
          do n =1,parts%dim
              do i =1,parts%ntotal+parts%nvirt
                  w(i,m,n) = 0
              enddo
          enddo
      enddo
      
      do k = 1,parts%niac
         i = parts%pair_i(k)
         j = parts%pair_j(k)
 !        if(i==5)then
 !            write(*,*)'k=',k,'i=',i,'j=',j
 !        endif
           
         !这里是不是dwdx如果对i用+，则对j就应该用-?
         w(i,1,1) = w(i,1,1) + (parts%x(1,j) - parts%x(1,i))*parts%dwdx(1,k)*parts%mass%r(j)/parts%rho%r(j)
         w(j,1,1) = w(j,1,1) - (parts%x(1,i) - parts%x(1,j))*parts%dwdx(1,k)*parts%mass%r(i)/parts%rho%r(i)
         w(i,1,2) = w(i,1,2) + (parts%x(1,j) - parts%x(1,i))*parts%dwdx(2,k)*parts%mass%r(j)/parts%rho%r(j)
         w(j,1,2) = w(j,1,2) - (parts%x(1,i) - parts%x(1,j))*parts%dwdx(2,k)*parts%mass%r(i)/parts%rho%r(i)
         w(i,2,1) = w(i,2,1) + (parts%x(2,j) - parts%x(2,i))*parts%dwdx(1,k)*parts%mass%r(j)/parts%rho%r(j)
         w(j,2,1) = w(j,2,1) - (parts%x(2,i) - parts%x(2,j))*parts%dwdx(1,k)*parts%mass%r(i)/parts%rho%r(i)
         w(i,2,2) = w(i,2,2) + (parts%x(2,j) - parts%x(2,i))*parts%dwdx(2,k)*parts%mass%r(j)/parts%rho%r(j)
         w(j,2,2) = w(j,2,2) - (parts%x(2,i) - parts%x(2,j))*parts%dwdx(2,k)*parts%mass%r(i)/parts%rho%r(i)
      enddo
      
      do i =1,parts%ntotal+parts%nvirt
          deter = w(i,1,1)*w(i,2,2)-w(i,1,2)*w(i,2,1)
          rew(i,1,1) =  w(i,2,2)/deter
          if(deter==0)then
          rew(i,1,1) =  0
          rew(i,1,2) =  0
          rew(i,2,1) =  0
          rew(i,2,2) =  0
          else
          rew(i,1,1) =  w(i,2,2)/deter
          rew(i,1,2) = -w(i,1,2)/deter
          rew(i,2,1) = -w(i,2,1)/deter
          rew(i,2,2) =  w(i,1,1)/deter
          endif
          do d =1,parts%dim
              d_rho(d,i) = 0
          enddo
      enddo
      
      
      do k = 1, parts%niac
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          d_rho(1,i) = d_rho(1,i) + (parts%rho%r(j) - parts%rho%r(i))*(rew(i,1,1)*parts%dwdx(1,k)+rew(i,1,2)*parts%dwdx(2,k))*parts%mass%r(j)/parts%rho%r(j)
          d_rho(2,i) = d_rho(2,i) + (parts%rho%r(j) - parts%rho%r(i))*(rew(i,2,1)*parts%dwdx(1,k)+rew(i,2,2)*parts%dwdx(2,k))*parts%mass%r(j)/parts%rho%r(j)
          d_rho(1,j) = d_rho(1,j) - (parts%rho%r(i) - parts%rho%r(j))*(rew(j,1,1)*parts%dwdx(1,k)+rew(j,1,2)*parts%dwdx(2,k))*parts%mass%r(i)/parts%rho%r(i)
          d_rho(2,j) = d_rho(2,j) - (parts%rho%r(i) - parts%rho%r(j))*(rew(j,2,1)*parts%dwdx(1,k)+rew(j,2,2)*parts%dwdx(2,k))*parts%mass%r(i)/parts%rho%r(i)
!          d_rho(1,j) = d_rho(1,j) + (parts%rho%r(i) - parts%rho%r(j))*(rew(j,1,1)*parts%dwdx(1,k)+rew(j,1,2)*parts%dwdx(2,k))*parts%mass%r(i)/parts%rho%r(i)
!          d_rho(2,j) = d_rho(2,j) + (parts%rho%r(i) - parts%rho%r(j))*(rew(j,2,1)*parts%dwdx(1,k)+rew(j,2,2)*parts%dwdx(2,k))*parts%mass%r(i)/parts%rho%r(i)
       enddo
      
        do k = 1,parts%niac
           i = parts%pair_i(k)
           j = parts%pair_j(k)
           rr = 0.e0
           do d=1,dim
              dx(d)  =  parts%x(d,i) -  parts%x(d,j)
              rr     = rr + dx(d)*dx(d)
           enddo
           muv = 2.0*(f%r(i)-f%r(j))/rr

         h = 0.d0
         do d=1,dim
            !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
!            h = h + dx(d)*muv*parts%dwdx(d,k)
           h = h + (dx(d)*muv-d_rho(d,i)-d_rho(d,j))*parts%dwdx(d,k)        
!           h = h + dx(d)*muv*parts%dwdx(d,k)      
         enddo
          df%r(i) = df%r(i) + delta*parts%hsml(i)*parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
          df%r(j) = df%r(j) - delta*parts%hsml(j)*parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)
!         local(i,it) = local(i,it) + delta*parts%hsml(i)*parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
!         local(j,it) = local(j,it) - delta*parts%hsml(j)*parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)
      enddo
!      WRITE(*,*) 'now is ',parts%itimestep
      
!      if(parts%itimestep==63789) then
!         write(*,*) 'sd'
!      endif

      
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
      type(p2r) vx_i(3), vx_j(3), av_i(3), av_j(3)
      integer i,j,k,d,ntotal,niac    

      ntotal = parts%ntotal
      niac   = parts%niac
      epsilon = parts%numeric%epsilon
      
      parts%av%x = 0.d0; parts%av%y = 0.d0
     
      do k=1,niac       
         i = parts%pair_i(k)
         j = parts%pair_j(k)       
         mrho = (parts%rho%r(i)+parts%rho%r(j))/2.0
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        av_i = parts%av%cmpt(i); av_j = parts%av%cmpt(j)
         do d=1,parts%dim
            !dvx(d) = parts%vx(d,i) - parts%vx(d,j)            
            dvx(d) = vx_i(d)%p - vx_j(d)%p            
            av_i(d)%p = av_i(d)%p - parts%mass%r(j)*dvx(d)/mrho*parts%w(k)
            av_j(d)%p = av_j(d)%p + parts%mass%r(i)*dvx(d)/mrho*parts%w(k)       
         enddo                    
      enddo  
        
!      do i = 1, ntotal
!         do d = 1, parts%dim
!            av_i(d)%p = epsilon * av_i(d)%p
!         enddo 
!      enddo

      parts%av%x = epsilon * parts%av%x; parts%av%y = epsilon * parts%av%y

      return
      end subroutine


!-------------------------------------------------------------------
      subroutine delta_sph_vx(parts)
!-------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: water
!      real(dp), dimension(:) :: f,df

      real(dp) dx(3),delta, muv, rr, h, dvx(3), delta_alpha
      type(p2r) vx_i(3), vx_j(3)
      integer i,j,k,d,ntotal,niac,dim

!      write(*,*) 'In art_density...'

      ntotal      =  parts%ntotal + parts%nvirt
      niac        =  parts%niac; dim = parts%dim            
      delta       = parts%numeric%delta
      delta_alpha = parts%numeric%delta_alpha
      water => parts%material
      
      do k=1,niac
         i = parts%pair_i(k)
         j = parts%pair_j(k)
         rr = 0.e0
         do d=1,dim
            dx(d)  =  parts%x(d,i) -  parts%x(d,j)
            rr     = rr + dx(d)*dx(d)
         enddo
         
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
           dvx(d) = vx_i(d)%p - vx_j(d)%p   
        enddo
!         muv = 2.0*(f(i)-f(j))/rr

         h = 0.d0
         if(rr/=0)then
         do d=1,dim
            !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
            h = h + dx(d)*dvx(d)/rr*parts%dwdx(1,k)
         enddo
         endif
         parts%dvx%x%r(i) = parts%dvx%x%r(i) + delta_alpha*parts%hsml(i)*parts%c%r(i)*water%rho0*parts%mass%r(j)*h/parts%rho%r(j)/parts%rho%r(i)
         parts%dvx%x%r(j) = parts%dvx%x%r(j) - delta_alpha*parts%hsml(j)*parts%c%r(j)*water%rho0*parts%mass%r(i)*h/parts%rho%r(i)/parts%rho%r(j)

         h = 0.d0
         if(rr/=0)then
         do d=1,dim
            !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
            h = h + dx(d)*dvx(d)/rr*parts%dwdx(2,k)
         enddo
         endif
         parts%dvx%y%r(i) = parts%dvx%y%r(i) + delta_alpha*parts%hsml(i)*parts%c%r(i)*water%rho0*parts%mass%r(j)*h/parts%rho%r(j)/parts%rho%r(i)
         parts%dvx%y%r(j) = parts%dvx%y%r(j) - delta_alpha*parts%hsml(j)*parts%c%r(j)*water%rho0*parts%mass%r(i)*h/parts%rho%r(i)/parts%rho%r(j)

      enddo
      

      return
      end subroutine

      
!----------------------------------------------------------------------      
      subroutine av_vel_omp(parts)
!----------------------------------------------------------------------   
      implicit none

      class(particles),target :: parts
      real(dp) vcc, dvx(3), epsilon, mrho
      type(p2r) vx_i(3), vx_j(3), av_i(3), av_j(3)
      real(dp),allocatable,dimension(:,:,:) :: local
      integer i,j,k,d,ntotal,niac,it,nthreads  

      ntotal = parts%ntotal
      niac   = parts%niac
      epsilon = parts%numeric%epsilon
      nthreads = parts%nthreads
      allocate(local(parts%dim,parts%ntotal+parts%nvirt,nthreads))
      call parts%get_niac_start_end

      parts%av%x = 0.d0; parts%av%y = 0.d0

!$omp parallel
!$omp do private(k,i,j,mrho,vx_i,vx_j,av_i,av_j,d,dvx) 
do it = 1,nthreads    
   do i = 1, ntotal
      do d= 1, parts%dim
         local(d,i,it)=0.d0
      enddo
   enddo
   do k=parts%niac_start(it),parts%niac_end(it)
         i = parts%pair_i(k)
         j = parts%pair_j(k)       
         mrho = (parts%rho%r(i)+parts%rho%r(j))/2.0
         vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!        av_i = parts%av%cmpt(i); av_j = parts%av%cmpt(j)
         do d=1,parts%dim
            !dvx(d) = parts%vx(d,i) - parts%vx(d,j)            
            dvx(d) = vx_i(d)%p - vx_j(d)%p        
            local(d,i,it) = local(d,i,it) - parts%mass%r(j)*dvx(d)/mrho*parts%w(k)
            local(d,j,it) = local(d,j,it) + parts%mass%r(i)*dvx(d)/mrho*parts%w(k)
!            av_i(d)%p = av_i(d)%p - parts%mass%r(j)*dvx(d)/mrho*parts%w(k)
!            av_j(d)%p = av_j(d)%p + parts%mass%r(i)*dvx(d)/mrho*parts%w(k)       
         enddo                    
   enddo  
enddo
!$omp end do
!$omp do private(it)
do i = 1, ntotal
    do it = 1,nthreads
       parts%av%x%r(i) = parts%av%x%r(i) + local(1,i,it)
       parts%av%y%r(i) = parts%av%y%r(i) + local(2,i,it)
    enddo   
enddo
!$omp end do
!$omp end parallel
!      do i = 1, ntotal
!         do d = 1, parts%dim
!            av_i(d)%p = epsilon * av_i(d)%p
!         enddo 
!      enddo

      parts%av%x = epsilon * parts%av%x; parts%av%y = epsilon * parts%av%y

      return
      end subroutine

      
!----------------------------------------------------------------------
      subroutine shear_strain_rate(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz
      type(p2r) vx_i(3), vx_j(3)
      integer i, j, k, d, dim, ntotal, niac

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac; dim = parts%dim
      
      parts%tab%x%r(1:ntotal) = 0.e0
      if(dim>=2)then
         parts%tab%y%r(1:ntotal) = 0.e0
         parts%tab%xy%r(1:ntotal) = 0.e0
      endif
      if(dim==3)then
         parts%tab%z%r(1:ntotal) = 0.e0
         parts%tab%xz%r(1:ntotal) = 0.e0
         parts%tab%yz%r(1:ntotal) = 0.e0
      endif
      
!     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

        do k=1,niac
          i = parts%pair_i(k)
          j = parts%pair_j(k)
          vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
          do d=1,dim
             !dvx(d) = parts%vx(d,j) - parts%vx(d,i)
             dvx(d) = vx_j(d)%p - vx_i(d)%p
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
             parts%tab%x%r(i) = parts%tab%x%r(i) + parts%mass%r(j)*hxx/parts%rho%r(j)
             parts%tab%x%r(j) = parts%tab%x%r(j) + parts%mass%r(i)*hxx/parts%rho%r(i)             
          else if (dim.eq.2) then           
             parts%tab%x%r(i) = parts%tab%x%r(i) + parts%mass%r(j)*hxx/parts%rho%r(j)
             parts%tab%x%r(j) = parts%tab%x%r(j) + parts%mass%r(i)*hxx/parts%rho%r(i)   
             parts%tab%xy%r(i) = parts%tab%xy%r(i) + parts%mass%r(j)*hxy/parts%rho%r(j)
             parts%tab%xy%r(j) = parts%tab%xy%r(j) + parts%mass%r(i)*hxy/parts%rho%r(i)            
             parts%tab%y%r(i) = parts%tab%y%r(i) + parts%mass%r(j)*hyy/parts%rho%r(j)
             parts%tab%y%r(j) = parts%tab%y%r(j) + parts%mass%r(i)*hyy/parts%rho%r(i)          
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

!---------------------------------------------------------------------
      subroutine Jaumann_rate(parts)
!----------------------------------------------------------------------
      implicit none

      class(particles) parts
      type(material), pointer :: soil
      real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz, G
      type(p2r) vx_i(3), vx_j(3) 
      integer i, j, k, d, dim, ntotal, niac

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac; dim = parts%dim
      soil => parts%material
      G = soil%e/(2.0*(1+soil%niu))

!     Hook's law

!$omp parallel do
      do i = 1, ntotal
         parts%dstr%x%r(i) = parts%dstr%x%r(i)+G*parts%tab%x%r(i)   ! No accumulation origionaly
         parts%dstr%xy%r(i) = parts%dstr%xy%r(i)+G*parts%tab%xy%r(i)
         parts%dstr%y%r(i) = parts%dstr%y%r(i)+G*parts%tab%y%r(i)
      enddo

!         if(parts%soil_pressure==2)then
!$omp parallel do
      do i = 1, ntotal
         parts%dp%r(i)   = parts%dp%r(i) - soil%k*parts%vcc%r(i)  !!! simultaneous pressure  
      enddo
!         endif         
     
!     spin tensor
      parts%wxy = 0.5*(parts%df4_omp(parts%vx%x,'y')-parts%df4_omp(parts%vx%y,'x'))

!   Jaumann rate

!$omp parallel do
      do i = 1, ntotal
         parts%dstr%x%r(i) = parts%dstr%x%r(i)+2.0*parts%str%xy%r(i)*parts%wxy%r(i)
         parts%dstr%xy%r(i) = parts%dstr%xy%r(i)-(parts%str%x%r(i)-parts%str%y%r(i))*parts%wxy%r(i)
         parts%dstr%y%r(i) = parts%dstr%y%r(i)-2.0*parts%str%xy%r(i)*parts%wxy%r(i)
      enddo         

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
         tmax  = sqrt(((soil%str%x%r(i)-soil%str%y%r(i))/2)**2+soil%str%xy%r(i)**2)
         if(tmax<1.e-6)cycle
         yield = cohesion*cos(phi)+soil%p%r(i)*sin(phi)

         if(yield<=0.)then    ! <
            yield=0.; soil%p%r(i)=-cohesion*tan(phi)**(-1.0)
            soil%rho%r(i) = property%rho0
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
            soil%str%x%r(i) = skale * soil%str%x%r(i)
            soil%str%xy%r(i) = skale * soil%str%xy%r(i)
            soil%str%y%r(i) = skale * soil%str%y%r(i)
         endif
         !sxx(i) = 0.; sxy(i) = 0.; syy(i) = 0.
      enddo
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
!$omp parallel do private(I1,J2,skale) reduction(+:k)
      do i = 1, ntotal
         !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         !if(tmax<1.e-6)cycle
         !yield = cohesion*cos(phi)+p(i)*sin(phi)


         if(soil%p%r(i)<0.)then    ! <
            !k = k + 1
            !soil%fail(i) = 1
            soil%p%r(i)=0.d0
            !soil%rho(i) = property%rho0
         endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         I1 = 3.*soil%p%r(i) 

         J2 = soil%str%x%r(i)**2.+2.*soil%str%xy%r(i)**2.+soil%str%y%r(i)**2.+(soil%str%x%r(i)+soil%str%y%r(i))**2.
         J2 = sqrt(J2/2.)+1.d-6
         !if(J2<1.e-6)cycle

         if(J2>alpha1*I1)then
!         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = alpha1*I1/J2
            soil%str%x%r(i) = skale * soil%str%x%r(i)
            soil%str%xy%r(i) = skale * soil%str%xy%r(i)
            soil%str%y%r(i) = skale * soil%str%y%r(i)
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

         I1 = 3.*soil%p%r(i) 

         J2 = soil%str%x%r(i)**2.+2.*soil%str%xy%r(i)**2.+soil%str%y%r(i)**2.+(soil%str%x%r(i)+soil%str%y%r(i))**2.
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

      dsxx => parts%dstr%x%r
      dsxy => parts%dstr%xy%r
      dsyy => parts%dstr%y%r
      vcc  => parts%vcc%r
      sxx => parts%str%x%r
      sxy => parts%str%xy%r
      syy => parts%str%y%r

      property => parts%material
      phi = property%phi
      k   = property%k; e = property%E; niu = property%niu
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      G = e/(2.0*(1+niu))
!$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
      do i = 1, ntotal
                            !if(parts%fail(i)==1)then

      exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%tab%xy%r(i)/2.
      eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

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
      parts%dp%r(i) = parts%dp%r(i) + 3.*k*alpha*dlambda

! Accumulative deviatoric strain
 
      parts%epsilon_p%r(i) = parts%epsilon_p%r(i)    &
                         + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            !endif ! Fail
      enddo
!$omp end parallel do
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

      dsxx2 => parts%dstr2%x%r
      dsxy2 => parts%dstr2%xy%r
      dsyy2 => parts%dstr2%y%r
      vcc  => parts%vcc%r
      sxx => parts%str%x%r
      sxy => parts%str%xy%r
      syy => parts%str%y%r

      property => parts%material
      phi = property%phi
      k   = property%k; e = property%E; niu = property%niu
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      G = e/(2.0*(1+niu))
!$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
      do i = 1, ntotal
                            !if(parts%fail(i)==1)then

      exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%tab%xy%r(i)/2.
      eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

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
      parts%dp2%r(i) =  3.*k*alpha*dlambda

! Accumulative deviatoric strain
 
!      parts%epsilon_p(i) = parts%epsilon_p(i) 
!     &                   + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            !endif ! Fail
      enddo
!$omp end parallel do
      
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

      dsxx => parts%dstr%x%r
      dsxy => parts%dstr%xy%r
      dsyy => parts%dstr%y%r
      vcc  => parts%vcc%r
      sxx => parts%str%x%r
      sxy => parts%str%xy%r
      syy => parts%str%y%r

      property => parts%material
      phi = property%phi
      k   = property%k; e = property%E; niu = property%niu
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      G = e/(2.0*(1+niu))
!$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
      do i = 1, ntotal
                            if(parts%fail(i)==1)then

      exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%tab%xy%r(i)/2.
      eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

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
      parts%dp%r(i) = parts%dp%r(i) + 3.*k*alpha*dlambda

! Accumulative deviatoric strain
 
      parts%epsilon_p%r(i) = parts%epsilon_p%r(i)       & 
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
      type(p2r) vx_i(3), vx_j(3), dvx_i(3), dvx_j(3)
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
          rrw = water%w(k)/(water%rho%r(i)*soil%rho%r(j))

! For staturated soil
!        if(volume_fraction) cf = water%vof(i)*water%rho(i)*(-gravity)/ks
!        if(water%volume_fraction) cf = water%vof(i)*soil%vof(j)*water%rho(i)*(-gravity)/ks
        cf = water%vof%r(i)*soil%vof%r(j)*water%rho%r(i)*(-gravity)/ks

          vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
          dvx_i = water%dvx%cmpt(i); dvx_j = soil%dvx%cmpt(j)
          do d=1,dim
             !sp = cf*(water%vx(d,i)-soil%vx(d,j))*rrw
             !water%dvx(d,i) = water%dvx(d,i) - soil%mass(j)*sp
             !soil%dvx(d,j)  = soil%dvx(d,j) + water%mass(i)*sp   
             sp = cf*(vx_i(d)%p-vx_j(d)%p)*rrw
             dvx_i(d)%p = dvx_i(d)%p - soil%mass%r(j)*sp
             dvx_j(d)%p = dvx_j(d)%p + water%mass%r(i)*sp   
          enddo

          !sp = cf*(water%vx%x%r(i)-soil%vx%x%r(j))*rrw
          !water%dvx%x%r(i) = water%dvx%x%r(i) - soil%mass(j)*sp
          !soil%dvx%x%r(j)  =  soil%dvx%x%r(j) + water%mass(i)*sp   
          !sp = cf*(water%vx%y%r(i)-soil%vx%y%r(j))*rrw
          !water%dvx%y%r(i) = water%dvx%y%r(i) - soil%mass(j)*sp
          !soil%dvx%y%r(j)  =  soil%dvx%y%r(j) + water%mass(i)*sp   

      enddo

      return
      end subroutine 

!---------------------------------------------------------------
      subroutine darcy_law_omp(water, soil)
!---------------------------------------------------------------
!      use param, only: volume_fraction
!      use m_particles
      implicit none

      type(particles) water, soil
      real(dp), allocatable, dimension(:,:,:) :: wlocal, slocal
      real(dp) dx(3), ks, ns, gw, cf, sp, rrw, gravity
      type(p2r) vx_i(3), vx_j(3), dvx_i(3), dvx_j(3)
      type(material), pointer :: h2o,sio2  
      type(numerical), pointer :: numeric
      integer i, j, k, d, dim, nthreads,it, wntotal, sntotal

      nthreads = water%nthreads
      wntotal = water%ntotal + water%nvirt
      sntotal = soil%ntotal + soil%nvirt

      h2o => water%material
      sio2=>  soil%material
      numeric => water%numeric
      gravity = numeric%gravity
      dim = water%dim

      gw = h2o%rho0*(-gravity); ns = sio2%porosity
      ks = sio2%permeability
      cf = gw*ns/ks
      !cf = 6.e6

      allocate(wlocal(dim,wntotal,nthreads))
      allocate(slocal(dim,sntotal,nthreads))
      call water%get_niac_start_end
    
!$omp parallel
!$omp do private(i,j,rrw,cf,vx_i,vx_j,d,sp,k) 
do it = 1,nthreads
       do i = 1, wntotal
          do d = 1,dim
             wlocal(d,i,it) = 0.d0
          enddo  
       enddo
       do j = 1, sntotal
          do d = 1,dim
             slocal(d,j,it) = 0.d0
          enddo  
       enddo          
       do k=water%niac_start(it),water%niac_end(it)
          i = water%pair_i(k)
          j = water%pair_j(k)  
          rrw = water%w(k)/(water%rho%r(i)*soil%rho%r(j))

! For staturated soil
!        if(volume_fraction) cf = water%vof(i)*water%rho(i)*(-gravity)/ks
!        if(water%volume_fraction) cf = water%vof(i)*soil%vof(j)*water%rho(i)*(-gravity)/ks
        cf = water%vof%r(i)*soil%vof%r(j)*water%rho%r(i)*(-gravity)/ks

          vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
!          dvx_i = water%dvx%cmpt(i); dvx_j = soil%dvx%cmpt(j)
          do d=1,dim
             !sp = cf*(water%vx(d,i)-soil%vx(d,j))*rrw
             !water%dvx(d,i) = water%dvx(d,i) - soil%mass(j)*sp
             !soil%dvx(d,j)  = soil%dvx(d,j) + water%mass(i)*sp   
             sp = cf*(vx_i(d)%p-vx_j(d)%p)*rrw
             !dvx_i(d)%p = dvx_i(d)%p - soil%mass%r(j)*sp
             !dvx_j(d)%p = dvx_j(d)%p + water%mass%r(i)*sp 
             wlocal(d,i,it) = wlocal(d,i,it) - soil%mass%r(j)*sp
             slocal(d,j,it) = slocal(d,j,it) + water%mass%r(i)*sp
          enddo
      enddo
enddo     
!$omp end do
!$omp barrier

!$omp do private(it)
do i = 1,wntotal
!   dvx_i = water%dvx%cmpt(i)
!      do d=1,dim
!      dvx_i(d)%p = local_water(d,i)+dvx_i(d)%p
!      enddo
   do it =1,nthreads
      water%dvx%x%r(i) = water%dvx%x%r(i) + wlocal(1,i,it)
      water%dvx%y%r(i) = water%dvx%y%r(i) + wlocal(2,i,it)
   enddo
enddo
!$omp end do

!$omp do private(it)
do j = 1,sntotal
!   dvx_j = soil%dvx%cmpt(j)
!      do d=1,dim
!      dvx_j(d)%p = local_soil(d,j)+dvx_j(d)%p
!      enddo
   do it = 1, nthreads
      soil%dvx%x%r(j) = soil%dvx%x%r(j) + slocal(1,j,it)
      soil%dvx%y%r(j) = soil%dvx%y%r(j) + slocal(2,j,it)      
   enddo 
enddo
!$omp end do
!$omp end parallel
      
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
         mprr = water%mass%r(i)*water%p%r(i)/(water%rho%r(i)*soil%rho%r(j))
!         mprr = water%mass(i)*(water%p(i)+soil%p(j))/       &      Bui2014
!                (water%rho(i)*soil%rho(j))
!         do d = 1, water%dim
!            soil%dvx(d,j) = soil%dvx(d,j) + mprr*water%dwdx(d,k)  &  !+
!                            *soil%vof%r(j)  
!         enddo

            soil%dvx%x%r(j) = soil%dvx%x%r(j) + mprr*water%dwdx(1,k)  &  !+
                            *soil%vof%r(j)  
            soil%dvx%y%r(j) = soil%dvx%y%r(j) + mprr*water%dwdx(2,k)  &  !+
                            *soil%vof%r(j)  

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

!-----------------------------------------------------
      subroutine pore_water_pressure_omp(water,soil)
!-----------------------------------------------------
!      use param
!      use m_particles
      implicit none
      real(dp), allocatable, dimension(:,:,:) :: local
      type(particles) water, soil
      real(dp) mprr
      integer i, j, k, d, ntotal, it, nthreads, dim

      nthreads = water%nthreads
      ntotal = soil%ntotal+soil%nvirt; dim = water%dim

      allocate(local(dim,ntotal,nthreads))
      call water%get_niac_start_end
         
!$omp parallel
!$omp do private(i,j,mprr,k) 
do it =1,nthreads
   do j = 1,ntotal
      local(1:dim,j,it) = 0.d0
   enddo 
   do k = water%niac_start(it),water%niac_end(it)
         i = water%pair_i(k)   ! water
         j = water%pair_j(k)   ! soil
         mprr = water%mass%r(i)*water%p%r(i)/(water%rho%r(i)*soil%rho%r(j))
!         mprr = water%mass(i)*(water%p(i)+soil%p(j))/       &      Bui2014
!                (water%rho(i)*soil%rho(j))
         do d = 1, dim
!            soil%dvx(d,j) = soil%dvx(d,j) + mprr*water%dwdx(d,k)  &  !+
!                            *soil%vof%r(j)  
             local(d,j,it) = local(d,j,it) + mprr*water%dwdx(d,k)  &  !+
                             *soil%vof%r(j)  
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
enddo
!$omp end do
!$omp barrier

!$omp do private(it)
do j = 1,soil%ntotal+soil%nvirt
   do it = 1,water%nthreads
      soil%dvx%x%r(j) = soil%dvx%x%r(j) + local(1,j,it)
      soil%dvx%y%r(j) = soil%dvx%y%r(j) + local(2,j,it)
   enddo  
enddo   
!$omp end do
!$omp end parallel
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
         water%vof2%r(i) = water%vof2%r(i)+soil%mass%r(j)*water%w(k)
      enddo

      do k = 1, water%ntotal+water%nvirt
         water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0
      enddo

      return
      end subroutine 

!-------------------------------------------------------------
      subroutine volume_fraction_water_omp(water, soil)
!-------------------------------------------------------------
!      use param
!      use m_particles
      implicit none

      type(particles) water, soil
      integer i,j,k,d, ntotal,nthreads,it
      real(dp), allocatable, dimension(:,:) :: local
      type(material), pointer :: sio2

      sio2 => soil%material
      ntotal = water%ntotal+water%nvirt
      nthreads =water%nthreads

      allocate(local(ntotal,nthreads))
      call water%get_niac_start_end

      water%vof2 = 0.d0

!$omp parallel
!$omp do private(i,j,k)
do it = 1, nthreads 
     do i = 1,ntotal
        local(i,it) = 0.d0
     enddo  
     do k = water%niac_start(it),water%niac_end(it)
         i = water%pair_i(k)
         j = water%pair_j(k)
!         water%vof2%r(i) = water%vof2%r(i)+soil%mass%r(j)*water%w(k)
         local(i,it) = local(i,it)+soil%mass%r(j)*water%w(k)
     enddo
enddo
!$omp end do

!$omp do private(it)
      do k = 1, water%ntotal+water%nvirt
!         water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0
         do it= 1, nthreads 
            water%vof2%r(k) = water%vof2%r(k) + local(k,it)
         enddo 
         water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0
      enddo
!$omp end do
!$omp end parallel
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
      type(p2r) vx_i(3), vx_j(3)

      sio2 => soil%material
      ntotal = water%ntotal+water%nvirt

!      water%dvof = 0.d0

      do k = 1, water%niac
         i = water%pair_i(k)
         j = water%pair_j(k)
         vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
         do d = 1, water%dim
            !dvx(d) = water%vx(d,i)-soil%vx(d,j)
            dvx(d) = vx_i(d)%p - vx_j(d)%p
         enddo 
         tmp = dvx(1)*water%dwdx(1,k)+dvx(2)*water%dwdx(2,k)
         water%dvof%r(i) = water%dvof%r(i)-soil%mass%r(j)*tmp/sio2%rho0
      enddo

      return
      end subroutine 

!-------------------------------------------------------------
      subroutine volume_fraction_water2_omp(water, soil)
!-------------------------------------------------------------
!      use param
!      use m_particles
      implicit none

      type(particles) water, soil
      integer i,j,k,d, ntotal,nthreads,it
      type(material), pointer :: sio2
      real(dp), allocatable, dimension(:,:) :: local
      double precision dvx(3),tmp
      type(p2r) vx_i(3), vx_j(3)

      sio2 => soil%material
      ntotal = water%ntotal+water%nvirt
      nthreads = water%nthreads

      allocate(local(ntotal,nthreads))
      call water%get_niac_start_end

!      water%dvof = 0.d0

!$omp parallel
!$omp do private(i,j,k,vx_i,vx_j,d,dvx,tmp)
do it = 1, nthreads
      do i = 1,ntotal
         local(i,it) = 0.d0
      enddo  
      do k = water%niac_start(it),water%niac_end(it)
         i = water%pair_i(k)
         j = water%pair_j(k)
         vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
         do d = 1, water%dim
            !dvx(d) = water%vx(d,i)-soil%vx(d,j)
            dvx(d) = vx_i(d)%p - vx_j(d)%p
         enddo 
         tmp = dvx(1)*water%dwdx(1,k)+dvx(2)*water%dwdx(2,k)
!         water%dvof%r(i) = water%dvof%r(i)-soil%mass%r(j)*tmp/sio2%rho0
         local(i,it) = local(i,it)-soil%mass%r(j)*tmp/sio2%rho0
      enddo
enddo
!$omp end do

!$omp do private(it)
do i = 1, ntotal
    do it = 1,nthreads
       water%dvof%r(i) = water%dvof%r(i) + local(i,it)
    enddo   
enddo
!$omp end do
!$omp end parallel

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
         parts%vof%r(i) = parts%rho%r(i)/sio2%rho0
      enddo

      return
      end subroutine
      
!--------------------------------------------------------------------
      subroutine volume_fraction_soil_omp(parts)
!--------------------------------------------------------------------
!      use param
!      use m_particles
      implicit none

      type(particles) parts
      integer i,j,k, ntotal
      type(material), pointer :: sio2

      sio2 => parts%material
      ntotal = parts%ntotal + parts%nvirt
      
!$omp parallel do
      do i = 1, ntotal
         parts%vof%r(i) = parts%rho%r(i)/sio2%rho0
      enddo
!$omp end parallel do 

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
      !water%dvx(d,i) = water%dvx(d,i) - water%vx(d,i)*cf/water%rho(i)
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

!--------------------------------------------------
subroutine LU_break(A,n,L,U)
!--------------------------------------------------
implicit none
integer k,i,j,t,n
real(dp) :: A(n,n),L(n,n),U(n,n)
real(dp) :: sum1,sum2
do k = 1,n
    sum1=0 ; sum2=0
    do j = k,n
        sum1=0
        do t=1,k-1
            sum1=sum1+L(k,t)*U(t,j)
        enddo
        U(k,j)=A(k,j)-sum1
    enddo
    do j = 1,k-1
        U(k,j)=0    
    enddo
    do i =k+1,n
        do t=1,k-1
            sum2=sum2+L(i,t)*U(t,k)
        enddo
!        if(U(k,k)==0) then
!            write(*,*) 'U(k,k)=0 at LU_break'
!            stop
!        endif
        L(i,k)=(A(i,k)-sum2)/U(k,k)    
        
    enddo
    do i = k+1,n
        L(k,i)=0
    enddo
enddo

do j = 1,n
    L(j,j)=1
enddo
return
end subroutine
    
!求矩阵转置
!--------------------------------------------------
subroutine transpose_matrix(A,n,m,B)
!--------------------------------------------------
implicit none
integer i,j,n,m
real(dp) :: A(n,m),B(m,n),t
do i = 1,m
    do j =1,n
        B(i,j)=A(j,i)
    enddo
enddo
return
end subroutine


!--------------------------------------------------
subroutine multiply_matrix(A,n,m,B,l,C)
!--------------------------------------------------
implicit none
integer i,j,k,l,n,m
real(dp) sum
real(dp) A(n,m),B(m,l),C(n,l)
do i =1,n
    do j =1,l
        sum=0
        do k =1,m
            sum=sum+A(i,k)*B(k,j)
        enddo
        C(i,j)=sum
    enddo
enddo
return
end subroutine
        

!--------------------------------------------------
subroutine inverse_uptri_matrix(A,n,inverse_A)
!--------------------------------------------------
implicit none
integer i,j,k,n
real(dp) sum
real(dp) A(n,n),inverse_A(n,n)

do i=1,n
    do j=1,n
    inverse_A(i,j)=0
    enddo
enddo

do i =1,n
    inverse_A(i,i)=1/A(i,i)
!    if(A(i,i)==0)then
!    write(*,*),'A(i,i)=0 at MLS'
!    stop
!    endif
enddo
do i =1,n
    do j =1,n
        if (i<j)then
            sum = 0 
            do k=i,j-1
                sum=sum+inverse_A(i,k)*A(k,j)
            enddo
            inverse_A(i,j)=-sum/A(j,j)
        endif
    enddo
enddo
return
end subroutine

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
