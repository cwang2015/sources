!---------------------------
    module m_particles
!---------------------------
use constants, ex_pi=>pi
implicit none

! Configuration of the water-jet problem; Problem specific
type geo

! size of each dimension
  double precision dl,dr,dh,ds,dj,D,L,H,immerse
! descrete intervals
  integer  mdl,mdr,mdh,mds,mdj,mD,mL,mH
  contains
    procedure :: set => set_geo_sub

end type

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

type zone
     character(len=32) name
     real(dp) phi_s, phi_f
end type

type particle
   integer itype
   real(dp) x, y, vol, hsml, mass
   contains
      procedure :: copy => copy_particle_sub
end type

type enumerator
   integer :: inside = 0, inlet = 1, boundary = 2
end type

type tensor
   real(dp) xx,xy,xz,yy,yz,zz
end type

type material

! For water

   real(dp) rho0,b,gamma,c,viscosity

! For soil

   real(dp) k,porosity,permeability,G,E,niu,cohesion,phi
   
end type

type numerical

! Gravitational acceleration
  real(dp) :: gravity = -9.8

! Artificial viscosity
! alpha: shear viscosity; beta: bulk viscosity; etq: parameter to avoid sigularities   
   real(dp) :: alpha=0.d0, beta=0.d0, etq=0.1d0

! Leonard_Johns repulsive force
   real(dp) :: rr0 = 0.d0, dd = 0.d0, p1 = 12, p2 = 4

! Delta-SPH
   real(dp) :: delta = 0.1d0
   
end type

type array
   integer :: ndim1
!   integer :: dim = 2, rank = 0 
   real(dp), pointer, dimension(:)   :: r  => null()
!   real(dp), pointer, dimension(:,:) :: rr => null()
!   real(dp), pointer, dimension(:) :: x,y,z,xy,xz,yz
!   type(array), pointer :: x => null(), y => null(), z => null()
!   type(array), pointer :: xy => null(), xz => null(), yz => null()
   contains
       procedure :: dealloc => array_dealloc
       procedure :: array_sum
       procedure :: array_sub
       procedure :: array_minus  ! reverse sign
       procedure :: array_mul 
       procedure :: array_equal
       generic :: operator(+) => array_sum
       generic :: operator(-) => array_sub, array_minus
       generic :: operator(*) => array_mul
       generic :: assignment(=) => array_equal
       final :: array_final
end type

! Particles in SPH method
type particles

   integer :: dim   = 2
   integer :: maxn  = 0
   integer :: niac  = 0
   integer :: max_interaction = 0
   integer :: ntotal = 0, nvirt = 0
   integer :: max_zone = 10

   character(len=32) :: imaterial
   class(*), pointer :: material => null()

! Background point
  type(block), pointer :: backpoint => null()
  integer nreal_zone, nvirtual_zone
  integer, pointer, dimension(:) :: real_zone => null(), virtual_zone => null()

! Numerical parameters
   class(*), pointer :: numeric => null()
   
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
       procedure :: take_virtual_points
       procedure :: setup_itype
!       procedure :: grad_scalar
!       procedure :: grad_tensor
!       generic :: grad => grad_scalar, grad_tensor
 
end type 

interface write_field
   module procedure :: write_scalar_field
   module procedure :: write_vector2d_field
end interface

!interface p2m
!   module procedure :: point_to_water
!end interface

interface operator(+)
   module procedure :: any_add 
   module procedure :: add_any 
end interface

!=======
contains
!=======

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
   subroutine copy_particle_sub(this,that)
!-------------------------------------------
implicit none
class(particle) this
class(particle)  that

this%x = that%x

return
end subroutine

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

!-----------------------------------------------
      subroutine take_virtual_points(this,tank)
!-----------------------------------------------
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
    recursive subroutine array_equal(a,b)
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
       subroutine array_dealloc(this)
!--------------------------------------------
implicit none
class(array) this

if(associated(this%r))deallocate(this%r)

return
end subroutine

!--------------------------------------------
       subroutine array_final(this)
!--------------------------------------------
implicit none
type(array) this

call this%dealloc

return
end subroutine

!--------------------------------------------
           function array_sum(a,b)
!--------------------------------------------
implicit none
class(array), intent(in) :: a, b
class(array), allocatable:: array_sum
integer ndim1

if(a%ndim1/=b%ndim1)stop 'Cannot add arrays!'
ndim1 = a%ndim1
allocate(array_sum)
allocate(array_sum%r(ndim1))
array_sum%ndim1 = ndim1

array_sum%r(1:ndim1) = a%r(1:ndim1) + b%r(1:ndim1)

end function

!--------------------------------------------
           function array_sub(a,b)
!--------------------------------------------
implicit none
class(array), intent(in) :: a, b
class(array), allocatable:: array_sub
integer ndim1

if(a%ndim1/=b%ndim1)stop 'Cannot add arrays!'
ndim1 = a%ndim1
allocate(array_sub)
allocate(array_sub%r(ndim1))
array_sub%ndim1 = ndim1

array_sub%r(1:ndim1) = a%r(1:ndim1) - b%r(1:ndim1)

end function

!--------------------------------------------
           function array_minus(a)
!--------------------------------------------
implicit none
class(array), intent(in) :: a
class(array), allocatable:: array_minus
integer ndim1

ndim1 = a%ndim1
allocate(array_minus)
allocate(array_minus%r(ndim1))
array_minus%ndim1 = ndim1

array_minus%r(1:ndim1) = -a%r(1:ndim1)

end function

!--------------------------------------------
           function array_mul(a,c)
!--------------------------------------------
implicit none
class(array), intent(in) :: a
class(array), allocatable:: array_mul
class(*),intent(in) :: c
integer ndim1

ndim1 = a%ndim1
allocate(array_mul)
allocate(array_mul%r(ndim1))
array_mul%ndim1 = ndim1

select type(c)
   type is (real(8))
      array_mul%r(1:ndim1) = c*a%r(1:ndim1)
   class default
      stop 'Errorwang: binary operation has not been defined!'
end select

end function

!------------------------------------------
!      subroutine grad_scalar(this,scalar)
!------------------------------------------
!implicit none
!class(particles) this
!type(array) scalar



!return
!end subroutine

!------------------------------------------
!      subroutine grad_tensor(this,scalar)
!------------------------------------------
!implicit none
!class(particles) this
!integer scalar



!return
!end subroutine

end module
