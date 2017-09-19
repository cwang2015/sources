!----------------------------
    module m_dao
!----------------------------
implicit none

type dao
  contains
    procedure :: get_size=>get_total_particles_number     
    procedure :: get_maxn=>get_maximum_particles_number
    procedure :: get_dim=>get_dimension
endtype
private :: get_total_particles_number, get_maximum_particles_number, &
           get_dimension
contains
!-----------------------------------------------------------        
   function get_total_particles_number(this) result(ntotal)
!-----------------------------------------------------------           
   implicit none
   class(dao) this
   integer ntotal
   ntotal = 0
   end function
!-----------------------------------------------------------        
   function get_maximum_particles_number(this) result(maxn)
!-----------------------------------------------------------           
   implicit none
   class(dao) this
   integer maxn
   maxn = 0
   end function
!-----------------------------------------------------------        
   function get_dimension(this) result(dim)
!-----------------------------------------------------------           
   implicit none
   class(dao) this
   integer dim
   dim = 0
   end function   
end module

!---------------------------
    module m_array
!---------------------------
use m_dao
implicit none
integer, parameter :: dp = kind(0.d0)

type p2r
     real(dp), pointer :: p => null()
end type 

type array
   character(len=16) :: name = ''
   integer :: ndim1
   integer :: rank = 0
   real(dp), pointer, dimension(:) :: r  => null()
   type(array), pointer :: x => null(), y => null(), z => null()
   type(array), pointer :: xy => null(), xz => null(), yz => null()
   type(array), pointer, dimension(:) :: v => null()
   class(dao), pointer :: parts => null()
   contains
       procedure :: cmpt
       procedure :: get_dim
       procedure :: get_size
       procedure :: get_maxn
       procedure :: max_alloc
       procedure :: alloc
       procedure :: set
       !procedure :: p2cmpt
       final :: array_final
end type

!PGI fortran doesn't support ::

interface assignment(=)
  !module procedure :: array_equal_array
  module procedure  array_equal_array
  module procedure  array_equal_double_real
  module procedure  array_equal_real
end interface

interface operator(+)
   module procedure  array_add_array 
   module procedure  array_add_double_real 
   module procedure  double_real_add_array 
end interface

interface operator(-)
   module procedure  array_sub_array
   module procedure  array_sub_real
   module procedure  array_sub_double_real
   module procedure  array_minus
end interface

interface operator(*)
   module procedure  array_mul_array
   module procedure  array_mul_real
   module procedure  array_mul_double_real
   module procedure  double_real_mul_array
   module procedure  real_mul_array
end interface

interface operator(/)
   module procedure  array_div_array
   module procedure  array_div_real
   module procedure  array_div_double_real
end interface

interface operator(**)
   module procedure  array_pow_real
   module procedure  array_pow_double_real
end interface

!=======
contains
!=======

!--------------------------------------------
    function get_dim(a) result(dim)
!--------------------------------------------            
implicit none
class(array) a
integer dim
dim = a%parts%get_dim()
return
end function

!--------------------------------------------
    function get_size(a) result(size)
!--------------------------------------------            
implicit none
class(array) a
integer size
if(.not.associated(a%parts))then
    write(*,*) a%name,  ' :parts not associated!'
endif    
size = a%parts%get_size()
return
end function

!--------------------------------------------
    function get_maxn(a) result(maxn)
!--------------------------------------------            
implicit none
class(array) a
integer maxn
maxn = a%parts%get_maxn()
return
end function

!--------------------------------------------
    subroutine set(this,name,rank,parts)
!--------------------------------------------
implicit none
class(array) this
class(dao), target, optional :: parts
integer, optional :: rank
character(len=*), optional :: name

if(present(rank)) this%rank = rank
if(present(parts)) this%parts => parts
if(present(name)) this%name = name
end subroutine

!--------------------------------------------
     subroutine max_alloc(this)
!--------------------------------------------
implicit none
class(array) this
integer dimrank, maxn
dimrank = this%get_dim()**this%rank
maxn = this%get_maxn()

if(dimrank==1)then
  allocate(this%r(maxn))
elseif(dimrank==2)then
  allocate(this%x,this%y)
  allocate(this%x%r(maxn),this%y%r(maxn))
  this%x%parts=>this%parts
  this%y%parts=>this%parts
  this%x%name = trim(this%name)//'.x'
  this%y%name = trim(this%name)//'.y'
elseif(dimrank==4)then
  allocate(this%x,this%y,this%xy,this%z)
  allocate(this%x%r(maxn),this%y%r(maxn),this%xy%r(maxn),this%z%r(maxn))
  this%x%parts=>this%parts
  this%y%parts=>this%parts
  this%xy%parts=>this%parts  
  this%z%parts=>this%parts
  this%x%name = trim(this%name)//'.x'
  this%y%name = trim(this%name)//'.y'
  this%xy%name = trim(this%name)//'.xy'
  this%z%name = trim(this%name)//'.z'
endif  

end subroutine

!--------------------------------------------
     subroutine alloc(this)
!--------------------------------------------
implicit none
class(array) this
integer dimrank, size
dimrank = this%get_dim()**this%rank
size = this%get_size()

if(dimrank==1)then
  allocate(this%r(size))
elseif(dimrank==2)then
  allocate(this%x,this%y)
  allocate(this%x%r(size),this%y%r(size))
  this%x%parts=>this%parts
  this%y%parts=>this%parts
  this%x%name = trim(this%name)//'.x'
  this%y%name = trim(this%name)//'.y'  
elseif(dimrank==4)then
  allocate(this%x,this%y,this%xy)
  allocate(this%x%r(size),this%y%r(size),this%xy%r(size))
  this%x%parts=>this%parts
  this%y%parts=>this%parts
  this%xy%parts=>this%parts    
  this%x%name = trim(this%name)//'.x'
  this%y%name = trim(this%name)//'.y'
  this%xy%name = trim(this%name)//'.xy'  
endif  

end subroutine

!--------------------------------------------
    subroutine array_equal_array(a,b)
!--------------------------------------------        
implicit none
type(array),intent(INOUT) :: a
type(array),intent(IN)    :: b
integer ndim1,i

!if(a%ndim1/=b%ndim1)stop 'Error: inequal length arrays!'
if(.not.associated(a%parts,b%parts))stop 'Error: inequal length arrays!'

!ndim1 = a%ndim1
ndim1 = a%get_size()

!$omp parallel do 
do i = 1, ndim1
   a%r(i) = b%r(i)
enddo
!$omp end parallel do

return
end subroutine

!--------------------------------------------
    subroutine array_equal_double_real(a,b)
!--------------------------------------------        
implicit none
type(array),intent(INOUT) :: a
real(dp),intent(IN)    :: b
integer ndim1,i

!ndim1 = a%ndim1
!if(ndim1==0) stop 'array_equal_real: ndim1=0'  !!!
ndim1 = a%get_size()

!$omp parallel do
do i = 1, ndim1
   a%r(i) = b
enddo
!$omp end parallel do

return
end subroutine
    
!--------------------------------------------
    subroutine array_equal_real(a,b)
!--------------------------------------------        
implicit none
type(array),intent(INOUT) :: a
real,intent(IN)    :: b
integer ndim1,i

!ndim1 = a%ndim1
!if(ndim1==0) stop 'array_equal_real: ndim1=0' !!!
ndim1 = a%get_size()

!$omp parallel do
do i = 1, ndim1
   a%r(i) = b
enddo
!$omp end parallel do

return
end subroutine
!--------------------------------------------
!      function cmpt(this,i) result(val)
!--------------------------------------------
!implicit none
!class(array) this
!integer i
!real(dp) val(3)

!val(1) = this%x%r(i)
!val(2) = this%y%r(i)
!val(3) = 0

!end function

!--------------------------------------------
     function cmpt(this,i) result(p)
!--------------------------------------------
implicit none
class(array) this
integer i
type(p2r) p(3)

p(1)%p => this%x%r(i)
p(2)%p => this%y%r(i)
p(3)%p => null()

end function

!--------------------------------------------
       subroutine array_final(this)
!--------------------------------------------
implicit none
type(array) :: this

if(associated(this%r))deallocate(this%r)

return
end subroutine

!--------------------------------------------
    function array_add_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a, b
type(array), allocatable :: c
integer ndim1,i

!if(a%ndim1/=b%ndim1)stop 'Cannot add arrays!'
if(.not.associated(a%parts,b%parts)) stop 'Cannot add arrays!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) + b%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_add_double_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: b
type(array), allocatable :: c
integer ndim1,i

!if(a%ndim1==0)stop 'array_add_double_real: ndim1=0!'
if(.not.associated(a%parts)) stop 'array_add_double_real: parts not associated!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) + b
enddo
!$omp end parallel do
c%parts => a%parts
end function
    
!--------------------------------------------
    function double_real_add_array(b,a) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: b
type(array), allocatable :: c
integer ndim1,i

!if(a%ndim1==0)stop 'double_real_add_array: ndim1=0!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) + b
enddo
!$omp end parallel do
c%parts => a%parts
end function    
    
!--------------------------------------------
    function array_sub_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a, b
type(array), allocatable:: c
integer ndim1,i

!if(a%ndim1/=b%ndim1)stop 'Cannot substract arrays!'
if(.not.associated(a%parts,b%parts))stop 'Cannot substract arrays!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) - b%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_sub_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: b
type(array), allocatable:: c
integer ndim1,i

!if(a%ndim1/=b%ndim1)stop 'Cannot substract arrays!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) - b
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_sub_double_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: b
type(array), allocatable:: c
integer ndim1,i

!if(a%ndim1/=b%ndim1)stop 'Cannot substract arrays!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) - b
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_minus(a) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
type(array), allocatable:: c
integer ndim1,i

!if(a%ndim1==0)stop 'array_minus: ndim1=0!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = -a%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_mul_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a,b
type(array), allocatable:: c
integer ndim1,i

!if(a%ndim1/=b%ndim1)stop 'Cannot multiply arrays!'
if(.not.associated(a%parts,b%parts))stop 'Cannot multiply arrays!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)*b%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_mul_real(a,r) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: r
type(array), allocatable:: c
integer ndim1,i

!if(a%ndim1==0)stop 'array_mul_real: ndim1=0!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------------
    function array_mul_double_real(a,r) result(c)
!--------------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: r
type(array), allocatable :: c
integer ndim1,i

!if(a%ndim1==0)stop 'array_mul_double_real: ndim1=0!'
!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------------
    function double_real_mul_array(r,a) result(c)
!--------------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: r
type(array), allocatable :: c
integer ndim1,i

!if(a%ndim1==0)stop 'double_real_mul_array: ndim1=0!'
!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function    

!--------------------------------------------------
    function real_mul_array(r,a) result(c)
!--------------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: r
type(array), allocatable :: c
integer ndim1,i

!if(a%ndim1==0)stop 'real_mul_array: ndim1=0!'
!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function    
    
!--------------------------------------------
    function array_div_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a,b
type(array), allocatable:: c
integer ndim1,i

!if(a%ndim1/=b%ndim1)stop 'Cannot divide arrays!'
if(.not.associated(a%parts,b%parts))stop 'Cannot divide arrays!'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)/b%r(i)
enddo
!$omp end parallel do
c%parts => a%parts
end function


!--------------------------------------------
    function array_div_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: b
type(array), allocatable:: c
integer ndim1,i

if(b==0)stop 'Cannot divide arrays! b=0'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)/b
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_div_double_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: b
type(array), allocatable:: c
integer ndim1,i

if(b==0)stop 'Cannot divide arrays! b=0'

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)/b
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_pow_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: b
type(array), allocatable:: c
integer ndim1,i

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)**b
enddo
!$omp end parallel do
c%parts => a%parts
end function

!--------------------------------------------
    function array_pow_double_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: b
type(array), allocatable:: c
integer ndim1,i

!ndim1 = a%ndim1
ndim1 = a%get_size()
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)**b
enddo
!$omp end parallel do
c%parts => a%parts
end function

end module
