!---------------------------
    module m_array
!---------------------------
implicit none
integer, parameter :: dp = kind(0.d0)

type p2r
     real(dp), pointer :: p => null()
end type 

type array
   integer :: ndim1
   real(dp), pointer, dimension(:) :: r  => null()
   type(array), pointer :: x => null(), y => null(), z => null()
   type(array), pointer :: xy => null(), xz => null(), yz => null()
   contains
       procedure :: cmpt
       !procedure :: p2cmpt
       final :: array_final
end type

interface assignment(=)
  module procedure :: array_equal_array
  module procedure :: array_equal_double_real
  module procedure :: array_equal_real
end interface

interface operator(+)
   module procedure :: array_add_array 
   module procedure :: array_add_double_real 
   module procedure :: double_real_add_array 
end interface

interface operator(-)
   module procedure :: array_sub_array
   module procedure :: array_minus
end interface

interface operator(*)
   module procedure :: array_mul_array
   module procedure :: array_mul_real
   module procedure :: array_mul_double_real
   module procedure :: double_real_mul_array
   module procedure :: real_mul_array
end interface

interface operator(/)
   module procedure :: array_div_array
   module procedure :: array_div_real
end interface

!=======
contains
!=======

!--------------------------------------------
    subroutine array_equal_array(a,b)
!--------------------------------------------        
implicit none
type(array),intent(INOUT) :: a
type(array),intent(IN)    :: b
integer ndim1,i

if(a%ndim1/=b%ndim1)stop 'Error: inequal length arrays!'

ndim1 = a%ndim1

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

ndim1 = a%ndim1
!if(ndim1==0) stop 'array_equal_real: ndim1=0'  !!!

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

ndim1 = a%ndim1
!if(ndim1==0) stop 'array_equal_real: ndim1=0' !!!

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
type(array) this

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

if(a%ndim1/=b%ndim1)stop 'Cannot add arrays!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) + b%r(i)
enddo
!$omp end parallel do

end function

!--------------------------------------------
    function array_add_double_real(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: b
type(array), allocatable :: c
integer ndim1,i

if(a%ndim1==0)stop 'array_add_double_real: ndim1=0!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) + b
enddo
!$omp end parallel do

end function
    
!--------------------------------------------
    function double_real_add_array(b,a) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: b
type(array), allocatable :: c
integer ndim1,i

if(a%ndim1==0)stop 'double_real_add_array: ndim1=0!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) + b
enddo
!$omp end parallel do

end function    
    
!--------------------------------------------
    function array_sub_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a, b
type(array), allocatable:: c
integer ndim1,i

if(a%ndim1/=b%ndim1)stop 'Cannot substract arrays!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i) - b%r(i)
enddo
!$omp end parallel do

end function

!--------------------------------------------
    function array_minus(a) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
type(array), allocatable:: c
integer ndim1,i

if(a%ndim1==0)stop 'array_minus: ndim1=0!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = -a%r(i)
enddo
!$omp end parallel do

end function

!--------------------------------------------
    function array_mul_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a,b
type(array), allocatable:: c
integer ndim1,i

if(a%ndim1/=b%ndim1)stop 'Cannot multiply arrays!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)*b%r(i)
enddo
!$omp end parallel do

end function

!--------------------------------------------
    function array_mul_real(a,r) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: r
type(array), allocatable:: c
integer ndim1,i

if(a%ndim1==0)stop 'array_mul_real: ndim1=0!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do

end function

!--------------------------------------------------
    function array_mul_double_real(a,r) result(c)
!--------------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: r
type(array), allocatable :: c
integer ndim1,i

if(a%ndim1==0)stop 'array_mul_double_real: ndim1=0!'
ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do

end function

!--------------------------------------------------
    function double_real_mul_array(r,a) result(c)
!--------------------------------------------------
implicit none
type(array), intent(in) :: a
real(dp), intent(in) :: r
type(array), allocatable :: c
integer ndim1,i

if(a%ndim1==0)stop 'double_real_mul_array: ndim1=0!'
ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do

end function    

!--------------------------------------------------
    function real_mul_array(r,a) result(c)
!--------------------------------------------------
implicit none
type(array), intent(in) :: a
real, intent(in) :: r
type(array), allocatable :: c
integer ndim1,i

if(a%ndim1==0)stop 'real_mul_array: ndim1=0!'
ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = r*a%r(i)
enddo
!$omp end parallel do

end function    
    
!--------------------------------------------
    function array_div_array(a,b) result(c)
!--------------------------------------------
implicit none
type(array), intent(in) :: a,b
type(array), allocatable:: c
integer ndim1,i

if(a%ndim1/=b%ndim1)stop 'Cannot divide arrays!'

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)/b%r(i)
enddo
!$omp end parallel do

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

ndim1 = a%ndim1
allocate(c)
allocate(c%r(ndim1))
c%ndim1 = ndim1

!$omp parallel do
do i = 1, ndim1
   c%r(i) = a%r(i)/b
enddo
!$omp end parallel do

end function

end module
