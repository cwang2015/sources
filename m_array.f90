!---------------------------
    module m_array
!---------------------------
implicit none
integer, parameter :: dp = kind(0.d0)

type array
   integer :: ndim1
   real(dp), pointer, dimension(:) :: r  => null()
   type(array), pointer :: x => null(), y => null(), z => null()
   type(array), pointer :: xy => null(), xz => null(), yz => null()
   contains
       final :: array_final
end type

interface assignment(=)
  module procedure :: array_equal_array
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

end module
