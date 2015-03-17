!----------------------------------------
             module m_anyArr
!----------------------------------------
use constants
use m_primitives
implicit none

private


!CONTAINERS
!==========

type, extends(dao), public :: array    !everything

character(len=32) :: nam=' '

!I'm a field or matrix?
!ima =1 field; =2 matrix
integer :: ima = 0

!数组参数
!type=1 int_array; =2 real_array;  
!=3 anyother array stored in d(:)
integer :: type = 2
integer :: ndim1, ndim2 = 1
integer :: maxdim1, maxdim2 = 1

integer, pointer, dimension(:) :: i => null()    ! 整型
real(dp), pointer, dimension(:) :: r => null()   ! 实型
class(dao), pointer, dimension(:) :: d => null() ! 其他类型
class(*), pointer, dimension(:) :: cmpt => null()

!场的参数

!loc='vol': volum filed; ='face': face field; 
!   ='vtx': vertex field
!   ='any': not related with a mesh
character(len=16) :: loc='any'
type(discrete_scheme), pointer :: scheme

!dim: 维数，2维还是3维；
!rank＝0: 标量场；＝1：矢量场；＝2：张量场
!取决于dim^rank 
integer :: dim = 2, rank = 0
type(array), pointer :: x => null(), y => null(), z => null()

!矩阵参数
!For matrix storage: Element By Element（EBE）or
!                    Edge Based Storage (EBS)
!kstor=0, any matrix; =1, EBE; =2, EBS
integer :: kstor = 1
!Matrix property: symmetry or not
integer :: ksymm = 0

!For EBE  
type(array), pointer :: diag=>null(), skew=>null()
!For EBS
type(array), pointer :: owner=>null(), neibor=>null()

!Surrogate
class(dao), pointer :: smesh
!与本场关联的几何对象
type(array), pointer :: geo

contains

  !常用方法：
  procedure :: write => array_write_sub
  procedure :: alloc => xallocate
  procedure :: dealloc => array_deallocate_sub
  procedure :: set => array_set_sub
  procedure :: info => array_get_information_sub 

  procedure :: cpstrc => copy_array_structure
  !procedure :: cmpts => components
  !procedure :: p2cmpt
  !procedure :: xcmpt
  !procedure :: AA  ! Array Algebra
  !procedure :: FA  ! Field Algebra
  procedure :: matvec

  procedure :: almat => matrix_allocate_sub
  procedure :: length => array_get_length_fn 

  !procedure :: one
  !procedure :: done
  procedure :: zero
  !procedure :: dzero
  !procedure :: vzero
  procedure :: one_to_N
  !procedure :: one_to_Nd0
  
  !procedure :: add_array => add_array_sub
  !procedure :: add_real  => add_real_sub 
  !generic,public :: add => add_array

  FINAL :: clean

end type

type, public :: equation
  type(array), pointer :: A,X,B
  contains
    procedure :: equation => construct_equation
end type

! Array constructor
interface array
  module procedure constructor
end interface

interface operator (+)

  !利用allocatable :: ar实现加法
  
  !module procedure array_add_fn
  !module procedure array_add_real 
  !module procedure real_add_array

  !利用临时工作数组w1,w2和指针pw实现加法
  !module procedure array_addition 
end interface

interface operator (-)
  !module procedure array_sub_fn
end interface

interface operator (*)
  !module procedure array_mul_fn
  !module procedure array_mul_real
  module procedure real_mul_array
  !module procedure array_multiply
end interface

interface dot_product
  !module procedure array_dot_product
end interface

interface assignment (=)
  module procedure array_equal_sub
  !module procedure array_equal_real
  !module procedure array_equal_integer
end interface

!public :: operator(+), operator(-), operator(*),  &
!          dot_product
public :: assignment(=), operator(*)

!type(array), target  :: w(3)
!type(array), pointer :: pw => null()
!integer :: iw = 0

contains

! Array constructor
!-----------------------------------------------------
           function constructor(this)
!-----------------------------------------------------
implicit none
class(array)this
type(array) constructor

end function

!-----------------------------------------------------
     recursive subroutine array_write_sub(this)
!-----------------------------------------------------
   implicit none
   class(array) this
   integer ndim1, ndim2, i, j, k
   logical :: debug = .false.
 
   if(debug) write(*,*) 'Debug array_write_sub ...'

   ndim1 = this%ndim1; ndim2 = this%ndim2

   if(this%type==1)then
      if(ndim2==1)then
         write(*,*) this%i
      else
         do i=1,ndim1
            k = ndim2*(i-1)
            write(*,*) (this%i(k+j),j=1,ndim2)
         enddo
      endif
   elseif(this%type==2)then
     write(*,*) this%r
   elseif(this%type==3)then
     call write_dao_array(this%d)
   else
     write(*,*) 'I dont know its type!'
   endif

   return
   end subroutine

!-------------------------------------------------
            subroutine xallocate(this,src)
!-------------------------------------------------
   implicit none
   class(array) this
   class(dao), optional :: src
   integer length
   logical :: debug = .false.

   if(debug) write(*,*) 'Entry: scalar_array_allocate_sub!'

   length = this%length()

   if(this%type==1)then
     allocate(this%i(length))
   elseif(this%type==2)then
     allocate(this%r(length))
   else
     if(.not.present(src)) Stop 'Source should present!'
     allocate(this%d(length), source=src)
     !allocate(this%d(length), mold=src)
   endif

   return
   end subroutine

!!-----------------------------------------------------------
!      recursive subroutine allocate_vector_array(this)
!!-----------------------------------------------------------
!implicit none
!class(array) this
!integer dimrank
!
!dimrank = this%dim**this%rank
!if(dimrank==1)then
!  return
!elseif(dimrank==2)then
! allocate(this%x,this%y)
! call this%x%cpstrc(this)
! call this%y%cpstrc(this)
! this%x%rank=0 
! this%y%rank=0
! call xalloc(this%x,this%y)
!elseif(dimrank==4)then
! call this%x%alvec
! call this%y%alvec
!endif
!
!return
!end subroutine

!----------------------------------------------------
       subroutine matrix_allocate_sub(this)
!----------------------------------------------------
implicit none
class(array) this

if(this%kstor==0)then
   return
elseif(this%kstor==1)then !EBE
   allocate(this%diag, this%skew)
elseif(this%kstor==2)then !EBS
   allocate(this%diag, this%owner, this%neibor)
endif

return
end subroutine


!!----------------------------------------------------
!   recursive subroutine array_allocate_sub(this,src)
!!----------------------------------------------------
!   implicit none
!   class(array) this
!   class(dao), optional :: src
!   integer length, dimrank
!   logical :: debug = .false.
!
!   if(debug) write(*,*) 'Entry: array_allocate_sub!'
!
!                 if(this%ima==0)then
!       
!   call xalloc(this,src)                  
!
!                 elseif(this%ima==1)then       ! I am a field!
!                                               ! alvec
!!   dimrank = this%dim**this%rank
!!   if(this%rank==0)then
!!      call xalloc(this,src)
!!   elseif(this%rank==1)then
!!      if(this%dim==2)then
!!         if(.not.associated(this%x)) allocate(this%x)
!!         if(.not.associated(this%y)) allocate(this%y)
!!         this%x%ndim1 = this%ndim1
!!         this%y%ndim1 = this%ndim1
!!         this%x%type = this%type
!!         this%y%type = this%type
!!         call this%x%alloc
!!         call this%y%alloc
!!      endif
!!   endif
!
!   call this%alvec
!
!                 elseif(this%ima==2) then   ! I am a matrix!
!                                            ! almat
!!  if(this%kstor==0)then
!!      call xalloc(this)
!!  elseif(this%kstor==1)then !EBS
!!      allocate(this%owner, this%neibor)
!!      this%owner%ndim1 = this%ndim1 
!!      this%neibor%ndim1 = this%ndim1
!!      this%owner%type = this%type
!!      this%neibor%type = this%type
!!      call this%owner%alloc
!!      if(this%ksymm==0)then
!!         call this%neibor%alloc
!!      else
!!         this%neibor=>this%owner
!!      endif
!!   elseif(this%kstor==2)then !EBE
!!      allocate(this%diag, this%skew)
!!   endif
!
!    call this%almat
!
!               endif
!
!   return
!   end subroutine

!----------------------------------------------------
       subroutine array_deallocate_sub(this)
!----------------------------------------------------
implicit none
class(array) this
integer istat  

!if(associated(this%cmpt)) deallocate(this%cmpt)
if(associated(this%d))then
     deallocate(this%d); this%d=>null()
endif
if(associated(this%i))then
    deallocate(this%i); this%i=>null()
endif
if(associated(this%r))then
    deallocate(this%r); this%r=>null()
endif
!if(associated(this%v2))then
!  deallocate(this%v2, stat=istat)
!  this%v2=>null()
!endif
if(associated(this%x))then
    deallocate(this%x); this%x=>null()
endif
if(associated(this%y))then
    deallocate(this%y); this%y=>null()
endif
if(associated(this%diag))then
    deallocate(this%diag); this%diag=>null()
endif
if(associated(this%skew))then
    deallocate(this%skew); this%skew=>null()
endif
if(associated(this%owner))then
    deallocate(this%owner); this%owner=>null()
endif
if(associated(this%neibor))then
    deallocate(this%neibor); this%neibor=>null()
endif

return
end subroutine

!----------------------------------------------------------------
   subroutine array_set_sub(this,ima, ndim1,ndim2,type,   &
              dim, rank, kstor, ksymm, smesh, nam, loc)
!----------------------------------------------------------------
   implicit none
   class(array) this
   class(dao), target, optional :: smesh
   integer, optional :: ima, ndim1,ndim2,type,dim,rank
   integer, optional :: kstor, ksymm
   character(len=*), optional :: nam, loc

   !write(*,*) 'In m_anyArr::array_setup_sub'
   if(present(nam)) this%nam=nam
   if(present(ima)) this%ima=ima 
   if(present(ndim1)) this%ndim1=ndim1
   if(present(ndim2)) this%ndim2=ndim2
   if(present(type)) this%type=type
   if(present(dim)) this%dim = dim
   if(present(rank)) this%rank = rank
   if(present(kstor)) this%kstor = kstor
   if(present(ksymm)) this%ksymm = ksymm
   if(present(smesh)) this%smesh => smesh
   if(present(loc)) this%loc=loc

   return
   end subroutine

!----------------------------------------------------------------
   subroutine copy_array_structure(this,arr)   
!----------------------------------------------------------------
   implicit none
   class(array) this
   type(array) arr

   !write(*,*) 'In m_anyArr::copy_array_structure'

   this%nam   = arr%nam
   this%ima   = arr%ima
   this%ndim1 = arr%ndim1
   this%ndim2 = arr%ndim2
   this%type  = arr%type
   this%dim   = arr%dim
   this%rank  = arr%rank
   this%kstor = arr%kstor
   this%ksymm = arr%ksymm
   this%smesh =>arr%smesh
   this%loc   = arr%loc

   return
   end subroutine

!-----------------------------------------------------
   function array_get_length_fn(this) result(length)
!-----------------------------------------------------
   implicit none 
   class(array) this
   integer length
   length= this%ndim1*this%ndim2
   return
   end function

!--------------------------------------------------------------
        subroutine array_get_information_sub(this)
!--------------------------------------------------------------
implicit none
class(array) this

write(*,*) 'My name is: ', this%nam
write(*,*) 'type=', this%type, ' ndim1=', this%ndim1
write(*,*) 'ndim2=', this%ndim2
write(*,*) 'ndim=', this%dim, ' rank=', this%rank

if(associated(this%d)) write(*,*) 'd associated!'
if(associated(this%r)) write(*,*) 'r associated!'
if(associated(this%x)) write(*,*) 'x associated!'
if(associated(this%y)) write(*,*) 'y assocated!'
if(associated(this%i)) write(*,*) 'i assocated!'

return
end subroutine

!-----------------------------------------------------------
               subroutine components(this)
!-----------------------------------------------------------
implicit none
class(array) this
type(vector2d), pointer, dimension(:) :: pv

if(.not.associated(this%x)) allocate(this%x)
if(.not.associated(this%y)) allocate(this%y)
call this%x%cpstrc(this); this%x%type=2
call this%y%cpstrc(this); this%y%type=2

call p2m(pv,this%d)
this%x%r=>pv%x
this%y%r=>pv%y

return
end subroutine

!-----------------------------------------------------------
               subroutine p2cmpt(this,p,cmpt) 
!-----------------------------------------------------------
implicit none
class(array) this
type(array)  p
character(len=*) cmpt
type(vector2d), pointer, dimension(:) :: pv

call this%cpstrc(p); this%type=2
call p2m(pv,p%d)
if(trim(cmpt)=='x')then
   this%r=>pv%x
elseif(trim(cmpt)=='y')then
   this%r=>pv%y
endif

return
end subroutine

!Better to use p2cmpt!!
!-----------------------------------------------------------
               function xcmpt(this) result(px)
!-----------------------------------------------------------
! Usage: p=>a%xcmpt()
!-----------------------------------------------------------
implicit none
class(array) this
type(array), pointer :: px
!type(array), target, save  :: x              ! Dangerous!
!type(array), target, allocatable  :: x
type(vector2d), pointer, dimension(:) :: pv

allocate(px)
call px%cpstrc(this); px%type=2
call p2m(pv,this%d)
px%r=>pv%x
!deallocate(px)                              ! Memory Leak!

!allocate(x)
!call x%cpstrc(this); x%type=2
!call p2m(pv,this%d)
!x%r=>pv%x
!px => x
!deallocate(x)

return
end function

!!=================Generic interface (=)====================
!-----------------------------------------------------------
       subroutine array_equal_sub(a,b)
!-----------------------------------------------------------
implicit none
type(array), intent(INOUT) :: a
type(array), intent(IN)    :: b

if(a%type==1)then
  a%i = b%i
elseif(a%type==2)then
  a%r = b%r
elseif(a%type==3.and.b%type==3)then
   if(.not.same_type_as(a%d,b%d)) stop 'Assignment not allowed!'
   select type (p1=>a%d)
      type is (vector2d)
   select type (p2=>b%d)
      type is (vector2d)
           p1=p2
   end select
   end select
endif

return
end subroutine

!-----------------------------------------------------------
       subroutine array_equal_real(a,r)
!-----------------------------------------------------------
implicit none
type(array), intent(INOUT) :: a
real(dp) r

if(a%type==2)then
  a%r = r
endif

return
end subroutine

!-----------------------------------------------------------
       subroutine array_equal_integer(a,i)
!-----------------------------------------------------------
implicit none
type(array), intent(INOUT) :: a
integer i

if(a%type==1)then
  a%i = i
endif

return
end subroutine

!DEC$IF(.FALSE.)
! Better than array_add_fn
!利用临时工作数组w1,w2和指针pw实现加法
!-----------------------------------------------------------
           function array_addition(a,b) result(add)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a, b
type(array), pointer :: add

!if(associated(pw,w1))then
!  pw=>w2
!  !write(*,*) 'associated w2'
!else
!  pw=>w1
!  !write(*,*) 'associated w1'
!endif

!iw=iw+1
!if(iw>3) iw=0
!pw=>w(iw)
!add => pw
!if(associated(add%r)) deallocate(add%r)
!allocate(add%r(size(a%r)))
!add%r = a%r + b%r

select type(pa=>a%cmpt)
  type is(real(dp))
     select type(pb=>b%cmpt)
       type is(real(dp))
          add%r = pa+pb
     end select
  type in(dao)
end select

end function

!-----------------------------------------------------------
           function array_multiply(a,b) result(mul)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a, b
type(array), pointer :: mul

!if(associated(pw,w1))then
!  pw=>w2
!  !write(*,*) 'associated w2'
!else
!  pw=>w1
!  !write(*,*) 'associated w1'
!endif

iw=iw+1
if(iw>3) iw=0
pw=>w(iw)
mul => pw
if(associated(mul%r)) deallocate(mul%r)
allocate(mul%r(size(a%r)))
mul%r = a%r * b%r

end function

!DEC$ENDIF
!!=================Generic bounding : add ==================
!-----------------------------------------------------------
           subroutine add_array_sub(this,a)
!-----------------------------------------------------------
implicit none
class(array) this
type(array) a

if(this%type==1)then
   this%i = this%i + a%i
elseif(this%type==2)then
   this%r = this%r + a%r
endif

return
end subroutine

!DEC$IF(.FALSE.)
!!=================Generic interface (+)====================
!-----------------------------------------------------------
           function array_add_fn(a,b) result(add)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a, b
type(array), target :: add

if(a%type==1)then
   allocate(add%wi(size(a%i)))
   add%i => add%wi
   add%i = a%i + b%i
elseif(a%type==2)then
   allocate(add%wr(size(a%r)))
   add%r => add%wr
   add%r = a%r + b%r
endif

return
end function

!-----------------------------------------------------------
           function array_add_real(a,r) result(add)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a
real(dp), intent(IN) :: r
type(array), target :: add

if(a%type==1)then
   allocate(add%wi(size(a%i)))
   add%i => add%wi
   add%i = a%i + r
elseif(a%type==2)then
   allocate(add%wr(size(a%r)))
   add%r => add%wr
   add%r = a%r + r
endif

return
end function

!-----------------------------------------------------------
           function real_add_array(r,a) result(add)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a
real(dp), intent(IN) :: r
type(array), target :: add

if(a%type==1)then
   allocate(add%wi(size(a%i)))
   add%i => add%wi
   add%i = a%i + r
elseif(a%type==2)then
   allocate(add%wr(size(a%r)))
   add%r => add%wr
   add%r = a%r + r
endif

return
end function

!!=================Generic interface (-)====================
!-----------------------------------------------------------
           function array_sub_fn(a,b) result(sub)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a, b
type(array), target :: sub
!type(array), target, allocatable :: sub

!allocate(sub)
if(a%type==1)then
   allocate(sub%wi(size(a%i)))
   sub%i => sub%wi
   sub%i = a%i - b%i
elseif(a%type==2)then
   allocate(sub%wr(size(a%r)))
   sub%r => sub%wr
   sub%r = a%r - b%r
endif

!call sub%set(nam='sub',type=1)
!allocate(sub%i(size(a%i)))
!sub%i = a%i-b%i
!write(*,*) 'sub%i= ', sub%i

return
end function

!!=================Generic interface (*)====================
!-----------------------------------------------------------
           function array_mul_fn(a,b) result(mul)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a, b
type(array), target :: mul
type(vector2d), pointer, dimension(:) :: pa, pb

if(a%type==1)then
   allocate(mul%wi(size(a%i)))
   mul%i => mul%wi
   mul%i = a%i * b%i
elseif(a%type==2.and.b%type==2)then
   allocate(mul%wr(size(a%r)))
   mul%r => mul%wr
   mul%r = a%r * b%r
elseif(a%type==3.and.b%type==3)then
   allocate(mul%wr(size(a%d)))
   call p2m(pa,a%d); call p2m(pb,b%d)
   mul%r => mul%wr
   mul%r = pa%x*pb%x+pa%y*pb%y
endif

return
end function

!-----------------------------------------------------------
           function array_mul_real(a,r) result(mul)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a
real(dp), intent(IN) :: r
type(array), target :: mul

!if(a%type==1)then
!   allocate(mul%ai(size(a%i)))
!   mul%i => mul%ai
!   mul%i = a%i * b%i
if(a%type==2)then
   allocate(mul%wr(size(a%r)))
   mul%r => mul%wr
   mul%r = a%r * r
endif

return
end function

!DEC$ENDIF
!-----------------------------------------------------------
           function real_mul_array(r,a) result(mul)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a
real(dp), intent(IN) :: r
type(array), target :: mul

!if(a%type==1)then
!   allocate(mul%ai(size(a%i)))
!   mul%i => mul%ai
!   mul%i = a%i * b%i
if(a%type==2)then
   mul%type=2
   allocate(mul%r(size(a%r)))
   mul%r = a%r * r
endif

return
end function

!DEC$IF(.FALSE.)
!!============Generic dot_product===========================
!-----------------------------------------------------------
           function array_dot_product(a,b) result(dot)
!-----------------------------------------------------------
implicit none
type(array), intent(IN) :: a, b
real(dp) dot

!if(a%type==1)then
!   allocate(mul%ai(size(a%i)))
!   mul%i => mul%ai
!   mul%i = a%i * b%i
!if(a%type==2)then
   dot = dot_product(a%r,b%r)
!endif

return
end function

!DEC$ENDIF

!-----------------------------------------------------------
    recursive subroutine AA(x, op, y, z, ic, rc, dc, mask)
!-----------------------------------------------------------
implicit none
class(array) x
character(len=*) op
type(array), optional :: y, z
integer, optional, dimension(:) :: mask
integer, optional :: ic
real(dp), optional :: rc
class(dao), optional :: dc

integer i
                               
if(trim(op)=='X=C')then
  if(x%type==1.and.present(ic))then
    x%i = ic
  elseif(x%type==2.and.present(rc))then
    x%r = rc
  else
    write(*,*) 'This case has not been implemented!'
  endif

elseif(trim(op)=='X=CX')then 
                            if(x%rank==0)then
  x%r = rc*x%r
                            elseif(x%rank==1)then
                              call aa(x%x,op,rc=rc)
                              call aa(x%y,op,rc=rc)
                            endif

elseif(trim(op)=='X=Y')then
                            if(x%rank==0)then
  if(x%type==2)then
     if(present(mask))then
          do i = 1, size(mask)   !从这里看，将bfaces等定义为array也是必要的
             x%r(mask(i))=y%r(i)
          enddo
      else
       x%r = y%r
     endif
  endif
                             elseif(x%rank==1)then
                                call aa(x%x,op,y%x,mask=mask)
                                call aa(x%y,op,y%y,mask=mask)
                             endif

elseif(op=='X=CY    ')then
  if(x%type==1.and.y%type==1.and.present(ic))then
    x%i = ic*y%i
  elseif(x%type==2.and.y%type==2.and.present(rc))then
    x%r = rc*y%r
  else
    write(*,*) 'This case has not been implemented!'
  endif

elseif(trim(op)=='X=XY')then
  x%r = x%r*y%r

elseif(op=='X=YZ    ')then
                              !if(x%rank==0)then
  if(x%type==1.and.y%type==1.and.z%type==1)then
     x%i = y%i*z%i
  elseif(x%type==1.and.y%type==1.and.z%type==2)then
     x%i = y%i*z%r
  elseif(x%type==2.and.y%type==2.and.z%type==1)then
     x%r = y%r*z%i
  elseif(x%type==2)then
     x%r = y%r*z%r
  else
    write(*,*) 'This case has not been implemented!'
  endif
                              !elseif(x%rank==1)then
                                !... 
                              !endif

elseif(op=='X+=YZ   ')then
                              !if(x%rank==0)then
  if(x%type==1.and.y%type==1.and.z%type==1)then
     x%i = x%i+y%i*z%i
  elseif(x%type==1.and.y%type==1.and.z%type==2)then
     x%i = x%i+y%i*z%r
  elseif(x%type==2.and.y%type==2.and.z%type==1)then
     x%r = x%r+y%r*z%i
  elseif(x%type==2)then
     x%r = x%r+y%r*z%r
  else
    write(*,*) 'This case has not been implemented!'
  endif
                              !elseif(x%rank==1)then
                                !... 
                              !endif

elseif(op=='X=CYZ    ')then

  if(x%type==1.and.y%type==1.and.z%type==1)then
     x%i = y%i*z%i
  elseif(x%type==1.and.y%type==1.and.z%type==2)then
     x%i = y%i*z%r
  elseif(x%type==2.and.y%type==2.and.z%type==1)then
     x%r = y%r*z%i
  elseif(x%type==2.and.y%type==2.and.z%type==2)then
     x%r = rc*y%r*z%r
  else
    write(*,*) 'This case has not been implemented!'
  endif

elseif(trim(op)=='X=X+Y')then
  
  x%r=x%r+y%r  

else
  write(*,*) 'This case has not been implemented!'
endif

return
end subroutine

!DEC$IF(.FALSE.)
!Field Algebra
!-----------------------------------------------------------
    recursive subroutine FA(x, op, y, z, ic, rc, dc, mask)
!-----------------------------------------------------------
implicit none
class(array) x
character(len=*) op
type(array), optional :: y, z
integer, optional, dimension(:) :: mask
integer, optional :: ic
real(dp), optional :: rc
class(dao), optional :: dc

integer i
                               
if(trim(op)=='X=C')then
  if(x%type==1.and.present(ic))then
    x%i = ic
  elseif(x%type==2.and.present(rc))then
    x%r = rc
  else
    write(*,*) 'This case has not been implemented!'
  endif

elseif(trim(op)=='X=CX')then 
                            if(x%rank==0)then
  x%r = rc*x%r
                            elseif(x%rank==1)then
                              call aa(x%x,op,rc=rc)
                              call aa(x%y,op,rc=rc)
                            endif

elseif(trim(op)=='X=Y')then
                            if(x%rank==0)then
  if(x%type==2)then
     if(present(mask))then
          do i = 1, size(mask)   !从这里看，将bfaces等定义为array也是必要的
             x%r(mask(i))=y%r(i)
          enddo
      else
       x%r = y%r
     endif
  endif
                             elseif(x%rank==1)then
                                call aa(x%x,op,y%x,mask=mask)
                                call aa(x%y,op,y%y,mask=mask)
                             endif

elseif(op=='X=CY    ')then
  if(x%type==1.and.y%type==1.and.present(ic))then
    x%i = ic*y%i
  elseif(x%type==2.and.y%type==2.and.present(rc))then
    x%r = rc*y%r
  else
    write(*,*) 'This case has not been implemented!'
  endif

elseif(trim(op)=='X=XY')then
  x%r = x%r*y%r

elseif(op=='X=YZ    ')then
                              !if(x%rank==0)then
  if(x%type==1.and.y%type==1.and.z%type==1)then
     x%i = y%i*z%i
  elseif(x%type==1.and.y%type==1.and.z%type==2)then
     x%i = y%i*z%r
  elseif(x%type==2.and.y%type==2.and.z%type==1)then
     x%r = y%r*z%i
  elseif(x%type==2)then
     x%r = y%r*z%r
  else
    write(*,*) 'This case has not been implemented!'
  endif
                              !elseif(x%rank==1)then
                                !... 
                              !endif
elseif(op=='X=CYZ    ')then
  if(x%rank==0.and.y%rank==1.and.z%rank==1)then
     call AA(x,'X=YZ    ',y%x,z%x)
     call AA(x,'X+=YZ   ',y%y,z%y)
     call AA(x,'X=CX    ',rc=rc)
  else
    write(*,*) 'This case has not been implemented!'
  endif

elseif(trim(op)=='X=X+Y')then
  
  x%r=x%r+y%r  

else
  write(*,*) 'This case has not been implemented!'
endif

return
end subroutine

!DEC$ENDIF

!--------------------------------------------------
    subroutine matvec(x,op,a,y,c,mesh,lego)
!--------------------------------------------------
implicit none
class(array) x,a,y
character(len=*) op
real(dp),optional :: c
type(array),optional :: mesh
logical,optional :: lego

real(dp), pointer, dimension(:) :: p
integer i,length,ndim1

ndim1 = a%ndim1
length = a%ndim1*a%ndim2

do i = 1, ndim1
   p => a%r(i:length:ndim1)   !下标三元组
   x%r(i) = dot_product(p,y%r)
enddo

return
end subroutine

!3 one
!---------------------------
    subroutine one(this)
!---------------------------
implicit none
class(array) this
  this%i = 1
return
end subroutine

!3 done
!---------------------------
    subroutine done(this)
!---------------------------
implicit none
class(array) this
 this%r = 1
return
end subroutine

!3 zero
!---------------------------
    subroutine zero(this)
!---------------------------
implicit none
class(array) this
  this%i = 0
return
end subroutine

!----------------------------
   subroutine dzero(this)
!----------------------------
implicit none
class(array) this
this%r = 0.d0
return
end subroutine

!------------------------------------
    subroutine vzero(this)
!------------------------------------
implicit none
class(array) this
type(vector2d), pointer, dimension(:) :: cmpt

if(debug==.true.)write(*,*) 'In zero---'
!this%cmpt%x = 0.d0
!this%cmpt%y = 0.d0

select type (p=>this%d)
 type is (vector2d)
  cmpt => p
end select

cmpt%x = 0.d0
cmpt%y = 0.d0

return
end subroutine

!---------------------------------
    subroutine one_to_N(this)
!---------------------------------
implicit none
class(array) this
integer i, length

length = this%length()
forall(i=1:length)
  this%i(i) = i
end forall

return
end subroutine

!3 one_to_Nd0
!---------------------------------
    subroutine one_to_Nd0(this)
!---------------------------------
implicit none
class(array) this
integer i

do i = 1, this%length()
   this%r(i) = i*1.d0
enddo

return
end subroutine

!---------------------------------
    subroutine clean(this)
!---------------------------------
implicit none
type(array) this   !! Can not be class!
integer ierr

!write(*,*) 'In cleaning... ', this%nam

!if(this%type==1.and.associated(this%i))then
if(associated(this%i))then
 !write(*,*) 'Cleanned ', this%nam
 deallocate(this%i,stat=ierr); this%i=>null()
endif
if(associated(this%r))then
 !write(*,*) 'Cleanned ', this%nam
 deallocate(this%r,stat=ierr); this%r=>null()
endif

return
end subroutine

!---------------------------------
!    subroutine clean1(this)
!---------------------------------
!implicit none
!type(array) this
!
!return
!end subroutine

!----------------------------------------------
     subroutine construct_equation(this,axb)
!----------------------------------------------
implicit none
class(equation) this, axb

!call move_alloc(axb%a%r,this%a%r)
this%a%r = axb%a%r
this%b%r = axb%b%r
this%x%r => axb%x%r

return
end subroutine

end module


