!----------------------------------------
             module m_anyArr
!----------------------------------------
use constants, ex_pi => pi
implicit none

private

type, public :: array

integer :: ndim1

integer :: dim = 2, rank = 0

real(dp), pointer, dimension(:) :: r => null()   

type(array), pointer :: x => null(), y => null(), z => null()
type(array), pointer :: xy => null(), xz => null(), yz => null()

class(*), pointer :: parts => null() 

contains

  procedure :: alloc => xallocate
  procedure :: dealloc => array_deallocate_sub
  procedure :: set => array_set_sub

  procedure :: cpstrc => copy_array_structure

  procedure :: almat => matrix_allocate_sub
  procedure :: length => array_get_length_fn 



  FINAL :: clean

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



end module


