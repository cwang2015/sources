!--------------------
subroutine test1(f)
!--------------------
use m_particles
implicit none

type(array) f
type(array), allocatable :: g

allocate(g)
g%r => f%r

return
end subroutine

 

