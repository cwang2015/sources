module m_sph_fo
contains

! Calculate partial derivatives of a field
!-------------------------------------------
          function df(f,x,parts)
!-------------------------------------------
use m_particles
implicit none

real(dp) f(:)
character(len=1) x
type(particles) parts
real(dp), allocatable :: df(:)

integer ntotal,niac
integer, pointer, dimension(:) :: pair_i, pair_j 
double precision, pointer, dimension(:) :: mass, rho, dwdx
real(dp) fwx
integer i, j, k

niac = parts%niac

pair_i   => parts%pair_i
pair_j   => parts%pair_j
mass     => parts%mass
rho      => parts%rho

allocate(df(size(f))); df = 0.

if(x=='x')dwdx=>parts%dwdx(1,:)
if(x=='y')dwdx=>parts%dwdx(2,:)

do k=1,niac
   i = pair_i(k)
   j = pair_j(k)
   fwx = (f(i)+f(j))*dwdx(k)
   df(i) = df(i) + mass(j)/rho(j)*fwx
   df(j) = df(j) - mass(i)/rho(i)*fwx
enddo

end function

end module              

 