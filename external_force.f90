!--------------------------------------------------------------------------
      subroutine repulsive_force(parts)
!--------------------------------------------------------------------------
! Subroutine to calculate the external forces, e.g. gravitational forces.      
! The forces from the interactions with boundary virtual particles 
! are also calculated here as external forces.
!--------------------------------------------------------------------------
!use param
use m_particles
implicit none

type(particles) parts

integer ntotal, niac
integer, pointer, dimension(:) :: pair_i, pair_j, itype
double precision, pointer, dimension(:) :: mass, hsml, x(:,:)
double precision, pointer, dimension(:,:) :: dvxdt

type(numerical), pointer :: numeric      
double precision dx(3), rr, f, rr0, dd, p1, p2     
integer i, j, k, d, dim

ntotal = parts%ntotal + parts%nvirt
niac   = parts%niac

pair_i   => parts%pair_i
pair_j   => parts%pair_j
itype    => parts%itype
mass     => parts%mass
x        => parts%x
hsml     => parts%hsml
dvxdt    => parts%dvx
numeric  => parts%numeric
dim      =  parts%dim
           
! Boundary particle force and penalty anti-penetration force. 

!rr0 =0.02  !1.25e-5
rr0 = numeric%rr0
!rr0 = hsml(1)
dd = numeric%dd !1.e-5 !1.e-2
p1 = numeric%p1
p2 = numeric%p2
      
do k=1,niac
   i = pair_i(k)
   j = pair_j(k)  
   if(itype(i).gt.0.and.itype(j).lt.0)then  
      rr = 0.      
      do d=1,dim
         dx(d) =  x(d,i) -  x(d,j)
         rr = rr + dx(d)*dx(d)
      enddo  
      rr = sqrt(rr)
      if(rr.lt.rr0)then
         f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
         do d = 1, dim
            dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
         enddo
      endif
   endif        
enddo   
       
return
end subroutine         
