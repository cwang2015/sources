subroutine input

!use param
use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k,stat,nvtx,ncell,c,tst
type(block) nozzle, tank
type(material), pointer :: property
double precision soil_surface,water_surface
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    
real :: p,p1,p2,p3
character(len=128) str
real,dimension(:):: s(50000)
integer,dimension(:) :: b1(50000),b2(50000),b3(50000)
stat = 0
! Set nozzle and tank geometry parameters

call tank%set(xl = 4.88d0,yl = 3.04d0,m=488,n=304)
open(unit=9,file='D:\LiuSPH\data\dambreak_wedge2.neu',status='old')!文件的名字是a.txt
parts%vol = 0.d0

do while(stat == 0)
 read(9,'(A)',iostat=stat) str
        if(index(str,'NUMMP')/=0.or.index(str,'NELEM')/=0)then
            read(9,*) nvtx, ncell
!write(*,*) 'sd',nvtx,ncell
            exit
        endif
enddo
!write(*,*) 'sd',nvtx,ncell

npoint = nvtx
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))

do while(stat==0)
     read(9,'(A)',iostat=stat) str
        if(index(str,'NODAL COORDINATES')/=0)then
            do i = 1, nvtx
                read(9,*) k, tank%x(k), tank%y(k)
            enddo
        elseif(index(str,'ELEMENTS/CELLS')/=0)then
            do i = 1, ncell
               read(9,*) c,j,j,b1(c),b2(c),b3(c)
            enddo
        endif
enddo

do i =1,ncell
    p1=sqrt((tank%x(b1(i)) - tank%x(b2(i)))**2 + (tank%y(b1(i)) - tank%y(b2(i)))**2)
    p2=sqrt((tank%x(b1(i)) - tank%x(b3(i)))**2 + (tank%y(b1(i)) - tank%y(b3(i)))**2)
    p3=sqrt((tank%x(b2(i)) - tank%x(b3(i)))**2 + (tank%y(b2(i)) - tank%y(b3(i)))**2)
    p = (p1+p2+p3)/2
    S(i) = sqrt(p*(p-p1)*(p-p2)*(p-p3))
    parts%vol_nm%r(b1(i)) = parts%vol_nm%r(b1(i)) + S(i)/3
    parts%vol_nm%r(b2(i)) = parts%vol_nm%r(b2(i)) + S(i)/3
    parts%vol_nm%r(b3(i)) = parts%vol_nm%r(b3(i)) + S(i)/3
enddo

write(*,*)'ssssssss'

! Zoning
tank%zone = 2
do i = 1, nvtx
   if(tank%x(i)<0.001) tank%zone(i) = 1
   if(tank%x(i)>2.20) tank%zone(i)=8
   if(tank%y(i)<0.001)  tank%zone(i)=5
   if(tank%x(i)>0.85.and.tank%y(i)<tank%x(i) - 0.849.and.tank%x(i)<1.03)tank%zone(i)= 5
   if(tank%x(i)>1.03.and.tank%y(i)<1.21 - tank%x(i).and.tank%x(i)<1.21)tank%zone(i)= 5 ! 
   if(tank%zone(i) == 2.and.tank%x(i)>0.5005) tank%zone(i) = 7
enddo

parts%nvtx = nvtx

call parts%take_real_new_mesh(tank,2)

call parts%take_virtual_new_mesh(tank,1)
call parts%sort_new_meshn(tank,1)
call parts%take_virtual_new_mesh(tank,5)
call parts%sort_new_meshl(tank,5)
call parts%take_virtual_new_mesh(tank,8)
call parts%sort_new_meshs(tank,8)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

tank%dx = 0.01  !special
parts%hsml = 2*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 1.01d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%x(2,i)>1.01)parts%p%r(i)=0.0
   if(parts%x(1,i)>0.501)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
do i = 1,parts%ntotal + parts%nvirt 
   if(parts%x(2,i)>0.496.and.parts%x(2,i)<0.504)then
        parts%vol%r(i) = parts%vol%r(i)*2.0d0
   endif
enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
end subroutine




