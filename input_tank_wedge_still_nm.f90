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
real,dimension(:):: s(10000)
integer,dimension(:) :: b1(20000),b2(20000),b3(20000)
stat = 0
! Set nozzle and tank geometry parameters


!call tank%set(xl=3.34d0,yl=1.86d0,m=334,n=186)
!call tank%set(xl=3.34d0,yl=1.86d0,m=167,n=93)
!call tank%set(xl=3.32d0,yl=1.86d0,m=166,n=93)
!call tank%set(xl=3.4d0,yl=1.9d0,m=34,n=19)
!call tank%set(xl=3.32d0,yl=1.86d0,m=664,n=372)
!call tank%set(xl=3.74d0,yl=2.06d0,m=374,n=206)
!call tank%set(xl=0.44d0,yl=0.22d0,m=352,n=176)

!1call tank%set(xl=2.04d0,yl=1.52d0,m=102,n=76)
!2call tank%set(xl=2.02d0,yl=1.51d0,m=202,n=151)
call tank%set(xl = 2.44d0,yl = 1.52d0,m=122,n=76)



!call tank%set(xl=3.26d0,yl=1.82d0,m=32,n=18)

open(unit=9,file='D:\LiuSPH\data\circle.neu',status='old')!文件的名字是a.txt
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

! Zoning
tank%zone = 2
do i = 1, nvtx
   if(tank%x(i)<0.01) tank%zone(i) = 1
!   if(tank%zone(i)==1.and.tank%x(i)<0.02.and.tank%y(i)>0.02) tank%zone(i)=3
!   if(tank%zone(i)==1.and.tank%x(i)>0.02.and.tank%x(i)<2.42) tank%zone(i)=4
   if(tank%x(i)>1.99) tank%zone(i)=8
   if(tank%y(i)<0.01)  tank%zone(i)=5
!   if(tank%x(i)>2.42.and.tank%y(i)<0.02)  tank%zone(i)=6
!   if(tank%zone(i)==2.and.tank%y(i)>0.42)tank%zone(i)= 8

   if(tank%x(i)>0.80.and.tank%x(i)<0.803.and.tank%y(i)<0.03)tank%zone(i)= 5
   if(tank%x(i)>0.805.and.tank%x(i)<0.807.and.tank%y(i)<0.045)tank%zone(i)= 5
   if(tank%x(i)>0.811.and.tank%x(i)<0.813.and.tank%y(i)<0.06)tank%zone(i)= 5
   if(tank%x(i)>0.82.and.tank%x(i)<0.823.and.tank%y(i)<0.08)tank%zone(i)= 5
   if(tank%x(i)>0.829.and.tank%x(i)<0.831.and.tank%y(i)<0.1)tank%zone(i)= 5
   if(tank%x(i)>0.84.and.tank%x(i)<0.845.and.tank%y(i)<0.12)tank%zone(i)= 5
   if(tank%x(i)>0.85.and.tank%x(i)<0.86.and.tank%y(i)<0.13)tank%zone(i)= 5
   if(tank%x(i)>0.865.and.tank%x(i)<0.87.and.tank%y(i)<0.15)tank%zone(i)= 5
   if(tank%x(i)>0.88.and.tank%x(i)<0.89.and.tank%y(i)<0.16)tank%zone(i)= 5
   if(tank%x(i)>0.9.and.tank%x(i)<0.905.and.tank%y(i)<0.17)tank%zone(i)= 5
   if(tank%x(i)>0.92.and.tank%x(i)<0.93.and.tank%y(i)<0.18)tank%zone(i)= 5
   if(tank%x(i)>0.94.and.tank%x(i)<0.95.and.tank%y(i)<0.19)tank%zone(i)= 5
   if(tank%x(i)>0.95.and.tank%x(i)<0.97.and.tank%y(i)<0.19)tank%zone(i)= 5
   if(tank%x(i)>0.97.and.tank%x(i)<0.98.and.tank%y(i)<0.19)tank%zone(i)= 5
   if(tank%x(i)>0.98.and.tank%x(i)<1.01.and.tank%y(i)<0.191)tank%zone(i)= 5
   if(tank%x(i)>1.02.and.tank%x(i)<1.03.and.tank%y(i)<0.190)tank%zone(i)= 5
   if(tank%x(i)>1.02.and.tank%x(i)<1.03.and.tank%y(i)<0.190)tank%zone(i)= 5
   if(tank%x(i)>1.04.and.tank%x(i)<1.05.and.tank%y(i)<0.187)tank%zone(i)= 5
   if(tank%x(i)>1.05.and.tank%x(i)<1.06.and.tank%y(i)<0.182)tank%zone(i)= 5
   if(tank%x(i)>1.07.and.tank%x(i)<1.08.and.tank%y(i)<0.175)tank%zone(i)= 5
   if(tank%x(i)>1.09.and.tank%x(i)<1.1.and.tank%y(i)<0.17)tank%zone(i)= 5
   if(tank%x(i)>1.11.and.tank%x(i)<1.12.and.tank%y(i)<0.16)tank%zone(i)= 5
   if(tank%x(i)>1.13.and.tank%x(i)<1.14.and.tank%y(i)<0.15)tank%zone(i)= 5
   if(tank%x(i)>1.14.and.tank%x(i)<1.15.and.tank%y(i)<0.13)tank%zone(i)= 5
   if(tank%x(i)>1.15.and.tank%x(i)<1.16.and.tank%y(i)<0.12)tank%zone(i)= 5
   if(tank%x(i)>1.16.and.tank%x(i)<1.17.and.tank%y(i)<0.1)tank%zone(i)= 5
   if(tank%x(i)>1.17.and.tank%x(i)<1.18.and.tank%y(i)<0.08)tank%zone(i)= 5
   if(tank%x(i)>1.18.and.tank%x(i)<1.19.and.tank%y(i)<0.06)tank%zone(i)= 5
   if(tank%x(i)>1.19.and.tank%x(i)<1.2.and.tank%y(i)<0.05)tank%zone(i)= 5
   if(tank%x(i)>1.19.and.tank%x(i)<1.2.and.tank%y(i)<0.03)tank%zone(i)= 5
!   if(tank%x(i)>0.86.and.tank%y(i)<tank%x(i) - 0.85.and.tank%x(i)<1.09)tank%zone(i)= 5
!   if(tank%x(i)>1.09.and.tank%y(i)<1.33 - tank%x(i).and.tank%x(i)<1.31)tank%zone(i)= 5
!   if(tank%zone(i)==2.and.tank%x(i)>2.22) tank%zone(i)= 8
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

tank%dx = 0.02  !special
parts%hsml = 2*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.5d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i)==1.and.parts%x(2,i)>0.51)parts%p%r(i)=0.0
   if(parts%zone(i)==8.and.parts%x(2,i)>0.51)parts%p%r(i)=0.0
!   if(parts%zone(i)/=2.and.parts%x(1,i)>0.11)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
do i = 1,parts%ntotal + parts%nvirt 
   if(parts%x(2,i)>0.49.and.parts%x(2,i)<0.59)then
        parts%vol%r(i) = parts%vol%r(i)*2.0d0
   endif
enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
end subroutine




