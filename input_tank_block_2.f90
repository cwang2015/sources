subroutine input

!use param
use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k,stat,nvtx,ncell,c,tst
type(block) nozzle, tank
type(material), pointer :: property
double precision soil_surface,water_surface
!real,dimension(:,:)::a(478,2)
logical :: dbg = .false.
double precision element_size, soil_submerged_depth   
real :: p,p1,p2,p3
character(len=128) str
real,dimension(:):: s(5000)
integer,dimension(:) :: b1(5000),b2(5000),b3(5000)
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
call tank%set(xl=1.d0,yl=1.d0,m=100,n=100)



!call tank%set(xl=3.26d0,yl=1.82d0,m=32,n=18)
!this%dx = this%xl/this%m
!this%dy = this%yl/this%n

!call tank%cell_center

open(unit=9,file='D:\LiuSPH\data\stillwater112.neu',status='old')!文件的名字是a.txt
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
allocate(tank%x(nvtx),tank%y(nvtx),tank%zone(nvtx))

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

!do i =1,nvtx
!    write(*,*) V(i)
!enddo
!read(9,*)((a(i,j),j=1,2),i=1,478)!把矩阵a读入了

!do i =1,478
!write(*,*)a(i,1),a(i,2)
!enddo

!k = 0
!do i = 1, npoint
! if(i<=478)then
!      k = k + 1
!      tank%x(k) = a(i,1)
!      tank%y(k) = a(i,2)
! endif
!enddo

! Zoning
tank%zone = 2
do i = 1, nvtx
   if(tank%x(i)<0.01) tank%zone(i) = 1
   if(tank%y(i)<0.01.and.tank%x(i)>0.01.and.tank%x(i)<0.39) tank%zone(i) = 3
   if(tank%x(i)<0.41.and.tank%x(i)>0.39.and.tank%y(i)<0.31) tank%zone(i)=4
   if(tank%x(i)<0.59.and.tank%x(i)>0.41.and.tank%y(i)<0.31.and.tank%y(i)>0.29) tank%zone(i)=5
   if(tank%x(i)<0.61.and.tank%x(i)>0.59.and.tank%y(i)<0.31) tank%zone(i)=6
   if(tank%x(i)>0.61.and.tank%y(i)<0.01)tank%zone(i)=7
   if(tank%x(i)>0.99) tank%zone(i)=8
!   if(tank%zone(i)==1.and.tank%x(i)>0.02.and.tank%x(i)<1.78) tank%zone(i)=4
!   if(tank%zone(i)==1.and.tank%x(i)>1.78) tank%zone(i)=5
!   if(tank%x(i)<0.02.and.tank%y(i)<0.02)  tank%zone(i)=6
!   if(tank%x(i)>1.78.and.tank%y(i)<0.02)  tank%zone(i)=6

   if(tank%zone(i)==2.and.tank%y(i)>0.52)tank%zone(i)= 80
!   if(tank%zone(i)==2.and.tank%y(i)<0.4.and.tank%x(i)<0.6.and.tank%x(i)>0.4)tank%zone(i)= 8
!   if(tank%zone(i)==2.and.tank%x(i)>0.98)tank%zone(i)= 8

enddo

parts%nvtx = nvtx

call parts%take_real_new_mesh(tank,2)

call parts%take_virtual_new_mesh(tank,1)
call parts%sort_new_meshn(tank,1)
call parts%take_virtual_new_mesh(tank,3)
call parts%take_virtual_new_mesh(tank,4)
call parts%sort_new_meshs(tank,4)
call parts%take_virtual_new_mesh(tank,5)
call parts%sort_new_meshl(tank,5)
call parts%take_virtual_new_mesh(tank,6)
call parts%sort_new_meshn(tank,6)
call parts%take_virtual_new_mesh(tank,7)
call parts%sort_new_meshl(tank,7)
call parts%take_virtual_new_mesh(tank,8)
call parts%sort_new_meshs(tank,8)
call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      

   
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

!parts%vol = tank%dx*tank%dy
parts%hsml = 2*tank%dx
parts%dspp = tank%dx

! itype is positive for real particles, negative for virtual particles.
 
call parts%setup_itype

! Set initial conditions of particles
! ...Velocity

parts%vx%x = 0.d0
parts%vx%y = 0.d0

! ...Pressure. You must define the free surface first.
water_surface = 0.52d0
property => parts%material
do i = 1,parts%ntotal+parts%nvirt
   parts%p%r(i) = property%rho0*gravity*(parts%x(2,i)-water_surface)
   if(parts%zone(i) ==1.and.parts%zone(i)>0.5) parts%p%r(i) = 0.d0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)

!do i = 1,parts%ntotal + parts%nvirt 
!    if(parts%zone(i)/=2)then 
!        parts%vol%r(i) = parts%vol%r(i)/2.0d0
!        if(i==1060) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!        if(i==1080) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!        if(i==1127) parts%vol%r(i)=parts%vol%r(i)/2.0d0
!        if(i==1147) parts%vol%r(i)=parts%vol%r(i)/2.0d0        
!        if(i==1099) parts%vol%r(i)=parts%vol%r(i)*3./2.0d0
!        if(i==1108) parts%vol%r(i)=parts%vol%r(i)*3./2.0d0
!    endif
!enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
end subroutine




