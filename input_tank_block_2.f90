subroutine input

!use param
use declarations_sph
implicit none     
      
integer  ntotal, bntotal, npoint
integer i, j, d, k
type(block) nozzle, tank
type(material), pointer :: property
double precision soil_surface,water_surface
real,dimension(:,:)::a(478,2)
logical :: dbg = .false.
double precision element_size, soil_submerged_depth    

! Set nozzle and tank geometry parameters

!call tank%set(xl=3.34d0,yl=1.86d0,m=334,n=186)
!call tank%set(xl=3.34d0,yl=1.86d0,m=167,n=93)
!call tank%set(xl=3.32d0,yl=1.86d0,m=166,n=93)
!call tank%set(xl=3.4d0,yl=1.9d0,m=34,n=19)
!call tank%set(xl=3.32d0,yl=1.86d0,m=664,n=372)
!call tank%set(xl=3.74d0,yl=2.06d0,m=374,n=206)
!call tank%set(xl=0.44d0,yl=0.22d0,m=352,n=176)

!1call tank%set(xl=2.04d0,yl=1.52d0,m=102,n=76)
call tank%set(xl=1.8d0,yl=1.d0,m=90,n=50)



!call tank%set(xl=3.26d0,yl=1.82d0,m=32,n=18)

npoint = tank%m*tank%n
allocate(tank%x(npoint),tank%y(npoint),tank%zone(npoint))
!1call tank%cell_center


open(unit=9,file='D:\LiuSPH\data\stillwater2.txt',status='old')!文件的名字是a.txt
!read(8,*)m,n
!allocate(a(m,n))
read(9,*)((a(i,j),j=1,2),i=1,478)!把矩阵a读入了

do i =1,478
write(*,*)a(i,1),a(i,2)
enddo

k = 0
do i = 1, npoint
 if(i<=478)then
      k = k + 1
      tank%x(k) = a(i,1)
      tank%y(k) = a(i,2)
 endif
enddo

!      write(*,*) 'x=', tank%x
!      write(*,*) 'y=', tank%y
!      just for test 1 

!tank%zone = 2
!do i = 1, tank%m*tank%n
!   if(tank%x(i)<0.1.or.tank%x(i)>3.1.or.tank%y(i)<0.1) tank%zone(i) = 1
!   if(tank%zone(i)==1.and.tank%x(i)<0.1.and.tank%y(i)>0.1) tank%zone(i)=3
!   if(tank%zone(i)==1.and.tank%x(i)>0.01.and.tank%x(i)<3.1) tank%zone(i)=4
!   if(tank%zone(i)==1.and.tank%x(i)>3.1) tank%zone(i)=5
!   if(tank%zone(i)==2.and.tank%x(i)>1.3)tank%zone(i)= 7
!   if(tank%zone(i)==2.and.tank%y(i)>0.7)tank%zone(i)= 7

!enddo


! Zoning
tank%zone = 2
do i = 1, 478
   if(tank%x(i)<0.01.or.tank%x(i)>0.99.or.tank%y(i)<0.01) tank%zone(i) = 1
   if(tank%x(i)<0.41.and.tank%x(i)>0.39.and.tank%y(i)<0.31) tank%zone(i)=1
   if(tank%x(i)<0.59.and.tank%x(i)>0.41.and.tank%y(i)<0.31.and.tank%y(i)>0.29) tank%zone(i)=1
!   if(tank%zone(i)==1.and.tank%x(i)>0.02.and.tank%x(i)<1.78) tank%zone(i)=4
!   if(tank%zone(i)==1.and.tank%x(i)>1.78) tank%zone(i)=5
!   if(tank%x(i)<0.02.and.tank%y(i)<0.02)  tank%zone(i)=6
!   if(tank%x(i)>1.78.and.tank%y(i)<0.02)  tank%zone(i)=6

   if(tank%zone(i)==2.and.tank%y(i)>0.52)tank%zone(i)= 8
!   if(tank%zone(i)==2.and.tank%y(i)<0.4.and.tank%x(i)<0.6.and.tank%x(i)>0.4)tank%zone(i)= 8
!   if(tank%zone(i)==2.and.tank%x(i)>0.98)tank%zone(i)= 8
   
   
enddo

call parts%take_real(tank,2)
call parts%take_virtual(tank,1)

!call parts%take_boundary_for_tank_block2(tank)

call parts%setup_ndim1

!      write(*,*) parts%ntotal, parts%nvirt
      
! Basic settings for particles (vol,hsml,itype)
! vol means the volume of a cell. We calculate the mass of each particle according to mass = rho*vol

parts%vol = tank%dx*tank%dy
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
   if(parts%zone(i)==3.and.parts%x(2,i)>0.52)parts%p%r(i)=0.0
   if(parts%zone(i)==5.and.parts%x(2,i)>0.52)parts%p%r(i)=0.0
enddo

parts%c = property%c

! Calulate desity from pressure according to Equation of State (EoS): p = B(rho/rho0-1)^gamma
  
call initial_density(parts)
      
!do i = 1,parts%ntotal +parts%nvirt 
!    if(tank%zone(i)==1) parts%vol%r(i) =  tank%dx*tank%dy/2.0d0
!    if(tank%zone(i)==3) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!    if(tank%zone(i)==4) parts%vol%r(i) =  tank%dx*tank%dy/4.0d0
!enddo

do i = 1,parts%ntotal + parts%nvirt 
    if(parts%zone(i)/=2)then 
        parts%vol%r(i) = parts%vol%r(i)/2.0d0
        if(i==1060) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(i==1080) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(i==1127) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(i==1147) parts%vol%r(i)=parts%vol%r(i)/2.0d0        
!        if(i==1149) parts%vol%r(i)=parts%vol%r(i)/2.0d0
        if(i==1099) parts%vol%r(i)=parts%vol%r(i)*3./2.0d0
        if(i==1108) parts%vol%r(i)=parts%vol%r(i)*3./2.0d0
!        if(parts%x(1,i)>0.98.and.parts%x(2,i)<0.02) parts%vol%r(i)=parts%vol%r(i)/2.0d0
    endif
enddo

parts%mass = parts%vol * parts%rho

call get_norm(parts)

return
end subroutine




