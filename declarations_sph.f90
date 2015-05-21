!-----------------------------
    module declarations_sph
!-----------------------------
use m_particles
use param
implicit none

! Control panel
logical :: self_gravity = .true.
!logical :: volume_fraction = .true.

! Paramters
double precision :: gravity = -9.8


! Logic unit of files

integer f_xv,f_state,f_other
integer xv_vp,state_vp,other_vp

! Files
character(len=32) results, res_water, res_soil, res_other

! Filenames
character(len=32) output_file_name(3)

! Boundary conditions
double precision inlet_velocity

! Geometry object
!type(geo) wasserjet
logical with_nozzle

! Particle size: large, small?
character(len=32) particle_size
character(len=32) operating_mode

logical single_phase
logical dry_soil

! Particle type : water, soil or mixture?
character(len=32) particle_type

type(material), target :: H2O, SiO2
type(parameters), target :: params
type(numerical), target :: numeric

! Particles 
type(particles), target :: parts, soil

! Working arrays
double precision, pointer, dimension(:) :: txx=>null(), tyy=>null(), txy=>null()

!Numerical paramters
double precision :: dt, time = 0.d0
integer :: maxtimestep = 0 , itimestep = 0

contains

!---------------------------
    subroutine open_files
!---------------------------
implicit none
logical :: dbg = .false.

if(dbg) write(*,*) 'In open_files...'

call get_unit(f_xv)
!open(f_xv,file="../f_xv.dat")
open(f_xv,file=trim(results))

call get_unit(f_state)
!open(f_state,file="../f_state.dat")
open(f_state,file=trim(res_soil))

call get_unit(f_other)
!open(f_other,file="../f_other.dat") 
open(f_other,file=trim(res_other)) 

!call get_unit(xv_vp)
!open(xv_vp,file="../xv_vp.dat")
!call get_unit(state_vp)
!open(state_vp,file="../state_vp.dat")
!call get_unit(other_vp)
!open(other_vp,file="../other_vp.dat")

return
end subroutine

!---------------------------
    subroutine close_files
!---------------------------
implicit none

close(f_xv)
close(f_state)
close(f_other)

return
end subroutine

!----------------------------------------------------------------------
      subroutine output
!----------------------------------------------------------------------           
      implicit none     
      
      integer, pointer, dimension(:) :: itype
      double precision, pointer, dimension(:,:) :: x
      double precision, pointer, dimension(:) :: mass,hsml,p

      integer ntotal, ntotal2
      integer i, d, npart, i1, i2, f1, f2, f3     

      ntotal =  parts%ntotal+parts%nvirt
      ntotal2=  soil%ntotal+soil%nvirt
     
      itype  => parts%itype
      x      => parts%x
      mass   => parts%mass
      hsml   => parts%hsml
      p      => parts%p

              if(trim(parts%imaterial)=='water')then

      write(f_xv,*) 'VARIABLES="X","Y","Pressure","VoF",'     
      write(f_xv,*) '"vx","vy","rho","zone","VoF2","mass" '
      write(f_xv,*) 'ZONE I=', ntotal, ' F=BLOCK'
      write(f_xv,*) x(1,1:ntotal) !, soil%x(1,1:ntotal2)
      write(f_xv,*) x(2,1:ntotal) !, soil%x(2,1:ntotal2)
!      write(f_xv,*) vx(1,1:ntotal)
!      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
!      write(f_xv,*) mass(1:ntotal)
      write(f_xv,*) p(1:ntotal) !, soil%p(1:ntotal2)
      write(f_xv,*) parts%vof(1:ntotal)
!      write(f_xv,*) u(1:ntotal)
!      write(f_xv,*) itype(1:ntotal)
!      write(f_xv,*) hsml(1:ntotal)                                        
      write(f_xv,*) parts%vx(1,1:ntotal)
      write(f_xv,*) parts%vx(2,1:ntotal)
      write(f_xv,*) parts%rho(1:ntotal)
      write(f_xv,*) parts%zone(1:ntotal)
      write(f_xv,*) parts%vof2(1:ntotal)
      write(f_xv,*) parts%mass(1:ntotal)

             elseif(trim(parts%imaterial)=='soil')then

      write(f_xv,*) 'VARIABLES="X","Y","Pressure", "VoF", "ep",' 
      write(f_xv,*) '"sxy", "sxx","syy","vx","vy", "rho", "mass", '
      write(f_xv,*) '"sigma_yy", "zone" '
      write(f_xv,*) 'ZONE I=', ntotal, ' F=BLOCK'
 
      write(f_xv,*)  parts%x(1,1:ntotal)
      write(f_xv,*)  parts%x(2,1:ntotal)
!      write(f_xv,*) vx(1,1:ntotal)
!      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
!      write(f_xv,*) mass(1:ntotal)
      write(f_xv,*)  parts%p(1:ntotal)
!      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_xv,*)  parts%vof(1:ntotal)
      write(f_xv,*)  parts%epsilon_p(1:ntotal)
      write(f_xv,*)  parts%sxy(1:ntotal)
      write(f_xv,*)  parts%sxx(1:ntotal)
      write(f_xv,*)  parts%syy(1:ntotal)
      write(f_xv,*)  parts%vx(1,1:ntotal)
      write(f_xv,*)  parts%vx(2,1:ntotal)
      write(f_xv,*)  parts%rho(1:ntotal)
      write(f_xv,*)  parts%mass(1:ntotal)
      write(f_xv,*)  -parts%p(1:ntotal) + parts%syy(1:ntotal)
      write(f_xv,*)  parts%zone(1:ntotal)
      write(f_other,*) time, -parts%p(395)+parts%syy(395)

             endif

      !return

      write(f_state,*) 'VARIABLES="X","Y","Pressure", "VoF", "ep", '
      write(f_state,*) '"sxy", "sxx","syy","vx","vy", "rho","mass" '
      write(f_state,*) '"sigma_yy", "zone" '
      write(f_state,*) 'ZONE I=', ntotal2, ' F=BLOCK'
 
      write(f_state,*)  soil%x(1,1:ntotal2)
      write(f_state,*)  soil%x(2,1:ntotal2)
!      write(f_xv,*) vx(1,1:ntotal)
!      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
!      write(f_xv,*) mass(1:ntotal)
      write(f_state,*)  soil%p(1:ntotal2)
!      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_state,*)  soil%vof(1:ntotal2)
      write(f_state,*)  soil%epsilon_p(1:ntotal2)
      write(f_state,*)  soil%sxy(1:ntotal2)
      write(f_state,*)  soil%sxx(1:ntotal2)
      write(f_state,*)  soil%syy(1:ntotal2)
      write(f_state,*)  soil%vx(1,1:ntotal2)
      write(f_state,*)  soil%vx(2,1:ntotal2)
      write(f_state,*)  soil%rho(1:ntotal2)
      write(f_state,*)  soil%mass(1:ntotal2)
      write(f_state,*)  -soil%p(1:ntotal2) + soil%syy(1:ntotal2)
      write(f_state,*)  soil%zone(1:ntotal2)
      !write(f_other,*) time, -soil%p(420)+soil%syy(420)
      !write(f_other,*) time, -soil%p(395)+soil%syy(395)
      !write(f_other,*) time, -parts%p(420)+parts%syy(420)

      return
      end subroutine

!--------------------------
   subroutine time_print
!--------------------------

    implicit none

   ! . local scalars.
   character ( len =  8 ) :: datstr
   character ( len = 10 ) :: timstr

   ! . Get the current date and time.
   call date_and_time ( datstr, timstr )

   ! . Write out the date and time.
   write ( * , "(/A)"  ) "                  Date = " // datstr(7:8) // "/" // &
                                          datstr(5:6) // "/" // &
                                          datstr(1:4)
   write ( * , "(A)"   ) "                  Time = " // timstr(1:2) // ":" // &
                                          timstr(3:4) // ":" // &
                                          timstr(5:10)
   write ( * , *)

   end subroutine

!------------------------------
   subroutine time_elapsed(s)
!------------------------------
    use ifport
    implicit none

    real(8) :: s

   s = rtc()

   end subroutine

end module
