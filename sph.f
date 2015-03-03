      program SPH
!     subroutine beham

c----------------------------------------------------------------------
c     This is a three dimensional SPH code. the followings are the 
c     basic parameters needed in this codeor calculated by this code
c----------------------------------------------------------------------
      use constants
      use param
      use declarations
      use declarations_sph
      implicit none     
c
      integer d, m, i, yesorno
      double precision s1, s2
c      character(len=32) casefileeee
c
      call test
      stop 'test'
      call time_print
      call time_elapsed(s1)

      call getarg(1,cas_file)
      write(*,*) 'cas_file:', cas_file
      call read_cas

      call open_files      

      call allocate_sph

      call input
      call output

      call parts%minimum_time_step

 1    write(*,*)'  ***************************************************' 
      write(*,*)'          Please input the maximal time steps '
      write(*,*)'  ***************************************************'
      read(*,*) maxtimestep
      !maxtimestep = 400000      

      call time_integration

      call time_elapsed(s2)

      write (*,*)'        Elapsed CPU time = ', s2-s1
      s1 = s2

      write(*,*)'  ***************************************************'
      write(*,*) 'Are you going to run more time steps ? (0=No, 1=yes)'
      write(*,*)'  ***************************************************'
      read (*,*) yesorno     
      if (yesorno.ne.0) go to 1

      call close_files 
  
      call time_print
      call time_elapsed(s2)

      write (*,*)'        Elapsed CPU time = ', s2-s1
                           
      end
