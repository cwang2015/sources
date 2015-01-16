	      program SPH
!              subroutine beham

c----------------------------------------------------------------------
c     This is a three dimensional SPH code. the followings are the 
c     basic parameters needed in this codeor calculated by this code

c     mass-- mass of particles                                      [in]
c     ntotal-- total particle number ues                            [in]
c     dt--- Time step used in the time integration                  [in]
c     itype-- types of particles                                    [in]
c     x-- coordinates of particles                              [in/out]
c     vx-- velocities of particles                              [in/out]
c     rho-- dnesities of particles                              [in/out]
c     p-- pressure  of particles                                [in/out]
c     u-- internal energy of particles                          [in/out]
c     hsml-- smoothing lengths of particles                     [in/out]
c     c-- sound velocity of particles                              [out]
c     s-- entropy of particles                                     [out]
c     e-- total energy of particles                                [out]

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
      write(*,*) 'dfjoajfoajfiopappppp'
      call time_print
      call time_elapsed(s1)

      call getarg(1,cas_file)
      write(*,*) 'cas_file:', cas_file
      call read_cas

      call open_files      

      call allocate_sph

      if (shocktube)   dt = 0.005
      if (shearcavity) dt = 5.e-5
      !if (waterjet)    dt = 0.0001    !single_phase
      !if (waterjet)    dt = 0.00001   !two_phase

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
