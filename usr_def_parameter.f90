!     ****************************
      subroutine usr_def_parameter      &
!     ****************************
          (pname, pvalu, ncoma)
! ----------------------------------------------------------
!  Version 1.0 June. 01, 2010, Chun WANG (SJTU) 34204490
! ----------------------------------------------------------
!  Function: assignment via internal file
! ----------------------------------------------------------
      use constants
      use declarations
      use param
      use declarations_sph
!
      implicit none
!
      character(len=*)   :: pname, pvalu
      integer(kind=ik_3)  :: ncoma, ndim, i
!
      pname = adjustl(pname)
      pvalu = adjustl(pvalu)
!     
 
      call make_upper(pname)
!            
      select case (trim(adjustl(pname)))
!
! This part can be modified by the users!
!
! For example:
!
      case('PRINT STEP')

          read(pvalu,*) print_step
          write(*,*) 'Print step = ', print_step

      case('SAVE STEP')

          read(pvalu,*) save_step
          write(*,*) 'Save step = ', save_step

      case('PARTICLE PAIRS FIND METHOD')

         read(pvalu,*) numeric%nnps
         write(*,*) 'Particle pairs find method = ', numeric%nnps 

      case('RESULTS FILE')
                
                  !ndim = 1
                  !if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  read(pvalu,*) results, res_soil, res_other
            write(*,*) 'Results file = ', trim(results), ', ',  &
                      trim(res_soil), ', ', trim(res_other)

      case('INLET VELOCITY')

         read(pvalu,*) inlet_velocity
         write(*,*) 'Inlet velocity = ', inlet_velocity

      case('SINGLE PHASE')
                
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  read(pvalu,*) single_phase
                  write(*,*) 'Single phase = ', single_phase

      case('FIRST PHASE')

                  read(pvalu,*) parts%imaterial
                  write(*,*) 'First phase = ', parts%imaterial

      case('SECOND PHASE')

                  read(pvalu,*) soil%imaterial
                  write(*,*) 'Second phase = ', soil%imaterial

      case('DRY SOIL')

                  read(pvalu,*) dry_soil
                  write(*,*) 'Dry soil = ', dry_soil

      case('PARTICLE SIZE')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  read(pvalu,*) particle_size
                  write(*,*) 'Particle size = ', particle_size

       case('OPERATING MODE')

                  read(pvalu,*) operating_mode
                  write(*,*) 'Operating mode = ', operating_mode

       case('IMMERSE DEPTH')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  !read(pvalu,*) wasserjet%immerse
                  !write(*,*) 'Immerse depth = ', wasserjet%immerse

       case('WITH NOZZLE')
                  
                  read(pvalu,*) with_nozzle
                  write(*,*) 'With nozzle = ', with_nozzle

      case('WATER PROPERTY')

          read(pvalu,*) h2o%rho0,h2o%b,h2o%gamma, h2o%c, h2o%viscosity
          write(*,*) 'Water property: rho0, b, gamma, c, viscosity' 
          write(*,*) h2o%rho0, h2o%b, h2o%gamma, h2o%c, h2o%viscosity

      case('SOIL PROPERTY')

          read(pvalu,*) sio2%rho0, sio2%k, sio2%E, sio2%niu, sio2%c
          write(*,*) 'Soil property: rho0, k, Youngs, Poisson, c' 
          write(*,*) sio2%rho0, sio2%k, sio2%E, sio2%niu, sio2%c

      case('DARCY PARAMETER')

          read(pvalu,*) sio2%porosity, sio2%permeability
          write(*,*) 'Darcy parameter: porosity, permeability '
          write(*,*) sio2%porosity, sio2%permeability

      case('INTEGRATION SCHEME')
          read(pvalu,*) integrate_scheme
          write(*,*) 'Integration scheme: ', integrate_scheme 

      case('SOIL FAILURE CRITERIA')
       
          read(pvalu,*) sio2%cohesion, SIO2%phi
          write(*,*) 'Soil failure criteria: cohesion, phi'
          write(*,*) sio2%cohesion, sio2%phi

      case('DYNAMIC VISCOSITY')
  
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) miu
!                  write(*,*) 'Dynamic viscosity =', miu

      case('MAXIMUM ITERATION NUMBER')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) cfg%ITMAX
!                  write(*,*) 'Maximum iteration number = ', cfg%ITMAX

      case('DENSITY')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) rho
!                  write(*,*) 'Density = ', rho
                
      case('TOLERANCE')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) cfg%tol
!                  write(*,*) 'Tolerance = ', cfg%tol

      case('RELAXATION FACTOR FOR VELOCITY')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) rfu
!                  write(*,*) 'Relaxation factor for velocity = ', rfu

      case('RELAXATION FACTOR FOR PRESSURE')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) rfp
!                  write(*,*) 'Relaxation factor for pressure = ',rfp

      case('STEADY')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) steady
!                  write(*,*) 'Steady = ',steady
      case('TIME STEP SIZE')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  read(pvalu,*) dt
                  write(*,*) 'Time step size = ', dt
      case('MAXIMUM TIME STEP')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) max_time_step
!                  write(*,*) 'Maximum time step = ', max_time_step


      case('BUOYANCY')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) buoyancy
!                  write(*,*) 'Buoyancy = ', buoyancy

      case('RELAXATION FACTOR FOR TEMPERATURE')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) rft
!                  write(*,*) 'Relaxation factor for temperature = ', rft

      case('PRESSURE REFERENCE POINT')
                  
!                  read(pvalu,*) ref_pnt%x, ref_pnt%y 
!                  write(*,*) 'Pressure reference point: ', ref_pnt%x, 
!     *                                                     ref_pnt%y
      case('SOLVER')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,'(A)') cfg%solver
!                  write(*,*) 'Solver = ', cfg%solver

      case('PRECONDITION')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,'(A)') cfg%precon
!                  write(*,*) 'Precondition = ', cfg%precon

      case('CONVECTION SCHEME')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,'(A)') convection_scheme
!                  write(*,*) 'Convection scheme = ', convection_scheme

       case('MATRIX STORAGE')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,'(A)') matrix_storage
!                  write(*,*) 'Matrix storage = ', matrix_storage

      case DEFAULT

                 write(*,*) "Error: No such keyword -- ", pname

      end select
!          

      return
      end


