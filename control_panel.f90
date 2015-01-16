!     ****************************
       subroutine control_panel
!     ****************************
!  Version 1.0 June. 01, 2010, Chun WANG (SJTU) 34204490
! ----------------------------------------------------------
!  Function: assignment via internal file
! ----------------------------------------------------------
      use param
      use declarations_sph
!
      implicit none

!      case('PRINT STEP')

          print_step = 1000
          write(*,*) 'Print step = ', print_step 

!      case('SAVE STEP')

          save_step = 1000
          write(*,*) 'Save step = ', save_step 

!      case('RESULTS FILE')
                
          results='res.dat'
          res_soil='res_soil.dat'
          write(*,*) 'Results file = ', trim(results), '  ',   &
                                        trim(res_soil)

!      case('INLET VELOCITY')

          inlet_velocity = -25.
          write(*,*) 'Inlet velocity = ', inlet_velocity

!      case('SINGLE PHASE')
                
          single_phase=.false.
          write(*,*) 'Single phase = ', single_phase
          
!      case('FIRST PHASE')

          parts%imaterial='water' 
          write(*,*) 'First phase = ', parts%imaterial

!      case('SECOND PHASE')

          soil%imaterial='soil'
          write(*,*) 'Second phase = ', soil%imaterial

!      case('DRY SOIL')

          dry_soil=.false.
          write(*,*) 'Dry soil = ', dry_soil

!      case('PARTICLE SIZE')
          
          particle_size='koshi'
          write(*,*) 'Particle size = ', particle_size

!       case('IMMERSE DEPTH')

           wasserjet%immerse = 0.
           write(*,*) 'Immerse depth = ', wasserjet%immerse

!       case('WITH NOZZLE')
                  
            with_nozzle=.false.
            write(*,*) 'With nozzle = ', with_nozzle

!      case('WATER PROPERTY')

          h2o%rho0=1000.; h2o%b=1.e5; h2o%gamma=7.; h2o%c=20.; h2o%viscosity=1.e-3
          write(*,*) 'Water property: rho0, b, gamma, c, viscosity' 
          write(*,*) h2o%rho0, h2o%b, h2o%gamma, h2o%c, h2o%viscosity

!      case('SOIL PROPERTY')

          !sio2%rho0=2700.; sio2%k=0.7e6; sio2%E=0.84e6; sio2%niu=0.3; sio2%c=30
          sio2%rho0=2700.; sio2%k=125.e6; sio2%E=150.e6; sio2%niu=0.3; sio2%c=215
          write(*,*) 'Soil property: rho0, k, Youngs, Poisson, c' 
          write(*,*) sio2%rho0, sio2%k, sio2%E, sio2%niu, sio2%c

!      case('DARCY PARAMETER')

          sio2%porosity=0.3; sio2%permeability=0.05
          write(*,*) 'Darcy parameter: porosity, permeability '
          write(*,*) sio2%porosity, sio2%permeability

!      case('SOIL FAILURE CRITERIA')
       
          sio2%cohesion =0.; SIO2%phi=0.04316
          write(*,*) 'Soil failure criteria: cohesion, phi'
          write(*,*) sio2%cohesion, sio2%phi

!      case('TIME STEP SIZE')

          dt = 0.000005
          write(*,*) 'Time step size = ', dt

      return
      end
