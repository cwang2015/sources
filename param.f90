!---------------------------------------------------------
!     Including file for parameters and constants used 
!     in the entire SPH software packages.
!---------------------------------------------------------
      module param
!---------------------------------------------------------
      implicit none

! dim : Dimension of the problem (1, 2 or 3)
integer :: dim = 2

!     maxn    : Maximum number of particles
!     max_interation : Maximum number of interaction pairs
integer :: maxn = 12000, max_interaction = 100 * 12000

!     Parameters for the computational geometry,  
!     x_maxgeom : Upper limit of allowed x-regime 
!     x_mingeom : Lower limit of allowed x-regime 
!     y_maxgeom : Upper limit of allowed y-regime 
!     y_mingeom : Lower limit of allowed y-regime 
!     z_maxgeom : Upper limit of allowed z-regime 
!     z_mingeom : Lower limit of allowed z-regime 
double precision :: x_maxgeom = 10.e0, x_mingeom = -10.e0,  &
                    y_maxgeom = 10.e0, y_mingeom = -10.e0,  &
                    z_maxgeom = 10.e0, z_mingeom = -10.e0
    
!     SPH algorithm for particle approximation (pa_sph)
!     pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
!              2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
integer :: pa_sph = 2 

!     Nearest neighbor particle searching (nnps) method
!     nnps = 1 : Simplest and direct searching
!            2 : Sorting grid linked list
!            3 : Tree algorithm
integer :: nnps = 1 

!     Smoothing length evolution (sle) algorithm
!     sle = 0 : Keep unchanged,
!           1 : h = fac * (m/rho)^(1/dim)
!           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
!           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) ) 

integer :: sle = 0 

!     Smoothing kernel function 
!     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!         = 2, Gauss kernel   (Gingold and Monaghan 1981) 
!         = 3, Quintic kernel (Morris 1997)
!         = 4, Wendland    
integer :: skf = 4 

!     Switches for different senarios

!     summation_density = .TRUE. : Use density summation model in the code, 
!                        .FALSE.: Use continuiity equation
!     average_velocity = .TRUE. : Monaghan treatment on average velocity,
!                       .FALSE.: No average treatment.
!     config_input = .TRUE. : Load initial configuration data,
!                   .FALSE.: Generate initial configuration.
!     virtual_part = .TRUE. : Use vritual particle,
!                   .FALSE.: No use of vritual particle.
!     vp_input = .TRUE. : Load virtual particle information,
!               .FALSE.: Generate virtual particle information.
!     visc = .true. : Consider viscosity,
!           .false.: No viscosity.
!     ex_force =.true. : Consider external force,
!               .false.: No external force.
!     visc_artificial = .true. : Consider artificial viscosity,
!                      .false.: No considering of artificial viscosity.
!     heat_artificial = .true. : Consider artificial heating,
!                      .false.: No considering of artificial heating.
!     self_gravity = .true. : Considering self_gravity,
!                    .false.: No considering of self_gravity
!     nor_density =  .true. : Density normalization by using CSPM,
!                    .false.: No normalization.

integer :: integrate_scheme = 1  ! =1, LF; =2, Verlet
logical :: summation_density  = .false.         
logical :: average_velocity  = .false.         
logical :: config_input  = .false. 
logical :: virtual_part  = .true. 
logical :: vp_input  = .false.  
logical :: visc  = .true.  
logical :: ex_force  = .true.
logical :: visc_artificial  = .false. 
logical :: heat_artificial  = .false. 
!logical :: self_gravity  = .true.      
logical :: nor_density  = .false.              

integer :: soil_pressure = 2  ! =1, eos; =2, mean trace
integer :: stress_integration = 1
integer :: yield_criterion = 2
integer :: plasticity = 3  ! =0 non; =1 Bui, =2 return mapping =3 Lopez
logical :: artificial_density = .true.                  
logical :: soil_artificial_stress = .false.

logical :: volume_fraction = .true.
logical :: water_artificial_volume = .true.
logical :: volume_fraction_renorm = .true.

! 0 ignor; 1 negative pressure to zero; 2 artficial stress
integer :: water_tension_instability = 0

! Symmetry of the problem
! nsym = 0 : no symmetry,
!      = 1 : axis symmetry,
!      = 2 : center symmetry.     
integer :: nsym = 0

! Control parameters for output 
! int_stat = .true. : Print statistics about SPH particle interactions.
!                     including virtual particle information.
! print_step: Print Timestep (On Screen)
! save_step : Save Timestep    (To Disk File)
! moni_particle: The particle number for information monitoring.
logical :: int_stat = .true.
integer :: print_step, save_step, moni_particle = 264
           
! double precision pi,gravity
double precision, parameter :: pi = 3.14159265358979323846 
!parameter ( gravity = -9.8 )
      
! Simulation cases
! shocktube = .true. : carry out shock tube simulation
! shearcavity = .true. : carry out shear cavity simulation
! logical shocktube, shearcavity, waterjet
logical :: shocktube  = .false. 
logical :: shearcavity  = .false. 
logical :: waterjet = .true.

!     Recorde time interval
integer :: save_step_from = 0, save_step_to = 100

end module
