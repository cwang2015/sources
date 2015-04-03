c---------------------------------------------------------
c     Including file for parameters and constants used 
c     in the entire SPH software packages.
c---------------------------------------------------------
      module param
c---------------------------------------------------------
      implicit none

c     dim : Dimension of the problem (1, 2 or 3)
      integer dim
      parameter ( dim = 2)

c     maxn    : Maximum number of particles
c     max_interation : Maximum number of interaction pairs
      integer maxn,max_interaction
      parameter ( maxn    = 12000    ,                  ! 12000
     &            max_interaction = 100 * maxn )

c     Parameters for the computational geometry,  
c     x_maxgeom : Upper limit of allowed x-regime 
c     x_mingeom : Lower limit of allowed x-regime 
c     y_maxgeom : Upper limit of allowed y-regime 
c     y_mingeom : Lower limit of allowed y-regime 
c     z_maxgeom : Upper limit of allowed z-regime 
c     z_mingeom : Lower limit of allowed z-regime 
      double precision x_maxgeom,x_mingeom,y_maxgeom,
     &                 y_mingeom,z_maxgeom,z_mingeom
      parameter ( x_maxgeom =  10.e0     ,
     &            x_mingeom = -10.e0     ,
     &            y_maxgeom =  10.e0     ,
     &            y_mingeom = -10.e0     ,
     &            z_maxgeom =  10.e0     ,
     &            z_mingeom = -10.e0     )
    
c     SPH algorithm for particle approximation (pa_sph)
c     pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
c              2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
      integer pa_sph 
      parameter(pa_sph = 2)

c     Nearest neighbor particle searching (nnps) method
c     nnps = 1 : Simplest and direct searching
c            2 : Sorting grid linked list
c            3 : Tree algorithm
      integer nnps 
      parameter(nnps = 1 )

c     Smoothing length evolution (sle) algorithm
c     sle = 0 : Keep unchanged,
c           1 : h = fac * (m/rho)^(1/dim)
c           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
c           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) ) 

      integer sle 
      parameter(sle = 0)

c     Smoothing kernel function 
c     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c         = 2, Gauss kernel   (Gingold and Monaghan 1981) 
c         = 3, Quintic kernel (Morris 1997)
c         = 4, Wendland    
      integer skf 
      parameter(skf = 4)

c     Switches for different senarios

c     summation_density = .TRUE. : Use density summation model in the code, 
c                        .FALSE.: Use continuiity equation
c     average_velocity = .TRUE. : Monaghan treatment on average velocity,
c                       .FALSE.: No average treatment.
c     config_input = .TRUE. : Load initial configuration data,
c                   .FALSE.: Generate initial configuration.
c     virtual_part = .TRUE. : Use vritual particle,
c                   .FALSE.: No use of vritual particle.
c     vp_input = .TRUE. : Load virtual particle information,
c               .FALSE.: Generate virtual particle information.
c     visc = .true. : Consider viscosity,
c           .false.: No viscosity.
c     ex_force =.true. : Consider external force,
c               .false.: No external force.
c     visc_artificial = .true. : Consider artificial viscosity,
c                      .false.: No considering of artificial viscosity.
c     heat_artificial = .true. : Consider artificial heating,
c                      .false.: No considering of artificial heating.
c     self_gravity = .true. : Considering self_gravity,
c                    .false.: No considering of self_gravity
c     nor_density =  .true. : Density normalization by using CSPM,
c                    .false.: No normalization.
      logical summation_density, average_velocity, config_input,
     &        virtual_part, vp_input, visc, ex_force, heat_artificial,
     &        visc_artificial, nor_density
      parameter ( summation_density  = .false. )       !
      parameter ( average_velocity  = .false. )       !
      parameter ( config_input  = .false. )
      parameter ( virtual_part  = .true. )
      parameter ( vp_input  = .false.  )
      parameter ( visc  = .true.  )
      parameter ( ex_force  = .true.)
      parameter ( visc_artificial  = .false. )         !
      parameter ( heat_artificial  = .false. )
      !parameter ( self_gravity  = .true. )      
      parameter ( nor_density  = .false. )             !

      integer soil_pressure ! =1, eos; =2, mean trace
      parameter(soil_pressure = 2)
      integer stress_integration
      parameter ( stress_integration = 1)
      integer yield_criterion
      parameter ( yield_criterion = 2)
      integer plasticity  ! =0 non; =1 Bui, =2 return mapping =3 Lopez
      parameter ( plasticity = 3 )
      logical artificial_density                  
      parameter (artificial_density = .true.)
      logical soil_artificial_stress
      parameter (soil_artificial_stress = .false.)

      logical volume_fraction
      parameter(volume_fraction=.true.)
      logical water_artificial_volume
      parameter (water_artificial_volume = .true.)
      logical, parameter :: volume_fraction_renorm = .true.

! 0 ignor; 1 negative pressure to zero; 2 artficial stress
      integer water_tension_instability
      parameter (water_tension_instability = 0)

c     Symmetry of the problem
c     nsym = 0 : no symmetry,
c          = 1 : axis symmetry,
c          = 2 : center symmetry.     
      integer    nsym
      parameter ( nsym = 0)

c     Control parameters for output 
c     int_stat = .true. : Print statistics about SPH particle interactions.
c                        including virtual particle information.
c     print_step: Print Timestep (On Screen)
c     save_step : Save Timestep    (To Disk File)
c     moni_particle: The particle number for information monitoring.
      logical int_stat
      parameter ( int_stat = .true. )
      integer print_step, save_step, moni_particle
      parameter ( !print_step = 20,
     &            !save_step = 20,
     &            moni_particle = 264   )
           
      !double precision pi,gravity
      double precision pi
      parameter ( pi = 3.14159265358979323846 )
      !parameter ( gravity = -9.8 )
      
c     Simulation cases
c     shocktube = .true. : carry out shock tube simulation
c     shearcavity = .true. : carry out shear cavity simulation
      logical shocktube, shearcavity, waterjet
      parameter ( shocktube  = .false. )
      parameter ( shearcavity  = .false. )
      parameter ( waterjet = .true.)

!     Recorde time interval
      integer save_step_from, save_step_to
      parameter (save_step_from = 0, save_step_to = 100)

      end module
