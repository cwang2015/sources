!-----------------------------
    module declarations_sph
!-----------------------------
use m_particles
!use param
implicit none

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
logical :: average_velocity  = .true.         
logical :: config_input  = .false. 
logical :: virtual_part  = .true. 
logical :: vp_input  = .false.  
logical :: visc  = .true.  
logical :: ex_force  = .true.
logical :: visc_artificial  = .true. 
logical :: heat_artificial  = .false. 
!logical :: self_gravity  = .true.      
logical :: nor_density  = .false.              

integer :: soil_pressure = 2  ! =1, eos; =2, mean trace
integer :: stress_integration = 1
integer :: yield_criterion = 2
integer :: plasticity = 3  ! =0 non; =1 Bui, =2 return mapping =3 Lopez
logical :: artificial_density = .true.                  
logical :: soil_artificial_stress = .true.

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
!double precision, parameter :: pi = 3.14159265358979323846 
!parameter ( gravity = -9.8 )
      

!     Recorde time interval
integer :: save_step_from = 0, save_step_to = 100

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
!type(parameters), target :: params
type(numerical), target :: numeric

! Particles 
type(particles), target :: parts, soil

! Working arrays
double precision, pointer, dimension(:) :: txx=>null(), tyy=>null(), txy=>null()

!Numerical paramters
double precision :: dt, time = 0.d0
integer :: maxtimestep = 0 , itimestep = 0

contains

!     *********************
      subroutine read_cas 
!     *********************
! Function: read the parameter file
! ----------------------------------------------------
      use constants
      use declarations
!
      implicit none
!
!     maxlr:  maximum length of a record 
!     maxlc:  maximum length of a command
!     cmd  :  command = parameter's name + it's value
!     filename : name of the steering file
!                                   
      character(len=maxlc ):: cmd
      character(len=maxlr) :: record
      integer              :: stat = 0, lenrc
      logical              :: alive 
! 
      inquire(file=cas_file,exist=alive)
      if(.not.alive)then
        write(*,*) "Steering file ", cas_file, " doesn't exist."
        stop
      endif
!
      open(unit=unprm,file=cas_file)
!
      cmd = ""
      do while(stat.eq.0)
         read(unprm, "(A)",iostat=stat) record
!         if(stat/=0)exit
!
! End of file: Treat the last unfinished command, and exit.
!
         record = adjustl(record)
         if(stat/=0 .or. record(1:1).eq.'&')then ! '&' may be used for namelist
            call treat_cmd(trim(cmd))      
            exit
         endif
!
         lenrc  = len_trim(record)
         if(lenrc.eq.0)cycle                   ! blank line
         if(record(1:1).eq."#")cycle           ! comments begins with #
!
         cmd=trim(cmd)//trim(record)
!
         if(record(lenrc:lenrc).eq.",")cycle   ! command not completed
!
         call treat_cmd(trim(cmd))             ! treat a complete command        
!
         cmd = ""                              ! new command begins
      enddo

!      if(cmd.eq."")return
!      write(*,*) 'cmd = ', cmd

!      write(*,*) 'Final parameter: ', reynb, contin, files


      close(unprm)
      return
      end subroutine

!     *************************
      subroutine treat_cmd(cmd)
!     *************************
! Function:
!     1. split the command into two parts: name + value
!     2. count how many ',' exist in the value parts
!     3. assign the value to the parameter
! ------------------------------------------------------
      use constants
!
      implicit none
!
!     cmd  :  command = name + value
!     pname:  parameter's name
!     pval :  parameter's value
!      
      character(len=*)        :: cmd
      character(len=len(cmd)) :: pname, pvalu
      integer(kind=ik_3)     :: lenc, idx
      integer(kind=ik_3)     :: ncoma, scount
      external                :: scount      
!       
      cmd  = trim(adjustl(cmd))
      lenc = len(cmd)
!
      if(lenc.eq.0)return
!
      if(cmd(lenc:lenc).eq.",") cmd(lenc:lenc)="" !remove the last ","
!      
      idx=index(cmd,'=')
!
      pname = cmd(1:idx-1)
      pvalu = cmd(idx+1:lenc)
!
      if(len_trim(adjustl(pname)).eq.0)return   
      if(len_trim(adjustl(pvalu)).eq.0)return
!
      ncoma = scount(pvalu,',')
      call usr_def_parameter(trim(adjustl(pname)),         & 
                        trim(adjustl(pvalu)), ncoma)

      return
      end subroutine

!     *****************************
      subroutine make_upper(string)
!     *****************************
! Function: lowercase -> uppercase
! ----------------------------------------------------
      implicit none
!     
      character(len=*) :: string
      integer :: len, i, code
      len = len_trim(string)
      do i = 1, len
         code = ichar(string(i:i))
         if(code>=ichar('a').and.code<=ichar('z'))then
            string(i:i)=char(code-32)
         endif
      enddo
!
      return
      end subroutine

!     ****************************
      subroutine usr_def_parameter      &
!     ****************************
          (pname, pvalu, ncoma)
!  Function: assignment via internal file
! ----------------------------------------------------------
      use constants
      use declarations
!      use param
!      use declarations_sph
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

      case DEFAULT

                 write(*,*) "Error: No such keyword -- ", pname

      end select
!          

      return
      end subroutine

!-------------------------------
     subroutine allocate_sph
!-------------------------------
implicit none
integer dim,maxn,max_interaction

!---------------------Real particles (water)---------------------
!parts%dim  = dim
!parts%maxn = maxn
!parts%max_interaction = max_interaction
dim = parts%dim
maxn = parts%maxn
max_interaction = parts%max_interaction

call allocate_particles(parts)

allocate(parts%zone(maxn)); parts%zone = 0
!allocate(parts%real_zone(parts%max_zone)); parts%real_zone = 0
!allocate(parts%virtual_zone(parts%max_zone)); parts%virtual_zone = 0

allocate(parts%w(max_interaction));             parts%w   = 0.d0
allocate(parts%dwdx(dim,max_interaction));      parts%dwdx= 0.d0
allocate(parts%pair_i(max_interaction));        parts%pair_i = 0
allocate(parts%pair_j(max_interaction));        parts%pair_j = 0
allocate(parts%countiac(maxn));                 parts%countiac=0
allocate(parts%grid(parts%maxngx,parts%maxngy,parts%maxngz)); parts%grid = 0
allocate(parts%xgcell(dim,maxn)); parts%xgcell = 0
allocate(parts%celldata(maxn)); parts%celldata = 0

! Fields variables

call allocate_particles_fields(parts)

! Working arrays outside particles type

allocate(txx(maxn), tyy(maxn), txy(maxn))
txx = 0.d0; tyy = 0.d0; txy = 0.d0

parts%txx => txx; parts%tyy => tyy; parts%txy => txy

if(trim(parts%imaterial)=='water')then
   allocate(parts%sxx,parts%sxy,parts%syy)
   parts%sxx%r => txx; parts%syy%r => tyy; parts%sxy%r => txy
elseif(trim(parts%imaterial)=='soil')then
   parts%wxy => txy
endif   

! Material assignment

if(trim(parts%imaterial)=='water') parts%material => H2O
if(trim(parts%imaterial)=='soil')  parts%material => SiO2

! Numerical parameters
parts%numeric => numeric

call parts%get_num_threads

!--------------- Soil particles --------------------------------

if(single_phase) return

soil%dim = dim
soil%maxn= maxn
soil%max_interaction = max_interaction

call allocate_particles(soil)

allocate(soil%zone(maxn)); soil%zone = 0
!allocate(soil%real_zone(soil%max_zone)); soil%real_zone = 0
!allocate(soil%virtual_zone(soil%max_zone)); soil%virtual_zone = 0

soil%w      => parts%w
soil%dwdx   => parts%dwdx
soil%pair_i => parts%pair_i
soil%pair_j => parts%pair_j
allocate(soil%countiac(maxn));    soil%countiac=0
soil%grid   => parts%grid; soil%xgcell => parts%xgcell
soil%celldata => parts%celldata

! Fields variables
call allocate_particles_fields(soil)

soil%txx => txx; soil%tyy => tyy; soil%txy => txy
soil%wxy => txy

soil%material => SiO2
soil%numeric => numeric

call soil%get_num_threads

return
end subroutine

!-------------------------------------------
     subroutine allocate_particles(parts)
!-------------------------------------------
implicit none

type(particles) parts
integer dim,maxn,max_interaction

dim  = parts%dim
maxn = parts%maxn
max_interaction = parts%max_interaction 

if(.not.associated(parts%itype)) allocate(parts%itype(maxn))
if(.not.associated(parts%x))     allocate(parts%x(dim,maxn))
if(.not.associated(parts%vol))   allocate(parts%vol(maxn))
if(.not.associated(parts%mass))  allocate(parts%mass(maxn))
if(.not.associated(parts%hsml))  allocate(parts%hsml(maxn))

return
end subroutine

!---------------------------------------------------------
     subroutine allocate_particles_fields(parts)
!---------------------------------------------------------
! Fields variables
!---------------------------------------------------------
implicit none

type(particles) parts
integer dim,maxn

dim  = parts%dim
maxn = parts%maxn

!allocate(parts%vx(dim,maxn));  parts%vx  = 0.d0
allocate(parts%vx); allocate(parts%vx%x,parts%vx%y)
allocate(parts%vx%x%r(maxn), parts%vx%y%r(maxn))
parts%vx%x%r = 0.d0; parts%vx%y%r = 0.d0
allocate(parts%rho); allocate(parts%rho%r(maxn)); parts%rho%r = 0.d0
allocate(parts%p); allocate(parts%p%r(maxn)); parts%p%r = 0.d0
allocate(parts%u(maxn));       parts%u   = 0.d0
allocate(parts%c(maxn));       parts%c   = 0.d0
!allocate(parts%eta(maxn));     parts%eta = 0.d0 

! Old value
allocate(parts%v_min(dim,maxn));   parts%v_min   = 0.d0
allocate(parts%rho_min(maxn));     parts%rho_min = 0.d0
allocate(parts%u_min(maxn));       parts%u_min   = 0.d0

! Accelerations
allocate(parts%drho(maxn));    parts%drho= 0.d0
allocate(parts%dvx); allocate(parts%dvx%x,parts%dvx%y)
allocate(parts%dvx%x%r(maxn)); parts%dvx%x%r = 0.d0  
allocate(parts%dvx%y%r(maxn)); parts%dvx%y%r = 0.d0  
allocate(parts%du(maxn));      parts%du  = 0.d0 
allocate(parts%av(dim,maxn));  parts%av  = 0.d0

! Volume Fraction
allocate(parts%vof); allocate(parts%vof%r(maxn))
parts%vof%r = 0.d0

if(trim(parts%imaterial)=='water')then
   allocate(parts%vof2(maxn))          ! phi_f = 1- phi_s
   parts%vof2 = 0.d0
   allocate(parts%dvof(maxn))
   parts%dvof = 0.d0
   allocate(parts%vof_min(maxn))
   parts%vof_min = 0.d0
endif

if(trim(parts%imaterial)=='soil')then
   allocate(parts%sxx,parts%sxy,parts%syy)
   allocate(parts%sxx%r(maxn), parts%syy%r(maxn), parts%sxy%r(maxn))
   parts%sxx%r = 0.d0; parts%sxy%r = 0.d0; parts%syy%r = 0.d0
   allocate(parts%sxx_min(maxn))
   allocate(parts%syy_min(maxn))
   allocate(parts%sxy_min(maxn))
   allocate(parts%p_min(maxn))
   allocate(parts%vcc(maxn))    
   !allocate(parts%szz(maxn))    !!! Intel Fortran Compiler is shit!!!
   parts%sxx_min = 0.d0; parts%sxy_min = 0.d0; parts%syy_min = 0.d0
   parts%p_min = 0.d0
   parts%vcc = 0.d0
   allocate(parts%dsxx(maxn),parts%dsyy(maxn),parts%dsxy(maxn))
   parts%dsxx = 0.d0; parts%dsxy = 0.d0; parts%dsyy = 0.d0
   allocate(parts%dp(maxn))
   parts%dp(maxn) = 0.d0
! For return mapping algorithm
   allocate(parts%dsxx2(maxn),parts%dsyy2(maxn),parts%dsxy2(maxn))
   parts%dsxx2 = 0.d0; parts%dsxy2 = 0.d0; parts%dsyy2 = 0.d0
   allocate(parts%dp2(maxn))
   parts%dp2(maxn) = 0.d0

   allocate(parts%fail(maxn)); parts%fail = 0

   allocate(parts%epsilon_p(maxn))
   parts%epsilon_p = 0.d0
endif

allocate(parts%drhodx(dim,maxn))
parts%drhodx = 0.d0

return
end subroutine

!----------------------------------------------------------------------
      subroutine time_integration
!----------------------------------------------------------------------
!      use param
!      use declarations_sph
      implicit none     
!
      integer :: i, j, k, d, ntotal, it
      type(particles), pointer :: pl
      
      call parts%setup_ndim1; call soil%setup_ndim1

      do it = 1, maxtimestep   
  
        itimestep = itimestep+1

        parts%itimestep = itimestep
        soil%itimestep  = itimestep

        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  
  
!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

        if (itimestep .ne. 1) then
           pl => parts
           call first_half
        
                             ! For soil
                             !if(.not.single_phase)then
                                pl => soil
                                call first_half
                                call drucker_prager_failure_criterion(pl)
                                !call volume_fraction_soil(soil)
                                !call volume_fraction_water(parts, soil)
                             !endif  

        endif

!---  Definition of variables out of the function vector:    

        call single_step

        if (itimestep .eq. 1) then
           pl => parts
           call first_step

                                !if(.not.single_phase)then
                                    pl => soil
                                    call first_step

                                    if(plasticity==2)then
                                       call plastic_or_not(pl)
                                       call return_mapping
                                    endif

                                    call drucker_prager_failure_criterion(pl)
                                    !call volume_fraction_soil(soil)
                                    !call volume_fraction_water(parts, soil)
                                !endif
        else
   
           pl => parts
           call second_half          

                                !if(.not.single_phase)then
                                   pl => soil
                                   call second_half

                                   if(plasticity==2)then
                                      call plastic_or_not(pl)
                                      call return_mapping
                                   endif

                                   call drucker_prager_failure_criterion(pl)
                                   !call volume_fraction_soil(soil)
                                   !call volume_fraction_water(parts, soil)
                                !endif

        endif 


! For saturated flow

!        if(.not.single_phase)then
!           call volume_fraction_soil(soil)
!           call volume_fraction_water(parts, soil)
!        endif        

!--------------------Velocity Inlet-----------------------


        time = time + dt

!	if (itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0) then
!          call output
!	endif 

!        if (mod(itimestep,print_step).eq.0) then
!          write(*,*)
!          write(*,101)'x','velocity', 'dvx'    
!          write(*,100)x(1,moni_particle), vx(1,moni_particle), 
!     &                dvx(1,moni_particle)    
!        endif
        
!101     format(1x,3(2x,a12))	 
!100     format(1x,3(2x,e13.6))

      enddo

      return

! -------------------
      contains
! -------------------

      subroutine first_half
! ---------------------------------------------------------------------
implicit none
type(p2r) vxi(3), dvxi(3)

      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            pl%rho_min(i) = pl%rho%r(i)
            pl%rho%r(i) = pl%rho%r(i) +(dt/2.)* pl%drho(i)
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof_min(i) = pl%vof%r(i)
            pl%vof%r(i) = pl%vof%r(i)+(dt/2.)*pl%dvof(i) 
         endif
 
         if(trim(pl%imaterial)=='soil')then
                  if(stress_integration==1)then
         !pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)   ! Nothing to do!
         !pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
         !pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
                  elseif(stress_integration==2)then
         pl%sxx_min(i) = pl%sxx%r(i)
         pl%sxy_min(i) = pl%sxy%r(i)
         pl%syy_min(i) = pl%syy%r(i)
         pl%p_min(i)   = pl%p%r(i)
         pl%sxx%r(i) = pl%sxx%r(i) + (dt/2.)*pl%dsxx(i)
         pl%sxy%r(i) = pl%sxy%r(i) + (dt/2.)*pl%dsxy(i)
         pl%syy%r(i) = pl%syy%r(i) + (dt/2.)*pl%dsyy(i)
         pl%p%r(i) = pl%p%r(i) + (dt/2.)*pl%dp(i)      !!!simultaneous pressure
                  endif
         endif
          
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%p2cmpt(i); dvxi = pl%dvx%p2cmpt(i)
         do d = 1, pl%dim
         !   pl%v_min(d, i) = pl%vx(d, i)
         !   pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
            pl%v_min(d, i) = vxi(d)%p
            vxi(d)%p = vxi(d)%p + (dt/2.)*dvxi(d)%p         
         enddo

            !pl%v_min(1, i) = pl%vx%x%r(i)
            !pl%vx%x%r(i) = pl%vx%x%r(1, i) + (dt/2.)*pl%dvx%x%r(i)
            !pl%v_min(2, i) = pl%vx%y%r(i)
            !pl%vx%y%r(i) = pl%vx%y%r(i) + (dt/2.)*pl%dvx%y%r(i)
      enddo

      return
      end subroutine

      subroutine first_step
! -------------------------------------------------------------------
      implicit none
      type(p2r) vxi(3), dvxi(3)

      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if (.not.summation_density ) then
            pl%rho%r(i) = pl%rho%r(i) + (dt/2.)* pl%drho(i)
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof%r(i) = pl%vof%r(i)+(dt/2.)*pl%dvof(i)
         endif

         if(trim(pl%imaterial)=='soil')then
                  if(stress_integration==1)then
         pl%sxx%r(i) = pl%sxx%r(i) + dt*pl%dsxx(i)
         pl%sxy%r(i) = pl%sxy%r(i) + dt*pl%dsxy(i)
         pl%syy%r(i) = pl%syy%r(i) + dt*pl%dsyy(i)
         pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%sxx%r(i) = pl%sxx%r(i) + (dt/2.)*pl%dsxx(i)
         pl%sxy%r(i) = pl%sxy%r(i) + (dt/2.)*pl%dsxy(i)
         pl%syy%r(i) = pl%syy%r(i) + (dt/2.)*pl%dsyy(i)
         pl%p%r(i)   = pl%p%r(i) + (dt/2.)*pl%dp(i)    !!! simultaneous pressure
                  endif
         endif
        
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%p2cmpt(i); dvxi = pl%dvx%p2cmpt(i)
         do d = 1, pl%dim        
            !pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = vxi(d)%p + (dt/2.) * dvxi(d)%p + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * vxi(d)%p            
         enddo           
            !pl%vx(1, i) = pl%vx(1, i) + (dt/2.) * pl%dvx%x%r(i)   &
            !            + pl%av(1, i)
            !pl%x(1, i) = pl%x(1, i) + dt * pl%vx(1, i)
            !pl%vx(2, i) = pl%vx(2, i) + (dt/2.) * pl%dvx%y%r(i)   &
            !            + pl%av(2, i)
            !pl%x(2, i) = pl%x(2, i) + dt * pl%vx(2, i)
      enddo 
      
      return
      end subroutine

      subroutine second_half
! ------------------------------------------------------------------
      implicit none
      type(p2r) vxi(3), dvxi(3)

      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if (.not.summation_density ) then 
            pl%rho%r(i) = pl%rho_min(i) + dt*pl%drho(i)
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof%r(i) = pl%vof_min(i)+dt*pl%dvof(i) 
         endif
                 
         if(trim(pl%imaterial)=='soil')then 
                  if(stress_integration==1)then 
         pl%sxx%r(i) = pl%sxx%r(i) + dt*pl%dsxx(i)
         pl%sxy%r(i) = pl%sxy%r(i) + dt*pl%dsxy(i)
         pl%syy%r(i) = pl%syy%r(i) + dt*pl%dsyy(i)
         pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%sxx%r(i) = pl%sxx_min(i) + dt*pl%dsxx(i)
         pl%sxy%r(i) = pl%sxy_min(i) + dt*pl%dsxy(i)
         pl%syy%r(i) = pl%syy_min(i) + dt*pl%dsyy(i)
         pl%p%r(i)   = pl%p_min(i) + dt*pl%dp(i)     !!! simultaneous pressure
                  endif
         endif

         if(pl%itype(i)<0)cycle
         vxi = pl%vx%p2cmpt(i); dvxi = pl%dvx%p2cmpt(i)
         do d = 1, pl%dim                   
            !pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = pl%v_min(d, i) + dt * dvxi(d)%p   &
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * vxi(d)%p
         enddo
            !pl%vx(1, i) = pl%v_min(1, i) + dt * pl%dvx%x%r(i)   &
            !            + pl%av(1, i)
            !pl%x(1, i) = pl%x(1, i) + dt * pl%vx(1, i)                  
            !pl%vx(2, i) = pl%v_min(2, i) + dt * pl%dvx%y%r(i)   &
            !            + pl%av(2, i)
            !pl%x(2, i) = pl%x(2, i) + dt * pl%vx(2, i)                              
      enddo

      return
      end subroutine

      subroutine return_mapping
! ------------------------------------------------------------------
      implicit none

      do i=1,pl%ntotal+pl%nvirt  ! origionally pl%ntotal            
                 
         !if(trim(pl%imaterial)=='soil')then
         if(pl%fail(i)==1)then 
                  if(stress_integration==1)then 
         pl%sxx%r(i) = pl%sxx%r(i) + dt*pl%dsxx2(i)
         pl%sxy%r(i) = pl%sxy%r(i) + dt*pl%dsxy2(i)
         pl%syy%r(i) = pl%syy%r(i) + dt*pl%dsyy2(i)
         pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp2(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%sxx%r(i) = pl%sxx_min(i) + dt*pl%dsxx(i)
         pl%sxy%r(i) = pl%sxy_min(i) + dt*pl%dsxy(i)
         pl%syy%r(i) = pl%syy_min(i) + dt*pl%dsyy(i)
         pl%p%r(i)   = pl%p_min(i) + dt*pl%dp(i)     !!! simultaneous pressure
                  endif
         endif
         !endif
         
      enddo

      return
      end subroutine

      end subroutine

!DEC$IF(.FALSE.)
!----------------------------------------------------------------------      
      subroutine time_integration_for_water
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it
type(particles), pointer :: pl
              
do it = 1, maxtimestep    
   itimestep = itimestep+1
   parts%itimestep = itimestep

!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

   if(itimestep .ne. 1)then
      pl => parts
      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            pl%rho_min(i) = pl%rho%r(i)
            pl%rho%r(i) = pl%rho%r(i) +(dt/2.)* pl%drho(i)
         endif
          
         if(pl%itype(i)<0)cycle

         !do d = 1, pl%dim
         !   pl%v_min(d, i) = pl%vx(d, i)
         !   pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
         !enddo
            pl%v_min(1, i) = pl%vx(1, i)
            pl%vx(1, i) = pl%vx(1, i) + (dt/2.)*pl%dvx%x%r(i)
            pl%v_min(2, i) = pl%vx(2, i)
            pl%vx(2, i) = pl%vx(2, i) + (dt/2.)*pl%dvx%y%r(i)
      enddo
           
   endif

   call single_step_for_water

   if(itimestep .eq. 1) then
      pl => parts
      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if(.not.summation_density )then
            pl%rho%r(i) = pl%rho%r(i) + (dt/2.)* pl%drho(i)
         endif

         if(pl%itype(i)<0)cycle

!         do d = 1, pl%dim        
!            pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   & 
!                        + pl%av(d, i)
!            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
!         enddo           
            pl%vx(1, i) = pl%vx(1, i) + (dt/2.) * pl%dvx%x%r(i)   & 
                        + pl%av(1, i)
            pl%x(1, i) = pl%x(1, i) + dt * pl%vx(1, i)
            pl%vx(2, i) = pl%vx(2, i) + (dt/2.) * pl%dvx%y%r(i)   & 
                        + pl%av(2, i)
            pl%x(2, i) = pl%x(2, i) + dt * pl%vx(2, i)
      enddo 
           
   else

       pl => parts
       do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if(.not.summation_density )then 
            pl%rho%r(i) = pl%rho_min(i) + dt*pl%drho(i)
         endif

         if(pl%itype(i)<0)cycle

!         do d = 1, pl%dim                   
!            pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)    & 
!                        + pl%av(d, i)
!            pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)                  
!         enddo
            pl%vx(1, i) = pl%v_min(1, i) + dt * pl%dvx%x%r(i)    & 
                        + pl%av(1, i)
            pl%x(1, i) = pl%x(1, i) + dt * pl%vx(1, i)                  
            pl%vx(2, i) = pl%v_min(2, i) + dt * pl%dvx%y%r(i)    & 
                        + pl%av(2, i)
            pl%x(2, i) = pl%x(2, i) + dt * pl%vx(2, i)                              
      enddo

   endif 

   time = time + dt

   if(mod(itimestep,print_step).eq.0)then
      write(*,*)'______________________________________________'
      write(*,*)'  current number of time step =',                &
                itimestep,'     current time=', real(time+dt)
      write(*,*)'______________________________________________'
   endif  

   if(itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0)then
      call output
   endif 

enddo

return
end subroutine      

!DEC$ENDIF
!DEC$IF(.FALSE.)
!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_verlet
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
type(particles), pointer :: pl
              
real(dp), allocatable, dimension(:,:) :: lastvx   
real(dp), allocatable, dimension(:,:) :: temp1    
real(dp), allocatable, dimension(:)   :: lastrho  
real(dp), allocatable, dimension(:)   :: temp2    

allocate(lastvx(2,parts%maxn))
allocate(temp1(2,parts%maxn))
allocate(lastrho(parts%maxn))
allocate(temp2(parts%maxn))
    pl => parts
    lastvx = pl%vx
    lastrho = pl%rho%r

do it = 1, maxtimestep 
    itimestep = itimestep+1
   call single_step_for_water

if(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal
       !do d = 1, pl%dim
       !   pl%x(d,i) = pl%x(d,i) + dt * pl%vx(d,i) +(dt**2./2.)*pl%dvx(d,i)
       !   lastvx(d,i) = lastvx(d,i) + 2.*dt*pl%dvx(d,i)
       
       !enddo
          pl%x(1,i) = pl%x(1,i) + dt * pl%vx(1,i) +(dt**2./2.)*pl%dvx%x%r(i)
          lastvx(1,i) = lastvx(1,i) + 2.*dt*pl%dvx%x%r(i)
          pl%x(2,i) = pl%x(2,i) + dt * pl%vx(2,i) +(dt**2./2.)*pl%dvx%y%r(i)
          lastvx(2,i) = lastvx(2,i) + 2.*dt*pl%dvx%y%r(i)       
    enddo
else
    do i = 1, pl%ntotal
       !do d = 1, pl%dim
       !   pl%x(d,i) = pl%x(d,i) + dt * pl%vx(d,i) +(dt**2./2.)*pl%dvx(d,i)
       !   lastvx(d,i) = pl%vx(d,i) + dt*pl%dvx(d,i)
       !enddo
          pl%x(1,i) = pl%x(1,i) + dt * pl%vx(1,i) +(dt**2./2.)*pl%dvx%x%r(i)
          lastvx(1,i) = pl%vx(1,i) + dt*pl%dvx%x%r(i)
          pl%x(2,i) = pl%x(2,i) + dt * pl%vx(2,i) +(dt**2./2.)*pl%dvx%y%r(i)
          lastvx(2,i) = pl%vx(2,i) + dt*pl%dvx%y%r(i)       
    enddo
endif
temp1 = pl%vx
pl%vx = lastvx
lastvx = temp1

if(itimestep .eq. 1) then
    do i = 1, pl%ntotal+pl%nvirt
          lastrho(i) = pl%rho%r(i) + dt * 2. * pl%drho(i)
    enddo
elseif(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal+pl%nvirt
          lastrho(i) = lastrho(i) + dt * 2. * pl%drho(i)
    enddo
else
   do i = 1, pl%ntotal+pl%nvirt
          lastrho(i) = pl%rho%r(i) + dt* pl%drho(i)
   enddo
endif
temp2 = pl%rho%r
pl%rho%r = lastrho
lastrho = temp2

   time = time + dt

   if(mod(itimestep,print_step).eq.0)then
      write(*,*)'______________________________________________'
      write(*,*)'  current number of time step =',                &
                itimestep,'     current time=', real(time+dt)
      write(*,*)'______________________________________________'
   endif  

   if(itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0)then
      call output
   endif 

enddo

deallocate(lastvx)
deallocate(lastrho)
deallocate(temp1)
deallocate(temp2)


return
end subroutine      
         
!DEC$ENDIF

!DEC$IF(.FALSE.)
!----------------------------------------------------------------------
             subroutine time_integration_for_soil
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it
type(particles), pointer :: pl
type(p2r) vxi(3), dvxi(3)
              
do it = 1, maxtimestep     
   itimestep = itimestep+1
   parts%itimestep = itimestep
  
! If not first time step, then update thermal energy, density and 
! velocity half a time step  

   if(itimestep .ne. 1)then

      pl => parts
!      call first_half
      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            pl%rho_min(i) = pl%rho%r(i)
            pl%rho%r(i) = pl%rho%r(i) +(dt/2.)* pl%drho(i)
         endif
 
         if(stress_integration==1)then
            !pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)   ! Nothing to do!
            !pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
            !pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
         elseif(stress_integration==2)then
            pl%sxx_min(i) = pl%sxx%r(i)
            pl%sxy_min(i) = pl%sxy%r(i)
            pl%syy_min(i) = pl%syy%r(i)
            pl%p_min(i)   = pl%p%r(i)
            pl%sxx%r(i) = pl%sxx%r(i) + (dt/2.)*pl%dsxx(i)
            pl%sxy%r(i) = pl%sxy%r(i) + (dt/2.)*pl%dsxy(i)
            pl%syy%r(i) = pl%syy%r(i) + (dt/2.)*pl%dsyy(i)
            pl%p%r(i) = pl%p%r(i) + (dt/2.)*pl%dp(i)      !!!simultaneous pressure
         endif
          
         if(pl%itype(i)<0)cycle
         !do d = 1, pl%dim
         !   pl%v_min(d, i) = pl%vx(d, i)
         !   pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
         !enddo
            !pl%v_min(1, i) = pl%vx(1, i)
            !pl%vx(1, i) = pl%vx(1, i) + (dt/2.)*pl%dvx%x%r(i)
            !pl%v_min(2, i) = pl%vx(2, i)
            !pl%vx(2, i) = pl%vx(2, i) + (dt/2.)*pl%dvx%y%r(i)
      enddo 

      call drucker_prager_failure_criterion(pl)

   endif

   call single_step_for_soil

   if(itimestep .eq. 1)then
      pl => parts
!           call first_step
      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if(.not.summation_density )then
            pl%rho%r(i) = pl%rho%r(i) + (dt/2.)* pl%drho(i)
         endif

         if(stress_integration==1)then
            pl%sxx%r(i) = pl%sxx%r(i) + dt*pl%dsxx(i)
            pl%sxy%r(i) = pl%sxy%r(i) + dt*pl%dsxy(i)
            pl%syy%r(i) = pl%syy%r(i) + dt*pl%dsyy(i)
            pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp(i)       !!! simultaneous pressure
         elseif(stress_integration==2)then
            pl%sxx%r(i) = pl%sxx%r(i) + (dt/2.)*pl%dsxx(i)
            pl%sxy%r(i) = pl%sxy%r(i) + (dt/2.)*pl%dsxy(i)
            pl%syy%r(i) = pl%syy%r(i) + (dt/2.)*pl%dsyy(i)
            pl%p%r(i)   = pl%p%r(i) + (dt/2.)*pl%dp(i)    !!! simultaneous pressure
         endif
        
         if(pl%itype(i)<0)cycle

         vxi = pl%vx%p2cmpt(i); dvxi = pl%dvx%p2cmpt(i)
         do d = 1, pl%dim        
            vxi%p(d) = vxi%p(d) + (dt/2.) * dvxi%p(d)   &
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * vxi%p(d)
         enddo           
            !pl%vx(1, i) = pl%vx(1, i) + (dt/2.) * pl%dvx%x%r(i)   &
            !            + pl%av(1, i)
            !pl%x(1, i) = pl%x(1, i) + dt * pl%vx(1, i)
            !pl%vx(2, i) = pl%vx(2, i) + (dt/2.) * pl%dvx%y%r(i)   &
            !            + pl%av(2, i)
            !pl%x(2, i) = pl%x(2, i) + dt * pl%vx(2, i)         

      enddo 

!              if(plasticity==2)then
!                 call plastic_or_not(pl)
!                 call return_mapping
!              endif

      call drucker_prager_failure_criterion(pl)

   else
   
      pl => parts
!           call second_half          
      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if(.not.summation_density )then 
            pl%rho%r(i) = pl%rho_min(i) + dt*pl%drho(i)
         endif
                 
         if(stress_integration==1)then 
            pl%sxx%r(i) = pl%sxx%r(i) + dt*pl%dsxx(i)
            pl%sxy%r(i) = pl%sxy%r(i) + dt*pl%dsxy(i)
            pl%syy%r(i) = pl%syy%r(i) + dt*pl%dsyy(i)
            pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp(i)       !!! simultaneous pressure
         elseif(stress_integration==2)then
            pl%sxx%r(i) = pl%sxx_min(i) + dt*pl%dsxx(i)
            pl%sxy%r(i) = pl%sxy_min(i) + dt*pl%dsxy(i)
            pl%syy%r(i) = pl%syy_min(i) + dt*pl%dsyy(i)
            pl%p%r(i)   = pl%p_min(i) + dt*pl%dp(i)     !!! simultaneous pressure
         endif

         if(pl%itype(i)<0)cycle
         vxi = pl%vx%p2cmpt(i); dvxi = pl%dvx%p2cmpt(i)
         do d = 1, pl%dim                   
            vxi%p(d) = pl%v_min(d, i) + dt * dvxi%p(d)   &
                        + pl%av(d, i)
            pl%x(d, i) = pl%x(d, i) + dt * vxi%p(d)                  
         enddo
            !pl%vx(1, i) = pl%v_min(1, i) + dt * pl%dvx%x%r(i)   &
            !            + pl%av(1, i)
            !pl%x(1, i) = pl%x(1, i) + dt * pl%vx(1, i)                 
            !pl%vx(2, i) = pl%v_min(2, i) + dt * pl%dvx%y%r(i)   &
            !            + pl%av(2, i)
            !pl%x(2, i) = pl%x(2, i) + dt * pl%vx(2, i)                              
      enddo

!              if(plasticity==2)then
!                 call plastic_or_not(pl)
!                 call return_mapping
!              endif

      call drucker_prager_failure_criterion(pl)

   endif 

   time = time + dt

   if(mod(itimestep,print_step).eq.0)then
      write(*,*)'______________________________________________'
      write(*,*)'  current number of time step =',              &
                itimestep,'     current time=', real(time+dt)
      write(*,*)'______________________________________________'
   endif  
        
   if(itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0) then
      call output
   endif 

enddo

return
end subroutine

!DEC$ENDIF

!----------------------------------------------------------------------      
                   subroutine single_step
!----------------------------------------------------------------------
!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
!use param 
!use declarations_sph
!use m_sph_fo
implicit none

integer  nphase, iphase
type(particles), pointer :: pl
logical :: dbg = .false.
integer i, ntotal
                 
if(dbg) write(*,*) 'In single_step...'

                                      nphase = 2
                                      if(single_phase) nphase=1

                                      do iphase = 1, nphase
                                         if(iphase==1) pl => parts  !!!
                                         if(iphase==2) pl => soil   !!!

pl%dvx%x%r = 0.d0; pl%dvx%y%r = 0.d0; pl%drho = 0.d0
if(trim(pl%imaterial)=='water')pl%dvof = 0.d0
if(trim(pl%imaterial)=='soil')then
   pl%dsxx = 0.d0; pl%dsxy = 0.d0; pl%dsyy = 0.d0; pl%dp = 0.d0
endif
 
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
if (pl%numeric%nnps.eq.1) then 
   call direct_find(pl)
else if (pl%numeric%nnps.eq.2) then
   call link_list(pl)     
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
!        call link_list(itimestep, parts%ntotal+parts%nvirt,
!     &       parts%hsml(1),parts%x,parts%niac,parts%pair_i,
!     &       parts%pair_j,parts%w,parts%dwdx,parts%countiac)
else if (pl%numeric%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
endif         

if(mod(itimestep,print_step).eq.0.and.int_stat) then
   call pl%interaction_statistics
endif   

!--- Added by Wang
!if(nor_density) call norm_density(pl)

!---  Density approximation or change rate
     
!if(summation_density)then      
!if(mod(itimestep,30)==0) call sum_density(pl)
!else             
    call con_density(pl)         
!endif
      
if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      call delta_sph(pl,pl%rho%r,pl%drho)
   !endif
endif

!---  Dynamic viscosity:

if(trim(pl%imaterial)=='water')then
   !if (visc) call viscosity(pl)
elseif(trim(pl%imaterial)=='soil')then
   !call shear_modulus(pl)
endif
       
!---  Internal forces:

call shear_strain_rate(pl)   
if(trim(pl%imaterial)=='soil')call velocity_divergence(pl)

call pressure(pl)

if(trim(pl%imaterial)=='water')then
   call newtonian_fluid(pl)
elseif(trim(pl%imaterial)=='soil')then

   if(yield_criterion == 1)then
      call mohr_coulomb_failure_criterion(pl)
   elseif(yield_criterion == 2)then
      call drucker_prager_failure_criterion(pl)
   endif

   if(mod(itimestep,print_step).eq.0)    &
            write(*,*) 'Failured points: ', pl%nfail
endif

!Calculate internal force for water phase !! -phi_f Grad(p)
if(pl%imaterial=='water')then
   if(pl%nthreads==1)then
   pl%dvx%x%r = -pl%vof%r*pl%df(pl%p%r,'x') + pl%df(pl%vof%r*pl%sxx%r,'x') + pl%df(pl%vof%r*pl%sxy%r,'y')
   pl%dvx%y%r = -pl%vof%r*pl%df(pl%p%r,'y') + pl%df(pl%vof%r*pl%sxy%r,'x') + pl%df(pl%vof%r*pl%syy%r,'y')
   else
   pl%dvx%x = -pl%vof*pl%df_omp(pl%p,'x') + pl%df_omp(pl%vof*pl%sxx,'x') + pl%df_omp(pl%vof*pl%sxy,'y')
   pl%dvx%y = -pl%vof*pl%df_omp(pl%p,'y') + pl%df_omp(pl%vof*pl%sxy,'x') + pl%df_omp(pl%vof*pl%syy,'y')           
   endif        

   !where (pl%rho.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
   !where (pl%rho.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho
else      
   !call int_force(pl)
   pl%dvx%x%r = -pl%df3_omp(pl%vof%r*pl%p%r,'x') + pl%df3_omp(pl%vof%r*pl%sxx%r,'x') + pl%df3_omp(pl%vof%r*pl%sxy%r,'y')
   pl%dvx%y%r = -pl%df3(pl%vof%r*pl%p%r,'y') + pl%df3(pl%vof%r*pl%sxy%r,'x') + pl%df3(pl%vof%r*pl%syy%r,'y')

   !where (pl%rho%r.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
   !where (pl%rho%r.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r 
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho 

endif      

!if(trim(pl%imaterial)=='water'.and.water_tension_instability==2) &
!   call tension_instability(pl) 

! --- Plasticity flow rule   ! This was done before Jaummann_rate, because we 
!     we need txx,tyy,tzz, which was destroyed in Jaumann_rate!

if(trim(pl%imaterial)=='soil')then   
   if(plasticity==1)then 
      call plastic_flow_rule(pl)
   elseif(plasticity==2)then     
      call plastic_flow_rule2(pl)
   elseif(plasticity==3)then
      call plastic_or_not(pl)
      call plastic_flow_rule3(pl)
   endif
endif

! --- Jaumann rate  !When???

if(trim(pl%imaterial)=='soil')call Jaumann_rate(pl)

!---  Artificial viscosity:

if (visc_artificial) call pl%art_visc

if(trim(pl%imaterial)=='soil'.and.soil_artificial_stress)then
        !call art_stress(pl)
   call pl%delta_sph(pl%p%r,pl%dp)
   call pl%delta_sph(pl%sxx%r,pl%dsxx)
   call pl%delta_sph(pl%sxy%r,pl%dsxy)
   call pl%delta_sph(pl%syy%r,pl%dsyy)
endif        
if(trim(pl%imaterial)=='water'.and.water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
        call pl%delta_sph(pl%vof%r,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
          call pl%repulsive_force
!      endif

pl%dvx%y%r = pl%dvx%y%r + gravity

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (average_velocity) call av_vel(pl) 

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

                                            enddo ! iphase

!call drag_force(parts,soil)   !!! Porous media

                if(.not.single_phase)then

!-------------------Water/soil interaction-------------------------------

!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
      !if (numeric%nnps.eq.1) then 
         call direct_find_2(parts,soil)
        !call direct_find(parts)
      !else if (numeric%nnps.eq.2) then
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      !else if (numeric%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      !endif         

      if(mod(itimestep,print_step).eq.0.and.int_stat) then
         call parts%interaction_statistics
      endif    

      call darcy_law(parts,soil)          
      call pore_water_pressure(parts,soil) 

      if(volume_fraction)then
         call volume_fraction_soil(soil)
         call volume_fraction_water2(parts,soil)
         call volume_fraction_water(parts,soil)  ! phi_f = 1- phi_s
         if(volume_fraction_renorm)then
            if(mod(itimestep,40).eq.0) then
               ntotal = parts%ntotal+parts%nvirt
               parts%rho%r(1:ntotal) = parts%rho%r(1:ntotal)/parts%vof%r(1:ntotal) 
               parts%vof%r = parts%vof2
               parts%rho%r(1:ntotal) = parts%rho%r(1:ntotal)*parts%vof%r(1:ntotal) 
            endif
         endif
      endif

!      call direct_find(parts) 
!      call art_volume_fraction_water2(parts)
     
             endif ! .not.single_phase

      if(itimestep>=save_step_from.and.   &
         mod(itimestep,save_step).eq.0)then
         call output
      endif 

return
end subroutine


!-------------------------------------------------
      subroutine single_step_for_water
!-------------------------------------------------
!use param
!use declarations_sph
!use m_sph_fo
implicit none
type(particles), pointer :: pl
type(material),pointer :: property

pl => parts
property => pl%material

pl%dvx%x%r = 0.d0; pl%dvx%y%r = 0.d0; pl%drho = 0.d0

!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length

call pl%find_pairs

if(mod(itimestep,print_step).eq.0.and.int_stat) then
   call pl%interaction_statistics
endif   

!---  Density approximation or change rate
     
!if(summation_density) call sum_density(pl)
call sum_density(pl)
pl%drho = -pl%rho%r*(pl%df2(pl%vx%x%r,'x')+pl%df2(pl%vx%y%r,'y'))
      
if(artificial_density)then
   !call renormalize_density_gradient(pl)
   !call art_density(pl)
   call delta_sph(pl,pl%rho%r,pl%drho)
endif
       
!---  Internal forces:

!Calculate pressure

where(pl%rho%r>0.0) pl%p%r = property%b*((pl%rho%r/property%rho0)**property%gamma-1)

!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

pl%txx = 2./3.*(2.0*pl%df(pl%vx%x%r,'x')-pl%df(pl%vx%y%r,'y'))
pl%txy = pl%df(pl%vx%x%r,'y')+pl%df(pl%vx%y%r,'x')
pl%tyy = 2./3.*(2.0*pl%df(pl%vx%y%r,'y')-pl%df(pl%vx%x%r,'x'))

!Newtonian fluid

pl%sxx%r = property%viscosity*pl%txx
pl%syy%r = property%viscosity*pl%tyy
pl%sxy%r = property%viscosity*pl%txy

!Calculate internal force

pl%dvx%x%r = - pl%df(pl%p%r,'x') + pl%df(pl%sxx%r,'x') + pl%df(pl%sxy%r,'y')
pl%dvx%y%r = - pl%df(pl%p%r,'y') + pl%df(pl%sxy%r,'x') + pl%df(pl%syy%r,'y')

where (pl%rho%r.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
where (pl%rho%r.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r
       
!if(water_tension_instability==2) call tension_instability(pl) 

!---  Artificial viscosity:

if (visc_artificial) call pl%art_visc
    
!---  External forces:

pl%dvx%y%r = pl%dvx%y%r + gravity

!call repulsive_force(pl)

! Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)
     
! Calculating average velocity of each partile for avoiding penetration

if (average_velocity) call av_vel(pl) 

if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

return
end subroutine


!----------------------------------------------------------------------      
               subroutine single_step_for_soil
!----------------------------------------------------------------------
!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
!use param 
!use declarations_sph
!use m_sph_fo
implicit none

type(particles), pointer :: pl
type(material), pointer :: property

pl => parts
property => pl%material

pl%dvx%x%r = 0.d0; pl%dvx%y%r = 0.d0; pl%drho = 0.d0
pl%dsxx = 0.d0; pl%dsxy = 0.d0; pl%dsyy = 0.d0
pl%dp = 0.d0
 
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
call pl%find_pairs

if(mod(itimestep,print_step).eq.0.and.int_stat) then
   call pl%interaction_statistics
endif   

!---  Density approximation or change rate
     
if (summation_density) then      
    call sum_density(pl)
else             
    call con_density(pl)         
endif
      
if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      call delta_sph(pl,pl%rho%r,pl%drho)
   !endif
endif

!call shear_modulus(pl)
       
!---  Internal forces:

call shear_strain_rate(pl)   
call velocity_divergence(pl)

call pressure(pl)

if(yield_criterion == 1)then
   call mohr_coulomb_failure_criterion(pl)
elseif(yield_criterion == 2)then
   call drucker_prager_failure_criterion(pl)
endif

!call int_force1(pl)

pl%dvx%x%r = - pl%df(pl%p%r,'x') + pl%df(pl%sxx%r,'x') + pl%df(pl%sxy%r,'y')
pl%dvx%y%r = - pl%df(pl%p%r,'y') + pl%df(pl%sxy%r,'x') + pl%df(pl%syy%r,'y')

where (pl%rho%r.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
where (pl%rho%r.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r       

! --- Plasticity flow rule   ! This was done before Jaummann_rate, because we 
!     we need txx,tyy,tzz, which was destroyed in Jaumann_rate!

if(plasticity==1)then 
   call plastic_flow_rule(pl)
elseif(plasticity==2)then     
   call plastic_flow_rule2(pl)
elseif(plasticity==3)then
   call plastic_or_not(pl)
   call plastic_flow_rule3(pl)
endif

! --- Jaumann rate  !When???

call Jaumann_rate(pl)

!---  Artificial viscosity:

if (visc_artificial) call pl%art_visc

!if(soil_artificial_stress) call art_stress(pl)
if(soil_artificial_stress)then
   call pl%delta_sph(pl%p%r,pl%dp)
   call pl%delta_sph(pl%sxx%r,pl%dsxx)
   call pl%delta_sph(pl%sxy%r,pl%dsxy)
   call pl%delta_sph(pl%syy%r,pl%dsyy)
endif   

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
!      endif

pl%dvx%y%r = pl%dvx%y%r + gravity

! Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)
     
! Calculating average velocity of each partile for avoiding penetration

if (average_velocity) call av_vel(pl) 

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

return
end subroutine




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
      p      => parts%p%r

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
      write(f_xv,*) parts%vof%r(1:ntotal)
!      write(f_xv,*) u(1:ntotal)
!      write(f_xv,*) itype(1:ntotal)
!      write(f_xv,*) hsml(1:ntotal)                                        
      write(f_xv,*) parts%vx%x%r(1:ntotal)
      write(f_xv,*) parts%vx%y%r(1:ntotal)
      write(f_xv,*) parts%rho%r(1:ntotal)
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
      write(f_xv,*)  parts%p%r(1:ntotal)
!      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_xv,*)  parts%vof%r(1:ntotal)
      write(f_xv,*)  parts%epsilon_p(1:ntotal)
      write(f_xv,*)  parts%sxy%r(1:ntotal)
      write(f_xv,*)  parts%sxx%r(1:ntotal)
      write(f_xv,*)  parts%syy%r(1:ntotal)
      write(f_xv,*)  parts%vx%x%r(1:ntotal)
      write(f_xv,*)  parts%vx%y%r(1:ntotal)
      write(f_xv,*)  parts%rho%r(1:ntotal)
      write(f_xv,*)  parts%mass(1:ntotal)
      write(f_xv,*)  -parts%p%r(1:ntotal) + parts%syy%r(1:ntotal)
      write(f_xv,*)  parts%zone(1:ntotal)
      write(f_other,*) time, -parts%p%r(395)+parts%syy%r(395)

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
      write(f_state,*)  soil%p%r(1:ntotal2)
!      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_state,*)  soil%vof%r(1:ntotal2)
      write(f_state,*)  soil%epsilon_p(1:ntotal2)
      write(f_state,*)  soil%sxy%r(1:ntotal2)
      write(f_state,*)  soil%sxx%r(1:ntotal2)
      write(f_state,*)  soil%syy%r(1:ntotal2)
      write(f_state,*)  soil%vx%x%r(1:ntotal2)
      write(f_state,*)  soil%vx%y%r(1:ntotal2)
      write(f_state,*)  soil%rho%r(1:ntotal2)
      write(f_state,*)  soil%mass(1:ntotal2)
      write(f_state,*)  -soil%p%r(1:ntotal2) + soil%syy%r(1:ntotal2)
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
