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
   type(array), pointer :: tab => null()
!double precision, pointer, dimension(:) :: txx=>null(), tyy=>null(), txy=>null()

!Numerical paramters
double precision :: dt, time = 0.d0
integer :: maxtimestep = 0 , itimestep = 0

!OpenMP Number of threads
integer :: nthreads = 1
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
      integer(kind=ik_3)     :: ncoma !, scount
!      external                :: scount      
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

!     ****************************************
      subroutine trim_pvalu(pvalu,delimt,ndim)
!     ****************************************
!  Version 1.0  Sept. 03, 08   Chun WANG (SJTU) 62933888
! ---------------------------------------------------------
! Funtion: Remove the substring after the ndim^th delimiter
! ---------------------------------------------------------
      use constants
!
      implicit none
!
      character(len=*)   :: pvalu
      character(len=1)   :: delimt
      integer(kind=ik_3)  :: ndim, scount, lenp, i
!      
      pvalu = adjustl(pvalu)
      lenp  = len_trim(pvalu)
!
      scount = 0
      do i = 1, lenp
         if(pvalu(i:i).eq.delimt)then
           scount = scount + 1
         endif
         if(scount.eq.ndim)then 
           pvalu(i:lenp) = ""
           exit
         endif
      enddo
      pvalu = trim(adjustl(pvalu))
!
      return
      end subroutine
         
!     ************************
      function scount(str,chr)
!     ************************
!  Version 1.0  Sept. 02, 08  Chun WANG (SJTU) 62933888
! -----------------------------------------------------
!  Function:  Find how many times chr appeared in str.
! -----------------------------------------------------
      use constants
!
      implicit none
!
      character(len=*) :: str
      character        :: chr
      integer(kind=ik_3)  :: scount, i, len_str
!      
      len_str = len_trim(str)
!
      scount = 0
      do i = 1, len_str
         if(str(i:i).eq.chr) scount = scount + 1
      enddo
!     
      return
      end function


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

      case('OPENMP NUMBER OF THREADS')

                  read(pvalu,*) nthreads
                  write(*,*) 'Openmp number of threads = ', nthreads

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
allocate(parts%ps(max_interaction));             parts%ps   = 0.d0
allocate(parts%rhos(max_interaction));             parts%rhos   = 0.d0
allocate(parts%gammaas(max_interaction));             parts%gammaas   = 0.d0
allocate(parts%dwdx(dim,max_interaction));      parts%dwdx= 0.d0
allocate(parts%n(dim,max_interaction));      parts%n= 0.d0
allocate(parts%mons(dim,max_interaction));      parts%mons= 0.d0
allocate(parts%dgu(dim,max_interaction));       parts%dgu= 0.d0
allocate(parts%pair_i(max_interaction));        parts%pair_i = 0
allocate(parts%pair_j(max_interaction));        parts%pair_j = 0
allocate(parts%countiac(maxn));                 parts%countiac=0
allocate(parts%grid(parts%maxngx,parts%maxngy,parts%maxngz)); parts%grid = 0
allocate(parts%xgcell(dim,maxn)); parts%xgcell = 0
allocate(parts%celldata(maxn)); parts%celldata = 0

! Fields variables

call allocate_particles_fields(parts)

! Working arrays outside particles type
allocate(tab)
allocate(tab%x, tab%y, tab%xy)
allocate(tab%x%r(maxn),tab%y%r(maxn),tab%xy%r(maxn))
tab%x = 0.d0; tab%y = 0.d0; tab%xy = 0.d0

parts%tab => tab

if(trim(parts%imaterial)=='water')then
!   allocate(parts%sxx,parts%sxy,parts%syy)
    allocate(parts%str);allocate(parts%str%x,parts%str%xy,parts%str%y)
    parts%str%x%r => tab%x%r; parts%str%y%r => tab%y%r; parts%str%xy%r => tab%xy%r
elseif(trim(parts%imaterial)=='soil')then
   parts%wxy%r => tab%xy%r
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
!soil%grid   => parts%grid; soil%xgcell => parts%xgcell
!soil%celldata => parts%celldata
allocate(soil%grid(soil%maxngx,soil%maxngy,soil%maxngz)); soil%grid = 0
allocate(soil%xgcell(dim,maxn)); soil%xgcell = 0
allocate(soil%celldata(maxn)); soil%celldata = 0

! Fields variables
call allocate_particles_fields(soil)

soil%tab => tab
soil%wxy => tab%xy

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
!if(.not.associated(parts%vol))   allocate(parts%vol(maxn))
if(.not.associated(parts%vol))then
    allocate(parts%vol); allocate(parts%vol%r(maxn))
endif
if(.not.associated(parts%vol_nm))then
    allocate(parts%vol_nm); allocate(parts%vol_nm%r(maxn))
endif
if(.not.associated(parts%mass)) then
    allocate(parts%mass); allocate(parts%mass%r(maxn))
endif
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
allocate(parts%mone); allocate(parts%mone%x,parts%mone%y)
allocate(parts%mone%x%r(maxn), parts%mone%y%r(maxn))
parts%mone%x%r = 0.d0; parts%mone%y%r = 0.d0
allocate(parts%gammaa); allocate(parts%gammaa%r(maxn))
parts%gammaa%r = 0.d0
!allocate(parts%mons); allocate(parts%mons%x,parts%mons%y)
!allocate(parts%mons%x%r(maxn), parts%mons%y%r(maxn))
!parts%mons%x%r = 0.d0; parts%mons%y%r = 0.d0
allocate(parts%rho); allocate(parts%rho%r(maxn)); parts%rho%r = 0.d0
!allocate(parts%rhos); allocate(parts%rhos%r(maxn)); parts%rhos%r = 0.d0
allocate(parts%divvx); allocate(parts%divvx%r(maxn));! parts%rhos%r = 0.d0
allocate(parts%p); allocate(parts%p%r(maxn)); parts%p%r = 0.d0
!allocate(parts%ps); allocate(parts%ps%r(maxn)); parts%ps%r = 0.d0
!allocate(parts%u(maxn));       parts%u   = 0.d0
allocate(parts%c); allocate(parts%c%r(maxn))
parts%c   = 0.d0
!allocate(parts%eta(maxn));     parts%eta = 0.d0 

! Old value
allocate(parts%v_min,parts%v_min%x,parts%v_min%y)
allocate(parts%v_min%x%r(maxn),parts%v_min%y%r(maxn))
parts%v_min%x = 0.d0; parts%v_min%y = 0.d0
allocate(parts%rho_min); allocate(parts%rho_min%r(maxn));  parts%rho_min = 0.d0
!allocate(parts%u_min); allocate(parts%u_min%r(maxn));  parts%u_min   = 0.d0

! Accelerations
allocate(parts%drho);allocate(parts%drho%r(maxn));    parts%drho= 0.d0
allocate(parts%dvx); allocate(parts%dvx%x,parts%dvx%y)
allocate(parts%dvx%x%r(maxn)); parts%dvx%x = 0.d0  
allocate(parts%dvx%y%r(maxn)); parts%dvx%y = 0.d0  
!allocate(parts%du(maxn));      parts%du  = 0.d0 
allocate(parts%av); allocate(parts%av%x,parts%av%y)
allocate(parts%av%x%r(maxn)); parts%av%x  = 0.d0
allocate(parts%av%y%r(maxn)); parts%av%y  = 0.d0


! Volume Fraction
allocate(parts%vof); allocate(parts%vof%r(maxn))
parts%vof%r = 0.d0

if(trim(parts%imaterial)=='water')then
   allocate(parts%vof2); allocate(parts%vof2%r(maxn))          ! phi_f = 1- phi_s
   parts%vof2 = 0.d0
   allocate(parts%dvof); allocate(parts%dvof%r(maxn))
   parts%dvof = 0.d0
   allocate(parts%vof_min); allocate(parts%vof_min%r(maxn))
   parts%vof_min = 0.d0
endif

if(trim(parts%imaterial)=='soil')then
    allocate(parts%str)
    allocate(parts%str%x,parts%str%xy,parts%str%y)
!   allocate(parts%sxx,parts%sxy,parts%syy)
   allocate(parts%str%x%r(maxn), parts%str%y%r(maxn), parts%str%xy%r(maxn))
   parts%str%x = 0.d0; parts%str%xy = 0.d0; parts%str%y = 0.d0
   allocate(parts%str_min)
   allocate(parts%str_min%x,parts%str_min%y,parts%str_min%xy)
   allocate(parts%str_min%x%r(maxn),parts%str_min%y%r(maxn),parts%str_min%xy%r(maxn))
   allocate(parts%p_min); allocate(parts%p_min%r(maxn))
   allocate(parts%vcc); allocate(parts%vcc%r(maxn))    
   !allocate(parts%szz(maxn))    !!! Intel Fortran Compiler is shit!!!
   parts%str_min%x = 0.d0; parts%str_min%xy = 0.d0; parts%str_min%y = 0.d0
   parts%p_min = 0.d0
   parts%vcc = 0.d0
   allocate(parts%dstr); allocate(parts%dstr%x, parts%dstr%y, parts%dstr%xy)
   allocate(parts%dstr%x%r(maxn),parts%dstr%y%r(maxn),parts%dstr%xy%r(maxn))
   parts%dstr%x = 0.d0; parts%dstr%xy = 0.d0; parts%dstr%y = 0.d0
   allocate(parts%dp); allocate(parts%dp%r(maxn))
   parts%dp%r(maxn) = 0.d0
! For return mapping algorithm
   allocate(parts%dstr2,parts%dstr2%x,parts%dstr2%y,parts%dstr2%xy)
   allocate(parts%dstr2%x%r(maxn),parts%dstr2%y%r(maxn),parts%dstr2%xy%r(maxn))
   parts%dstr2%x = 0.d0; parts%dstr2%xy = 0.d0; parts%dstr2%y = 0.d0
   allocate(parts%dp2); allocate(parts%dp2%r(maxn))
   parts%dp2%r(maxn) = 0.d0

   allocate(parts%fail(maxn)); parts%fail = 0

   allocate(parts%epsilon_p); allocate(parts%epsilon_p%r(maxn))
   parts%epsilon_p = 0.d0
endif

!allocate(parts%drhodx); allocate(parts%drhodx%r(maxn))
!parts%drhodx = 0.d0

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
type(p2r) vxi(3), dvxi(3), v_mini(3)

!      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            pl%rho_min = pl%rho
            pl%rho = pl%rho +(dt/2.)* pl%drho
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof_min = pl%vof
            pl%vof = pl%vof+(dt/2.)*pl%dvof 
         endif
 
         if(trim(pl%imaterial)=='soil')then
                  if(stress_integration==1)then
         !pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)   ! Nothing to do!
         !pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
         !pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
                  elseif(stress_integration==2)then
         pl%str_min%x = pl%str%x
         pl%str_min%xy = pl%str%xy
         pl%str_min%y = pl%str%y
         pl%p_min   = pl%p
         pl%str%x = pl%str%x + (dt/2.)*pl%dstr%x
         pl%str%xy = pl%str%xy + (dt/2.)*pl%dstr%xy
         pl%str%y = pl%str%y + (dt/2.)*pl%dstr%y
         !pl%str = pl%str + (dt/2.) * pl%dstr
         pl%p = pl%p + (dt/2.)*pl%dp      !!!simultaneous pressure
                  endif
         endif
          
      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i); v_mini = pl%v_min%cmpt(i)
         do d = 1, pl%dim
         !   pl%v_min(d, i) = pl%vx(d, i)
         !   pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
            v_mini(d)%p = vxi(d)%p
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
      type(p2r) vxi(3), dvxi(3), v_mini(3), avi(3)

!      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if (.not.summation_density ) then
            pl%rho = pl%rho + (dt/2.)* pl%drho
         endif

         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof = pl%vof+(dt/2.)*pl%dvof
         endif

         if(trim(pl%imaterial)=='soil')then
                  if(stress_integration==1)then
         pl%str%x = pl%str%x + dt*pl%dstr%x
         pl%str%xy = pl%str%xy + dt*pl%dstr%xy
         pl%str%y = pl%str%y + dt*pl%dstr%y
         pl%p   = pl%p   + dt*pl%dp       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%str%x = pl%str%x + (dt/2.)*pl%dstr%x
         pl%str%xy = pl%str%xy + (dt/2.)*pl%dstr%xy
         pl%str%y = pl%str%y + (dt/2.)*pl%dstr%y
         pl%p   = pl%p + (dt/2.)*pl%dp    !!! simultaneous pressure
                  endif
         endif
        
      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i); avi = pl%av%cmpt(i)
         do d = 1, pl%dim        
            !pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = vxi(d)%p + (dt/2.) * dvxi(d)%p + avi(d)%p
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
      type(p2r) vxi(3), dvxi(3), v_mini(3),avi(3)
      !$omp parallel do private(d,vxi,dvxi,v_mini,avi)
      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if (.not.summation_density ) then 
            pl%rho%r(i) = pl%rho_min%r(i) + dt*pl%drho%r(i)
         endif
   
         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
            pl%vof%r(i) = pl%vof_min%r(i)+dt*pl%dvof%r(i) 
         endif
                 
         if(trim(pl%imaterial)=='soil')then 
                  if(stress_integration==1)then 
         pl%str%x%r(i) = pl%str%x%r(i) + dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp%r(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%str%x%r(i) = pl%str_min%x%r(i) + dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str_min%xy%r(i) + dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str_min%y%r(i) + dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p_min%r(i) + dt*pl%dp%r(i)     !!! simultaneous pressure
                  endif
         endif
         
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i) ; v_mini = pl%v_min%cmpt(i)
         avi = pl%av%cmpt(i)
         do d = 1, pl%dim                   
            !pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = v_mini(d)%p + dt * dvxi(d)%p   &
                        + avi(d)%p
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
         pl%str%x%r(i) = pl%str%x%r(i) + dt*pl%dstr2%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + dt*pl%dstr2%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + dt*pl%dstr2%y%r(i)
         pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp2%r(i)       !!! simultaneous pressure
                  elseif(stress_integration==2)then
         pl%str%x%r(i) = pl%str_min%x%r(i) + dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str_min%xy%r(i) + dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str_min%y%r(i) + dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p_min%r(i) + dt*pl%dp%r(i)     !!! simultaneous pressure
                  endif
         endif
         !endif
         
      enddo

      return
      end subroutine

      end subroutine

!----------------------------------------------------------------------      
      subroutine time_integration_for_water
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it
type(particles), pointer :: pl
type(p2r) vxi(3), dvxi(3), v_mini(3)
              
call parts%setup_ndim1
pl => parts
do it = 1, maxtimestep    
   itimestep = itimestep+1
   parts%itimestep = itimestep


        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  
   
!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

   if(itimestep .ne. 1)then
pl => parts
      call first_half
!      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
                       
   endif
pl => parts
   call single_step_for_water
!    call single_step

   if(itimestep .eq. 1) then
pl => parts
      call first_step

   else
pl => parts
       call second_half 
       
       

   endif 

   time = time + dt

!  if(mod(itimestep,print_step).eq.0)then
!      write(*,*)'______________________________________________'
!      write(*,*)'  current number of time step =',                &
!                itimestep,'     current time=', real(time+dt)
!      write(*,*)'______________________________________________'
!   endif  

!   if(itimestep>=save_step_from.and.mod(itimestep,save_step).eq.0)then
!      call output
!   endif 

enddo

return

! -------------------
      contains
! -------------------

      subroutine first_half
! ---------------------------------------------------------------------
implicit none
integer i,d
type(p2r) vxi(3), dvxi(3), v_mini(3)

!      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.summation_density) then    
            parts%rho_min = parts%rho
            parts%rho = parts%rho +(dt/2.)* parts%drho
         endif

!         if(trim(pl%imaterial)=='water'.and.volume_fraction)then
!            pl%vof_min = pl%vof
!            pl%vof = pl%vof+(dt/2.)*pl%dvof 
!         endif
          
      do i = 1, parts%ntotal +parts%nvirt    ! originally only pl%ntotal       
         if(parts%itype(i)<0)cycle
        
         vxi = parts%vx%cmpt(i); dvxi = parts%dvx%cmpt(i); v_mini = parts%v_min%cmpt(i)
         do d = 1, parts%dim
         !   pl%v_min(d, i) = pl%vx(d, i)
         !   pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
            v_mini(d)%p = vxi(d)%p
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
      type(p2r) vxi(3), dvxi(3), v_mini(3), avi(3)

!      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if (.not.summation_density ) then
            parts%rho = parts%rho + (dt/2.)* parts%drho
         endif

      do i=1,parts%ntotal +parts%nvirt     ! origionally pl%ntotal
         if(parts%itype(i)<0)cycle
         vxi = parts%vx%cmpt(i); dvxi = parts%dvx%cmpt(i); avi = parts%av%cmpt(i)
         do d = 1, parts%dim        
            !pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = vxi(d)%p + (dt/2.) * dvxi(d)%p + avi(d)%p
            parts%x(d, i) = parts%x(d, i) + dt * vxi(d)%p            
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
      type(p2r) vxi(3), dvxi(3), v_mini(3),avi(3)
      !$omp parallel do private(d,vxi,dvxi,v_mini,avi)
      do i=1,parts%ntotal +parts%nvirt  ! origionally pl%ntotal            
            
         if (.not.summation_density ) then 
            parts%rho%r(i) = parts%rho_min%r(i) + dt*parts%drho%r(i)
         endif

 !        if(trim(pl%imaterial)=='water'.and.volume_fraction)then
 !           pl%vof%r(i) = pl%vof_min%r(i)+dt*pl%dvof%r(i) 
 !        endif
                 


         if(parts%itype(i)<0)cycle
         vxi = parts%vx%cmpt(i); dvxi = parts%dvx%cmpt(i) ; v_mini = parts%v_min%cmpt(i)
         avi = parts%av%cmpt(i)
         do d = 1, parts%dim                   
            !pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = v_mini(d)%p + dt * dvxi(d)%p   &
                        + avi(d)%p
            parts%x(d, i) = parts%x(d, i) + dt * vxi(d)%p
         enddo
            !pl%vx(1, i) = pl%v_min(1, i) + dt * pl%dvx%x%r(i)   &
            !            + pl%av(1, i)
            !pl%x(1, i) = pl%x(1, i) + dt * pl%vx(1, i)                  
            !pl%vx(2, i) = pl%v_min(2, i) + dt * pl%dvx%y%r(i)   &
            !            + pl%av(2, i)
            !pl%x(2, i) = pl%x(2, i) + dt * pl%vx(2, i)                              
      enddo
      !$omp end parallel do
      return
      end subroutine
      end subroutine
      
!在verlet算法中添加了XSPH,但不知道对不对。如何保证只用在密度和运动方程而不用在动量方程呢？
!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_verlet2
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
type(particles), pointer :: pl
type(array) :: lastvx_x,lastvx_y,temp1_x,temp1_y,lastrho,temp2

lastvx_x%ndim1 = parts%ntotal+parts%nvirt
lastvx_y%ndim1 = parts%ntotal+parts%nvirt
temp1_x%ndim1 = parts%ntotal+parts%nvirt
temp1_y%ndim1 = parts%ntotal+parts%nvirt
lastrho%ndim1 = parts%ntotal+parts%nvirt
temp2%ndim1 = parts%ntotal+parts%nvirt
allocate(lastvx_x%r(parts%ntotal+parts%nvirt))
allocate(lastvx_y%r(parts%ntotal+parts%nvirt))
allocate(temp1_x%r(parts%ntotal+parts%nvirt))
allocate(temp1_y%r(parts%ntotal+parts%nvirt))
allocate(lastrho%r(parts%ntotal+parts%nvirt))
allocate(temp2%r(parts%ntotal+parts%nvirt))
    pl => parts
    call parts%setup_ndim1
    lastvx_x = pl%vx%x
    lastvx_y = pl%vx%y
    temp1_x = pl%vx%x
    temp1_y = pl%vx%y
if(itimestep .eq. 1) lastrho = pl%rho

do it = 1, maxtimestep 
    itimestep = itimestep+1
   call single_step_for_water

if(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal  
           pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i) + (dt**2./2.)*pl%dvx%x%r(i)
           lastvx_x%r(i) = lastvx_x%r(i) + parts%av%x%r(i) + 2.*dt*pl%dvx%x%r(i)
           pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i) + (dt**2./2.)*pl%dvx%y%r(i)
           lastvx_y%r(i) = lastvx_y%r(i) + parts%av%y%r(i) + 2.*dt*pl%dvx%y%r(i)
    enddo
else
    do i = 1, pl%ntotal
           pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i) + (dt**2./2.)*pl%dvx%x%r(i)
           lastvx_x%r(i) = pl%vx%x%r(i) + parts%av%x%r(i) + dt*pl%dvx%x%r(i)
           pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i) + (dt**2./2.)*pl%dvx%y%r(i)
           lastvx_y%r(i) = pl%vx%y%r(i) + parts%av%y%r(i) + dt*pl%dvx%y%r(i)
    enddo
endif
temp1_x = pl%vx%x;temp1_y = pl%vx%y
pl%vx%x = lastvx_x;pl%vx%y = lastvx_y
lastvx_x = temp1_x;lastvx_y = temp1_y

if(itimestep .eq. 1) then
    do i = 1, pl%ntotal+pl%nvirt
          lastrho%r(i) = pl%rho%r(i) + dt * 2. * pl%drho%r(i)
    enddo
elseif(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal+pl%nvirt
          lastrho%r(i) = lastrho%r(i) + dt * 2. * pl%drho%r(i)
    enddo
else
   do i = 1, pl%ntotal+pl%nvirt
          lastrho%r(i) = pl%rho%r(i) + dt* pl%drho%r(i)
   enddo
endif
temp2 = pl%rho
pl%rho = lastrho
lastrho = temp2

   time = time + dt

enddo



return
end subroutine      

      
!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_verlet
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
type(particles), pointer :: pl
type(array) :: lastvx_x,lastvx_y,temp1_x,temp1_y,lastrho,temp2

lastvx_x%ndim1 = parts%ntotal+parts%nvirt
lastvx_y%ndim1 = parts%ntotal+parts%nvirt
temp1_x%ndim1 = parts%ntotal+parts%nvirt
temp1_y%ndim1 = parts%ntotal+parts%nvirt
lastrho%ndim1 = parts%ntotal+parts%nvirt
temp2%ndim1 = parts%ntotal+parts%nvirt
allocate(lastvx_x%r(parts%ntotal+parts%nvirt))
allocate(lastvx_y%r(parts%ntotal+parts%nvirt))
allocate(temp1_x%r(parts%ntotal+parts%nvirt))
allocate(temp1_y%r(parts%ntotal+parts%nvirt))
allocate(lastrho%r(parts%ntotal+parts%nvirt))
allocate(temp2%r(parts%ntotal+parts%nvirt))
    pl => parts
    call parts%setup_ndim1
    lastvx_x = pl%vx%x
    lastvx_y = pl%vx%y
    temp1_x = pl%vx%x
    temp1_y = pl%vx%y
if(itimestep .eq. 1) lastrho = pl%rho

do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water

if(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal
       !do d = 1, pl%dim
       !   pl%x(d,i) = pl%x(d,i) + dt * pl%vx(d,i) +(dt**2./2.)*pl%dvx(d,i)
       !   lastvx(d,i) = lastvx(d,i) + 2.*dt*pl%dvx(d,i)
       
       !enddo
!          pl%x(1,i) = pl%x(1,i) + dt * pl%vx(1,i) +(dt**2./2.)*pl%dvx%x%r(i)
!          lastvx(1,i) = lastvx(1,i) + 2.*dt*pl%dvx%x%r(i)
!          pl%x(2,i) = pl%x(2,i) + dt * pl%vx(2,i) +(dt**2./2.)*pl%dvx%y%r(i)
!          lastvx(2,i) = lastvx(2,i) + 2.*dt*pl%dvx%y%r(i)       
           pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i) + (dt**2./2.)*pl%dvx%x%r(i)
           lastvx_x%r(i) = lastvx_x%r(i) + 2.*dt*pl%dvx%x%r(i)
           pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i) + (dt**2./2.)*pl%dvx%y%r(i)
           lastvx_y%r(i) = lastvx_y%r(i) + 2.*dt*pl%dvx%y%r(i)
    enddo
else
    do i = 1, pl%ntotal
       !do d = 1, pl%dim
       !   pl%x(d,i) = pl%x(d,i) + dt * pl%vx(d,i) +(dt**2./2.)*pl%dvx(d,i)
       !   lastvx(d,i) = pl%vx(d,i) + dt*pl%dvx(d,i)
       !enddo
!          pl%x(1,i) = pl%x(1,i) + dt * pl%vx(1,i) +(dt**2./2.)*pl%dvx%x%r(i)
!          lastvx(1,i) = pl%vx(1,i) + dt*pl%dvx%x%r(i)
!          pl%x(2,i) = pl%x(2,i) + dt * pl%vx(2,i) +(dt**2./2.)*pl%dvx%y%r(i)
!          lastvx(2,i) = pl%vx(2,i) + dt*pl%dvx%y%r(i)       
           pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i) + (dt**2./2.)*pl%dvx%x%r(i)
           lastvx_x%r(i) = pl%vx%x%r(i) + dt*pl%dvx%x%r(i)
           pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i) + (dt**2./2.)*pl%dvx%y%r(i)
           lastvx_y%r(i) = pl%vx%y%r(i) + dt*pl%dvx%y%r(i)
    enddo
endif

temp1_x = pl%vx%x;temp1_y = pl%vx%y
pl%vx%x = lastvx_x;pl%vx%y = lastvx_y
lastvx_x = temp1_x;lastvx_y = temp1_y

if(itimestep .eq. 1) then
    do i = 1, pl%ntotal+pl%nvirt
          lastrho%r(i) = pl%rho%r(i) + dt * 2. * pl%drho%r(i)
    enddo
elseif(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal+pl%nvirt
          lastrho%r(i) = lastrho%r(i) + dt * 2. * pl%drho%r(i)
    enddo
else
   do i = 1, pl%ntotal+pl%nvirt
          lastrho%r(i) = pl%rho%r(i) + dt* pl%drho%r(i)
   enddo
endif
temp2 = pl%rho
pl%rho = lastrho
lastrho = temp2

   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

return
      end subroutine      


!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_unified
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
real(dp) dvx(3) 
type(particles), pointer :: pl
type(p2r) vx_i(3), vx_j(3) 



    pl => parts
    call parts%setup_ndim1
do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water
do i = 1,pl%ntotal
pl%vx%x%r(i) = pl%vx%x%r(i) + dt * pl%dvx%x%r(i)
pl%vx%y%r(i) = pl%vx%y%r(i) + dt * pl%dvx%y%r(i)
pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
enddo

do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
    do d=1,parts%dim
      dvx(d) = vx_i(d)%p - vx_j(d)%p
    enddo
    parts%drho%r(i) = parts%drho%r(i) + parts%mass%r(j)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))!书上的（4.34）
    parts%drho%r(j) = parts%drho%r(j) + parts%mass%r(i)*(dvx(1)*parts%dwdx(1,k)+dvx(2)*parts%dwdx(2,k))
    if(parts%itype(i)*parts%itype(j)>0)cycle
    parts%drho%r(i) = parts%drho%r(i) - parts%rho%r(i)*(parts%dgu(1,k) * vx_i(1)%p + parts%dgu(2,k) *vx_i(2)%p)
enddo

do i = 1,pl%ntotal
    pl%rho%r(i) = pl%rho%r(i) + dt * (pl%drho%r(i)/parts%gammaa%r(i)+10**-8)
enddo


   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

return
end subroutine      


!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_unified2
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab,rr0,dd,p1,p2,f,ii
type(particles), pointer :: pl
type(p2r) vx_i(3), vx_j(3) 

    pl => parts
    call parts%setup_ndim1
do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water2
do i = 1,pl%ntotal
pl%vx%x%r(i) = pl%vx%x%r(i) + dt * pl%dvx%x%r(i)
pl%vx%y%r(i) = pl%vx%y%r(i) + dt * pl%dvx%y%r(i)
!pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
!pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
enddo

do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    rr = 0.
    c0 = 10.* sqrt(9.81*0.5)  !因CASE不同而不同·
    ca = c0 * (parts%rho%r(i)/1000.)**3
    cb = c0 * (parts%rho%r(j)/1000.)**3
    cab = max(ca,cb)
    rhoab = parts%rho%r(i) - parts%rho%r(j) - 9.81*1000/c0**2*(parts%x(2,j)-parts%x(2,i))!这里很有可能是加号.
    vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
    do d=1,parts%dim
      dvx(d) = vx_i(d)%p - vx_j(d)%p
      dx(d) = parts%x(d,i) - parts%x(d,j)
      rr = rr + dx(d)**2
    enddo
    parts%drho%r(i) = parts%drho%r(i) + parts%rho%r(i)*parts%mass%r(j)/parts%rho%r(j)*((dvx(1) + cab*dx(1)/sqrt(rr)*rhoab/parts%rho%r(i))  &
                     * parts%dwdx(1,k)+ (dvx(2) + cab*dx(2)/sqrt(rr)*rhoab/parts%rho%r(i))*parts%dwdx(2,k))
    parts%drho%r(j) = parts%drho%r(j) + parts%rho%r(j)*parts%mass%r(i)/parts%rho%r(i)*((-dvx(1) + cab*(-dx(1))/sqrt(rr)*(-rhoab)/parts%rho%r(j))  &
                     * (-parts%dwdx(1,k))+((-dvx(2)) + cab*(-dx(2))/sqrt(rr)*(-rhoab)/parts%rho%r(j))*(-parts%dwdx(2,k)))
    if(parts%itype(i)*parts%itype(j)>0)cycle
    parts%drho%r(i) = parts%drho%r(i) - parts%rho%r(i)*(parts%dgu(1,k) * vx_i(1)%p + parts%dgu(2,k) *vx_i(2)%p)
 !下面是为了给靠近边界的粒子一个镜像，优点是不会穿透出去，也不会违反能量守恒，但效果不规律，压强不规律
    if(j==1160.or.j==1169.or.j==1178.or.j==1218)then
     rr0 = parts%dspp; dd = parts%numeric%dd
     p1 = parts%numeric%p1; p2 = parts%numeric%p2
!    if(sqrt(rr)<(parts%hsml(1)*0.45))then
!        pl%vx%x%r(i) = -pl%vx%x%r(i)
!        pl%vx%y%r(i) = -pl%vx%y%r(i)
!    endif
!    endif
!    if(sqrt(rr)<(parts%hsml(1)*0.45))then
!      if(parts%n(1,j)==0) pl%vx%y%r(i) = -pl%vx%y%r(i)
!      if(parts%n(2,j)==0) pl%vx%x%r(i) = -pl%vx%x%r(i)
!    endif
      rr = sqrt(rr)
      !if(rr.lt.rr0)then
      if(rr.gt.rr0)cycle
         f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
         ii = i
         if(parts%itype(i)<0)ii=j 
         parts%dvx%x%r(ii) = parts%dvx%x%r(ii) + dd * dx(1)*f
         parts%dvx%y%r(ii) = parts%dvx%y%r(ii) + dd * dx(2)*f
endif


enddo

do i = 1,pl%ntotal
    pl%rho%r(i) = pl%rho%r(i) + dt * (pl%drho%r(i)/(parts%gammaa%r(i)+10**-8))
    pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
    pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
enddo

   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

return
end subroutine      


!-------------------------------------------------------------------        ---      
      subroutine time_integration_for_water_by_unified2s  !也在这里把连续性方程求解改成能量守恒的形式
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab
type(particles), pointer :: pl
type(p2r) vx_i(3), vx_j(3) 

    pl => parts
    call parts%setup_ndim1
do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water2
do i = 1,pl%ntotal
if(parts%zone(i)==7) then
pl%vx%x%r(i) = 0.
pl%dvx%x%r(i) = 0.
else
pl%vx%x%r(i) = pl%vx%x%r(i) + dt * pl%dvx%x%r(i)
pl%vx%y%r(i) = pl%vx%y%r(i) + dt * pl%dvx%y%r(i)
endif
pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
enddo


do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    rr = 0.
    c0 = 25.* sqrt(9.81*0.5)
    ca = c0 * (parts%rho%r(i)/1000.)**3
    cb = c0 * (parts%rho%r(j)/1000.)**3
    cab = max(ca,cb)
    rhoab = parts%rho%r(i) - parts%rho%r(j) - 9.81*1000/c0**2*(parts%x(2,j)-parts%x(2,i))!这里很有可能是加号.
    vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
    do d=1,parts%dim
      dvx(d) = vx_i(d)%p - vx_j(d)%p
      dx(d) = parts%x(d,i) - parts%x(d,j)
      rr = rr + dx(d)**2
    enddo
    parts%drho%r(i) = parts%drho%r(i) + parts%mass%r(j)*((dvx(1) + cab*dx(1)/sqrt(rr)*rhoab/parts%rho%r(i))  &
                     * parts%dwdx(1,k)+ (dvx(2) + cab*dx(2)/sqrt(rr)*rhoab/parts%rho%r(i))*parts%dwdx(2,k))
    parts%drho%r(j) = parts%drho%r(j) + parts%mass%r(i)*((-dvx(1) + cab*(-dx(1))/sqrt(rr)*(-rhoab)/parts%rho%r(j))  &
                     * (-parts%dwdx(1,k))+((-dvx(2)) + cab*(-dx(2))/sqrt(rr)*(-rhoab)/parts%rho%r(j))*(-parts%dwdx(2,k)))
    if(parts%itype(i)*parts%itype(j)>0)cycle
    parts%drho%r(i) = parts%drho%r(i) - parts%rhos(k)*(parts%dgu(1,k) * vx_i(1)%p + parts%dgu(2,k) *vx_i(2)%p)
enddo

do i = 1,pl%ntotal
    pl%rho%r(i) = pl%rho%r(i) + dt * (pl%drho%r(i)/parts%gammaa%r(i))
enddo

   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

return
end subroutine      

!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_unified3
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab
type(particles), pointer :: pl
type(array) :: temp
type(p2r) vx_i(3), vx_j(3) 

temp%ndim1 = parts%ntotal+parts%nvirt
allocate(temp%r(parts%ntotal+parts%nvirt))
    pl => parts
    call parts%setup_ndim1
do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water3
   
!do i = 1,pl%ntotal
!    if(pl%x(1,i)<0.014) then
!        parts%dvx%x%r(i)= 0.
!        parts%vx%x%r(i)= 0.
!        parts%dvx%y%r(i)=0.
!        parts%vx%y%r(i)=0.
!    endif
!enddo
   
do i = 1,pl%ntotal
pl%vx%x%r(i) = pl%vx%x%r(i) + dt * pl%dvx%x%r(i)
pl%vx%y%r(i) = pl%vx%y%r(i) + dt * pl%dvx%y%r(i)
pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
!if(pl%x(2,i)>1.5) then
!    pl%x(2,i)=1.5
!endif
enddo
!temp = 0.
!do k=1,parts%niac
!    i = parts%pair_i(k)
!    j = parts%pair_j(k)
!    temp%r(i) = temp%r(i) + parts%mass%r(j) * parts%w(k)
!    temp%r(j) = temp%r(j) + parts%mass%r(i) * parts%w(k)
!enddo

!do i = 1,pl%ntotal
!    pl%rho%r(i) = pl%rho%r(i) * pl%gammaa%r(i) - temp%r(i)
!enddo


   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

close(1000)

return
      end subroutine      
      

!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_unified_verlet
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
type(particles), pointer :: pl
type(array) :: lastvx_x,lastvx_y,temp1_x,temp1_y,lastrho,temp2

lastvx_x%ndim1 = parts%ntotal+parts%nvirt
lastvx_y%ndim1 = parts%ntotal+parts%nvirt
temp1_x%ndim1 = parts%ntotal+parts%nvirt
temp1_y%ndim1 = parts%ntotal+parts%nvirt
temp2%ndim1 = parts%ntotal+parts%nvirt
allocate(lastvx_x%r(parts%ntotal+parts%nvirt))
allocate(lastvx_y%r(parts%ntotal+parts%nvirt))
allocate(temp1_x%r(parts%ntotal+parts%nvirt))
allocate(temp1_y%r(parts%ntotal+parts%nvirt))
allocate(temp2%r(parts%ntotal+parts%nvirt))
    pl => parts
    call parts%setup_ndim1
    lastvx_x = pl%vx%x
    lastvx_y = pl%vx%y
    temp1_x = pl%vx%x
    temp1_y = pl%vx%y

do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water3

if(mod(itimestep,50) .ne. 0) then
    do i = 1, pl%ntotal  
           pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i) + (dt**2./2.)*pl%dvx%x%r(i)
           lastvx_x%r(i) = lastvx_x%r(i) + parts%av%x%r(i) + 2.*dt*pl%dvx%x%r(i)
           pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i) + (dt**2./2.)*pl%dvx%y%r(i)
           lastvx_y%r(i) = lastvx_y%r(i) + parts%av%y%r(i) + 2.*dt*pl%dvx%y%r(i)
    enddo
else
    do i = 1, pl%ntotal
           pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i) + (dt**2./2.)*pl%dvx%x%r(i)
           lastvx_x%r(i) = pl%vx%x%r(i) + parts%av%x%r(i) + dt*pl%dvx%x%r(i)
           pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i) + (dt**2./2.)*pl%dvx%y%r(i)
           lastvx_y%r(i) = pl%vx%y%r(i) + parts%av%y%r(i) + dt*pl%dvx%y%r(i)
    enddo
endif

temp1_x = pl%vx%x;temp1_y = pl%vx%y
pl%vx%x = lastvx_x;pl%vx%y = lastvx_y
lastvx_x = temp1_x;lastvx_y = temp1_y


   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

return
      end subroutine      

      
      
!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_unified_Couette
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab
type(particles), pointer :: pl
type(array) :: temp
type(p2r) vx_i(3), vx_j(3) 

temp%ndim1 = parts%ntotal+parts%nvirt
allocate(temp%r(parts%ntotal+parts%nvirt))
    pl => parts
    call parts%setup_ndim1
    
do i = 1,pl%ntotal + pl%nvirt
    if(parts%zone(i) == 1) pl%vx%x%r(i) = 0.000025
enddo

do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water_Couette
   
do i = 1,pl%ntotal
pl%vx%x%r(i) = pl%vx%x%r(i) + dt * pl%dvx%x%r(i)
!pl%vx%y%r(i) = pl%vx%y%r(i) + dt * pl%dvx%y%r(i)
pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
!pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
!if(pl%x(1,i)>0.001) pl%x(1,i)=pl%x(1,i) - 0.001
if(pl%x(1,i)>0.001025) pl%x(1,i)=pl%x(1,i) - 0.001025
enddo
!do i = pl%ntotal + 1,pl%ntotal + pl%nvirt
!  pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
!enddo


!temp = 0.
!do k=1,parts%niac
!    i = parts%pair_i(k)
!    j = parts%pair_j(k)
!    temp%r(i) = temp%r(i) + parts%mass%r(j) * parts%w(k)
!    temp%r(j) = temp%r(j) + parts%mass%r(i) * parts%w(k)
!enddo

!do i = 1,pl%ntotal
!    pl%rho%r(i) = pl%rho%r(i) * pl%gammaa%r(i) - temp%r(i)
!enddo


   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

close(1000)
close(1001)
close(1002)
return
      end subroutine     
      
!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_unified4  !把连续性方程改造一下，可以算是原创
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab
type(particles), pointer :: pl
type(p2r) vx_i(3), vx_j(3) 

    pl => parts
    call parts%setup_ndim1
do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water2
do i = 1,pl%ntotal
pl%vx%x%r(i) = pl%vx%x%r(i) + dt * pl%dvx%x%r(i)
pl%vx%y%r(i) = pl%vx%y%r(i) + dt * pl%dvx%y%r(i)
pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
enddo

do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    rr = 0.
    c0 = 10.* sqrt(9.81*0.6)
    ca = c0 * (parts%rho%r(i)/1000.)**3
    cb = c0 * (parts%rho%r(j)/1000.)**3
    cab = max(ca,cb)
    rhoab = parts%rho%r(i) - parts%rho%r(j) - 9.81*1000/c0**2*(parts%x(2,j)-parts%x(2,i))!这里很有可能是加号.
    vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
    do d=1,parts%dim
      dvx(d) = vx_i(d)%p - vx_j(d)%p
      dx(d) = parts%x(d,i) - parts%x(d,j)
      rr = rr + dx(d)**2
    enddo
    
    parts%drho%r(i) = parts%drho%r(i) + parts%mass%r(j)*((dvx(1) + cab*dx(1)/sqrt(rr)*rhoab/parts%rho%r(i))  &
                     * parts%dwdx(1,k)+ (dvx(2) + cab*dx(2)/sqrt(rr)*rhoab/parts%rho%r(i))*parts%dwdx(2,k))
    parts%drho%r(j) = parts%drho%r(j) + parts%mass%r(i)*((-dvx(1) + cab*(-dx(1))/sqrt(rr)*(-rhoab)/parts%rho%r(j))  &
                     * (-parts%dwdx(1,k))+((-dvx(2)) + cab*(-dx(2))/sqrt(rr)*(-rhoab)/parts%rho%r(j))*(-parts%dwdx(2,k)))
    if(parts%itype(i)*parts%itype(j)>0)cycle
    parts%drho%r(i) = parts%drho%r(i) - parts%rho%r(j)*(parts%dgu(1,k) * vx_i(1)%p + parts%dgu(2,k) *vx_i(2)%p)
enddo

do i = 1,pl%ntotal
    pl%rho%r(i) = pl%rho%r(i) + dt * (pl%drho%r(i)/(parts%gammaa%r(i)+10**-8))
enddo


   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

return
      end subroutine      

!----------------------------------------------------------------------      
      subroutine time_integration_for_water_by_unified5 !完全按照第一篇文章的来
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it 
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab
type(particles), pointer :: pl
type(array) :: temp
type(p2r) vx_i(3), vx_j(3) 

temp%ndim1 = parts%ntotal+parts%nvirt
allocate(temp%r(parts%ntotal+parts%nvirt))
    pl => parts
    call parts%setup_ndim1
do it = 1, maxtimestep 
    itimestep = itimestep+1
    parts%itimestep = itimestep
   call single_step_for_water5
do i = 1,pl%ntotal
pl%vx%x%r(i) = pl%vx%x%r(i) + dt * pl%dvx%x%r(i)
pl%vx%y%r(i) = pl%vx%y%r(i) + dt * pl%dvx%y%r(i)
pl%x(1,i) = pl%x(1,i) + dt * pl%vx%x%r(i)
pl%x(2,i) = pl%x(2,i) + dt * pl%vx%y%r(i)
enddo
!temp = 0.
!do k=1,parts%niac
!    i = parts%pair_i(k)
!    j = parts%pair_j(k)
!    temp%r(i) = temp%r(i) + parts%mass%r(j) * parts%w(k)
!    temp%r(j) = temp%r(j) + parts%mass%r(i) * parts%w(k)
!enddo

!do i = 1,pl%ntotal
!    pl%rho%r(i) = pl%rho%r(i) * pl%gammaa%r(i) - temp%r(i)
!enddo


   time = time + dt
   
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif  

   
enddo

return
      end subroutine      
      
      
!DEC$IF(.FALSE.)
!----------------------------------------------------------------------
             subroutine time_integration_for_soil
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it
type(particles), pointer :: pl
type(p2r) vxi(3), dvxi(3), v_mini(3)
              
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
            pl%rho_min%r(i) = pl%rho%r(i)
            pl%rho%r(i) = pl%rho%r(i) +(dt/2.)* pl%drho%r(i)
         endif
 
         if(stress_integration==1)then
            !pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)   ! Nothing to do!
            !pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
            !pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
         elseif(stress_integration==2)then
            pl%str_min%x%r(i) = pl%str%x%r(i)
            pl%str_min%xy%r(i) = pl%str%xy%r(i)
            pl%str_min%y%r(i) = pl%str%y%r(i)
            pl%p_min%r(i)   = pl%p%r(i)
            pl%str%x%r(i) = pl%str%x%r(i) + (dt/2.)*pl%dstr%x%r(i)
            pl%str%xy%r(i) = pl%str%xy%r(i) + (dt/2.)*pl%dstr%xy%r(i)
            pl%str%y%r(i) = pl%str%y%r(i) + (dt/2.)*pl%dstr%y%r(i)
            pl%p%r(i) = pl%p%r(i) + (dt/2.)*pl%dp%r(i)      !!!simultaneous pressure
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
            pl%rho%r(i) = pl%rho%r(i) + (dt/2.)* pl%drho%r(i)
         endif

         if(stress_integration==1)then
            pl%str%x%r(i) = pl%str%x%r(i) + dt*pl%dstr%x%r(i)
            pl%str%xy%r(i) = pl%str%xy%r(i) + dt*pl%dstr%xy%r(i)
            pl%str%y%r(i) = pl%str%y%r(i) + dt*pl%dstr%y%r(i)
            pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp%r(i)       !!! simultaneous pressure
         elseif(stress_integration==2)then
            pl%str%x%r(i) = pl%str%x%r(i) + (dt/2.)*pl%dstr%x%r(i)
            pl%str%xy%r(i) = pl%str%xy%r(i) + (dt/2.)*pl%dstr%xy%r(i)
            pl%str%y%r(i) = pl%str%y%r(i) + (dt/2.)*pl%dstr%y%r(i)
            pl%p%r(i)   = pl%p%r(i) + (dt/2.)*pl%dp%r(i)    !!! simultaneous pressure
         endif
        
         if(pl%itype(i)<0)cycle

         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i)
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
            pl%rho%r(i) = pl%rho_min%r(i) + dt*pl%drho%r(i)
         endif
                 
         if(stress_integration==1)then 
            pl%str%x%r(i) = pl%str%x%r(i) + dt*pl%dstr%x%r(i)
            pl%str%xy%r(i) = pl%str%xy%r(i) + dt*pl%dstr%xy%r(i)
            pl%str%y%r(i) = pl%str%y%r(i) + dt*pl%dstr%y%r(i)
            pl%p%r(i)   = pl%p%r(i)   + dt*pl%dp%r(i)       !!! simultaneous pressure
         elseif(stress_integration==2)then
            pl%str%x%r(i) = pl%str_min%x%r(i) + dt*pl%dstr%x%r(i)
            pl%str%xy%r(i) = pl%str%xy_min(i) + dt*pl%dstr%xy%r(i)
            pl%str%y%r(i) = pl%str%yy_min(i) + dt*pl%dstr%y%r(i)
            pl%p%r(i)   = pl%p_min%r(i) + dt*pl%dp%r(i)     !!! simultaneous pressure
         endif

         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i); v_mini = pl%vx%cmpt(i)
         do d = 1, pl%dim                   
            vxi%p(d) = pl%v_mini%p(d) + dt * dvxi%p(d)   &
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

pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
if(trim(pl%imaterial)=='water')pl%dvof = 0.d0
if(trim(pl%imaterial)=='soil')then
   pl%dstr%x = 0.d0; pl%dstr%xy = 0.d0; pl%dstr%y = 0.d0; pl%dp = 0.d0
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
    !call con_density(pl)         
    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
      
if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      call delta_sph_omp(pl,pl%rho,pl%drho)
   !endif
endif

!---  Dynamic viscosity:

if(trim(pl%imaterial)=='water')then
   !if (visc) call viscosity(pl)
elseif(trim(pl%imaterial)=='soil')then
   !call shear_modulus(pl)
endif
       
!---  Internal forces:

!call shear_strain_rate(pl)
pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

if(trim(pl%imaterial)=='soil') pl%vcc = pl%div_omp(pl%vx)
!call velocity_divergence(pl)

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
   !if(pl%nthreads==1)then
!   pl%dvx%x%r = -pl%vof%r*pl%df(pl%p%r,'x') + pl%df(pl%vof%r*pl%str%x%r,'x') + pl%df(pl%vof%r*pl%str%xy%r,'y')
!   pl%dvx%y%r = -pl%vof%r*pl%df(pl%p%r,'y') + pl%df(pl%vof%r*pl%str%xy%r,'x') + pl%df(pl%vof%r*pl%str%y%r,'y')
   pl%dvx%x = -pl%vof*pl%df(pl%p,'x') + pl%df(pl%vof*pl%str%x,'x') + pl%df(pl%vof*pl%str%xy,'y')
   pl%dvx%y = -pl%vof*pl%df(pl%p,'y') + pl%df(pl%vof*pl%str%xy,'x') + pl%df(pl%vof*pl%str%y,'y')   
   !else
   !pl%dvx%x = -pl%vof*pl%df_omp(pl%p,'x') + pl%df_omp(pl%vof*pl%str%x,'x') + pl%df_omp(pl%vof*pl%str%xy,'y')
   !pl%dvx%y = -pl%vof*pl%df_omp(pl%p,'y') + pl%df_omp(pl%vof*pl%str%xy,'x') + pl%df_omp(pl%vof*pl%str%y,'y')  
   !endif        
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
   !where (pl%rho%r.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
   !where (pl%rho%r.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
else      
   !call int_force(pl)
   pl%dvx%x = -pl%df3_omp(pl%vof*pl%p,'x') + pl%df3_omp(pl%vof*pl%str%x,'x') + pl%df3_omp(pl%vof*pl%str%xy,'y')
   pl%dvx%y = -pl%df3_omp(pl%vof*pl%p,'y') + pl%df3_omp(pl%vof*pl%str%xy,'x') + pl%df3_omp(pl%vof*pl%str%y,'y')
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

if (visc_artificial) call pl%art_visc_omp

if(trim(pl%imaterial)=='soil'.and.soil_artificial_stress)then
        !call art_stress(pl)
   call pl%delta_sph_omp(pl%p,pl%dp)
   call pl%delta_sph_omp(pl%str%x,pl%dstr%x)
   call pl%delta_sph_omp(pl%str%xy,pl%dstr%xy)
   call pl%delta_sph_omp(pl%str%y,pl%dstr%y)
endif        
if(trim(pl%imaterial)=='water'.and.water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
        call pl%delta_sph_omp(pl%vof,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
          call pl%repulsive_force_omp
!      endif

pl%dvx%y = pl%dvx%y + gravity

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (average_velocity) call av_vel_omp(pl) 

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
!         call direct_find_2(parts,soil)
          call link_list2_omp(parts,soil)
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

      call darcy_law_omp(parts,soil)          
      call pore_water_pressure_omp(parts,soil) 

      if(volume_fraction)then
         call volume_fraction_soil_omp(soil)
         call volume_fraction_water2_omp(parts,soil)
         call volume_fraction_water_omp(parts,soil)  ! phi_f = 1- phi_s
         if(volume_fraction_renorm)then
            if(mod(itimestep,40).eq.0) then
               ntotal = parts%ntotal+parts%nvirt
               parts%rho = parts%rho/parts%vof 
               parts%vof = parts%vof2
               parts%rho = parts%rho*parts%vof 
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

!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
!use param 
!use declarations_sph
!use m_sph_fo
implicit none

integer  nphase
type(particles), pointer :: pl
logical :: dbg = .false.
integer i, ntotal
type(material),pointer :: water
water => parts%material
pl => parts        
if(dbg) write(*,*) 'In single_step...'
call pl%setup_ndim1
pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
!pl%dvof = 1.d0
 
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
!if(itimestep<62000) then
if(mod(itimestep,25)==0) then
call sum_density(pl)
else             
!    call sum_density(pl)         
    
    pl%drho = -pl.rho*pl.div2(pl.vx)
endif
!else
!if(mod(itimestep,25)==0) then
!call sum_density(pl)
!else             
!    call sum_density(pl)         
    
!    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
!endif
    
    
if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!      call delta_sph_omp(pl,pl%rho,pl%drho)
       call delta_rho(pl,pl%rho,pl%drho)
   !endif
endif

!---  Dynamic viscosity:
     
water => parts%material
parts%p = water%b*((parts%rho/(water%rho0))**water%gamma-1.d0)

!第二种状态方程
!parts%p = water%c**2*(parts%rho-water%rho0)

!call freesurface(pl) 
!parts%c%r(1:ntotal) = water%c*(parts%rho%r(1:ntotal)/(water%rho0))**3.0    



!---  Internal forces:

!call shear_strain_rate(pl)
pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!call velocity_divergence(pl)

!call pressure(pl)

!   call newtonian_fluid(pl)

      parts%str%x = water%viscosity*parts%tab%x
      parts%str%y = water%viscosity*parts%tab%y
      parts%str%xy = water%viscosity*parts%tab%xy
! write(*,*) parts%tab%x%r(862),parts%tab%y%r(861),parts%tab%y%r(882),parts%tab%y%r(883)
! write(*,*) parts%str%x%r(862),parts%str%y%r(861),parts%str%y%r(882),parts%str%y%r(883)
!Calculate internal force for water phase !! -phi_f Grad(p)

   pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
   pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
! write(*,*) parts%dvx%x%r(862),parts%dvx%y%r(862),parts%dvx%y%r(861),parts%dvx%y%r(861)
!if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!       call delta_sph_vx(pl)
   !endif
!endif
!---  Artificial viscosity:

!if (visc_artificial) call pl%art_visc
       
!if(trim(pl%imaterial)=='water'.and.water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
 !       call pl%delta_sph_omp(pl%vof,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
call pl%repulsive_force_omp                ! can be tried
!      endif

pl%dvx%y = pl%dvx%y + gravity
!write(*,*) parts%dvx%x%r(862),parts%dvx%y%r(862),parts%dvx%y%r(861),parts%dvx%y%r(861),parts%dvx%y%r(860),parts%dvx%y%r(860)
!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (average_velocity) call av_vel(pl) 


!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

      if(itimestep>=save_step_from.and.   &
         mod(itimestep,save_step).eq.0)then
         call output
         endif 

!      if(itimestep>=8650.and.itimestep<=8777)then
!         call output
!      endif 


         
return
end subroutine
!-------------------------------------------------
      subroutine single_step_for_water2
!-------------------------------------------------

!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
!use param 
!use declarations_sph
!use m_sph_fo
implicit none

integer  nphase
type(particles), pointer :: pl
real(dp)  dguaa(parts%dim,parts%ntotal)
logical :: dbg = .false.
integer i, ntotal,d,k,j
type(material),pointer :: water
water => parts%material
pl => parts        
if(dbg) write(*,*) 'In single_step...'
call pl%setup_ndim1
pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
!pl%dvof = 1.d0
 
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
!if(mod(itimestep,2)==0) then
!   call get_gammaa4(pl)
!else 
!if(itimestep==1) then
   call get_gammaa2(pl)
!endif
   !else
!do i = 1,parts%ntotal
!    do d = 1, parts%dim
!        dguaa(d,i) = 0.
!    enddo
!enddo
!do k = 1, parts%niac
!   i = parts%pair_i(k)
 !  j = parts%pair_j(k)
!   if(parts%itype(i)*parts%itype(j)>0) cycle   
!   dguaa(d,i) = dguaa(d,i) + parts%dgu(d,k)
!enddo
!!do i = 1,parts%ntotal
!   parts%gammaa%r(i) = parts%gammaa%r(i) + dt * (dguaa(1,i) * parts%vx%x%r(i) + dguaa(2,j) * parts%vx%y%r(j))
!enddo
!endif
!------unified get rho of nvirt particles
call nvirt_density_unified2(pl)       

!call pressure_nvirt(pl)
!______Analytical value of delta gamma
!if(mod(itimestep,2)==0) then
!  call delta_gamma_unified6(pl)
!else 
  call delta_gamma_unified3(pl)
!endif

call real_density_unified2(pl)

!---  Density approximation or change rate
!if(summation_density)then   
!if(itimestep<62000) then
!if(mod(itimestep,25)==0) then
!call sum_density(pl)
!else             
!    call sum_density(pl)         
    
!    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
!else
!if(mod(itimestep,25)==0) then
!call sum_density(pl)
!else             
!    call sum_density(pl)         
    
!    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
!endif
!do i = parts%ntotal +1,parts%ntotal+parts%nvirt
!    pl%drho%r(i) = 0
!enddo

     
if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      call delta_sph_omp(pl,pl%rho,pl%drho)
!       call delta_rho(pl,pl%rho,pl%drho)
   !endif
endif

!---  Dynamic viscosity:
     
!water => parts%material
!parts%p = water%b*((parts%rho/(water%rho0))**water%gamma-1.d0)

!down is for unified condition
water => parts%material
do i = 1,parts%ntotal 
parts%p%r(i) = water%b*((parts%rho%r(i)/(water%rho0))**water%gamma-1.d0)
enddo
!call momentum_equation_unified6(pl)  

!up is for unified condition
 
!call pressure_nvirt(pl)
!第二种状态方程
!parts%p = water%c**2*(parts%rho-water%rho0)

!call freesurface(pl) 
!parts%c%r(1:ntotal) = water%c*(parts%rho%r(1:ntotal)/(water%rho0))**3.0    


!---  Internal forces:

!call shear_strain_rate(pl)
!3~~~~~~~~~down is before
!3pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!3pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!3pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!3pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!3pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!call velocity_divergence(pl)

!call pressure(pl)

!   call newtonian_fluid(pl)

!3      parts%str%x = water%viscosity*parts%tab%x
!3      parts%str%y = water%viscosity*parts%tab%y
!3      parts%str%xy = water%viscosity*parts%tab%xy

!Calculate internal force for water phase !! -phi_f Grad(p)


!3   pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!3   pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   

!down is for dealt-sph
!   pl%dvx%x = pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!   pl%dvx%y = pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y') 
!   call pressure_nvirt(pl)
!   pl%dvx%x = pl%dvx%x - pl%df(pl%p,'x')
!   pl%dvx%y = pl%dvx%y - pl%df(pl%p,'y')
!up is for dealt-sph

!3   pl%dvx%x = pl%dvx%x/pl%rho
!3   pl%dvx%y = pl%dvx%y/pl%rho
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
!3~~~~~~~~~up is before
!if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!       call delta_sph_vx(pl)
   !endif
!endif
!---  Artificial viscosity:

!4~~~~not need now
!if (visc_artificial) call pl%art_visc
!4~~~~not need now       

!if(trim(pl%imaterial)=='water'.and.water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
 !       call pl%delta_sph_omp(pl%vof,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)

!5~~~~not need now
!call pl%repulsive_force_omp                ! can be tried
!      endif
!5~~~~not need now     

!6~~~~not need now
!pl%dvx%y = pl%dvx%y + gravity
!6~~~~not need now

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

!7~~~~not need now
!if (average_velocity) call av_vel(pl) 
!7~~~~not need now

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

      if(itimestep>=save_step_from.and.   &
         mod(itimestep,save_step).eq.0)then
         call output
         endif 

!      if(itimestep>=8650.and.itimestep<=8777)then
!         call output
!      endif 


         
return
      end subroutine
      
      
!-------------------------------------------------
      subroutine single_step_for_water3
!-------------------------------------------------
!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
!use param 
!use declarations_sph
!use m_sph_fo
implicit none

integer  nphase
type(particles), pointer :: pl
logical :: dbg = .false.
integer i, ntotal,k,j,d
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab,hv(3),r,selfdens
type(p2r) vx_i(3), vx_j(3) 
type(array) temp
type(material),pointer :: water

temp%ndim1 = parts%ntotal+parts%nvirt
allocate(temp%r(parts%ntotal+parts%nvirt))
water => parts%material
pl => parts        
if(dbg) write(*,*) 'In single_step...'
call pl%setup_ndim1
pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
!pl%dvof = 1.d0
 
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length

        
      hv = 0.d0
      r=0.d0


if(itimestep/=1)then
 do i=1,parts%ntotal
   call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
   temp%r(i)=-selfdens*parts%mass%r(i)  !由于rho不能要自加，这里也要自加
 enddo
do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    if(parts%itype(i)>0)    temp%r(i) = temp%r(i) - parts%mass%r(j) * parts%w(k)
    if(parts%itype(j)>0)    temp%r(j) = temp%r(j) - parts%mass%r(i) * parts%w(k)
enddo
endif

do i =1,parts%ntotal
    temp%r(i) = temp%r(i) + parts%gammaa%r(i) * parts%rho%r(i)
enddo

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

call get_gammaa2(pl)
if(itimestep/=1)then
    
 do i=1,parts%ntotal
   call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
   temp%r(i)=temp%r(i) + selfdens*parts%mass%r(i)  !由于rho不能要自加，这里也要自加
 enddo
 
do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    rr = 0.
    c0 = 30!.*sqrt(9.8*0.5)!因case不同而不同
    ca = c0 * (parts%rho%r(i)/1000.)**3
    cb = c0 * (parts%rho%r(j)/1000.)**3
    cab = max(ca,cb)
    rhoab = parts%rho%r(i) - parts%rho%r(j) - 9.81*1000/c0**2*(parts%x(2,j)-parts%x(2,i))!这里很有可能是加号.
    vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
    do d=1,parts%dim
      dvx(d) = vx_i(d)%p - vx_j(d)%p
      dx(d) = parts%x(d,i) - parts%x(d,j)
      rr = rr + dx(d)**2
    enddo

    if(parts%itype(i)>0)  temp%r(i) = temp%r(i) + parts%mass%r(j) * parts%w(k)
    if(parts%itype(j)>0)  temp%r(j) = temp%r(j) + parts%mass%r(i) * parts%w(k)
    if(parts%itype(i)>0.and.parts%itype(j)>0) then
        temp%r(i) = temp%r(i) + dt*(parts%mass%r(j)*cab/sqrt(rr)*rhoab/parts%rho%r(j)*(dx(1)*parts%dwdx(1,k) + dx(2)*parts%dwdx(2,k)))
        temp%r(j) = temp%r(j) + dt*(parts%mass%r(i)*cab/sqrt(rr)*(-rhoab)/parts%rho%r(i)*((-dx(1))*(-parts%dwdx(1,k)) + (-dx(2))*(-parts%dwdx(2,k))))        
    endif
enddo

do i = 1,pl%ntotal
    pl%rho%r(i) = temp%r(i)/parts%gammaa%r(i)
enddo
endif

!do i = 1,pl%ntotal
!    if(pl%x(1,i)<0.01) pl%rho%r(i) =1000
!enddo

    
!------unified get rho of nvirt particles


!call shear_strain_rate(pl)
!pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!      parts%str%x = water%viscosity*parts%tab%x
!      parts%str%y = water%viscosity*parts%tab%y
!      parts%str%xy = water%viscosity*parts%tab%xy


!if(itimestep>1)then
!    call nvirt_density_unified201(pl)
!else
!call nvirt_density_unified2(pl)    
!endif!这样就用了deltaP的方式来求了
!call pressure_nvirt_unified(pl)
!______Analytical value of delta gamma
call delta_gamma_unified2(pl)

call real_density_unified2(pl)
call nvirt_density_unified3(pl)    

!---  Density approximation or change rate
!if(summation_density)then   
!if(itimestep<62000) then
!if(mod(itimestep,25)==0) then
!call sum_density(pl)
!else             
!    call sum_density(pl)         
    
!    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
!else
!if(mod(itimestep,25)==0) then
!call sum_density(pl)
!else             
!    call sum_density(pl)         
    
!    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
!endif
!do i = parts%ntotal +1,parts%ntotal+parts%nvirt
!    pl%drho%r(i) = 0
!enddo

     
!`````````````````if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!`````````````````      call delta_sph_omp(pl,pl%rho,pl%drho)
!       call delta_rho(pl,pl%rho,pl%drho)
   !endif
!````````````````endif

!---  Dynamic viscosity:
     
!water => parts%material
!parts%p = water%b*((parts%rho/(water%rho0))**water%gamma-1.d0)

!down is for unified condition
water => parts%material
do i = 1,parts%ntotal 
parts%p%r(i) = water%b*((parts%rho%r(i)/(water%rho0))**water%gamma-1.d0)
enddo
!if(itimestep>1)then
!call momentum_equation_unified401(pl) 
!else
call momentum_equation_unified4half(pl)      
!endif!这样就用了deltaP的方式来求了


!pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!call velocity_divergence(pl)

!call pressure(pl)

!   call newtonian_fluid(pl)

!      parts%str%x = water%viscosity*parts%tab%x
!      parts%str%y = water%viscosity*parts%tab%y
!      parts%str%xy = water%viscosity*parts%tab%xy
! write(*,*) parts%tab%x%r(862),parts%tab%y%r(861),parts%tab%y%r(882),parts%tab%y%r(883)
! write(*,*) parts%str%x%r(862),parts%str%y%r(861),parts%str%y%r(882),parts%str%y%r(883)
!Calculate internal force for water phase !! -phi_f Grad(p)

!  pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!  pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   
!  pl%dvx%x = pl%dvx%x/pl%rho
!  pl%dvx%y = pl%dvx%y/pl%rho

!        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!          if(2187<i.and.i<2219)write(*,*) parts%p%r(i),i
!        enddo
!call freesurface2(pl)  !freesurface2通过使用divr来求解。通过比较后发现，这种形式在相应算例中优于freesurface

!     do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        parts%p%r(i)=parts%p%r(i)/wi(i)
!        if(2187<i.and.i<2219)write(*,*) parts%p%r(i),i
!     enddo
        

!up is for unified condition
 
!call pressure_nvirt(pl)
!第二种状态方程
!parts%p = water%c**2*(parts%rho-water%rho0)

!call freesurface(pl) 
!parts%c%r(1:ntotal) = water%c*(parts%rho%r(1:ntotal)/(water%rho0))**3.0    


!---  Internal forces:

!call shear_strain_rate(pl)
!3~~~~~~~~~down is before
!3pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!3pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!3pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!3pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!3pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!call velocity_divergence(pl)

!call pressure(pl)

!   call newtonian_fluid(pl)

!3      parts%str%x = water%viscosity*parts%tab%x
!3      parts%str%y = water%viscosity*parts%tab%y
!3      parts%str%xy = water%viscosity*parts%tab%xy

!Calculate internal force for water phase !! -phi_f Grad(p)


!3   pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!3   pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   

!down is for dealt-sph
!   pl%dvx%x = pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!   pl%dvx%y = pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y') 
!   call pressure_nvirt(pl)
!   pl%dvx%x = pl%dvx%x - pl%df(pl%p,'x')
!   pl%dvx%y = pl%dvx%y - pl%df(pl%p,'y')
!up is for dealt-sph

!3   pl%dvx%x = pl%dvx%x/pl%rho
!3   pl%dvx%y = pl%dvx%y/pl%rho
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
!3~~~~~~~~~up is before
!if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!       call delta_sph_vx(pl)
   !endif
!endif
!---  Artificial viscosity:

!4~~~~not need now
!if (visc_artificial) call pl%art_visc
!4~~~~not need now       

!if(trim(pl%imaterial)=='water'.and.water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
 !       call pl%delta_sph_omp(pl%vof,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)

!5~~~~not need now
!call pl%repulsive_force_omp                ! can be tried
!      endif
!5~~~~not need now     

!6~~~~not need now
!pl%dvx%y = pl%dvx%y + gravity
!6~~~~not need now

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration
!7
!7~~~~not need now
!if (average_velocity) call av_vel(pl) 
!7~~~~not need now
!parts%res_F = parts%p%r(5475) + parts%p%r(5476) + parts%p%r(5477) + parts%p%r(5478) + parts%p%r(5479) + parts%p%r(5480) + parts%p%r(5481) + parts%p%r(5482) +  &
!              parts%p%r(5483) + parts%p%r(5484) + parts%p%r(5485) + parts%p%r(5486) + parts%p%r(5487) + parts%p%r(5488) + parts%p%r(5489) + parts%p%r(5490) +  &
!              parts%p%r(5491) + parts%p%r(5492) + parts%p%r(5493) + parts%p%r(5494) + parts%p%r(5495) + parts%p%r(5496) + parts%p%r(5497) + parts%p%r(5498) +  &
!              parts%p%r(5499) + parts%p%r(5500) + parts%p%r(5501) + parts%p%r(5502) + parts%p%r(5503) + parts%p%r(5504) + parts%p%r(5505) + parts%p%r(5506) +  &
!              parts%p%r(5507) + parts%p%r(5508) + parts%p%r(5509) + parts%p%r(5510) + parts%p%r(5511) + parts%p%r(5512) + parts%p%r(5513) + parts%p%r(5514) +  &
!              parts%p%r(5515) + parts%p%r(5516) + parts%p%r(5517) + parts%p%r(5518) + parts%p%r(5519) + parts%p%r(5520) + parts%p%r(5521) + parts%p%r(5522) +  & 
!              parts%p%r(5523)
!2parts%res_F = parts%p%r(5904) + parts%p%r(5905) + parts%p%r(5906) + parts%p%r(5907) + parts%p%r(5908) + parts%p%r(5909) + parts%p%r(5910) + parts%p%r(5911) +  &
!2              parts%p%r(5912) + parts%p%r(5913) + parts%p%r(5914) + parts%p%r(5915) + parts%p%r(5916) + parts%p%r(5917) + parts%p%r(5918) + parts%p%r(5919) +  &
!2              parts%p%r(5920) + parts%p%r(5921) + parts%p%r(5922) + parts%p%r(5924) + parts%p%r(5925) + parts%p%r(5926) + parts%p%r(5927) 

!---  Convert velocity, force, and energy to f and dfdt  

!2      if(itimestep>=save_step_from.and.   &
!2         mod(itimestep,save_step).eq.0)then
!2open(1000,file = 'res_f2224half.dat')
!2write(1000,*) parts%res_F
!2      endif 


if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

      if(itimestep>=save_step_from.and.   &
         mod(itimestep,save_step).eq.0)then
         call output
         endif 

!      if(itimestep>=8650.and.itimestep<=8777)then
!         call output
!      endif 


         
return
      end subroutine
    



!-------------------------------------------------
      subroutine single_step_for_water_Couette
!-------------------------------------------------
!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
!use param 
!use declarations_sph
!use m_sph_fo
implicit none

integer  nphase
type(particles), pointer :: pl
logical :: dbg = .false.
integer i, ntotal,k,j,d
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab,hv(3),r,selfdens
type(p2r) vx_i(3), vx_j(3) 
type(array) temp
type(material),pointer :: water

temp%ndim1 = parts%ntotal+parts%nvirt
allocate(temp%r(parts%ntotal+parts%nvirt))
water => parts%material
pl => parts        
if(dbg) write(*,*) 'In single_step...'
call pl%setup_ndim1
pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
!pl%dvof = 1.d0
 
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length

        
      hv = 0.d0
      r=0.d0



if (pl%numeric%nnps.eq.1) then 
   call direct_find_Couette(pl)
!   call direct_find(pl)
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

! write(*,*) 'niac=',parts%niac


      
!--- Added by Wang
!if(nor_density) call norm_density(pl)

call get_gammaa2(pl)
!      write(*,*)parts%gammaa%r(1),parts%gammaa%r(41),parts%gammaa%r(81),parts%gammaa%r(361)
      
!if(itimestep/=1)then
    
! do i=1,parts%ntotal
!   call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
!   temp%r(i)=temp%r(i) + selfdens*parts%mass%r(i)  !由于rho不能要自加，这里也要自加
! enddo
  
!do k=1,parts%niac
!    i = parts%pair_i(k)
!    j = parts%pair_j(k)
!    rr = 0.
!    c0 = 30!.*sqrt(9.8*0.5)!因case不同而不同
!    ca = c0 * (parts%rho%r(i)/1000.)**3
!    cb = c0 * (parts%rho%r(j)/1000.)**3
!    cab = max(ca,cb)
!    rhoab = parts%rho%r(i) - parts%rho%r(j) - 9.81*1000/c0**2*(parts%x(2,j)-parts%x(2,i))!这里很有可能是加号.
!    vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
!    do d=1,parts%dim
!      dvx(d) = vx_i(d)%p - vx_j(d)%p
!      dx(d) = parts%x(d,i) - parts%x(d,j)
!      rr = rr + dx(d)**2
!    enddo

!    if(parts%itype(i)>0)  temp%r(i) = temp%r(i) + parts%mass%r(j) * parts%w(k)
!    if(parts%itype(j)>0)  temp%r(j) = temp%r(j) + parts%mass%r(i) * parts%w(k)
!    if(parts%itype(i)>0.and.parts%itype(j)>0) then
!        temp%r(i) = temp%r(i) + dt*(parts%mass%r(j)*cab/sqrt(rr)*rhoab/parts%rho%r(j)*(dx(1)*parts%dwdx(1,k) + dx(2)*parts%dwdx(2,k)))
!        temp%r(j) = temp%r(j) + dt*(parts%mass%r(i)*cab/sqrt(rr)*(-rhoab)/parts%rho%r(i)*((-dx(1))*(-parts%dwdx(1,k)) + (-dx(2))*(-parts%dwdx(2,k))))        
!    endif
!enddo

!do i = 1,pl%ntotal
!    pl%rho%r(i) = temp%r(i)/parts%gammaa%r(i)
!enddo
!endif

!------unified get rho of nvirt particles


!call shear_strain_rate(pl)
!pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!      parts%str%x = water%viscosity*parts%tab%x
!      parts%str%y = water%viscosity*parts%tab%y
!      parts%str%xy = water%viscosity*parts%tab%xy


!if(itimestep>1)then
!    call nvirt_density_unified201(pl)
!else
!call nvirt_density_unified2(pl)    
!endif!这样就用了deltaP的方式来求了
!call pressure_nvirt_unified(pl)
!______Analytical value of delta gamma
!call delta_gamma_unified2(pl)

call real_density_unified2(pl)
call nvirt_density_unified2_Couette(pl)    

call delta_gamma_unified3(pl)
do i =1,parts%ntotal + parts%nvirt
    parts%rho%r(i) = 1000
    parts%p%r(i) = 0.d0
enddo

!---  Density approximation or change rate
!if(summation_density)then   
!if(itimestep<62000) then
!if(mod(itimestep,25)==0) then
!call sum_density(pl)
!else             
!    call sum_density(pl)         
    
!    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
!else
!if(mod(itimestep,25)==0) then
!call sum_density(pl)
!else             
!    call sum_density(pl)         
    
!    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
!endif
!do i = parts%ntotal +1,parts%ntotal+parts%nvirt
!    pl%drho%r(i) = 0
!enddo

     
!`````````````````if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!`````````````````      call delta_sph_omp(pl,pl%rho,pl%drho)
!       call delta_rho(pl,pl%rho,pl%drho)
   !endif
!````````````````endif

!---  Dynamic viscosity:
     
!water => parts%material
!parts%p = water%b*((parts%rho/(water%rho0))**water%gamma-1.d0)

!down is for unified condition
water => parts%material
!do i = 1,parts%ntotal 
!parts%p%r(i) = water%b*((parts%rho%r(i)/(water%rho0))**water%gamma-1.d0)
!enddo
!if(itimestep>1)then
!call momentum_equation_unified401(pl) 
!else
call momentum_equation_unified_Couette(pl)      
!endif!这样就用了deltaP的方式来求了

!pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!call velocity_divergence(pl)

!call pressure(pl)

!   call newtonian_fluid(pl)

!      parts%str%x = water%viscosity*parts%tab%x
!      parts%str%y = water%viscosity*parts%tab%y
!      parts%str%xy = water%viscosity*parts%tab%xy
! write(*,*) parts%tab%x%r(862),parts%tab%y%r(861),parts%tab%y%r(882),parts%tab%y%r(883)
! write(*,*) parts%str%x%r(862),parts%str%y%r(861),parts%str%y%r(882),parts%str%y%r(883)
!Calculate internal force for water phase !! -phi_f Grad(p)

!  pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!  pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   
!  pl%dvx%x = pl%dvx%x/pl%rho
!  pl%dvx%y = pl%dvx%y/pl%rho

!        do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!          if(2187<i.and.i<2219)write(*,*) parts%p%r(i),i
!        enddo
!call freesurface2(pl)  !freesurface2通过使用divr来求解。通过比较后发现，这种形式在相应算例中优于freesurface

!     do i=parts%ntotal +1,parts%ntotal + parts%nvirt
!        parts%p%r(i)=parts%p%r(i)/wi(i)
!        if(2187<i.and.i<2219)write(*,*) parts%p%r(i),i
!     enddo
        

!up is for unified condition
 
!call pressure_nvirt(pl)
!第二种状态方程
!parts%p = water%c**2*(parts%rho-water%rho0)

!call freesurface(pl) 
!parts%c%r(1:ntotal) = water%c*(parts%rho%r(1:ntotal)/(water%rho0))**3.0    


!---  Internal forces:

!call shear_strain_rate(pl)
!3~~~~~~~~~down is before
!3pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!3pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!3pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!3pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!3pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!call velocity_divergence(pl)

!call pressure(pl)

!   call newtonian_fluid(pl)

!3      parts%str%x = water%viscosity*parts%tab%x
!3      parts%str%y = water%viscosity*parts%tab%y
!3      parts%str%xy = water%viscosity*parts%tab%xy

!Calculate internal force for water phase !! -phi_f Grad(p)


!3   pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!3   pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   

!down is for dealt-sph
!   pl%dvx%x = pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!   pl%dvx%y = pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y') 
!   call pressure_nvirt(pl)
!   pl%dvx%x = pl%dvx%x - pl%df(pl%p,'x')
!   pl%dvx%y = pl%dvx%y - pl%df(pl%p,'y')
!up is for dealt-sph

!3   pl%dvx%x = pl%dvx%x/pl%rho
!3   pl%dvx%y = pl%dvx%y/pl%rho
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
!3~~~~~~~~~up is before
!if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!       call delta_sph_vx(pl)
   !endif
!endif
!---  Artificial viscosity:

!4~~~~not need now
!if (visc_artificial) call pl%art_visc
!4~~~~not need now       

!if(trim(pl%imaterial)=='water'.and.water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
 !       call pl%delta_sph_omp(pl%vof,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)

!5~~~~not need now
!call pl%repulsive_force_omp                ! can be tried
!      endif
!5~~~~not need now     

!6~~~~not need now
!pl%dvx%y = pl%dvx%y + gravity
!6~~~~not need now

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

!7~~~~not need now
!if (average_velocity) call av_vel(pl) 
!7~~~~not need now

      if(itimestep.eq.100)then
open(1000,file = 'res_f_Couette001.dat')
write(1000,*) parts%vx%x%r(1),parts%vx%x%r(2),parts%vx%x%r(3),parts%vx%x%r(4),parts%vx%x%r(5),parts%vx%x%r(6),parts%vx%x%r(7),parts%vx%x%r(8), &
              parts%vx%x%r(9),parts%vx%x%r(10),parts%vx%x%r(11),parts%vx%x%r(12),parts%vx%x%r(13),parts%vx%x%r(14),parts%vx%x%r(15),parts%vx%x%r(16), &
              parts%vx%x%r(17),parts%vx%x%r(18),parts%vx%x%r(19),parts%vx%x%r(20),parts%vx%x%r(21),parts%vx%x%r(22),parts%vx%x%r(23),parts%vx%x%r(24), &
              parts%vx%x%r(25),parts%vx%x%r(26),parts%vx%x%r(27),parts%vx%x%r(28),parts%vx%x%r(29),parts%vx%x%r(30),parts%vx%x%r(31), &
              parts%vx%x%r(32),parts%vx%x%r(33),parts%vx%x%r(34),parts%vx%x%r(35),parts%vx%x%r(36),parts%vx%x%r(37),parts%vx%x%r(38),parts%vx%x%r(39)
      endif 
      
if(itimestep.eq.1000)then
open(1001,file = 'res_f_Couette01.dat')
write(1001,*) parts%vx%x%r(1),parts%vx%x%r(2),parts%vx%x%r(3),parts%vx%x%r(4),parts%vx%x%r(5),parts%vx%x%r(6),parts%vx%x%r(7),parts%vx%x%r(8), &
              parts%vx%x%r(9),parts%vx%x%r(10),parts%vx%x%r(11),parts%vx%x%r(12),parts%vx%x%r(13),parts%vx%x%r(14),parts%vx%x%r(15),parts%vx%x%r(16), &
              parts%vx%x%r(17),parts%vx%x%r(18),parts%vx%x%r(19),parts%vx%x%r(20),parts%vx%x%r(21),parts%vx%x%r(22),parts%vx%x%r(23),parts%vx%x%r(24), &
              parts%vx%x%r(25),parts%vx%x%r(26),parts%vx%x%r(27),parts%vx%x%r(28),parts%vx%x%r(29),parts%vx%x%r(30),parts%vx%x%r(31), &
              parts%vx%x%r(32),parts%vx%x%r(33),parts%vx%x%r(34),parts%vx%x%r(35),parts%vx%x%r(36),parts%vx%x%r(37),parts%vx%x%r(38),parts%vx%x%r(39)
endif 
      
if(itimestep.eq.10000)then
open(1002,file = 'res_f_Couette1.dat')
write(1002,*) parts%vx%x%r(1),parts%vx%x%r(2),parts%vx%x%r(3),parts%vx%x%r(4),parts%vx%x%r(5),parts%vx%x%r(6),parts%vx%x%r(7),parts%vx%x%r(8), &
              parts%vx%x%r(9),parts%vx%x%r(10),parts%vx%x%r(11),parts%vx%x%r(12),parts%vx%x%r(13),parts%vx%x%r(14),parts%vx%x%r(15),parts%vx%x%r(16), &
              parts%vx%x%r(17),parts%vx%x%r(18),parts%vx%x%r(19),parts%vx%x%r(20),parts%vx%x%r(21),parts%vx%x%r(22),parts%vx%x%r(23),parts%vx%x%r(24), &
              parts%vx%x%r(25),parts%vx%x%r(26),parts%vx%x%r(27),parts%vx%x%r(28),parts%vx%x%r(29),parts%vx%x%r(30),parts%vx%x%r(31), &
              parts%vx%x%r(32),parts%vx%x%r(33),parts%vx%x%r(34),parts%vx%x%r(35),parts%vx%x%r(36),parts%vx%x%r(37),parts%vx%x%r(38),parts%vx%x%r(39)
     endif 

!      if(itimestep.eq.1000)then
!open(1000,file = 'res_f_Couette001.dat')
!write(1000,*) parts%vx%x%r(1320),parts%vx%x%r(1319),parts%vx%x%r(1318),parts%vx%x%r(1317),parts%vx%x%r(1316),parts%vx%x%r(1315),parts%vx%x%r(1314),parts%vx%x%r(1313), &
!              parts%vx%x%r(1312),parts%vx%x%r(1311),parts%vx%x%r(1310),parts%vx%x%r(1309),parts%vx%x%r(1308),parts%vx%x%r(1307),parts%vx%x%r(1306),parts%vx%x%r(1305), &
!              parts%vx%x%r(1304),parts%vx%x%r(1303),parts%vx%x%r(1302),parts%vx%x%r(1301),parts%vx%x%r(1300),parts%vx%x%r(1299),parts%vx%x%r(1298),parts%vx%x%r(1297), &
!              parts%vx%x%r(1296),parts%vx%x%r(1295),parts%vx%x%r(1294),parts%vx%x%r(1293),parts%vx%x%r(1292),parts%vx%x%r(1291),parts%vx%x%r(1290),parts%vx%x%r(1289), &
!             parts%vx%x%r(1288),parts%vx%x%r(1287),parts%vx%x%r(1286),parts%vx%x%r(1285),parts%vx%x%r(1284),parts%vx%x%r(1283),parts%vx%x%r(1282),parts%vx%x%r(1281)
!      endif 

!      if(itimestep.eq.10000)then
!open(1001,file = 'res_f_Couette01.dat')
!write(1001,*) parts%vx%x%r(1320),parts%vx%x%r(1319),parts%vx%x%r(1318),parts%vx%x%r(1317),parts%vx%x%r(1316),parts%vx%x%r(1315),parts%vx%x%r(1314),parts%vx%x%r(1313), &
!              parts%vx%x%r(1312),parts%vx%x%r(1311),parts%vx%x%r(1310),parts%vx%x%r(1309),parts%vx%x%r(1308),parts%vx%x%r(1307),parts%vx%x%r(1306),parts%vx%x%r(1305), &
!              parts%vx%x%r(1304),parts%vx%x%r(1303),parts%vx%x%r(1302),parts%vx%x%r(1301),parts%vx%x%r(1300),parts%vx%x%r(1299),parts%vx%x%r(1298),parts%vx%x%r(1297), &
!              parts%vx%x%r(1296),parts%vx%x%r(1295),parts%vx%x%r(1294),parts%vx%x%r(1293),parts%vx%x%r(1292),parts%vx%x%r(1291),parts%vx%x%r(1290),parts%vx%x%r(1289), &
!             parts%vx%x%r(1288),parts%vx%x%r(1287),parts%vx%x%r(1286),parts%vx%x%r(1285),parts%vx%x%r(1284),parts%vx%x%r(1283),parts%vx%x%r(1282),parts%vx%x%r(1281)
!      endif 
      
!      if(itimestep.eq.100000)then
!open(1002,file = 'res_f_Couette1.dat')
!write(1002,*) parts%vx%x%r(1320),parts%vx%x%r(1319),parts%vx%x%r(1318),parts%vx%x%r(1317),parts%vx%x%r(1316),parts%vx%x%r(1315),parts%vx%x%r(1314),parts%vx%x%r(1313), &
!              parts%vx%x%r(1312),parts%vx%x%r(1311),parts%vx%x%r(1310),parts%vx%x%r(1309),parts%vx%x%r(1308),parts%vx%x%r(1307),parts%vx%x%r(1306),parts%vx%x%r(1305), &
!              parts%vx%x%r(1304),parts%vx%x%r(1303),parts%vx%x%r(1302),parts%vx%x%r(1301),parts%vx%x%r(1300),parts%vx%x%r(1299),parts%vx%x%r(1298),parts%vx%x%r(1297), &
!              parts%vx%x%r(1296),parts%vx%x%r(1295),parts%vx%x%r(1294),parts%vx%x%r(1293),parts%vx%x%r(1292),parts%vx%x%r(1291),parts%vx%x%r(1290),parts%vx%x%r(1289), &
!              parts%vx%x%r(1288),parts%vx%x%r(1287),parts%vx%x%r(1286),parts%vx%x%r(1285),parts%vx%x%r(1284),parts%vx%x%r(1283),parts%vx%x%r(1282),parts%vx%x%r(1281)
!      endif 
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

      if(itimestep>=save_step_from.and.   &
         mod(itimestep,save_step).eq.0)then
         call output
         endif 

!      if(itimestep>=8650.and.itimestep<=8777)then
!         call output
!      endif 


         
return
end subroutine
    
      
!-------------------------------------------------
      subroutine single_step_for_water5 !完全按照第一篇文章的来
!-------------------------------------------------

!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------
!use param 
!use declarations_sph
!use m_sph_fo
implicit none

integer  nphase
type(particles), pointer :: pl
logical :: dbg = .false.
integer i, ntotal,k,j,d
real(dp) dvx(3),dx(3) ,rr,ca,cb,cab,c0,rhoab
type(p2r) vx_i(3), vx_j(3) 
type(array) temp
type(material),pointer :: water

temp%ndim1 = parts%ntotal+parts%nvirt
allocate(temp%r(parts%ntotal+parts%nvirt))
water => parts%material
pl => parts        
if(dbg) write(*,*) 'In single_step...'
call pl%setup_ndim1
pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
!pl%dvof = 1.d0
 
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length

if(itimestep/=1)then
  temp = 0.
do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    if(parts%itype(i)>0)    temp%r(i) = temp%r(i) - parts%mass%r(j) * parts%w(k)
    if(parts%itype(j)>0)    temp%r(j) = temp%r(j) - parts%mass%r(i) * parts%w(k)
enddo
endif

do i =1,parts%ntotal
    temp%r(i) = temp%r(i) + parts%gammaa%r(i) * parts%rho%r(i)
enddo

if (pl%numeric%nnps.eq.1) then 
   call direct_find(pl)
else if (pl%numeric%nnps.eq.2) then
   call link_list(pl)     

else if (pl%numeric%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
endif         

if(mod(itimestep,print_step).eq.0.and.int_stat) then
   call pl%interaction_statistics
endif   

!--- Added by Wang
!if(nor_density) call norm_density(pl)

call get_gammaa2(pl)
if(itimestep/=1)then

do k=1,parts%niac
    i = parts%pair_i(k)
    j = parts%pair_j(k)
    if(parts%itype(i)>0)  temp%r(i) = temp%r(i) + parts%mass%r(j) * parts%w(k)
    if(parts%itype(j)>0)  temp%r(j) = temp%r(j) + parts%mass%r(i) * parts%w(k)
enddo

do i = 1,pl%ntotal
    pl%rho%r(i) = temp%r(i)/parts%gammaa%r(i)
enddo
endif

!------unified get rho of nvirt particles

call nvirt_density_unified2(pl)       

!call pressure_nvirt(pl)
!______Analytical value of delta gamma
call delta_gamma_unified3(pl)

call real_density_unified2(pl)

water => parts%material
do i = 1,parts%ntotal 
parts%p%r(i) = water%b*((parts%rho%r(i)/(water%rho0))**water%gamma-1.d0)
enddo
call momentum_equation_unified(pl)  

!up is for unified condition
 
!call pressure_nvirt(pl)
!第二种状态方程
!parts%p = water%c**2*(parts%rho-water%rho0)

!call freesurface(pl) 
!parts%c%r(1:ntotal) = water%c*(parts%rho%r(1:ntotal)/(water%rho0))**3.0    


!---  Internal forces:

!call shear_strain_rate(pl)
!3~~~~~~~~~down is before
!3pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!3pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
!Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

!3pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
!3pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
!3pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

!call velocity_divergence(pl)

!call pressure(pl)

!   call newtonian_fluid(pl)

!3      parts%str%x = water%viscosity*parts%tab%x
!3      parts%str%y = water%viscosity*parts%tab%y
!3      parts%str%xy = water%viscosity*parts%tab%xy

!Calculate internal force for water phase !! -phi_f Grad(p)


!3   pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!3   pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   

!down is for dealt-sph
!   pl%dvx%x = pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
!   pl%dvx%y = pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y') 
!   call pressure_nvirt(pl)
!   pl%dvx%x = pl%dvx%x - pl%df(pl%p,'x')
!   pl%dvx%y = pl%dvx%y - pl%df(pl%p,'y')
!up is for dealt-sph

!3   pl%dvx%x = pl%dvx%x/pl%rho
!3   pl%dvx%y = pl%dvx%y/pl%rho
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
!3~~~~~~~~~up is before
!if(artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
!       call delta_sph_vx(pl)
   !endif
!endif
!---  Artificial viscosity:

!4~~~~not need now
!if (visc_artificial) call pl%art_visc
!4~~~~not need now       

!if(trim(pl%imaterial)=='water'.and.water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
 !       call pl%delta_sph_omp(pl%vof,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)

!5~~~~not need now
!call pl%repulsive_force_omp                ! can be tried
!      endif
!5~~~~not need now     

!6~~~~not need now
!pl%dvx%y = pl%dvx%y + gravity
!6~~~~not need now

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

!7~~~~not need now
!if (average_velocity) call av_vel(pl) 
!7~~~~not need now

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(itimestep,print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

      if(itimestep>=save_step_from.and.   &
         mod(itimestep,save_step).eq.0)then
         call output
         endif 

!      if(itimestep>=8650.and.itimestep<=8777)then
!         call output
!      endif 


         
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
pl%dstr%x = 0.d0; pl%dstr%xy = 0.d0; pl%dstr%y = 0.d0
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
      call delta_sph(pl,pl%rho%r,pl%drho%r)
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

pl%dvx%x = - pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
pl%dvx%y = - pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')

pl%dvx%x = pl%dvx%x/pl%rho
pl%dvx%y = pl%dvx%y/pl%rho       

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
   call pl%delta_sph(pl%p%r,pl%dp%r)
   call pl%delta_sph(pl%str%x%r,pl%dstr%x%r)
   call pl%delta_sph(pl%str%xy%r,pl%dstr%xy%r)
   call pl%delta_sph(pl%str%y%r,pl%dstr%y%r)
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
      mass   => parts%mass%r
      hsml   => parts%hsml
      p      => parts%p%r

              if(trim(parts%imaterial)=='water')then

      write(f_xv,*) 'VARIABLES="X","Y","Pressure","VoF",'     
      write(f_xv,*) '"vx","vy","rho","zone","VoF2","mass" '
      write(f_xv,*) 'ZONE I=', ntotal, ' F=BLOCK'
      write(f_xv,*) (x(1,i),i=1,ntotal) !, soil%x(1,1:ntotal2)
      write(f_xv,*) (x(2,i),i=1,ntotal) !, soil%x(2,1:ntotal2)
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
      write(f_xv,*) parts%vof2%r(1:ntotal)
      write(f_xv,*) parts%mass%r(1:ntotal)

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
      write(f_xv,*)  parts%epsilon_p%r(1:ntotal)
      write(f_xv,*)  parts%str%xy%r(1:ntotal)
      write(f_xv,*)  parts%str%x%r(1:ntotal)
      write(f_xv,*)  parts%str%y%r(1:ntotal)
      write(f_xv,*)  parts%vx%x%r(1:ntotal)
      write(f_xv,*)  parts%vx%y%r(1:ntotal)
      write(f_xv,*)  parts%rho%r(1:ntotal)
      write(f_xv,*)  parts%mass%r(1:ntotal)
      write(f_xv,*)  -parts%p%r(1:ntotal) + parts%str%y%r(1:ntotal)
      write(f_xv,*)  parts%zone(1:ntotal)
      write(f_other,*) time, -parts%p%r(395)+parts%str%y%r(395)

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
      write(f_state,*)  soil%epsilon_p%r(1:ntotal2)
      write(f_state,*)  soil%str%xy%r(1:ntotal2)
      write(f_state,*)  soil%str%x%r(1:ntotal2)
      write(f_state,*)  soil%str%y%r(1:ntotal2)
      write(f_state,*)  soil%vx%x%r(1:ntotal2)
      write(f_state,*)  soil%vx%y%r(1:ntotal2)
      write(f_state,*)  soil%rho%r(1:ntotal2)
      write(f_state,*)  soil%mass%r(1:ntotal2)
      write(f_state,*)  -soil%p%r(1:ntotal2) + soil%str%y%r(1:ntotal2)
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
