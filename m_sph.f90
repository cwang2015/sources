!-----------------------------
    module declarations_sph
!-----------------------------
use m_particles
implicit none

! Logic unit of files

integer f_xv,f_state,f_other,f_tmp
integer xv_vp,state_vp,other_vp

! Files
character(len=32) results, res_water, res_soil, res_other

! Boundary conditions
double precision inlet_velocity

type(material), target :: H2O, SiO2

! Particles 
type(particles), target :: parts, soil


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

          read(pvalu,*) parts%print_step
          write(*,*) 'Print step = ', parts%print_step
          soil%print_step = parts%print_step

      case('SAVE STEP')

          read(pvalu,*) parts%save_step
          write(*,*) 'Save step = ', parts%save_step
          soil%save_step = parts%save_step

      case('PARTICLE PAIRS FIND METHOD')
         read(pvalu,*) parts%nnps
         write(*,*) 'Particle pairs find method = ', parts%nnps 
         soil%nnps = parts%nnps

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
                  read(pvalu,*) parts%single_phase
                  write(*,*) 'Single phase = ', parts%single_phase

      case('FIRST PHASE')

                  read(pvalu,*) parts%imaterial
                  write(*,*) 'First phase = ', parts%imaterial

      case('SECOND PHASE')

                  read(pvalu,*) soil%imaterial
                  write(*,*) 'Second phase = ', soil%imaterial

!      case('DRY SOIL')

!                  read(pvalu,*) dry_soil
!                  write(*,*) 'Dry soil = ', dry_soil

      case('PARTICLE SIZE')
!                  ndim = 1
!                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  read(pvalu,*) parts%particle_size
                  write(*,*) 'Particle size = ', parts%particle_size

!       case('OPERATING MODE')
!
!                  read(pvalu,*) operating_mode
!                  write(*,*) 'Operating mode = ', operating_mode

       case('IMMERSE DEPTH')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  !read(pvalu,*) wasserjet%immerse
                  !write(*,*) 'Immerse depth = ', wasserjet%immerse

!       case('WITH NOZZLE')
!                  
!                  read(pvalu,*) with_nozzle
!                  write(*,*) 'With nozzle = ', with_nozzle

      case('WATER PROPERTY')

          read(pvalu,*) h2o%rho0,h2o%b,h2o%gamma, h2o%c, h2o%viscosity
          write(*,*) 'Water property: rho0, b, gamma, c, viscosity' 
          write(*,*) h2o%rho0, h2o%b, h2o%gamma, h2o%c, h2o%viscosity

      case('WATER PRESSURE')
          read(pvalu,*) parts%water_pressure
          write(*,*) 'Water pressure: 0=hydro-static, 1=dynamic'
          write(*,*) parts%water_pressure    

      case('SOIL PROPERTY')

          read(pvalu,*) sio2%rho0, sio2%k, sio2%E, sio2%niu, sio2%c
          write(*,*) 'Soil property: rho0, k, Youngs, Poisson, c' 
          write(*,*) sio2%rho0, sio2%k, sio2%E, sio2%niu, sio2%c

      case('DARCY PARAMETER')

          read(pvalu,*) sio2%porosity, sio2%permeability
          write(*,*) 'Darcy parameter: porosity, permeability '
          write(*,*) sio2%porosity, sio2%permeability

      case('INTEGRATION SCHEME')
          read(pvalu,*) parts%integrate_scheme
          write(*,*) 'Integration scheme: ', parts%integrate_scheme 

      case('SOIL FAILURE CRITERIA')
       
          read(pvalu,*) sio2%cohesion, SIO2%phi, sio2%psi
          write(*,*) 'Soil failure criteria: cohesion, phi, psi'
          write(*,*) sio2%cohesion, sio2%phi, sio2%psi

      case('FLOW RULE')

          read(pvalu,*) parts%flow_rule
          write(*,*) 'Flow rule: 1, associated, 2, non-associated'
          write(*,*) parts%flow_rule

      case('CRITICAL STATE')
          read(pvalu,*) parts%critical_state
          write(*,*) 'Critical state: '
          write(*,*) parts%critical_state
          soil%critical_state = parts%critical_state

      case('TIME STEP SIZE')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
                  read(pvalu,*) parts%dt
                  write(*,*) 'Time step size = ', parts%dt
      case('MAXIMUM TIME STEP')
                  ndim = 1
                  if(ncoma>=ndim) call trim_pvalu(trim(pvalu),',',ndim)
!                  read(pvalu,*) max_time_step
!                  write(*,*) 'Maximum time step = ', max_time_step

      case('OPENMP NUMBER OF THREADS')

                  read(pvalu,*) parts%nthreads
                  write(*,*) 'Openmp number of threads = ', parts%nthreads
                  soil%nthreads = parts%nthreads

      case('BACKGROUND GRID RANGE')
                  read(pvalu,*) parts%x_mingeom,parts%x_maxgeom,  & 
                                parts%y_mingeom,parts%y_maxgeom
          write(*,*) 'Background grid grange: x_min, x_max, y_min, y_max'
          write(*,*) parts%x_mingeom,parts%x_maxgeom,  & 
                     parts%y_mingeom,parts%y_maxgeom

      case('USAW')
                  read(pvalu,*) parts%usaw
          write(*,*) 'unified semi-analytical wall boundary condition'
          write(*,*) parts%usaw
          soil%usaw = parts%usaw

      case('PARTICLE SHIFT')
          read(pvalu,*) parts%particle_shift
          write(*,*) 'Particle shift = ', parts%particle_shift
          !soil%usaw = parts%usaw          

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

allocate(parts%w(max_interaction));             parts%w   = 0.d0
allocate(parts%dwdx(dim,max_interaction));      parts%dwdx= 0.d0
allocate(parts%pair_i(max_interaction));        parts%pair_i = 0
allocate(parts%pair_j(max_interaction));        parts%pair_j = 0
allocate(parts%countiac(maxn));                 parts%countiac=0
allocate(parts%grid(parts%maxngx,parts%maxngy,parts%maxngz)); parts%grid = 0
allocate(parts%xgcell(dim,maxn)); parts%xgcell = 0
allocate(parts%celldata(maxn)); parts%celldata = 0
allocate(parts%norm(2,10000))

! Sherpard Filter
allocate(parts%wi); call parts%wi%set(rank=0,parts=parts,name = ' ')
call parts%wi%max_alloc
parts%wi%r = 0.d0

! Boundary information
allocate(parts%bn); call parts%bn%set(rank=1,parts=parts,name='bn')
call parts%bn%max_alloc
parts%bn = 0.d0

allocate(parts%bs); call parts%bs%set(rank=0,parts=parts, name='bs')
call parts%bs%max_alloc
parts%bs = 0.d0

! Fields variables

call allocate_particles_fields(parts)

! Working arrays outside particles type
allocate(parts%tab); call parts%tab%set(rank=2,parts=parts,name = 'tab')
!allocate(parts%tab%x, parts%tab%y, parts%tab%xy)
!allocate(parts%tab%x%r(maxn),parts%tab%y%r(maxn),parts%tab%xy%r(maxn))
call parts%tab%max_alloc
parts%tab%x%r = 0.d0; parts%tab%y%r = 0.d0; parts%tab%xy%r = 0.d0

!parts%tab => tab

if(trim(parts%imaterial)=='water')then
!   allocate(parts%sxx,parts%sxy,parts%syy)
!    allocate(parts%str);allocate(parts%str%x,parts%str%xy,parts%str%y)
!    parts%str%x%r => parts%tab%x%r
!    parts%str%y%r => parts%tab%y%r
!    parts%str%xy%r => parts%tab%xy%r
   parts%str => parts%tab
elseif(trim(parts%imaterial)=='soil')then
   parts%wxy => parts%tab%xy
endif   

! Material assignment

if(trim(parts%imaterial)=='water') parts%material => H2O
if(trim(parts%imaterial)=='soil')  parts%material => SiO2

! Numerical parameters
!parts%numeric => numeric
if(.not.associated(parts%numeric))allocate(parts%numeric)

call parts%get_num_threads
allocate(parts%niac_start(parts%nthreads),parts%niac_end(parts%nthreads))

!--------------- Soil particles --------------------------------

if(parts%single_phase) return

soil%dim = dim
soil%maxn= maxn
soil%max_interaction = max_interaction

call allocate_particles(soil)

allocate(soil%zone(maxn)); soil%zone = 0

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

allocate(soil%wi); call soil%wi%set(rank=0,parts=soil,name = ' ')
call soil%wi%max_alloc
soil%wi%r = 0.d0

! Fields variables
call allocate_particles_fields(soil)

!soil%tab => parts%tab
allocate(soil%tab); call soil%tab%set(rank=2, parts=soil, name='tab')
call soil%tab%max_alloc
soil%tab%x%r = 0.d0; soil%tab%y%r = 0.d0; soil%tab%xy%r = 0.d0

soil%wxy => soil%tab%xy

soil%material => SiO2
soil%numeric => parts%numeric

soil%flow_rule = parts%flow_rule

soil%x_maxgeom = parts%x_maxgeom; soil%x_mingeom = parts%x_mingeom
soil%y_maxgeom = parts%y_maxgeom; soil%y_mingeom = parts%y_mingeom

call soil%get_num_threads
allocate(soil%niac_start(soil%nthreads),soil%niac_end(soil%nthreads))

if(soil%critical_state==1)then
   allocate(soil%dvof); call soil%dvof%set(rank=0,parts=soil,name='dvof')
   call soil%dvof%max_alloc; soil%dvof%r = 0.d0   
   allocate(soil%vof_min)
   call soil%vof_min%set(rank=0,parts=soil,name='vof_min')
   call soil%vof_min%max_alloc; soil%vof_min%r = 0.d0   

   allocate(soil%psi); call soil%psi%set(rank=0,parts=soil,name='psi')
   call soil%psi%max_alloc; soil%psi%r = 0.d0  
   allocate(soil%rho0); call soil%rho0%set(rank=0,parts=soil,name='rho0')
   call soil%rho0%max_alloc; soil%rho0%r = 0.d0      
endif

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
    allocate(parts%vol); 
    !allocate(parts%vol%r(maxn))
    call parts%vol%set(rank=0,parts=parts,name='vol')
    call parts%vol%max_alloc
endif
if(.not.associated(parts%mass)) then
    allocate(parts%mass); call parts%mass%set(rank=0,parts=parts,name='mass')
    !allocate(parts%mass%r(maxn))
    call parts%mass%max_alloc
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

allocate(parts%vx); call parts%vx%set(rank=1,parts=parts,name='vx')
call parts%vx%max_alloc
parts%vx%x%r = 0.d0; parts%vx%y%r = 0.d0
allocate(parts%rho); call parts%rho%set(rank=0,parts=parts,name='rho')
call parts%rho%max_alloc
parts%rho%r = 0.d0
allocate(parts%p); call parts%p%set(rank=0,parts=parts,name='p')
call parts%p%max_alloc
parts%p%r = 0.d0
allocate(parts%c); call parts%c%set(rank=0,parts=parts,name='c')
call parts%c%max_alloc
parts%c%r = 0.d0

! Old value
allocate(parts%v_min); call parts%v_min%set(rank=1,parts=parts,name='v_min')
call parts%v_min%max_alloc
parts%v_min%x%r = 0.d0; parts%v_min%y%r = 0.d0

allocate(parts%rho_min); call parts%rho_min%set(rank=0,parts=parts,name='rho_min')
call parts%rho_min%max_alloc
parts%rho_min%r = 0.d0

! Accelerations
allocate(parts%drho);call parts%drho%set(rank=0,parts=parts,name='drho')
call parts%drho%max_alloc; parts%drho%r = 0.d0
allocate(parts%dvx); call parts%dvx%set(rank=1,parts=parts,name='dvx')
call parts%dvx%max_alloc
parts%dvx%x%r = 0.d0; parts%dvx%y%r = 0.d0  
allocate(parts%av); call parts%av%set(rank=1,parts=parts,name='av')
call parts%av%max_alloc
parts%av%x%r  = 0.d0; parts%av%y%r  = 0.d0

! For particle shifting
allocate(parts%r0);call parts%r0%set(rank=0,parts=parts,name='r0')
call parts%r0%max_alloc; parts%r0%r = 0.d0
allocate(parts%mt);call parts%mt%set(rank=0,parts=parts,name='mt')
call parts%mt%max_alloc; parts%mt%r = 0.d0
allocate(parts%ps); call parts%ps%set(rank=1,parts=parts,name='ps')
call parts%ps%max_alloc
parts%ps%x%r  = 0.d0; parts%ps%y%r  = 0.d0

! Volume Fraction
allocate(parts%vof); call parts%vof%set(rank=0,parts=parts,name='vof')
call parts%vof%max_alloc; parts%vof%r = 0.d0

if(trim(parts%imaterial)=='water')then
   ! phi_f = 1- phi_s
   allocate(parts%vof2); call parts%vof2%set(rank=0,parts=parts,name='vof2')
   call parts%vof2%max_alloc; parts%vof2%r = 0.d0
   allocate(parts%dvof); call parts%dvof%set(rank=0,parts=parts,name='dvof')
   call parts%dvof%max_alloc; parts%dvof%r = 0.d0   
   allocate(parts%vof_min)
   call parts%vof_min%set(rank=0,parts=parts,name='vof_min')
   call parts%vof_min%max_alloc; parts%vof_min%r = 0.d0   
endif

if(trim(parts%imaterial)=='soil')then
    allocate(parts%str)
!    allocate(parts%str%x,parts%str%xy,parts%str%y)
!   allocate(parts%sxx,parts%sxy,parts%syy)
!   allocate(parts%str%x%r(maxn), parts%str%y%r(maxn), parts%str%xy%r(maxn))
   call parts%str%set(rank=2,parts=parts,name='str'); call parts%str%max_alloc
   parts%str%x%r = 0.d0; parts%str%xy%r = 0.d0; parts%str%y%r = 0.d0   !!! %r
   
   
    !artificial stress
    allocate(parts%strp); 
    call parts%strp%set(rank=2,parts=parts,name='strp')
    call parts%strp%max_alloc
    parts%strp%x%r = 0.d0; parts%strp%y%r = 0.d0; parts%strp%xy%r = 0.d0
    
    allocate(parts%ast); call parts%ast%set(rank=2,parts=parts,name='ast')
    call parts%ast%max_alloc
    parts%ast%x%r = 0.d0; parts%ast%y%r = 0.d0; parts%ast%xy%r = 0.d0
    
    allocate(parts%astp); call parts%astp%set(rank=2,parts=parts,name='astp')
    call parts%astp%max_alloc
    parts%astp%x%r = 0.d0; parts%astp%y%r = 0.d0; parts%astp%xy%r = 0.d0
    
    allocate(parts%thetai); call parts%thetai%set(rank=0,parts=parts,name='thetai')
    call parts%thetai%max_alloc
    parts%thetai%r = 0.d0
   

   allocate(parts%str_min)
   !allocate(parts%str_min%x,parts%str_min%y,parts%str_min%xy)
   !allocate(parts%str_min%x%r(maxn),parts%str_min%y%r(maxn),parts%str_min%xy%r(maxn))
   call parts%str_min%set(rank=2,parts=parts,name='str_min')
   call parts%str_min%max_alloc
   parts%str_min%x%r = 0.d0; parts%str_min%xy%r = 0.d0; parts%str_min%y%r = 0.d0 !!!%r

   allocate(parts%p_min); call parts%p_min%set(rank=0,parts=parts,name='p_min')
   call parts%p_min%max_alloc

   allocate(parts%vcc); call parts%vcc%set(rank=0,parts=parts,name='vcc')
   call parts%vcc%max_alloc  
   !allocate(parts%szz(maxn))    !!! Intel Fortran Compiler is shit!!!

   parts%p_min%r = 0.d0
   parts%vcc%r = 0.d0

   allocate(parts%dstr)
   !allocate(parts%dstr%x, parts%dstr%y, parts%dstr%xy)
   !allocate(parts%dstr%x%r(maxn),parts%dstr%y%r(maxn),parts%dstr%xy%r(maxn))
   call parts%dstr%set(rank=2,parts=parts,name='dstr'); call parts%dstr%max_alloc
   parts%dstr%x%r = 0.d0; parts%dstr%xy%r = 0.d0; parts%dstr%y%r = 0.d0  !!!%r

   allocate(parts%dp); call parts%dp%set(rank=0,parts=parts,name='dp')
   call parts%dp%max_alloc; parts%dp%r = 0.d0 !!!%r

! For return mapping algorithm
   allocate(parts%dstr2)
   call parts%dstr2%set(rank=2,parts=parts,name='dstr2'); call parts%dstr2%max_alloc
   !,parts%dstr2%x,parts%dstr2%y,parts%dstr2%xy)
   !allocate(parts%dstr2%x%r(maxn),parts%dstr2%y%r(maxn),parts%dstr2%xy%r(maxn))
   parts%dstr2%x%r = 0.d0; parts%dstr2%xy%r = 0.d0; parts%dstr2%y%r = 0.d0 !!!%r
   allocate(parts%dp2); call parts%dp2%set(rank=0,parts=parts,name='dp2')
   !allocate(parts%dp2%r(maxn))
   call parts%dp2%max_alloc
   parts%dp2%r = 0.d0

   allocate(parts%fail(maxn)); parts%fail = 0

   allocate(parts%epsilon_p)
   call parts%epsilon_p%set(rank=0,parts=parts,name='ep')
   !allocate(parts%epsilon_p%r(maxn))
   call parts%epsilon_p%max_alloc
   parts%epsilon_p%r = 0.d0  !!!%r
   
   allocate(parts%spp)   !!! VOF
   call parts%spp%set(rank=0,parts=parts,name='spp'); call parts%spp%max_alloc
   parts%spp%r = 0.d0

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
      
      !call parts%setup_ndim1; call soil%setup_ndim1

      do it = 1, parts%maxtimestep   
  
        parts%itimestep = parts%itimestep+1

        !parts%itimestep = itimestep
        soil%itimestep  = parts%itimestep

        if (mod(parts%itimestep,parts%print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   parts%itimestep,'     current time=', real(parts%time+parts%dt)
         write(*,*)'______________________________________________'
        endif  
        

!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

        if (parts%itimestep .ne. 1) then
           pl => parts
           call first_half
        
                             ! For soil
                             !if(.not.single_phase)then
                                pl => soil
                                call first_half

                                if(pl%critical_state==0)then
                                call drucker_prager_failure_criterion(pl)
                                else
                        call  drucker_prager_failure_criterion_critical_state(pl)
                                endif

                                !call volume_fraction_soil(soil)
                                !call volume_fraction_water(parts, soil)
                             !endif  

        endif

!---  Definition of variables out of the function vector:    

!        call single_step_VOF
        call single_step

!****************************************************************        
!        do i = 1, parts%ntotal + parts%nvirt
!           if(parts%zone(i)==4)then
!              parts%dvx%x%r(i)=0.d0; parts%dvx%y%r(i)=0.d0
!              parts%av%x%r(i)=0.d0; parts%av%y%r(i)=0.d0
!              parts%dvof%r(i)=0.d0
!              parts%drho%r(i)=0.d0
!           endif   
!        enddo
!***************************************************************

        if (parts%itimestep .eq. 1) then
           pl => parts
           call first_step

                                !if(.not.single_phase)then
                                    pl => soil
                                    call first_step

                                    if(pl%plasticity==2)then
                                       call plastic_or_not(pl)
                                       call return_mapping
                                    endif
                                    
                                    if(pl%critical_state==0)then
                                    call drucker_prager_failure_criterion(pl)
                                    else
                   call drucker_prager_failure_criterion_critical_state(pl)
                                    endif

                                    !call volume_fraction_soil(soil)
                                    !call volume_fraction_water(parts, soil)
                                !endif
        else
   
           pl => parts
           call second_half          

                                !if(.not.single_phase)then
                                   pl => soil
                                   call second_half

                                   if(pl%plasticity==2)then
                                      call plastic_or_not(pl)
                                      call return_mapping
                                   endif

                                    if(pl%critical_state==0)then
                                    call drucker_prager_failure_criterion(pl)
                                    else
                   call drucker_prager_failure_criterion_critical_state(pl)
                                    endif

                                   
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

!        call vof_sph
!        cycle
  
        parts%time = parts%time + parts%dt
        soil%time = parts%time

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
            
         if (.not.pl%summation_density) then    
            pl%rho_min = pl%rho
            pl%rho = pl%rho +(pl%dt/2.)* pl%drho
         endif

         if(trim(pl%imaterial)=='water'.and.pl%volume_fraction)then
            pl%vof_min = pl%vof
            pl%vof = pl%vof+(pl%dt/2.)*pl%dvof 
         endif

         if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
            pl%vof_min = pl%vof
            pl%vof = pl%vof+(pl%dt/2.)*pl%dvof 
         endif


         if(trim(pl%imaterial)=='soil')then
                  if(pl%stress_integration==1)then
         !pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)   ! Nothing to do!
         !pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
         !pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
                  elseif(pl%stress_integration==2)then
         pl%str_min%x = pl%str%x
         pl%str_min%xy = pl%str%xy
         pl%str_min%y = pl%str%y
         pl%p_min   = pl%p
         pl%str%x = pl%str%x + (pl%dt/2.)*pl%dstr%x
         pl%str%xy = pl%str%xy + (pl%dt/2.)*pl%dstr%xy
         pl%str%y = pl%str%y + (pl%dt/2.)*pl%dstr%y
         !pl%str = pl%str + (dt/2.) * pl%dstr
         pl%p = pl%p + (pl%dt/2.)*pl%dp      !!!simultaneous pressure
                  endif
         endif
      !write(*,*) pl%ntotal + pl%nvirt    
      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i); v_mini = pl%v_min%cmpt(i)
         do d = 1, pl%dim
         !   pl%v_min(d, i) = pl%vx(d, i)
         !   pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
            v_mini(d)%p = vxi(d)%p
            vxi(d)%p = vxi(d)%p + (pl%dt/2.)*dvxi(d)%p         
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
      type(p2r) vxi(3), dvxi(3), v_mini(3), avi(3),psi(3)

!      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if (.not.pl%summation_density ) then
            pl%rho = pl%rho + (pl%dt/2.)* pl%drho
         endif

         if(trim(pl%imaterial)=='water'.and.pl%volume_fraction)then
            pl%vof = pl%vof+(pl%dt/2.)*pl%dvof
         endif

         if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
            pl%vof = pl%vof+(pl%dt/2.)*pl%dvof
         endif


         if(trim(pl%imaterial)=='soil')then
                  if(pl%stress_integration==1)then
         pl%str%x = pl%str%x + pl%dt*pl%dstr%x
         pl%str%xy = pl%str%xy + pl%dt*pl%dstr%xy
         pl%str%y = pl%str%y + pl%dt*pl%dstr%y
         pl%p   = pl%p   + pl%dt*pl%dp       !!! simultaneous pressure
                  elseif(pl%stress_integration==2)then
         pl%str%x = pl%str%x + (pl%dt/2.)*pl%dstr%x
         pl%str%xy = pl%str%xy + (pl%dt/2.)*pl%dstr%xy
         pl%str%y = pl%str%y + (pl%dt/2.)*pl%dstr%y
         pl%p   = pl%p + (pl%dt/2.)*pl%dp    !!! simultaneous pressure
                  endif
         endif
        
      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i); avi = pl%av%cmpt(i)
         psi = pl%ps%cmpt(i)
         do d = 1, pl%dim        
            !pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = vxi(d)%p + (pl%dt/2.) * dvxi(d)%p + avi(d)%p
            !pl%x(d, i) = pl%x(d, i) + pl%dt * vxi(d)%p            
            pl%x(d, i) = pl%x(d, i) + pl%dt * (vxi(d)%p+psi(d)%p)   !!! Shifting            
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
      type(p2r) vxi(3), dvxi(3), v_mini(3),avi(3), psi(3)
      !$omp parallel do private(d,vxi,dvxi,v_mini,avi)
      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if (.not.pl%summation_density ) then 
            pl%rho%r(i) = pl%rho_min%r(i) + pl%dt*pl%drho%r(i)
         endif
   
         if(trim(pl%imaterial)=='water'.and.pl%volume_fraction)then
            pl%vof%r(i) = pl%vof_min%r(i)+pl%dt*pl%dvof%r(i) 
         endif
                 
         if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
            pl%vof%r(i) = pl%vof_min%r(i)+pl%dt*pl%dvof%r(i) 
         endif


         if(trim(pl%imaterial)=='soil')then 
                  if(pl%stress_integration==1)then 
         pl%str%x%r(i) = pl%str%x%r(i) + pl%dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + pl%dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + pl%dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p%r(i)   + pl%dt*pl%dp%r(i)       !!! simultaneous pressure
                  elseif(pl%stress_integration==2)then
         pl%str%x%r(i) = pl%str_min%x%r(i) + pl%dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str_min%xy%r(i) + pl%dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str_min%y%r(i) + pl%dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p_min%r(i) + pl%dt*pl%dp%r(i)     !!! simultaneous pressure
                  endif
         endif
        
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i) ; v_mini = pl%v_min%cmpt(i)
         avi = pl%av%cmpt(i); psi = pl%ps%cmpt(i)
         do d = 1, pl%dim                   
            !pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = v_mini(d)%p + pl%dt * dvxi(d)%p   &
                        + avi(d)%p
            !pl%x(d, i) = pl%x(d, i) + pl%dt * vxi(d)%p
             pl%x(d, i) = pl%x(d, i) + pl%dt * (vxi(d)%p + psi(d)%p) !!!Shifting
         !if(i>4160)write(*,*) i,vxi(2)%p
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
                  if(pl%stress_integration==1)then 
         pl%str%x%r(i) = pl%str%x%r(i) + pl%dt*pl%dstr2%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + pl%dt*pl%dstr2%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + pl%dt*pl%dstr2%y%r(i)
         pl%p%r(i)   = pl%p%r(i)   + pl%dt*pl%dp2%r(i)       !!! simultaneous pressure
                  elseif(pl%stress_integration==2)then
         pl%str%x%r(i) = pl%str_min%x%r(i) + pl%dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str_min%xy%r(i) + pl%dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str_min%y%r(i) + pl%dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p_min%r(i) + pl%dt*pl%dp%r(i)     !!! simultaneous pressure
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
do it = 1, parts%maxtimestep    
   parts%itimestep = parts%itimestep+1
   !parts%itimestep = itimestep


        if (mod(parts%itimestep,parts%print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',              &
                   parts%itimestep,'     current time=', real(parts%time+parts%dt)
         write(*,*)'______________________________________________'
        endif  
   
!     If not first time step, then update thermal energy, density and 
!     velocity half a time step  

   if(parts%itimestep .ne. 1)then
pl => parts
      call first_half
!      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
                       
   endif
pl => parts
   call single_step_for_water
!    call single_step

   if(parts%itimestep .eq. 1) then
pl => parts
      call first_step

   else
pl => parts
       call second_half 
       
       

   endif 

   parts%time = parts%time + parts%dt

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
            
         if (.not.parts%summation_density) then    
            parts%rho_min = parts%rho
            parts%rho = parts%rho +(parts%dt/2.)* parts%drho
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
            vxi(d)%p = vxi(d)%p + (parts%dt/2.)*dvxi(d)%p         
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
         
         if (.not.parts%summation_density ) then
            parts%rho = parts%rho + (parts%dt/2.)* parts%drho
         endif

      do i=1,parts%ntotal +parts%nvirt     ! origionally pl%ntotal
         if(parts%itype(i)<0)cycle
         vxi = parts%vx%cmpt(i); dvxi = parts%dvx%cmpt(i); avi = parts%av%cmpt(i)
         do d = 1, parts%dim        
            !pl%vx(d, i) = pl%vx(d, i) + (dt/2.) * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = vxi(d)%p + (parts%dt/2.) * dvxi(d)%p + avi(d)%p
            parts%x(d, i) = parts%x(d, i) + parts%dt * vxi(d)%p            
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
            
         if (.not.parts%summation_density ) then 
            parts%rho%r(i) = parts%rho_min%r(i) + parts%dt*parts%drho%r(i)
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
            vxi(d)%p = v_mini(d)%p + parts%dt * dvxi(d)%p   &
                        + avi(d)%p
            parts%x(d, i) = parts%x(d, i) + parts%dt * vxi(d)%p
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

#ifdef REMOVE      
!!DEC$IF(.FALSE.)
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
         
!!DEC$ENDIF
#endif

!----------------------------------------------------------------------
             subroutine time_integration_for_soil
!----------------------------------------------------------------------
!use param
!use declarations_sph
implicit none     

integer :: i, j, k, d, ntotal, it
type(particles), pointer :: pl
type(p2r) vxi(3), dvxi(3), v_mini(3), avi(3)
              
do it = 1, parts%maxtimestep     
   parts%itimestep = parts%itimestep+1
   !parts%itimestep = itimestep
  
! If not first time step, then update thermal energy, density and 
! velocity half a time step  

   if(parts%itimestep .ne. 1)then

      pl => parts
!      call first_half
!      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
            
         if (.not.pl%summation_density) then    
            pl%rho_min = pl%rho
            pl%rho = pl%rho +(pl%dt/2.)* pl%drho
         endif
 
         if(pl%stress_integration==1)then
            !pl%sxx(i) = pl%sxx(i) + dt*pl%dsxx(i)   ! Nothing to do!
            !pl%sxy(i) = pl%sxy(i) + dt*pl%dsxy(i)
            !pl%syy(i) = pl%syy(i) + dt*pl%dsyy(i)
         elseif(pl%stress_integration==2)then
            pl%str_min%x = pl%str%x
            pl%str_min%xy = pl%str%xy
            pl%str_min%y = pl%str%y
            pl%p_min   = pl%p
            pl%str%x = pl%str%x + (pl%dt/2.)*pl%dstr%x
            pl%str%xy = pl%str%xy + (pl%dt/2.)*pl%dstr%xy
            pl%str%y = pl%str%y + (pl%dt/2.)*pl%dstr%y
            pl%p = pl%p + (pl%dt/2.)*pl%dp      !!!simultaneous pressure
         endif
          
      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i); v_mini = pl%v_min%cmpt(i)
         do d = 1, pl%dim
         !   pl%v_min(d, i) = pl%vx(d, i)
         !   pl%vx(d, i) = pl%vx(d, i) + (dt/2.)*pl%dvx(d, i)
            v_mini(d)%p = vxi(d)%p
            vxi(d)%p = vxi(d)%p + (pl%dt/2.)*dvxi(d)%p                  
         enddo
            !pl%v_min(1, i) = pl%vx(1, i)
            !pl%vx(1, i) = pl%vx(1, i) + (dt/2.)*pl%dvx%x%r(i)
            !pl%v_min(2, i) = pl%vx(2, i)
            !pl%vx(2, i) = pl%vx(2, i) + (dt/2.)*pl%dvx%y%r(i)
      enddo 

      call drucker_prager_failure_criterion(pl)

   endif

   call single_step_for_soil

   if(parts%itimestep .eq. 1)then
      pl => parts
!           call first_step
!      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         
         if(.not.pl%summation_density )then
            pl%rho = pl%rho + (pl%dt/2.)* pl%drho
         endif

         if(pl%stress_integration==1)then
            pl%str%x = pl%str%x + pl%dt*pl%dstr%x
            pl%str%xy = pl%str%xy + pl%dt*pl%dstr%xy
            pl%str%y = pl%str%y + pl%dt*pl%dstr%y
            pl%p   = pl%p   + pl%dt*pl%dp       !!! simultaneous pressure
         elseif(pl%stress_integration==2)then
            pl%str%x = pl%str%x + (pl%dt/2.)*pl%dstr%x
            pl%str%xy = pl%str%xy + (pl%dt/2.)*pl%dstr%xy
            pl%str%y = pl%str%y + (pl%dt/2.)*pl%dstr%y
            pl%p   = pl%p + (pl%dt/2.)*pl%dp    !!! simultaneous pressure
         endif
        
      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i); avi = pl%av%cmpt(i) 
         do d = 1, pl%dim        
            vxi(d)%p = vxi(d)%p + (pl%dt/2.) * dvxi(d)%p + avi(d)%p
            pl%x(d, i) = pl%x(d, i) + pl%dt * vxi(d)%p
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
!      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
            
         if(.not.pl%summation_density )then 
            pl%rho = pl%rho_min + pl%dt*pl%drho
         endif
                 
         if(pl%stress_integration==1)then 
            pl%str%x = pl%str%x + pl%dt*pl%dstr%x
            pl%str%xy = pl%str%xy + pl%dt*pl%dstr%xy
            pl%str%y = pl%str%y + pl%dt*pl%dstr%y
            pl%p  = pl%p + pl%dt*pl%dp       !!! simultaneous pressure
         elseif(pl%stress_integration==2)then
            pl%str%x = pl%str_min%x + pl%dt*pl%dstr%x
            pl%str%xy = pl%str_min%xy + pl%dt*pl%dstr%xy
            pl%str%y = pl%str_min%y + pl%dt*pl%dstr%y
            pl%p   = pl%p_min + pl%dt*pl%dp     !!! simultaneous pressure
         endif

      do i=1,pl%ntotal +pl%nvirt  ! origionally pl%ntotal            
         if(pl%itype(i)<0)cycle
         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i)
         !v_mini = pl%vx%cmpt(i)
         v_mini = pl%v_min%cmpt(i)
         avi = pl%av%cmpt(i)
         do d = 1, pl%dim                   
            vxi(d)%p = v_mini(d)%p + pl%dt * dvxi(d)%p   &
                        + avi(d)%p
            pl%x(d, i) = pl%x(d, i) + pl%dt * vxi(d)%p                  
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

   parts%time = parts%time + parts%dt

   if(mod(parts%itimestep,parts%print_step).eq.0)then
      write(*,*)'______________________________________________'
      write(*,*)'  current number of time step =',              &
                parts%itimestep,'     current time=', real(parts%time+parts%dt)
      write(*,*)'______________________________________________'
   endif  
        
   if(parts%itimestep>=parts%save_step_from.and.mod(parts%itimestep,parts%save_step).eq.0) then
      call output
   endif 



enddo

return
end subroutine


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
type(material), pointer :: mat
logical :: dbg = .false.
integer i, ntotal
                 
if(dbg) write(*,*) 'In single_step...'

                                      nphase = 2
                                      if(parts%single_phase) nphase=1

                                      do iphase = 1, nphase
                                         if(iphase==1) pl => parts  !!!
                                         if(iphase==2) pl => soil   !!!

pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
if(trim(pl%imaterial)=='water')pl%dvof = 0.d0
if(trim(pl%imaterial)=='soil')then
   pl%dstr%x = 0.d0; pl%dstr%xy = 0.d0; pl%dstr%y = 0.d0; pl%dp = 0.d0
   if(pl%critical_state==1) pl%dvof = 0.d0
endif

if(trim(pl%imaterial)=='water') call inlet_boundary2
call symmetry_boundary(pl)   !!! subjetting
!if(iphase==2) call slip_boundary(pl)   !!! submerged column collapse
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
if (pl%nnps.eq.1) then 
   call direct_find(pl)
else if (pl%nnps.eq.2) then
   call link_list(pl)     
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
!        call link_list(itimestep, parts%ntotal+parts%nvirt,
!     &       parts%hsml(1),parts%x,parts%niac,parts%pair_i,
!     &       parts%pair_j,parts%w,parts%dwdx,parts%countiac)
else if (pl%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
endif         

call count_interaction(pl)   ! Added for particle shifting

if(mod(pl%itimestep,pl%print_step).eq.0.and.pl%int_stat) then
   call pl%interaction_statistics
endif   
!write(*,*) associated(pl%wi)
call Sherpard_filter(pl)

!--- Added by Wang
!if(nor_density) call norm_density(pl)

!---  Density approximation or change rate
     
!if(summation_density)then      
!if(mod(pl%itimestep,30)==0) call sum_density(pl)
!call sum_density(pl)
!else             
    !call con_density(pl)         
    pl%drho = -pl%rho*pl%div2(pl%vx)
!endif
    

if(pl%artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      call delta_sph_omp(pl,pl%rho,pl%drho)
      !call output_field(pl,pl%rho)
      !call output_field(pl,pl%drho)
      !pause
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
!pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
   call volume_fraction_soil_critical_state(pl)
   !write(*,*) 'sdfa'
   !write(*,*) pl%vof%r(1:10)
   !write(*,*) pl%psi%r(1:10)
   !write(*,*) pl%dvof%r(1:10)
endif

if(trim(pl%imaterial)=='soil') pl%vcc = pl%div_omp(pl%vx)
!call velocity_divergence(pl)

call pressure(pl)

if(trim(pl%imaterial)=='water')then
   call newtonian_fluid(pl)
elseif(trim(pl%imaterial)=='soil')then

   if(pl%yield_criterion == 1)then
      call mohr_coulomb_failure_criterion(pl)
   elseif(pl%yield_criterion == 2)then
           if(pl%critical_state==0)then
      call drucker_prager_failure_criterion(pl)
           else
      call drucker_prager_failure_criterion_critical_state(pl)
           endif
   endif

   if(mod(pl%itimestep,pl%print_step).eq.0)    &
            write(*,*) 'Failured points: ', pl%nfail
endif

!Calculate internal force for water phase !! -phi_f Grad(p)
if(pl%imaterial=='water')then
   !if(pl%nthreads==1)then
!   pl%dvx%x%r = -pl%vof%r*pl%df(pl%p%r,'x') + pl%df(pl%vof%r*pl%str%x%r,'x') + pl%df(pl%vof%r*pl%str%xy%r,'y')
!   pl%dvx%y%r = -pl%vof%r*pl%df(pl%p%r,'y') + pl%df(pl%vof%r*pl%str%xy%r,'x') + pl%df(pl%vof%r*pl%str%y%r,'y')

!!! Replaced by Laplace
!   pl%dvx%x = -pl%vof*pl%df(pl%p,'x') + pl%df(pl%vof*pl%str%x,'x') + pl%df(pl%vof*pl%str%xy,'y')
!   pl%dvx%y = -pl%vof*pl%df(pl%p,'y') + pl%df(pl%vof*pl%str%xy,'x') + pl%df(pl%vof*pl%str%y,'y')   
   mat => pl%material
   pl%dvx%x = -pl%vof*pl%df(pl%p,'x') + pl%lap_omp(mat%viscosity,pl%vof*pl%vx%x)
   pl%dvx%y = -pl%vof*pl%df(pl%p,'y') + pl%lap_omp(mat%viscosity,pl%vof*pl%vx%y)      

   !else
   !pl%dvx%x = -pl%vof*pl%df_omp(pl%p,'x') + pl%df_omp(pl%vof*pl%str%x,'x') + pl%df_omp(pl%vof*pl%str%xy,'y')
   !pl%dvx%y = -pl%vof*pl%df_omp(pl%p,'y') + pl%df_omp(pl%vof*pl%str%xy,'x') + pl%df_omp(pl%vof*pl%str%y,'y')  
   !endif        
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
   !where (pl%rho%r.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
   !where (pl%rho%r.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho

!   call output_field(pl,pl%df(pl%p,'x'))
!   call output_field(pl,pl%df(pl%p,'y'))
!   call output_field(pl,pl%lap_omp(mat%viscosity,pl%vof*pl%vx%x))
!   call output_field(pl,pl%lap_omp(mat%viscosity,pl%vof*pl%vx%y))
!   call output_field(pl,pl%dvx%x)
!   call output_field(pl,pl%dvx%y)
!   pause

   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)
else      
   !call int_force(pl)
   pl%dvx%x = -pl%df3_omp(pl%vof*pl%p,'x') + pl%df3_omp(pl%vof*pl%str%x,'x') + pl%df3_omp(pl%vof*pl%str%xy,'y')
   pl%dvx%y = -pl%df3_omp(pl%vof*pl%p,'y') + pl%df3_omp(pl%vof*pl%str%xy,'x') + pl%df3_omp(pl%vof*pl%str%y,'y')
   !where (pl%rho%r.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
   !where (pl%rho%r.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r 
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho 
!   call output_field(pl,pl%dvx%x)
!   call output_field(pl,pl%dvx%y)
!   pause
endif      

!if(trim(pl%imaterial)=='water'.and.water_tension_instability==2) &
!   call tension_instability(pl) 

! --- Plasticity flow rule   ! This was done before Jaummann_rate, because we 
!     we need txx,tyy,tzz, which was destroyed in Jaumann_rate!

if(trim(pl%imaterial)=='soil')then   
   if(pl%plasticity==1)then 
      call plastic_flow_rule(pl)
   elseif(pl%plasticity==2)then     
      call plastic_flow_rule2(pl)
   elseif(pl%plasticity==3)then

           if(pl%critical_state==0)then
      call plastic_or_not(pl)
           else
      call plastic_or_not_critical_state(pl)
           endif

      if(pl%flow_rule==1) call plastic_flow_rule3(pl)
      if(pl%flow_rule==2)then
         if(pl%critical_state==0)then
              call non_associated_plastic_flow_rule3(pl)
         else
              call non_associated_flow_rule_critical_state(pl)
         endif
      endif        
   endif
endif

! --- Jaumann rate  !When???

if(trim(pl%imaterial)=='soil')then
   call Jaumann_rate(pl)
!   call output_field(pl,pl%dstr%x)
!   call output_field(pl,pl%dstr%y)
!   call output_field(pl,pl%dstr%xy)
!   call output_field(pl,pl%dp)
!   call output_field(pl,pl%vcc)
!   call output_field(pl,pl%wxy)
!   pause
endif

!---  Artificial viscosity:

if (pl%visc_artificial) call pl%art_visc_omp

if(trim(pl%imaterial)=='soil'.and.pl%soil_artificial_stress)then
        !call art_stress(pl)
   call pl%delta_sph_omp(pl%p,pl%dp)             !! Jump
   call pl%delta_sph_omp(pl%str%x,pl%dstr%x)     !! Jump
   call pl%delta_sph_omp(pl%str%xy,pl%dstr%xy)   !! No jump
   call pl%delta_sph_omp(pl%str%y,pl%dstr%y)     !! Jump
   !call output_field(pl,pl%dstr%x)
   !call output_field(pl,pl%dstr%y)
   !call output_field(pl,pl%dstr%xy)
   !call output_field(pl,pl%dp)
   !call output_field(pl,pl%wi)
   !pause   
endif        
if(trim(pl%imaterial)=='water'.and.pl%water_artificial_volume)then
        !call art_volume_fraction_water2(pl)
        !call output_field(pl,pl%vof)
        call pl%delta_sph_omp(pl%vof,pl%dvof)
!        call output_field(pl,pl%dvof)
!        pause
endif        

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
          call pl%repulsive_force_omp           !!!!!!!#########
!      endif

pl%dvx%y = pl%dvx%y + pl%gravity

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (pl%average_velocity) call av_vel_omp(pl) !!!!!!!##########
if(iphase==1.and.pl%particle_shift)then
   call particle_shifting(pl)
endif   

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(pl%itimestep,pl%print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

                                            enddo ! iphase

!call drag_force(parts,soil)   !!! Porous media

                if(.not.parts%single_phase)then

!-------------------Water/soil interaction-------------------------------

!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
      if (parts%nnps.eq.1) then 
         call direct_find_2(parts,soil)
      else if (parts%nnps.eq.2) then
          call link_list2_omp(parts,soil)
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      !else if (numeric%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      endif         

      if(mod(parts%itimestep,parts%print_step).eq.0.and.parts%int_stat) then
         call parts%interaction_statistics
      endif    

      call darcy_law_omp(parts,soil)          
      if(parts%water_pressure==1)then 
         call interpolation(parts,parts%p,soil,soil%spp)
         !call output_field(soil,soil%spp)
         !pause
         call pore_water_pressure_omp(parts,soil)
      else  !!Hydrostatic pressure
         soil%dvx%y = soil%dvx%y - pl%gravity*(1./2.7)
      endif

      if(parts%volume_fraction)then
              if(parts%critical_state==0)then
         call volume_fraction_soil_omp(soil)
              else
         !call volume_fraction_soil_critical_state(soil)   !!! Put forward! 
              endif
         call volume_fraction_water2_omp(parts,soil)  !! water%dvof
         call volume_fraction_water_omp(parts,soil)  ! phi_f = 1- phi_s, vof2
         !call volume_fraction_water_omp_with_filter(parts,soil)  ! phi_f = 1- phi_s
         if(parts%volume_fraction_renorm)then
            if(mod(parts%itimestep,40).eq.0) then
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

      if(parts%itimestep>=parts%save_step_from.and.   &
         mod(parts%itimestep,parts%save_step).eq.0)then
         call output
      endif 

      parts%nvirt = parts%nvirt - parts%nsymm
      soil%nvirt = soil%nvirt - soil%nsymm

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
  
if (pl%nnps.eq.1) then 
   call direct_find(pl)
else if (pl%nnps.eq.2) then
   call link_list(pl)     
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
!        call link_list(itimestep, parts%ntotal+parts%nvirt,
!     &       parts%hsml(1),parts%x,parts%niac,parts%pair_i,
!     &       parts%pair_j,parts%w,parts%dwdx,parts%countiac)
else if (pl%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
endif         

if(mod(pl%itimestep,pl%print_step).eq.0.and.pl%int_stat) then
   call pl%interaction_statistics
endif   

!--- Added by Wang
!if(nor_density) call norm_density(pl)

!---  Density approximation or change rate
     
!if(summation_density)then      
!if(mod(pl%itimestep,30)==0) call sum_density(pl)
!else             
!    call con_density(pl)         
    pl%drho = -pl%rho*pl%div2(pl%vx)
!endif
      
if(pl%artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      call delta_sph_omp(pl,pl%rho,pl%drho)
   !endif
endif

!---  Dynamic viscosity:
     


water => parts%material
parts%p = water%b*((parts%rho/(water%rho0))**water%gamma-1.d0)

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

!Calculate internal force for water phase !! -phi_f Grad(p)

   pl%dvx%x = -pl%df(pl%p,'x') + pl%df(pl%str%x,'x') + pl%df(pl%str%xy,'y')
   pl%dvx%y = -pl%df(pl%p,'y') + pl%df(pl%str%xy,'x') + pl%df(pl%str%y,'y')   
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho
   !write(*,*) pl%dvx%x%r(1:50),pl%dvx%y%r(1:50)


!---  Artificial viscosity:

if (pl%visc_artificial) call pl%art_visc_omp
       
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

pl%dvx%y = pl%dvx%y + pl%gravity

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (pl%average_velocity) call av_vel(pl) 

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(pl%itimestep,pl%print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

      if(pl%itimestep>=pl%save_step_from.and.   &
         mod(pl%itimestep,pl%save_step).eq.0)then
         call output
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
integer i

pl => parts
property => pl%material

pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
pl%dstr%x = 0.d0; pl%dstr%xy = 0.d0; pl%dstr%y = 0.d0
pl%dp = 0.d0
!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length

call pl%find_pairs

if(mod(pl%itimestep,pl%print_step).eq.0.and.pl%int_stat) then
   call pl%interaction_statistics
endif  

! For USAW
if(pl%usaw)then
   call Sherpard_filter(pl)
   call get_boundary_field(pl,pl%rho)
endif

!---  Density approximation or change rate
     
if (pl%summation_density) then      
    call sum_density(pl)
else             
    !call con_density(pl)
    pl%drho = -pl%rho*pl%div2_USAW(pl%vx)    
endif



      
if(pl%artificial_density)then
   !if(trim(pl%imaterial)=='water')then
      !!call renormalize_density_gradient(pl)
      !call art_density(pl)
      !if(.not.pl%usaw)then
      call delta_sph_omp(pl,pl%rho,pl%drho)
!      call diffusive_sph_omp(pl,pl%rho,pl%drho)
      !endif
   !endif
endif

!call shear_modulus(pl)
       
!---  Internal forces:

!call shear_strain_rate(pl)   
!call velocity_divergence(pl)

pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
!write(*,*) pl%tab%x%ndim1,pl%vx%x%ndim1
pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4_omp(pl%vx%x,'x')-pl%df4_omp(pl%vx%y,'y'))
pl%tab%xy = pl%df4_omp(pl%vx%x,'y')+pl%df4_omp(pl%vx%y,'x')
pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4_omp(pl%vx%y,'y')-pl%df4_omp(pl%vx%x,'x'))

pl%vcc = pl%div_omp(pl%vx)

if(pl%usaw)then
   call get_boundary_field(pl,pl%tab%x)
   call get_boundary_field(pl,pl%tab%xy)
   call get_boundary_field(pl,pl%tab%y)
   call get_boundary_field(pl,pl%vcc)
endif


call pressure(pl)

if(pl%usaw)then
   call get_boundary_field(pl,pl%str%x)
   call get_boundary_field(pl,pl%str%xy)
   call get_boundary_field(pl,pl%str%y)
   call get_boundary_field(pl,pl%p)
endif


if(pl%yield_criterion == 1)then
   call mohr_coulomb_failure_criterion(pl)
elseif(pl%yield_criterion == 2)then
   call drucker_prager_failure_criterion(pl)
endif

!call int_force1(pl)

pl%dvx%x = -pl%df3_omp(pl%p,'x') + pl%df3_omp(pl%str%x,'x') + pl%df3_omp(pl%str%xy,'y')
pl%dvx%y = -pl%df3_omp(pl%p,'y') + pl%df3_omp(pl%str%xy,'x') + pl%df3_omp(pl%str%y,'y')

     if(pl%usaw)then

do i = 1, pl%ntotal + pl%nvirt
   if(pl%itype(i)<0)cycle
   pl%dvx%x%r(i) = pl%dvx%x%r(i)/pl%rho%r(i)
   pl%dvx%y%r(i) = pl%dvx%y%r(i)/pl%rho%r(i)
enddo

     else

pl%dvx%x = pl%dvx%x/pl%rho
pl%dvx%y = pl%dvx%y/pl%rho     

     endif

! --- Plasticity flow rule   ! This was done before Jaummann_rate, because we 
!     we need txx,tyy,tzz, which was destroyed in Jaumann_rate!

if(pl%plasticity==1)then 
   call plastic_flow_rule(pl)
elseif(pl%plasticity==2)then     
   call plastic_flow_rule2(pl)
elseif(pl%plasticity==3)then
   call plastic_or_not(pl)
!   call plastic_flow_rule3(pl)
   call non_associated_plastic_flow_rule3(pl)
endif

! --- Jaumann rate  !When???

call Jaumann_rate(pl)

!                                 if(.not.pl%usaw)then

!---  Artificial viscosity:

if (pl%visc_artificial) call pl%art_visc_omp


!if(soil_artificial_stress) call art_stress(pl)

if(pl%soil_artificial_stress) call artificial_stress(pl)   ! LU Shi Yang

!                                 endif

if(pl%soil_artificial_stress)then
   call pl%delta_sph_omp(pl%p,pl%dp)
   call pl%delta_sph_omp(pl%str%x,pl%dstr%x)
   call pl%delta_sph_omp(pl%str%xy,pl%dstr%xy)
   call pl%delta_sph_omp(pl%str%y,pl%dstr%y)
!   call pl%diffusive_sph_omp(pl%p,pl%dp)
!   call pl%diffusive_sph_omp(pl%str%x,pl%dstr%x)
!   call pl%diffusive_sph_omp(pl%str%xy,pl%dstr%xy)
!   call pl%diffusive_sph_omp(pl%str%y,pl%dstr%y)
endif   

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
!      endif

pl%dvx%y = pl%dvx%y + pl%gravity
! Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)
     
! Calculating average velocity of each partile for avoiding penetration

!                                                if(.not.pl%usaw)then
if (pl%average_velocity) call av_vel_omp(pl) 
!                                                endif

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(pl%itimestep,pl%print_step).eq.0) then     
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

call get_unit(f_tmp)
open(f_tmp,file="tmp.dat")

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
close(f_tmp)

return
end subroutine

!-----------------------------------
    subroutine output_parameters
!-----------------------------------
implicit none
type(numerical), pointer :: numeric
type(material), pointer :: matt
numeric => parts%numeric
matt => parts%material

write(f_other,*) 'Numerical parameters:'
write(f_other,*) 'nnps = ', parts%nnps 
write(f_other,*) 'gravity = ', parts%gravity 
write(f_other,*) 'alpha, beta, etq = ', numeric%alpha, numeric%beta, numeric%etq
write(f_other,*) 'dd, p1, p2 = ', numeric%dd, numeric%p1, numeric%p2
write(f_other,*) 'delta = ', numeric%delta
write(f_other,*) 'epsilon = ', numeric%epsilon

write(f_other,*) 'Background grid range:'
write(f_other,*) 'x_maxgeom, x_mingeom = ', parts%x_maxgeom, parts%x_mingeom
write(f_other,*) 'y_maxgeom, y_mingeom = ', parts%y_maxgeom, parts%y_mingeom
write(f_other,*) 'z_maxgeom, z_mingeom = ', parts%z_maxgeom, parts%z_mingeom
write(f_other,*) 'maxngx, maxngy, maxngz = '
write(f_other,*) parts%maxngx, parts%maxngy, parts%maxngz

write(f_other,*) ' '
write(f_other,*) 'maxn, max_interaction = ', parts%maxn, parts%max_interaction
write(f_other,*) 'pa_sph = ', parts%pa_sph
write(f_other,*) 'sle = ', parts%sle
write(f_other,*) 'skf = ', parts%skf

write(f_other,*) 'dt, print step, save step: '
write(f_other,*) parts%dt, parts%print_step, parts%save_step

if(parts%imaterial=='water')then
   write(f_other,*) 'rho0, b, gamma, c, viscosity = '
   write(f_other,*) matt%rho0, matt%b, matt%gamma, matt%c, matt%viscosity
elseif(parts%imaterial=='soil')then
   write(f_other,*) 'rho0,k,Youngs,Poissons,c = '
   write(f_other,*) matt%rho0,matt%k,matt%E,matt%niu,matt%c
   write(f_other,*) 'porosity, permeability = '
   write(f_other,*) matt%porosity, matt%permeability
   write(f_other,*) 'Cohesion, phi = '
   write(f_other,*) matt%cohesion, matt%phi
endif

write(f_other,*) 'single phase = ', parts%single_phase

if(parts%single_phase) return
matt => soil%material

   write(f_other,*) 'rho0,k,Youngs,Poissons,c = '
   write(f_other,*) matt%rho0,matt%k,matt%E,matt%niu,matt%c
   write(f_other,*) 'porosity, permeability = '
   write(f_other,*) matt%porosity, matt%permeability
   write(f_other,*) 'Cohesion, phi = '
   write(f_other,*) matt%cohesion, matt%phi, matt%psi

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
!      write(f_xv,*) 'ZONE I=', ntotal, ' F=BLOCK'
      write(f_xv,*) 'ZONE T="time ', parts%time, '"', ' I=', ntotal, ' F=BLOCK'
      
      write(f_xv,*) (x(1,i),i=1,ntotal) !, soil%x(1,1:ntotal2)
      write(f_xv,*) (x(2,i),i=1,ntotal) !, soil%x(2,1:ntotal2)
!      write(f_xv,*) vx(1,1:ntotal)
!      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
!      write(f_xv,*) mass(1:ntotal)
      write(f_xv,*) (p(i),i=1,ntotal) !, soil%p(1:ntotal2)
      write(f_xv,*) (parts%vof%r(i),i=1,ntotal)
!      write(f_xv,*) u(1:ntotal)
!      write(f_xv,*) itype(1:ntotal)
!      write(f_xv,*) hsml(1:ntotal)                                        
      write(f_xv,*) (parts%vx%x%r(i),i=1,ntotal)
      write(f_xv,*) (parts%vx%y%r(i),i=1,ntotal)
      write(f_xv,*) (parts%rho%r(i),i=1,ntotal)
      write(f_xv,*) (parts%zone(i),i=1,ntotal)
      write(f_xv,*) (parts%vof2%r(i),i=1,ntotal)
      write(f_xv,*) (parts%mass%r(i),i=1,ntotal)

             elseif(trim(parts%imaterial)=='soil')then

      write(f_xv,*) 'VARIABLES="X","Y","Pressure", "VoF", "ep",' 
      write(f_xv,*) '"sxy", "sxx","syy","vx","vy", "rho", "mass", '
      write(f_xv,*) '"sigma_yy", "zone","filter " '
!      write(f_xv,*) 'ZONE I=', ntotal, ' F=BLOCK'
      write(f_xv,*) 'ZONE T="time ', parts%time, '"', ' I=', ntotal, ' F=BLOCK'
 
      write(f_xv,*)  (parts%x(1,i),i=1,ntotal)
      write(f_xv,*)  (parts%x(2,i),i=1,ntotal)
!      write(f_xv,*) vx(1,1:ntotal)
!      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
!      write(f_xv,*) mass(1:ntotal)
      write(f_xv,*)  (parts%p%r(i),i=1,ntotal)
!      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_xv,*)  (parts%vof%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%epsilon_p%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%str%xy%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%str%x%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%str%y%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%vx%x%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%vx%y%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%rho%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%mass%r(i),i=1,ntotal)
      write(f_xv,*)  (-parts%p%r(i) + parts%str%y%r(i),i=1,ntotal)
      write(f_xv,*)  (parts%zone(i),i=1,ntotal)
      write(f_xv,*)  (parts%wi%r(i),i=1,ntotal)
      write(f_other,*) parts%time, -parts%p%r(395)+parts%str%y%r(395)

             endif

      if(parts%single_phase) return

      write(f_state,*) 'VARIABLES="X","Y","Pressure", "VoF", "ep", '
      write(f_state,*) '"sxy", "sxx","syy","vx","vy", "rho","mass" '
      write(f_state,*) '"sigma_yy", "zone" '
      if(soil%critical_state==1) write(f_state,*) ' "psi" '
!      write(f_state,*) 'ZONE I=', ntotal2, ' F=BLOCK'
      write(f_state,*) 'ZONE T="time ', soil%time, '"', ' I=', ntotal2, ' F=BLOCK'
 
      write(f_state,*)  (soil%x(1,i),i=1,ntotal2)
      write(f_state,*)  (soil%x(2,i),i=1,ntotal2)
!      write(f_xv,*) vx(1,1:ntotal)
!      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
!      write(f_xv,*) mass(1:ntotal)
      write(f_state,*)  (soil%p%r(i),i=1,ntotal2)
!      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_state,*)  (soil%vof%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%epsilon_p%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%str%xy%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%str%x%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%str%y%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%vx%x%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%vx%y%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%rho%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%mass%r(i),i=1,ntotal2)
      write(f_state,*)  (-soil%p%r(i) + soil%str%y%r(i),i=1,ntotal2)
      write(f_state,*)  (soil%zone(i),i=1,ntotal2)
      !write(f_other,*) time, -soil%p(420)+soil%syy(420)
      !write(f_other,*) time, -soil%p(395)+soil%syy(395)
      !write(f_other,*) time, -parts%p(420)+parts%syy(420)

      if(soil%critical_state==1) write(f_state,*) (soil%psi%r(i), i= 1, ntotal2)
      !write(*,*) soil%psi%r(1:10)
      return
      end subroutine

!--------------------------------------------------
    subroutine output_field(this,field)
!--------------------------------------------------
implicit none
type(particles) this
type(array) field
integer i, ntotal

ntotal = this%ntotal+this%nvirt

      write(f_tmp,*) 'VARIABLES="X","Y","field" '     
      write(f_tmp,*) 'ZONE I=', ntotal, ' F=BLOCK'
      write(f_tmp,*) (this%x(1,i),i=1,ntotal) 
      write(f_tmp,*) (this%x(2,i),i=1,ntotal)
      write(f_tmp,*) (field%r(i),i=1,ntotal) 


return
end subroutine

!---------------------------------------------------------------------
!      subroutine output_cas
!---------------------------------------------------------------------
!implicit none

!write(f_cas,*) ''single_phase

!return
!end subroutine

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
!   subroutine time_elapsed(s)
!------------------------------
!    use ifport
!    implicit none

!    real(8) :: s

!   s = rtc()

!   end subroutine

end module
    
#ifdef REMOVE    
!!DEC$IF(.FALSE.)
!----------------------------------------------------------------------      
                   subroutine single_step_VOF
!----------------------------------------------------------------------
!   Subroutine to determine the right hand side of a differential 
!   equation in a single step for performing time integration 
!----------------------------------------------------------------------

                  
!use param 
use declarations_sph
!use m_sph_fo
!implicit none
PARAMETER (IBAR2=400,JBAR2=320,NPRTS=1,MESHX=1,MESHY=1,NVOR=10)
      common /two_fluid/ phif(ibar2,jbar2), phifn(ibar2,jbar2)
      COMMON /FV/ ACOM(1), UN(IBAR2,JBAR2), VN(IBAR2,JBAR2),   &
       PN(IBAR2,JBAR2),                                      &
       FN(IBAR2,JBAR2), U(IBAR2,JBAR2), V(IBAR2,JBAR2),      &
       P(IBAR2,JBAR2),                                       &
       F(IBAR2,JBAR2), PETA(IBAR2,JBAR2), BETA(IBAR2,JBAR2), &
       NF(IBAR2,JBAR2), TANTH(IBAR2,JBAR2), PS(IBAR2,JBAR2)  
target :: P, phif
integer  nphase, iphase
type(particles), pointer :: pl
logical :: dbg = .false.
integer i, ntotal
type(array) spp  ! soil particle pressure
!real(dp) p(:,:) !,phif(:,:)
type(material), pointer :: property
type(realarray) pre,phiff
                 
if(dbg) write(*,*) 'In single_step...'

                                      nphase = 2
                                      if(parts%single_phase) nphase=1

                                      do iphase = 1, nphase
                                         if(iphase==1)then
                                             pl => parts  !!!
                                             cycle 
                                         endif    
                                         if(iphase==2) pl => soil   !!!

pl%dvx%x = 0.d0; pl%dvx%y = 0.d0; pl%drho = 0.d0
if(trim(pl%imaterial)=='water')pl%dvof = 0.d0
if(trim(pl%imaterial)=='soil')then
   pl%dstr%x = 0.d0; pl%dstr%xy = 0.d0; pl%dstr%y = 0.d0; pl%dp = 0.d0
endif

!if(trim(pl%imaterial)=='water') call inlet_boundary2

!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
if (pl%nnps.eq.1) then 
   call direct_find(pl)
else if (pl%nnps.eq.2) then
   call link_list(pl)     
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
!        call link_list(itimestep, parts%ntotal+parts%nvirt,
!     &       parts%hsml(1),parts%x,parts%niac,parts%pair_i,
!     &       parts%pair_j,parts%w,parts%dwdx,parts%countiac)
else if (pl%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
endif         

if(mod(pl%itimestep,pl%print_step).eq.0.and.pl%int_stat) then
   call pl%interaction_statistics
endif   

! GET water-pore-pressure for each soil particles (For VOF method)
!call spp%set(name='spp',parts=soil)

!call spp%max_alloc
pre%r=>p
phiff%r=>phif
!write(*,*) pre%r
call cell2parts(soil,pre,soil%spp)
!write(*,*) 'aa',soil%spp%r(1:10)
phif = 0.0
call parts2cell(soil,soil%vof,phiff)
phif = 1.0-phif
!write(*,*) 'bb', phif(8,2),phifn(8,2)

!--- Added by Wang
!if(nor_density) call norm_density(pl)

!---  Density approximation or change rate
     
!if(summation_density)then      
!if(mod(itimestep,30)==0) call sum_density(pl)
!else             
    !call con_density(pl)         
    pl%drho = -pl.rho*pl.div2(pl.vx)
!endif
      
if(pl%artificial_density)then
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
!pl%tab%x%ndim1 = pl%ntotal+pl%nvirt
!pl%tab%xy%ndim1 = pl%tab%x%ndim1; pl%tab%y%ndim1 = pl%tab%x%ndim1
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

   if(pl%yield_criterion == 1)then
      call mohr_coulomb_failure_criterion(pl)
   elseif(pl%yield_criterion == 2)then
      call drucker_prager_failure_criterion(pl)
   endif

   if(mod(pl%itimestep,pl%print_step).eq.0)    &
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
   if(pl%plasticity==1)then 
      call plastic_flow_rule(pl)
   elseif(pl%plasticity==2)then     
      call plastic_flow_rule2(pl)
   elseif(pl%plasticity==3)then
      call plastic_or_not(pl)
      call plastic_flow_rule3(pl)
   endif
endif

! --- Jaumann rate  !When???

if(trim(pl%imaterial)=='soil')call Jaumann_rate(pl)

!---  Artificial viscosity:

if (pl%visc_artificial) call pl%art_visc_omp
!write(*,*) '0002',soil%itimestep, soil%dvx%x%r(7),soil%dvx%y%r(7)
if(trim(pl%imaterial)=='soil'.and.pl%soil_artificial_stress)then
        !call art_stress(pl)
   call pl%delta_sph_omp(pl%p,pl%dp)
   call pl%delta_sph_omp(pl%str%x,pl%dstr%x)
   call pl%delta_sph_omp(pl%str%xy,pl%dstr%xy)
   call pl%delta_sph_omp(pl%str%y,pl%dstr%y)
endif        
if(trim(pl%imaterial)=='water'.and.pl%water_artificial_volume)  &
        !call art_volume_fraction_water2(pl)
        call pl%delta_sph_omp(pl%vof,pl%dvof)

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
          call pl%repulsive_force_omp           !!!!!!!#########
!      endif

pl%dvx%y = pl%dvx%y + pl%gravity

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (pl%average_velocity) call av_vel_omp(pl) !!!!!!!##########

!---  Convert velocity, force, and energy to f and dfdt  
      
if(mod(pl%itimestep,pl%print_step).eq.0) then     
!  call pl%particle_monitor
   call pl%minimum_time_step  
endif

                                            enddo ! iphase

!call drag_force(parts,soil)   !!! Porous media

                if(.not.parts%single_phase)then

!-------------------Water/soil interaction-------------------------------

!---  Interaction parameters, calculating neighboring particles
!     and optimzing smoothing length
  
      if (parts%nnps.eq.1) then 
         !call direct_find_2(parts,soil)
      else if (parts%nnps.eq.2) then
          !call link_list2_omp(parts,soil)
!        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      !else if (numeric%nnps.eq.3) then 
!        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
!     &       pair_j,w,dwdx,ns)
      endif         

      if(mod(parts%itimestep,parts%print_step).eq.0.and.parts%int_stat) then
         call parts%interaction_statistics
      endif    

!      write(*,*) '0001',soil%itimestep, soil%dvx%x%r(7),soil%dvx%y%r(7)
      !call darcy_law_omp(parts,soil)          
      !call pore_water_pressure_omp(parts,soil)
      property => soil%material
      soil%dvx%x = soil%dvx%x - soil%df(soil%spp,'x')/property%rho0
      soil%dvx%y = soil%dvx%y - soil%df(soil%spp,'y')/property%rho0
!      write(*,*) '0000',soil%itimestep, soil%dvx%x%r(7),soil%dvx%y%r(7)
      
      

      !if(parts%volume_fraction)then
         call volume_fraction_soil_omp(soil)
      !   call volume_fraction_water2_omp(parts,soil)
      !   call volume_fraction_water_omp(parts,soil)  ! phi_f = 1- phi_s
      !   if(parts%volume_fraction_renorm)then
      !      if(mod(parts%itimestep,40).eq.0) then
      !         ntotal = parts%ntotal+parts%nvirt
      !         parts%rho = parts%rho/parts%vof 
      !         parts%vof = parts%vof2
      !         parts%rho = parts%rho*parts%vof 
      !      endif
      !   endif
      !endif

!      call direct_find(parts) 
!      call art_volume_fraction_water2(parts)
     
             endif ! .not.single_phase

      if(parts%itimestep>=parts%save_step_from.and.   &
         mod(parts%itimestep,parts%save_step).eq.0)then
         call output
         endif 
         
         
return
end subroutine

!!DEC$ENDIF
#endif
