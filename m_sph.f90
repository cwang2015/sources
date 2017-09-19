!-----------------------------
    module declarations_sph
!-----------------------------
#ifdef _OPENACC
  use openacc
#endif

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
!$acc declare create(parts,soil)

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

      case('SOIL MODEL')
          read(pvalu,*) parts%soil_model
          write(*,*) 'Soil model: '
          write(*,*) parts%soil_model
          soil%soil_model = parts%soil_model

      case('CAMBRIDGE MODEL1')
          !read(pvalu,*) sio2%zramda,sio2%zkapa,sio2%eps0,sio2%poi,sio2%zm, &
          !              sio2%zms,sio2%zbr,sio2%zmb,sio2%zmf,sio2%sig00,    &
          !              sio2%ocr,sio2%rs0,sio2%zeta0
          read(pvalu,*) sio2%zramda,sio2%zkapa,sio2%eps0,sio2%poi,sio2%zm, &
                        sio2%zms

          write(*,*) 'Cambridge model1: '
          write(*,*) sio2%zramda,sio2%zkapa,sio2%eps0,sio2%poi,sio2%zm, &
                        sio2%zms

      case('CAMBRIDGE MODEL2')
          read(pvalu,*) sio2%zbr,sio2%zmb,sio2%zmf,sio2%sig00,    &
                        sio2%ocr,sio2%rs0,sio2%zeta0

          write(*,*) 'Cambridge model2: '
          write(*,*) sio2%zbr,sio2%zmb,sio2%zmf,sio2%sig00,    &
                        sio2%ocr,sio2%rs0,sio2%zeta0          


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
!!$acc enter data create(parts,soil)
!---------------------Real particles (water)---------------------
!parts%dim  = dim
!parts%maxn = maxn
!parts%max_interaction = max_interaction
dim = parts%dim
maxn = parts%maxn
max_interaction = parts%max_interaction

call allocate_particles(parts)

allocate(parts%zone(maxn)); parts%zone = 0
!$acc enter data create(parts%zone)

allocate(parts%w(max_interaction));             parts%w   = 0.d0
allocate(parts%dwdx(dim,max_interaction));      parts%dwdx= 0.d0
allocate(parts%pair_i(max_interaction));        parts%pair_i = 0
allocate(parts%pair_j(max_interaction));        parts%pair_j = 0
allocate(parts%countiac(maxn));                 parts%countiac=0
allocate(parts%grid(parts%maxngx,parts%maxngy,parts%maxngz)); parts%grid = 0
allocate(parts%xgcell(dim,maxn)); parts%xgcell = 0
allocate(parts%celldata(maxn)); parts%celldata = 0
allocate(parts%norm(2,10000))

!$acc enter data create(parts%grid(parts%maxngx,parts%maxngy,parts%maxngz))
!$acc enter data create(parts%celldata(maxn))

! Sherpard Filter
allocate(parts%wi); call parts%wi%set(rank=0,parts=parts,name = ' ')
call parts%wi%max_alloc
parts%wi%r = 0.d0
!$acc enter data create(parts%wi,parts%wi%r)

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
call parts%tab%max_alloc
parts%tab%x%r = 0.d0; parts%tab%y%r = 0.d0; parts%tab%xy%r = 0.d0
parts%tab%z%r = 0.d0

!$acc enter data create(parts%tab,parts%tab%x, parts%tab%y, parts%tab%xy)
!$acc enter data create(parts%tab%x%r(maxn),parts%tab%y%r(maxn),parts%tab%xy%r(maxn))

if(trim(parts%imaterial)=='water')then
!   allocate(parts%sxx,parts%sxy,parts%syy)
!    allocate(parts%str);allocate(parts%str%x,parts%str%xy,parts%str%y)
!    parts%str%x%r => parts%tab%x%r
!    parts%str%y%r => parts%tab%y%r
!    parts%str%xy%r => parts%tab%xy%r
   parts%str => parts%tab
elseif(trim(parts%imaterial)=='soil')then

!Changed according to Wu Hao:
!-----------------------------
!parts%wxy => parts%tab%xy

allocate(parts%wxy); call parts%wxy%set(rank=0,parts=parts,name = 'wxy ')
call parts%wxy%max_alloc
parts%wxy%r = 0.d0
!------------------------------------------------------------------        

!$acc enter data create(parts%wxy,parts%wxy%r(maxn))

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
!$acc enter data create(soil%zone)

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

!$acc enter data create(soil%grid(soil%maxngx,soil%maxngy,soil%maxngz))
!$acc enter data create(soil%celldata(maxn))

allocate(soil%wi); call soil%wi%set(rank=0,parts=soil,name = ' ')
call soil%wi%max_alloc
soil%wi%r = 0.d0
!$acc enter data create(soil%wi,soil%wi%r)

! Fields variables
call allocate_particles_fields(soil)

!soil%tab => parts%tab
allocate(soil%tab); call soil%tab%set(rank=2, parts=soil, name='tab')
call soil%tab%max_alloc
soil%tab%x%r = 0.d0; soil%tab%y%r = 0.d0; soil%tab%xy%r = 0.d0
soil%tab%z%r = 0.d0

!$acc enter data create(soil%tab,soil%tab%x,soil%tab%y,soil%tab%xy,soil%tab%z)
!$acc enter data create(soil%tab%x%r(maxn),soil%tab%y%r(maxn),soil%tab%xy%r(maxn),soil%tab%z%r(maxn))

!Changed according to Wu Hao:
!-----------------------------
!soil%wxy => soil%tab%xy

allocate(soil%wxy); call soil%wxy%set(rank=0,parts=soil,name = 'wxy ')
call soil%wxy%max_alloc
soil%wxy%r = 0.d0
!------------------------------------------------------------------
!$acc enter data create(soil%wxy,soil%wxy%r(maxn))

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
!$acc enter data create(soil%dvof,soil%vof_min)
!$acc enter data create(soil%dvof%r(:),soil%vof_min%r(:))


   allocate(soil%psi); call soil%psi%set(rank=0,parts=soil,name='psi')
   call soil%psi%max_alloc; soil%psi%r = 0.d0  
   allocate(soil%rho0); call soil%rho0%set(rank=0,parts=soil,name='rho0')
   call soil%rho0%max_alloc; soil%rho0%r = 0.d0      
!$acc enter data create(soil%psi,soil%rho0)
!$acc enter data create(soil%psi%r(:),soil%rho0%r(:))

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
!$acc enter data create(parts%itype(:))

if(.not.associated(parts%x))     allocate(parts%x(dim,maxn))
!$acc enter data create(parts%x(dim,maxn))

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
!$acc enter data create(parts%mass)
!$acc enter data create(parts%mass%r(:))
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
!$acc enter data create(parts%vx,parts%vx%x,parts%vx%y)
!$acc enter data create(parts%vx%x%r(:),parts%vx%y%r(:))


allocate(parts%rho); call parts%rho%set(rank=0,parts=parts,name='rho')
call parts%rho%max_alloc
parts%rho%r = 0.d0
!$acc enter data create(parts%rho, parts%rho%r(:))

allocate(parts%p); call parts%p%set(rank=0,parts=parts,name='p')
call parts%p%max_alloc
parts%p%r = 0.d0
!$acc enter data create(parts%p,parts%p%r(:))

allocate(parts%c); call parts%c%set(rank=0,parts=parts,name='c')
call parts%c%max_alloc
parts%c%r = 0.d0

! Old value
allocate(parts%v_min); call parts%v_min%set(rank=1,parts=parts,name='v_min')
call parts%v_min%max_alloc
parts%v_min%x%r = 0.d0; parts%v_min%y%r = 0.d0
!$acc enter data create(parts%v_min,parts%v_min%x,parts%v_min%y)
!$acc enter data create(parts%v_min%x%r(:),parts%v_min%y%r(:))

allocate(parts%rho_min); call parts%rho_min%set(rank=0,parts=parts,name='rho_min')
call parts%rho_min%max_alloc
parts%rho_min%r = 0.d0
!$acc enter data create(parts%rho_min)
!$acc enter data create(parts%rho_min%r(:))


! Accelerations
allocate(parts%drho);call parts%drho%set(rank=0,parts=parts,name='drho')
call parts%drho%max_alloc; parts%drho%r = 0.d0
!$acc enter data create(parts%drho,parts%drho%r(:))

allocate(parts%dvx); call parts%dvx%set(rank=1,parts=parts,name='dvx')
call parts%dvx%max_alloc
parts%dvx%x%r = 0.d0; parts%dvx%y%r = 0.d0  
!$acc enter data create(parts%dvx,parts%dvx%x,parts%dvx%y)
!$acc enter data create(parts%dvx%x%r(:),parts%dvx%y%r(:))

allocate(parts%av); call parts%av%set(rank=1,parts=parts,name='av')
call parts%av%max_alloc
parts%av%x%r  = 0.d0; parts%av%y%r  = 0.d0

!$acc enter data create(parts%av,parts%av%x,parts%av%y)
!$acc enter data create(parts%av%x%r(:),parts%av%y%r(:))

! For particle shifting
allocate(parts%r0);call parts%r0%set(rank=0,parts=parts,name='r0')
call parts%r0%max_alloc; parts%r0%r = 0.d0
allocate(parts%mt);call parts%mt%set(rank=0,parts=parts,name='mt')
call parts%mt%max_alloc; parts%mt%r = 0.d0
allocate(parts%ps); call parts%ps%set(rank=1,parts=parts,name='ps')
call parts%ps%max_alloc
parts%ps%x%r  = 0.d0; parts%ps%y%r  = 0.d0

!$acc enter data create(parts%ps,parts%ps%x,parts%ps%y)
!$acc enter data create(parts%ps%x%r(:),parts%ps%y%r(:))

! Volume Fraction
allocate(parts%vof); call parts%vof%set(rank=0,parts=parts,name='vof')
call parts%vof%max_alloc; parts%vof%r = 0.d0
!$acc enter data create(parts%vof)
!$acc enter data create(parts%vof%r(:))

if(trim(parts%imaterial)=='water')then
   ! phi_f = 1- phi_s
   allocate(parts%vof2); call parts%vof2%set(rank=0,parts=parts,name='vof2')
   call parts%vof2%max_alloc; parts%vof2%r = 0.d0
   allocate(parts%dvof); call parts%dvof%set(rank=0,parts=parts,name='dvof')
   call parts%dvof%max_alloc; parts%dvof%r = 0.d0   
   allocate(parts%vof_min)
   call parts%vof_min%set(rank=0,parts=parts,name='vof_min')
   call parts%vof_min%max_alloc; parts%vof_min%r = 0.d0
!$acc enter data create(parts%vof2,parts%vof2%r)   
!$acc enter data create(parts%dvof,parts%vof_min)
!$acc enter data create(parts%dvof%r(:),parts%vof_min%r(:))
endif

if(trim(parts%imaterial)=='soil')then
    allocate(parts%str)
!    allocate(parts%str%x,parts%str%xy,parts%str%y)
!   allocate(parts%sxx,parts%sxy,parts%syy)
!   allocate(parts%str%x%r(maxn), parts%str%y%r(maxn), parts%str%xy%r(maxn))
   call parts%str%set(rank=2,parts=parts,name='str'); call parts%str%max_alloc
   parts%str%x%r = 0.d0; parts%str%xy%r = 0.d0; parts%str%y%r = 0.d0   !!! %r
   parts%str%z%r = 0.d0

!$acc enter data create(parts%str)
!$acc enter data create(parts%str%x,parts%str%xy,parts%str%y,parts%str%z)
!$acc enter data create(parts%str%x%r(:),parts%str%xy%r,parts%str%y%r,       &
!$acc parts%str%z%r)
   
   
    !artificial stress
    allocate(parts%strp); 
    call parts%strp%set(rank=2,parts=parts,name='strp')
    call parts%strp%max_alloc
    parts%strp%x%r = 0.d0; parts%strp%y%r = 0.d0; parts%strp%xy%r = 0.d0
    parts%strp%z%r = 0.d0
    
    allocate(parts%ast); call parts%ast%set(rank=2,parts=parts,name='ast')
    call parts%ast%max_alloc
    parts%ast%x%r = 0.d0; parts%ast%y%r = 0.d0; parts%ast%xy%r = 0.d0
    parts%ast%z%r = 0.d0
    
    allocate(parts%astp); call parts%astp%set(rank=2,parts=parts,name='astp')
    call parts%astp%max_alloc
    parts%astp%x%r = 0.d0; parts%astp%y%r = 0.d0; parts%astp%xy%r = 0.d0
    parts%astp%z%r = 0.d0
    
    allocate(parts%thetai); call parts%thetai%set(rank=0,parts=parts,name='thetai')
    call parts%thetai%max_alloc
    parts%thetai%r = 0.d0
   

   allocate(parts%str_min)
   call parts%str_min%set(rank=2,parts=parts,name='str_min')
   call parts%str_min%max_alloc
   parts%str_min%x%r = 0.d0; parts%str_min%xy%r = 0.d0; parts%str_min%y%r = 0.d0 !!!%r
   parts%str_min%z%r = 0.d0
 !$acc enter data create(parts%str_min)  
 !$acc enter data create(parts%str_min%x,parts%str_min%y,parts%str_min%xy)
 !$acc enter data create(parts%str_min%x%r(:),parts%str_min%y%r(:),parts%str_min%xy%r(:))

   allocate(parts%p_min); call parts%p_min%set(rank=0,parts=parts,name='p_min')
   call parts%p_min%max_alloc
   parts%p_min%r = 0.d0
!$acc enter data create(parts%p_min,parts%p_min%r(:))

   allocate(parts%vcc); call parts%vcc%set(rank=0,parts=parts,name='vcc')
   call parts%vcc%max_alloc 
   parts%vcc%r = 0.d0
!$acc enter data create(parts%vcc,parts%vcc%r(:))

   !allocate(parts%szz(maxn))    !!! Intel Fortran Compiler is shit!!!

   allocate(parts%dstr)
   call parts%dstr%set(rank=2,parts=parts,name='dstr'); call parts%dstr%max_alloc
   parts%dstr%x%r = 0.d0; parts%dstr%xy%r = 0.d0; parts%dstr%y%r = 0.d0  !!!%r
   parts%dstr%z%r = 0.d0

!$acc enter data create(parts%dstr)   
!$acc enter data create(parts%dstr%x, parts%dstr%y, parts%dstr%xy, parts%dstr%z)
!$acc enter data create(parts%dstr%x%r(:),parts%dstr%y%r(:),   &
!$acc parts%dstr%xy%r(:),parts%dstr%z%r(:))

   allocate(parts%dp); call parts%dp%set(rank=0,parts=parts,name='dp')
   call parts%dp%max_alloc; parts%dp%r = 0.d0 !!!%r
!$acc enter data create(parts%dp,parts%dp%r(:))

! For return mapping algorithm
   allocate(parts%dstr2)
   call parts%dstr2%set(rank=2,parts=parts,name='dstr2'); call parts%dstr2%max_alloc
   !,parts%dstr2%x,parts%dstr2%y,parts%dstr2%xy)
   !allocate(parts%dstr2%x%r(maxn),parts%dstr2%y%r(maxn),parts%dstr2%xy%r(maxn))
   parts%dstr2%x%r = 0.d0; parts%dstr2%xy%r = 0.d0; parts%dstr2%y%r = 0.d0 !!!%r
   parts%dstr2%z%r = 0.d0
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
!$acc enter data create(parts%epsilon_p,parts%epsilon_p%r(:))   
  
   allocate(parts%epsilon_p_min)
   call parts%epsilon_p_min%set(rank=0,parts=parts,name='ep_min')
   !allocate(parts%epsilon_p_min%r(maxn))
   call parts%epsilon_p_min%max_alloc
   parts%epsilon_p_min%r = 0.d0  !!!%r
!$acc enter data create(parts%epsilon_p_min,parts%epsilon_p_min%r(:))   

   allocate(parts%depsilon_p)
   call parts%depsilon_p%set(rank=0,parts=parts,name='dep')
   !allocate(parts%depsilon_p%r(maxn))
   call parts%depsilon_p%max_alloc
   parts%depsilon_p%r = 0.d0  !!!%r
!$acc enter data create(parts%depsilon_p,parts%depsilon_p%r(:))   

   allocate(parts%spp)   !!! VOF
   call parts%spp%set(rank=0,parts=parts,name='spp'); call parts%spp%max_alloc
   parts%spp%r = 0.d0
!$acc enter data create(parts%spp,parts%spp%r(:))   

   allocate(parts%spp0)   !!! VOF
   call parts%spp0%set(rank=0,parts=parts,name='spp0'); call parts%spp0%max_alloc
   parts%spp0%r = 0.d0
!$acc enter data create(parts%spp0,parts%spp0%r(:))   


endif

!allocate(parts%drhodx); allocate(parts%drhodx%r(maxn))
!parts%drhodx = 0.d0

if(trim(parts%imaterial)=='soil'.and.parts%soil_model==1)then
   allocate(parts%zr,parts%zrs,parts%zzeta,parts%zharbeta,parts%tnone,parts%vp0,parts%evi,parts%vp)
   call parts%zr%set(rank=0,parts=parts,name='zr '); call parts%zr%max_alloc
   parts%zr%r = 0.d0
   call parts%zrs%set(rank=0,parts=parts,name='zrs '); call parts%zrs%max_alloc
   parts%zrs%r = 0.d0
   call parts%zzeta%set(rank=0,parts=parts,name='zzeta '); call parts%zzeta%max_alloc
   parts%zzeta%r = 0.d0
   call parts%zharbeta%set(rank=2,parts=parts,name='zharbeta '); call parts%zharbeta%max_alloc
   parts%zharbeta%x%r = 0.d0;    parts%zharbeta%y%r =0.d0
   parts%zharbeta%xy%r =  0.d0;  parts%zharbeta%z%r =0.d0
   call parts%tnone%set(rank=0,parts=parts,name='tnone '); call parts%tnone%max_alloc
   parts%tnone%r = 0.d0   
   call parts%vp0%set(rank=0,parts=parts,name='vp0 '); call parts%vp0%max_alloc
   parts%vp0%r = 0.d0
   call parts%vp%set(rank=0,parts=parts,name='vp '); call parts%vp%max_alloc
   parts%vp%r = 0.d0   
   call parts%evi%set(rank=0,parts=parts,name='evi '); call parts%evi%max_alloc
   parts%evi%r = 0.d0


   allocate(parts%DRS,parts%DVP)
   allocate(parts%dharbeta1,parts%dharbeta2,parts%dharbeta3,parts%dharbeta4)

   call parts%DRS%set(rank=0,parts=parts,name='DRS ')
   call parts%DRS%max_alloc; parts%DRS%r = 0.d0

   call parts%dharbeta1%set(rank=0,parts=parts,name='dharbeta1 ')
   call parts%dharbeta1%max_alloc; parts%dharbeta1%r = 0.d0

   call parts%dharbeta2%set(rank=0,parts=parts,name='dharbeta2 ')
   call parts%dharbeta2%max_alloc; parts%dharbeta2%r = 0.d0

   call parts%dharbeta3%set(rank=0,parts=parts,name='dharbeta3 ')
   call parts%dharbeta3%max_alloc; parts%dharbeta3%r = 0.d0

   call parts%dharbeta4%set(rank=0,parts=parts,name='dharbeta4 ')
   call parts%dharbeta4%max_alloc; parts%dharbeta4%r = 0.d0   

   call parts%DVP%set(rank=0,parts=parts,name='DVP '); call parts%DVP%max_alloc
   parts%DVP%r = 0.d0   


   allocate(parts%ZRS_min,parts%zharbeta_min,parts%VP_min)

   call parts%ZRS_min%set(rank=0,parts=parts,name='ZRS_min ')
   call parts%ZRS_min%max_alloc; parts%ZRS_min%r = 0.d0

   call parts%zharbeta_min%set(rank=2,parts=parts,name='zharbeta_min ')
   call parts%zharbeta_min%max_alloc
   parts%zharbeta_min%x%r = 0.d0;    parts%zharbeta_min%y%r =0.d0
   parts%zharbeta_min%xy%r =  0.d0;  parts%zharbeta_min%z%r =0.d0   

   call parts%vp_min%set(rank=0,parts=parts,name='vp_min '); call parts%vp_min%max_alloc
   parts%vp_min%r = 0.d0   


endif

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
                                pl => soil
                                call first_half
                                if(pl%soil_model==0)call pl%scaling_back
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

                                    pl => soil
                                    call first_step
        else
   
           pl => parts
           call second_half          

                                   pl => soil
                                   call second_half
                                   if(pl%soil_model==0)call pl%scaling_back
        endif 

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

!Cambridge model         
         if(trim(pl%imaterial)=='soil'.and.pl%soil_model==1)then
            pl%ZRS_min = pl%ZRS
            pl%ZRS = pl%ZRS+(pl%dt/2.)*pl%DRS

            pl%zharbeta_min%x = pl%zharbeta%x
            pl%zharbeta_min%y = pl%zharbeta%y
            pl%zharbeta_min%z = pl%zharbeta%z
            pl%zharbeta_min%xy = pl%zharbeta%xy

            pl%zharbeta%x = pl%zharbeta%x + pl%dt/2.*pl%dharbeta1
            pl%zharbeta%y = pl%zharbeta%y + pl%dt/2.*pl%dharbeta2
            pl%zharbeta%z = pl%zharbeta%z + pl%dt/2.*pl%dharbeta3
            pl%zharbeta%xy = pl%zharbeta%xy + pl%dt/2.*pl%dharbeta4

            pl%VP_min = pl%VP
            pl%VP = pl%VP + pl%dt/2.*pl%DVP

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
         pl%epsilon_p_min = pl%epsilon_p

         pl%str%x = pl%str%x + (pl%dt/2.)*pl%dstr%x
         pl%str%xy = pl%str%xy + (pl%dt/2.)*pl%dstr%xy
         pl%str%y = pl%str%y + (pl%dt/2.)*pl%dstr%y
         !pl%str = pl%str + (dt/2.) * pl%dstr
         pl%p = pl%p + (pl%dt/2.)*pl%dp      !!!simultaneous pressure
         pl%epsilon_p = pl%epsilon_p + (pl%dt/2.)*pl%depsilon_p

                  endif
         endif
      !write(*,*) pl%ntotal + pl%nvirt    

!!!!! Here we need parallization!      
!$omp parallel do private(d,vxi,dvxi,v_mini)
      do i = 1, pl%ntotal +pl%nvirt    ! originally only pl%ntotal       
         if(pl%itype(i)<0)cycle
!         if(pl%x(2,i)>0.3.and.(pl%zone(i)==6.or.pl%zone(i)==100))cycle           
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
!$omp end parallel do

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

         if(trim(pl%imaterial)=='soil'.and.pl%soil_model==1)then
            pl%ZRS = pl%ZRS+(pl%dt/2.)*pl%DRS

            pl%zharbeta%x = pl%zharbeta%x + pl%dt/2.*pl%dharbeta1
            pl%zharbeta%y = pl%zharbeta%y + pl%dt/2.*pl%dharbeta2
            pl%zharbeta%z = pl%zharbeta%z + pl%dt/2.*pl%dharbeta3
            pl%zharbeta%xy = pl%zharbeta%xy + pl%dt/2.*pl%dharbeta4

            pl%VP = pl%VP + pl%dt/2.*pl%DVP            
         endif


         if(trim(pl%imaterial)=='soil')then
                  if(pl%stress_integration==1)then
         pl%str%x = pl%str%x + pl%dt*pl%dstr%x
         pl%str%xy = pl%str%xy + pl%dt*pl%dstr%xy
         pl%str%y = pl%str%y + pl%dt*pl%dstr%y
         pl%p   = pl%p   + pl%dt*pl%dp       !!! simultaneous pressure
         pl%epsilon_p = pl%epsilon_p + pl%dt*pl%depsilon_p
                  elseif(pl%stress_integration==2)then
         pl%str%x = pl%str%x + (pl%dt/2.)*pl%dstr%x
         pl%str%xy = pl%str%xy + (pl%dt/2.)*pl%dstr%xy
         pl%str%y = pl%str%y + (pl%dt/2.)*pl%dstr%y
         pl%p   = pl%p + (pl%dt/2.)*pl%dp    !!! simultaneous pressure
         pl%epsilon_p = pl%epsilon_p + (pl%dt/2.)*pl%depsilon_p
            !write(*,*) 'first step, soil%str'
            !write(*,*) pl%str%x%r(1:10)                             
            !write(*,*) pl%str%xy%r(1:10)                             
            !write(*,*) pl%str%y%r(1:10)                             
            !write(*,*) pl%p%r(1:10)                             
                  endif
         endif
        
!!!!! Here we need parallization!      
!$omp parallel do private(d,vxi,dvxi,avi,psi)
      do i=1,pl%ntotal +pl%nvirt     ! origionally pl%ntotal
         if(pl%itype(i)<0)cycle
!         if(pl%x(2,i)>0.3.and.(pl%zone(i)==6.or.pl%zone(i)==100))cycle           
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

            !if(i<=10)write(*,*) i, vxi(1)%p,vxi(2)%p
            !if(i<=10)write(*,*) i, pl%x(1,i),pl%x(2,i)
      enddo 
      !$omp end parallel do     
      return
      end subroutine

      subroutine second_half
! ------------------------------------------------------------------
      implicit none
      type(p2r) vxi(3), dvxi(3), v_mini(3),avi(3), psi(3)


!$omp parallel do private(d,vxi,dvxi,v_mini,avi,psi)
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

         if(trim(pl%imaterial)=='soil'.and.pl%soil_model==1)then
            pl%ZRS%r(i) = pl%ZRS%r(i)+(pl%dt)*pl%DRS%r(i)

            pl%zharbeta%x%r(i) = pl%zharbeta%x%r(i) + pl%dt*pl%dharbeta1%r(i)
            pl%zharbeta%y%r(i) = pl%zharbeta%y%r(i) + pl%dt*pl%dharbeta2%r(i)
            pl%zharbeta%z%r(i) = pl%zharbeta%z%r(i) + pl%dt*pl%dharbeta3%r(i)
            pl%zharbeta%xy%r(i) = pl%zharbeta%xy%r(i) + pl%dt*pl%dharbeta4%r(i)

            pl%VP%r(i) = pl%VP%r(i) + pl%dt*pl%DVP%r(i)                        
         endif

         if(trim(pl%imaterial)=='soil')then 
                  if(pl%stress_integration==1)then 
         pl%str%x%r(i) = pl%str%x%r(i) + pl%dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + pl%dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + pl%dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p%r(i)   + pl%dt*pl%dp%r(i)       !!! simultaneous pressure
         pl%epsilon_p%r(i) = pl%epsilon_p%r(i) + pl%dt*pl%depsilon_p%r(i)
                  elseif(pl%stress_integration==2)then
         pl%str%x%r(i) = pl%str_min%x%r(i) + pl%dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str_min%xy%r(i) + pl%dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str_min%y%r(i) + pl%dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p_min%r(i) + pl%dt*pl%dp%r(i)     !!! simultaneous pressure
         pl%epsilon_p%r(i) = pl%epsilon_p_min%r(i) + pl%dt*pl%depsilon_p%r(i)
                  endif
         endif
        
         if(pl%itype(i)<0)cycle
!         if(pl%x(2,i)>0.3.and.(pl%zone(i)==6.or.pl%zone(i)==100))cycle           

         vxi = pl%vx%cmpt(i); dvxi = pl%dvx%cmpt(i) ; v_mini = pl%v_min%cmpt(i)
         avi = pl%av%cmpt(i); psi = pl%ps%cmpt(i)
         do d = 1, pl%dim                   
            !pl%vx(d, i) = pl%v_min(d, i) + dt * pl%dvx(d, i)   &
            !            + pl%av(d, i)
            !pl%x(d, i) = pl%x(d, i) + dt * pl%vx(d, i)
            vxi(d)%p = v_mini(d)%p + pl%dt * dvxi(d)%p + avi(d)%p
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

      end subroutine


!----------------------------------------------------------------------
      subroutine time_integration_GPU
!----------------------------------------------------------------------
      implicit none     
!
      integer :: i, j, k, d, ntotal, it
     
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
           call first_half_GPU(parts)
           call first_half_GPU(soil)
           call scaling_back(soil)
        endif

!---  Definition of variables out of the function vector:    
         call single_step_GPU
        !call single_step_water_GPU
        !call single_step_soil_GPU

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
           call first_step_GPU(parts)
           call first_step_GPU(soil)
        else
           call second_half_GPU(parts)
           call second_half_GPU(soil)
           call scaling_back(soil)
        endif 

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

      subroutine first_half_GPU(pl)
! ---------------------------------------------------------------------
implicit none
type(particles) pl
integer i, ntotal
real(dp) dt

ntotal = pl%ntotal + pl%nvirt
dt = pl%dt

         if (.not.pl%summation_density) then
!$omp parallel do                 
!$acc kernels loop independent vector
            do i = 1, ntotal 
               pl%rho_min%r(i) = pl%rho%r(i) 
               pl%rho%r(i) = pl%rho%r(i)+(dt/2.)*pl%drho%r(i)
            enddo
!$omp end parallel do
!$omp barrier
         endif

         if(trim(pl%imaterial)=='water'.and.pl%volume_fraction)then
!$omp parallel do                
!$acc kernels loop independent vector
             do i = 1, ntotal
                pl%vof_min%r(i) = pl%vof%r(i)
                pl%vof%r(i)  = pl%vof%r(i) + dt/2.*pl%dvof%r(i)
             enddo
!$omp end parallel do
!$omp barrier
         endif

         if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
!$omp parallel do                 
!$acc kernels loop independent vector
             do i = 1, ntotal 
                pl%vof_min%r(i) = pl%vof%r(i)
                pl%vof%r(i)  = pl%vof%r(i) + dt/2.*pl%dvof%r(i)
             enddo   
!$omp end parallel do
!$omp barrier

         endif


         if(trim(pl%imaterial)=='soil')then
!$omp parallel do                 
!$acc kernels loop independent vector
         do i = 1, ntotal                   
         pl%str_min%x%r(i) = pl%str%x%r(i)
         pl%str_min%xy%r(i) = pl%str%xy%r(i)
         pl%str_min%y%r(i) = pl%str%y%r(i)
         pl%p_min%r(i)   = pl%p%r(i)
         pl%epsilon_p_min%r(i) = pl%epsilon_p%r(i)

         pl%str%x%r(i) = pl%str%x%r(i) + (dt/2.)*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + (dt/2.)*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + (dt/2.)*pl%dstr%y%r(i)
         !pl%str = pl%str + (dt/2.) * pl%dstr
         pl%p%r(i) = pl%p%r(i) + (dt/2.)*pl%dp%r(i)      
         pl%epsilon_p%r(i) = pl%epsilon_p%r(i) + (dt/2.)*pl%depsilon_p%r(i)
         enddo
!$omp end parallel do
!$omp barrier
         endif

!$omp parallel do                 
!$acc kernels loop independent vector present(pl%itype)
      do i = 1, ntotal     
         if(pl%itype(i)<0)cycle
!         if(pl%x(2,i)>0.3.and.(pl%zone(i)==6.or.pl%zone(i)==100))cycle           

            pl%v_min%x%r(i) =  pl%vx%x%r(i)         
            pl%v_min%y%r(i) =  pl%vx%y%r(i)         
            pl%vx%x%r(i) = pl%vx%x%r(i) + dt/2.d0*pl%dvx%x%r(i)         
            pl%vx%y%r(i) = pl%vx%y%r(i) + dt/2.d0*pl%dvx%y%r(i)         

      enddo
!$omp end parallel do
!$omp barrier
      return
      end subroutine


      subroutine first_step_GPU(pl)
! -------------------------------------------------------------------
      implicit none
      type(particles) pl
      integer i,ntotal
      real(dp) dt

      ntotal = pl%ntotal+pl%nvirt
      dt = pl%dt

         if (.not.pl%summation_density ) then
!$omp parallel do                 
!$acc kernels loop independent vector  
            do i = 1, ntotal
               pl%rho%r(i) = pl%rho%r(i) + (dt/2.)* pl%drho%r(i)
            enddo
!$omp end parallel do
!$omp barrier
         endif

         if(trim(pl%imaterial)=='water'.and.pl%volume_fraction)then
!$omp parallel do                 
!$acc kernels loop independent 
            do i = 1, ntotal
               pl%vof%r(i) = pl%vof%r(i)+(dt/2.)*pl%dvof%r(i)
            enddo
!$omp end parallel do
!$omp barrier
         endif

         if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
!$omp parallel do                 
!$acc kernels loop independent vector
            do i = 1, ntotal
               pl%vof%r(i) = pl%vof%r(i)+(dt/2.)*pl%dvof%r(i)
            enddo   
!$omp end parallel do
!$omp barrier
         endif


         if(trim(pl%imaterial)=='soil')then
                  if(pl%stress_integration==1)then
         pl%str%x = pl%str%x + pl%dt*pl%dstr%x
         pl%str%xy = pl%str%xy + pl%dt*pl%dstr%xy
         pl%str%y = pl%str%y + pl%dt*pl%dstr%y
         pl%p   = pl%p   + pl%dt*pl%dp       !!! simultaneous pressure
         pl%epsilon_p = pl%epsilon_p + pl%dt*pl%depsilon_p
                  elseif(pl%stress_integration==2)then

!$omp parallel do                 
!$acc kernels loop independent vector 
         do i = 1, ntotal
         pl%str%x%r(i) = pl%str%x%r(i) + (dt/2.)*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + (dt/2.)*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + (dt/2.)*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p%r(i) + (dt/2.)*pl%dp%r(i)    !!! simultaneous pressure
         pl%epsilon_p%r(i) = pl%epsilon_p%r(i) + (dt/2.)*pl%depsilon_p%r(i)
         enddo
!$omp end parallel do
!$omp barrier
                  endif
         endif

!$omp parallel do                 
!$acc kernels loop independent vector present(pl%itype,pl%x) 
      do i=1,ntotal    
         if(pl%itype(i)<0)cycle
!         if(pl%x(2,i)>0.3.and.(pl%zone(i)==6.or.pl%zone(i)==100))cycle           

          pl%vx%x%r(i) = pl%vx%x%r(i) + (dt/2.) * pl%dvx%x%r(i)   &
                       + pl%av%x%r(i)
          pl%vx%y%r(i) = pl%vx%y%r(i) + (dt/2.) * pl%dvx%y%r(i)   &
                       + pl%av%y%r(i)            
          pl%x(1, i) = pl%x(1, i) + dt * (pl%vx%x%r(i)+pl%ps%x%r(i)) 
          pl%x(2, i) = pl%x(2, i) + dt * (pl%vx%y%r(i)+pl%ps%y%r(i))
      enddo 
!$omp end parallel do
!$omp barrier

      return
      end subroutine

      subroutine second_half_GPU(pl)
! ------------------------------------------------------------------
      implicit none
      type(particles) pl
      integer i,ntotal
      real(dp) dt

      ntotal = pl%ntotal+pl%nvirt
      dt = pl%dt

            
         if (.not.pl%summation_density ) then 
                 
!$omp parallel do                 
!$acc kernels loop independent vector
            do i=1,ntotal              
               pl%rho%r(i) = pl%rho_min%r(i) + dt*pl%drho%r(i)
            enddo   
!$omp end parallel do
!$omp barrier
         endif
   
         if(trim(pl%imaterial)=='water'.and.pl%volume_fraction)then
!$omp parallel do                 
!$acc kernels loop independent vector
            do i=1,ntotal             
               pl%vof%r(i) = pl%vof_min%r(i)+dt*pl%dvof%r(i) 
            enddo   
!$omp end parallel do
!$omp barrier
         endif
                 
         if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
!$omp parallel do                 
!$acc kernels loop independent vector
            do i=1,ntotal
               pl%vof%r(i) = pl%vof_min%r(i)+dt*pl%dvof%r(i) 
            enddo   
!$omp end parallel do
!$omp barrier
         endif


         if(trim(pl%imaterial)=='soil')then 
                  if(pl%stress_integration==1)then 
      do i=1,ntotal              
         pl%str%x%r(i) = pl%str%x%r(i) + pl%dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str%xy%r(i) + pl%dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str%y%r(i) + pl%dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p%r(i)   + pl%dt*pl%dp%r(i)       !!! simultaneous pressure
         pl%epsilon_p%r(i) = pl%epsilon_p%r(i) + pl%dt*pl%depsilon_p%r(i)
      enddo   
                  elseif(pl%stress_integration==2)then
!$omp parallel do                 
!$acc kernels loop independent vector
      do i=1,ntotal              
         pl%str%x%r(i) = pl%str_min%x%r(i) + dt*pl%dstr%x%r(i)
         pl%str%xy%r(i) = pl%str_min%xy%r(i) + dt*pl%dstr%xy%r(i)
         pl%str%y%r(i) = pl%str_min%y%r(i) + dt*pl%dstr%y%r(i)
         pl%p%r(i)   = pl%p_min%r(i) + dt*pl%dp%r(i)     !!! simultaneous pressure
         pl%epsilon_p%r(i) = pl%epsilon_p_min%r(i) + dt*pl%depsilon_p%r(i)
      enddo    
!$omp end parallel do
!$omp barrier
                  endif
         endif

!$omp parallel do                 
!$acc kernels loop independent vector present(pl%itype,pl%x)
      do i=1,ntotal              
         if(pl%itype(i)<0)cycle
!         if(pl%x(2,i)>0.3.and.(pl%zone(i)==6.or.pl%zone(i)==100))cycle           

          pl%vx%x%r(i) = pl%v_min%x%r(i) + dt * pl%dvx%x%r(i)   &
                       + pl%av%x%r(i)
          pl%vx%y%r(i) = pl%v_min%y%r(i) + dt * pl%dvx%y%r(i)   &
                       + pl%av%y%r(i)
            
          pl%x(1, i) = pl%x(1, i) + dt * (pl%vx%x%r(i)+pl%ps%x%r(i))                  
          pl%x(2, i) = pl%x(2, i) + dt * (pl%vx%y%r(i)+pl%ps%y%r(i))                              
      enddo
!$omp end parallel do
!$omp barrier
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
   pl%depsilon_p = 0.d0; pl%dstr%z =0.d0
   if(pl%critical_state==1) pl%dvof = 0.d0
endif

!if(trim(pl%imaterial)=='water') call inlet_boundary2
!call over_flow(pl)  ! Outlet boundary !Doesn't work!
!call symmetry_boundary(pl)   !!! subjetting
!if(iphase==2) call slip_boundary(pl)   !!! submerged column collapse, take_virtual 
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
pl%tab%x = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%x,'x')-pl%df4(pl%vx%y,'y'))
pl%tab%xy = pl%df4(pl%vx%x,'y')+pl%df4(pl%vx%y,'x')
pl%tab%y = 2.d0/3.d0*(2.d0*pl%df4(pl%vx%y,'y')-pl%df4(pl%vx%x,'x'))

! Calculate dilatancy angle 'psi' and volume fraction rate of change 'dvof' with critical
! state theory
if(trim(pl%imaterial)=='soil'.and.pl%critical_state==1)then
   call soil_critical_state(pl)
endif

if(trim(pl%imaterial)=='soil') pl%vcc = pl%div_omp(pl%vx)
!call velocity_divergence(pl)

call pressure(pl)

if(trim(pl%imaterial)=='water')then
   call newtonian_fluid(pl)
endif   

if(trim(pl%imaterial)=='soil')then

!Modified by Wu Hao:--------------------
    if(pl%soil_model==0)call stress_rate(pl)
    if(pl%soil_model==1) call KISUPERCAM(pl) 
!-----------------------------------------------------------------   

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
   !where (pl%rho%r.gt.0.0) pl%dvx%x%r = pl%dvx%x%r/pl%rho%r
   !where (pl%rho%r.gt.0.0) pl%dvx%y%r = pl%dvx%y%r/pl%rho%r
   pl%dvx%x = pl%dvx%x/pl%rho
   pl%dvx%y = pl%dvx%y/pl%rho

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


!---  Artificial viscosity:
if (pl%visc_artificial) call pl%art_visc_omp

if(trim(pl%imaterial)=='soil'.and.pl%soil_artificial_stress)then
        !call art_stress(pl)
        ! call artificial_stress(pl)             ! Added by Lu Shiyang
    call pl%delta_sph_omp(pl%p,pl%dp)             !! Jump
    call pl%delta_sph_omp(pl%str%x,pl%dstr%x)     !! Jump
    call pl%delta_sph_omp(pl%str%xy,pl%dstr%xy)   !! No jump
    call pl%delta_sph_omp(pl%str%y,pl%dstr%y)     !! Jump
endif        
if(trim(pl%imaterial)=='water'.and.pl%water_artificial_volume)then
        !call art_volume_fraction_water2(pl)
         call pl%delta_sph_omp(pl%vof,pl%dvof)   !!! Important!
endif        

!--- Damping
!       if(trim(pl%imaterial)=='soil') call damping_stress(pl)
    
!---  External forces:

      !if (ex_force) call ext_force(pl)
!      if (ex_force)then
!          if(self_gravity) call gravity_force(pl)
!          call repulsive_force(pl)
          call pl%repulsive_force_omp           !!!!!!!#### Not important!
!      endif

pl%dvx%y = pl%dvx%y + pl%gravity

!     Calculating the neighboring particles and undating HSML
      
!if (sle.ne.0) call h_upgrade(pl)

!     Calculating average velocity of each partile for avoiding penetration

if (pl%average_velocity) call av_vel_omp(pl) !!!!!!!#### Not important!
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
         call interpolation(parts,parts%p,soil,soil%spp)  !!!important
         call pore_water_pressure_omp(parts,soil)         !!!for instability
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
            if(mod(parts%itimestep,40).eq.0) then  !!without this, no
               ntotal = parts%ntotal+parts%nvirt   !!transitional zone
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
                 
!----------------------------------------------------------------------      
                   subroutine single_step_test
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
integer i, j,k,ntotal,it
                 
if(dbg) write(*,*) 'In single_step...'

                                      nphase = 2
                                      it = parts%itimestep
     if(mod(it,10000)==0)then

                                      write(*,*) 'it = ', it
                                      if(parts%single_phase) nphase=1
                                      

                                      do iphase = 1, nphase
                                         if(iphase==1) pl => parts  !!!
                                         if(iphase==2) pl => soil   !!!


call direct_find(pl)     

if(iphase==1)then
        write(*,*) 'sss', parts%get_scale_k()*parts%hsml(1)
do k=1,pl%niac
   i = pl%pair_i(k)
   j = pl%pair_j(k)
   if(i==1485) write(*,*) i,j,sqrt((pl%x(1,i)-pl%x(1,j))**2+(pl%x(2,i)-pl%x(2,j))**2) 
   if(j==1485) write(*,*) j,i,sqrt((pl%x(1,i)-pl%x(1,j))**2+(pl%x(2,i)-pl%x(2,j))**2) 
enddo
endif
write(*,*) '------------------'


if(iphase==2)then
do k=1,pl%niac
   i = pl%pair_i(k)
   j = pl%pair_j(k)
   if(i==703) write(*,*) i,j,sqrt((pl%x(1,i)-pl%x(1,j))**2+(pl%x(2,i)-pl%x(2,j))**2) 
   if(j==703) write(*,*) j,i,sqrt((pl%x(1,i)-pl%x(1,j))**2+(pl%x(2,i)-pl%x(2,j))**2) 
enddo
endif
write(*,*) '---------------'

                                            enddo ! iphase


                if(.not.parts%single_phase)then

         call direct_find_2(parts,soil)

do k=1,parts%niac
   i = parts%pair_i(k)
   j = parts%pair_j(k)
   if(i==1485) write(*,*) i,j 
   if(j==703) write(*,*) 's',j,i 
enddo
write(*,*) '------------------'
     
             endif ! .not.single_phase

    endif !mod

return
end subroutine


!----------------------------------------------------------------------      
      subroutine single_step_GPU
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by using a sorting grid linked list  
!----------------------------------------------------------------------
implicit none

integer i, j, d, dim,scale_k, ntotal, niac, ntotal2    
real(dp) hsml,dr,r,dx(3),tdwdx(3),x(3),dvx(3),tmp(5)
real(dp) :: t1,t2,t3, tw, small = 1.0d-8, miu
integer, parameter :: max_neibor = 30
integer :: nneibor, neibor_j, neibor(max_neibor)
real(dp) :: dudx, dudy, dvdx, dvdy
real(dp) :: rrw,cf,sp,mprr,gravity
real(dp) :: delta, water_c, soil_c
real(dp) :: alpha,beta,etq,epsilon
real(dp), parameter :: pi = 3.14159265358979323846
real(dp) :: factor,q
real(dp) :: phi_eq,J2
real(dp) :: drho,dvx_x,dvx_y,dvof

      integer xgcell(3),xcell, ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)

integer itmod
itmod = mod(parts%itimestep,10000)

dim = parts%dim
scale_k = parts%get_scale_k()
hsml = parts%hsml(1)             !!! hsml is same for all particles!
gravity = parts%gravity
ntotal = parts%ntotal + parts%nvirt
ntotal2= soil%ntotal+soil%nvirt
delta = parts%numeric%delta
water_c  = parts%material%c
soil_c = soil%material%c
alpha = parts%numeric%alpha; beta = parts%numeric%beta; etq=parts%numeric%etq
epsilon = parts%numeric%epsilon

phi_eq = 0.58d0

!$acc update self(parts%x(1:2,ntotal))
!$acc update self(soil%x(1:2,ntotal2))

call create_link_list(parts)
call create_link_list(soil)

!$acc update device(parts%grid,parts%celldata(1:ntotal),parts%maxgridx,parts%mingridx) 
!$acc update device(parts%ngridx,parts%dgeomx,parts%ghsmlx,parts%dim)
!$acc update device(soil%grid, soil%celldata(1:ntotal2),soil%maxgridx,soil%mingridx)
!$acc update device(soil%ngridx,soil%dgeomx,soil%ghsmlx,soil%dim)

ntotal  = parts%ntotal + parts%nvirt

!$omp parallel do private(i) shared(ntotal,H2O)

!$acc kernels loop independent vector
do i = 1, ntotal
  parts%p%r(i) = H2O%b*((parts%rho%r(i)/(H2O%rho0*parts%vof%r(i)))**H2O%gamma-1.d0)
enddo
!$omp end parallel do
!$omp barrier

!call cpu_time(t2)
!write(*,*) 'time:', t2-t1

!$omp parallel do &
!$omp private(i,d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp) &
!$omp private(factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,dim,scale_k,hsml,gravity,H2O,alpha,beta,etq,parts,epsilon)

!$acc kernels loop independent vector &
!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp) &
!$acc private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) 
do i=1,ntotal

parts%drho%r(i) = 0.d0 
parts%dvx%x%r(i) = 0.d0
parts%dvx%y%r(i) = 0.d0
parts%dvof%r(i) = 0.d0
parts%av%x%r(i) = 0.d0
parts%av%y%r(i) = 0.d0

call find_neibors(parts,i,scale_k,hsml,nneibor,neibor,max_neibor)

do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = parts%x(d,i) - parts%x(d,j)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = parts%vx%x%r(i) - parts%vx%x%r(j)
   dvx(2) = parts%vx%y%r(i) - parts%vx%y%r(j)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)  
 
   parts%drho%r(i)=parts%drho%r(i) + (dvx(1)*tdwdx(1)+dvx(2)*tdwdx(2))*parts%mass%r(j)

!---Delta SPH

   tmp(1)=(dx(1)*tdwdx(1)+dx(2)*tdwdx(2))/(dr**2.d0+small)

!Without delta-sph, OK, but particles disordered. With delta-sph, particles well organized,
!e.g. line by line.
if(parts%vof%r(i)+parts%vof%r(j)>1.8.or.parts%vof%r(i)+parts%vof%r(j)<1.2)then
   parts%drho%r(i) = parts%drho%r(i) + 2.d0*delta*hsml*water_c*   &
                    (parts%rho%r(i)-parts%rho%r(j))*tmp(1)*parts%mass%r(j)/parts%rho%r(j)
   parts%dvof%r(i) = parts%dvof%r(i) + 2.d0*delta*hsml*water_c*   &
                    (parts%vof%r(i)-parts%vof%r(j))*tmp(1)*parts%mass%r(j)/parts%rho%r(j)
endif
   
   tmp(2) = -parts%vof%r(i)*parts%mass%r(j)*(parts%p%r(i)+parts%p%r(j))/   &
            (parts%rho%r(i)*parts%rho%r(j))
   parts%dvx%x%r(i) = parts%dvx%x%r(i) + tmp(2) * tdwdx(1)
   parts%dvx%y%r(i) = parts%dvx%y%r(i) + tmp(2) * tdwdx(2)

   tmp(2) = 2.0*H2O%viscosity*tmp(1)*parts%mass%r(j)/(parts%rho%r(i)*parts%rho%r(j))
   parts%dvx%x%r(i) = parts%dvx%x%r(i) + tmp(2) *                                &
         (parts%vof%r(i)*parts%vx%x%r(i)-parts%vof%r(j)*parts%vx%x%r(j))
   parts%dvx%y%r(i) = parts%dvx%y%r(i) + tmp(2) *                                &
         (parts%vof%r(i)*parts%vx%y%r(i)-parts%vof%r(j)*parts%vx%y%r(j))

!---artifical viscosity----
   tmp(3) = dx(1)*dvx(1) + dx(2)*dvx(2)
   if(tmp(3)<0.d0)then
      tmp(3) = hsml*tmp(3)/(dr**2.d0+hsml*hsml*etq*etq)
      tmp(3) = 2.d0*(beta*tmp(3)-alpha*water_c)*tmp(3)*parts%mass%r(j)/        &
              (parts%rho%r(i)+parts%rho%r(j))
      parts%dvx%x%r(i) = parts%dvx%x%r(i) - tmp(3)*tdwdx(1)
      parts%dvx%y%r(i) = parts%dvx%y%r(i) - tmp(3)*tdwdx(2)

!------Average velocity---------

      parts%av%x%r(i) = parts%av%x%r(i) - epsilon*2.0*parts%mass%r(j)*dvx(1)*tw/ &
                        (parts%rho%r(i)+parts%rho%r(j))
      parts%av%y%r(i) = parts%av%y%r(i) - epsilon*2.0*parts%mass%r(j)*dvx(2)*tw/ &
                        (parts%rho%r(i)+parts%rho%r(j))

   endif 
enddo ! neibor
   
   parts%dvx%y%r(i) = parts%dvx%y%r(i) + gravity

enddo !i
!$omp end parallel do
!$omp barrier

!call cpu_time(t1)
!write(*,*) 'time2: ', t1-t2

!--------------------------internal force for soil----------------------------------

ntotal = soil%ntotal + soil%nvirt

!$omp parallel do &
!$omp private(d,i,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$omp private(dudx,dudy,dvdx,dvdy) &
!$omp private(factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,dim,scale_k,hsml,alpha,beta,etq,soil,epsilon) &
!$omp private(J2)

!$acc kernels loop independent vector &
!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$acc private(dudx,dudy,dvdx,dvdy,J2)
do i = 1, ntotal

   soil%drho%r(i) = 0.d0
   soil%dvx%x%r(i) = 0.d0
   soil%dvx%y%r(i) = 0.d0
   soil%dp%r(i) = 0.d0
   soil%dstr%x%r(i) = 0.d0
   soil%dstr%xy%r(i) = 0.d0
   soil%dstr%y%r(i) = 0.d0
   soil%av%x%r(i) = 0.d0
   soil%av%y%r(i) = 0.d0
   dudx = 0.d0; dudy = 0.d0; dvdx = 0.d0; dvdy = 0.d0

   call find_neibors(soil,i,scale_k,hsml,nneibor,neibor,max_neibor)

do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = soil%x(d,i) - soil%x(d,j)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = soil%vx%x%r(i) - soil%vx%x%r(j)
   dvx(2) = soil%vx%y%r(i) - soil%vx%y%r(j)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)

!---continuous equation

    soil%drho%r(i)=soil%drho%r(i) + (dvx(1)*tdwdx(1)+dvx(2)*tdwdx(2))*soil%mass%r(j)

!---Delta SPH

   tmp(1)=(dx(1)*tdwdx(1)+dx(2)*tdwdx(2))/(dr**2.d0+small)

if(soil%vof%r(i)+soil%vof%r(j)>1.8d0.or.soil%vof%r(i)+soil%vof%r(j)<1.2d0)then

   soil%drho%r(i) = soil%drho%r(i) + 2.d0*delta*hsml*soil_c*   &
                    (soil%rho%r(i)-soil%rho%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    
   soil%dp%r(i) = soil%dp%r(i) + 2.d0*delta*hsml*soil_c*   &
                    (soil%p%r(i)-soil%p%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)   
   soil%dstr%x%r(i) = soil%dstr%x%r(i) + 2.d0*delta*hsml*soil_c*   &
              (soil%str%x%r(i)-soil%str%x%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    
   soil%dstr%xy%r(i) = soil%dstr%xy%r(i) + 2.d0*delta*hsml*soil_c*   &
              (soil%str%xy%r(i)-soil%str%xy%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    
   soil%dstr%y%r(i) = soil%dstr%y%r(i) + 2.d0*delta*hsml*soil_c*   &
              (soil%str%y%r(i)-soil%str%y%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    

endif            

!---momentum equation

    tmp(1) = soil%vof%r(i)*soil%p%r(i)/soil%rho%r(i)**2                  &
           + soil%vof%r(j)*soil%p%r(j)/soil%rho%r(j)**2

    tmp(2) = soil%vof%r(i)*soil%str%x%r(i)/soil%rho%r(i)**2              &
           + soil%vof%r(j)*soil%str%x%r(j)/soil%rho%r(j)**2   

    tmp(3) = soil%vof%r(i)*soil%str%xy%r(i)/soil%rho%r(i)**2             &
           + soil%vof%r(j)*soil%str%xy%r(j)/soil%rho%r(j)**2  
    
    tmp(4) = soil%vof%r(i)*soil%str%y%r(i)/soil%rho%r(i)**2              &
           + soil%vof%r(j)*soil%str%y%r(j)/soil%rho%r(j)**2

           
    soil%dvx%x%r(i) = soil%dvx%x%r(i) + soil%mass%r(j)*                  &
                    (-tmp(1)*tdwdx(1)+tmp(2)*tdwdx(1)+tmp(3)*tdwdx(2)) 

    soil%dvx%y%r(i) = soil%dvx%y%r(i) + soil%mass%r(j)*                  &
                    (-tmp(1)*tdwdx(2)+tmp(3)*tdwdx(1)+tmp(4)*tdwdx(2))

!---artifical viscosity----

   tmp(3) = dx(1)*dvx(1) + dx(2)*dvx(2)
   if(tmp(3)<0.d0)then
      tmp(3) = hsml*tmp(3)/(dr**2.d0+hsml*hsml*etq*etq)
      tmp(3) = 2.d0*(beta*tmp(3)-alpha*soil_c)*tmp(3)*soil%mass%r(j)/        &
              (soil%rho%r(i)+soil%rho%r(j))
      soil%dvx%x%r(i) = soil%dvx%x%r(i) - tmp(3)*tdwdx(1)
      soil%dvx%y%r(i) = soil%dvx%y%r(i) - tmp(3)*tdwdx(2)
   endif        

!----Average velocity--------

      soil%av%x%r(i) = soil%av%x%r(i) - epsilon*2.0*soil%mass%r(j)*dvx(1)*tw/ &
                        (soil%rho%r(i)+soil%rho%r(j))
      soil%av%y%r(i) = soil%av%y%r(i) - epsilon*2.0*soil%mass%r(j)*dvx(2)*tw/ &
                        (soil%rho%r(i)+soil%rho%r(j))

!---velocity gradient

    dudx = dudx - soil%mass%r(j)*dvx(1)*tdwdx(1)/soil%rho%r(j)
    dudy = dudy - soil%mass%r(j)*dvx(1)*tdwdx(2)/soil%rho%r(j)
    dvdx = dvdx - soil%mass%r(j)*dvx(2)*tdwdx(1)/soil%rho%r(j)
    dvdy = dvdy - soil%mass%r(j)*dvx(2)*tdwdx(2)/soil%rho%r(j)
enddo !neibor_j

    soil%dvx%y%r(i) = soil%dvx%y%r(i) + gravity
    
    soil%tab%x%r(i) = 2.d0/3.d0*(2.d0*dudx-dvdy)
    soil%tab%xy%r(i)= dudy + dvdx
    soil%tab%y%r(i) = 2.d0/3.d0*(2.d0*dvdy-dudx)

    soil%vcc%r(i)   = dudx + dvdy
    soil%wxy%r(i)   = 0.5d0*(dudy-dvdx)

! No critical state    
!    soil%vof%r(i) = soil%rho%r(i)/sio2%rho0

!Critical state
    J2 = soil%tab%x%r(i)**2.+2.*soil%tab%xy%r(i)**2.+soil%tab%y%r(i)**2.+  &
         (soil%tab%x%r(i)+soil%tab%y%r(i))**2.
    J2 = sqrt(J2/2)

         !if(parts%p%r(i)>1.d-3)then
         !   phi_eq = 0.58 - 25.d0*0.012*J2/parts%p%r(i)         !!! 0.60
         !else
            !phi_eq = 0.58
         !endif

!parts%psi%r(i) = atan(0.5109d0*(parts%vof%r(i)-phi_eq))  !!!K3=4.09,Pouliquen
soil%psi%r(i) = atan(4.09d0*(soil%vof%r(i)-phi_eq))  !!!K3=4.09,Pouliquen
soil%dvof%r(i) = -soil%vof%r(i)*tan(soil%psi%r(i))*J2

enddo  ! i
!$omp end parallel do
!$omp barrier

call cpu_time(t1)
!write(*,*) 'time3: ', t1-t2

call stress_rate_GPU(soil)

!-------------interaction for water---------------------
!call cpu_time(t1)
ntotal = parts%ntotal+parts%nvirt

!$omp parallel do &
!$omp private(d,i,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$omp private(rrw,cf,sp)  &
!$omp private(factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,scale_k,hsml,parts,soil)

!$acc kernels loop independent vector &
!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$acc private(rrw,cf,sp)
do i = 1, ntotal 

call find_neibors2(parts,i,soil,scale_k,hsml,nneibor,neibor,max_neibor)

do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = parts%x(d,i) - soil%x(d,j)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = parts%vx%x%r(i) - soil%vx%x%r(j)
   dvx(2) = parts%vx%y%r(i) - soil%vx%y%r(j)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)  !!!!!!!!


   rrw = tw/(parts%rho%r(i)*soil%rho%r(j))
   cf = 150.*1.e-3*(1-parts%vof%r(i))**2/(parts%vof%r(i)*0.225e-3**2)

   sp = cf*dvx(1)*rrw
   parts%dvx%x%r(i) = parts%dvx%x%r(i)-soil%mass%r(j)*sp
   sp = cf*dvx(2)*rrw
   parts%dvx%y%r(i) = parts%dvx%y%r(i)-soil%mass%r(j)*sp


   parts%dvof%r(i) = parts%dvof%r(i) -    &
   soil%mass%r(j)*soil%vof%r(j)*(dvx(1)*tdwdx(1)+dvx(2)*tdwdx(2))/soil%rho%r(j)

   enddo !neibor_j   
enddo ! i
!$omp end parallel do
!$omp barrier

!call cpu_time(t2)
!write(*,*) 'time4: ', t2-t1

!--------------darcy_law for soil-------------------------------

ntotal = soil%ntotal + soil%nvirt
!$omp parallel do &
!$omp private(d,i,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$omp private(rrw,cf,sp,mprr,factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,dim,scale_k,parts,soil)

!$acc kernels loop independent vector &
!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$acc private(rrw,cf,sp,mprr)

do i = 1, ntotal

   call find_neibors2(soil,i,parts,scale_k,hsml,nneibor,neibor,max_neibor)

soil%spp%r(i) = 0.d0
do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = parts%x(d,j) - soil%x(d,i)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = parts%vx%x%r(j) - soil%vx%x%r(i)
   dvx(2) = parts%vx%y%r(j) - soil%vx%y%r(i)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)  !!!!!!!!

   rrw = tw/(parts%rho%r(j)*soil%rho%r(i))
   cf = 150.*1.e-3*(1-parts%vof%r(j))**2/(parts%vof%r(j)*0.225e-3**2)

   sp = cf*dvx(1)*rrw
   soil%dvx%x%r(i) = soil%dvx%x%r(i)+parts%mass%r(j)*sp
   sp = cf*dvx(2)*rrw
   soil%dvx%y%r(i) = soil%dvx%y%r(i)+parts%mass%r(j)*sp

!--- Pore water pressure

   mprr = parts%mass%r(j)*(parts%p%r(j)-soil%spp0%r(i))/(parts%rho%r(j)*soil%rho%r(i))*soil%vof%r(i)
   soil%dvx%x%r(i) = soil%dvx%x%r(i) + mprr * tdwdx(1)
   soil%dvx%y%r(i) = soil%dvx%y%r(i) + mprr * tdwdx(2)

!--- Pore water pressure interpolation

    soil%spp%r(i) = soil%spp%r(i) +  &
                    parts%mass%r(j)*parts%p%r(j)*tw/parts%rho%r(j)

   enddo !neibor_j   
enddo ! i
!$omp end parallel do
!$omp barrier

!$omp parallel do private(i) shared(ntotal)
!$acc kernels loop independent vector
do i = 1, ntotal
   soil%spp0%r(i) = soil%spp%r(i)
enddo
!$omp end parallel do
!$omp barrier

      if(parts%itimestep>=parts%save_step_from.and.   &
         mod(parts%itimestep,parts%save_step).eq.0)then
         call output
      endif 

end subroutine


!----------------------------------------------------------------------      
      subroutine single_step_water_GPU
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by using a sorting grid linked list  
!----------------------------------------------------------------------
implicit none

integer i, j, d, dim,scale_k, ntotal, niac, ntotal2    
real(dp) hsml,dr,r,dx(3),tdwdx(3),x(3),dvx(3),tmp(5)
real(dp) :: t1,t2,t3, tw, small = 1.0d-8, miu
integer, parameter :: max_neibor = 30
integer :: nneibor, neibor_j, neibor(max_neibor)
real(dp) :: dudx, dudy, dvdx, dvdy
real(dp) :: rrw,cf,sp,mprr,gravity
real(dp) :: delta, water_c, soil_c
real(dp) :: alpha,beta,etq
real(dp), parameter :: pi = 3.14159265358979323846
real(dp) :: factor,q

      integer xgcell(3),xcell, ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)

integer itmod
itmod = mod(parts%itimestep,10000)

dim = parts%dim
scale_k = parts%get_scale_k()
hsml = parts%hsml(1)             !!! hsml is same for all particles!
gravity = parts%gravity
ntotal = parts%ntotal + parts%nvirt
ntotal2= soil%ntotal+soil%nvirt
delta = parts%numeric%delta
water_c  = parts%material%c
soil_c = soil%material%c
alpha = parts%numeric%alpha; beta = parts%numeric%beta; etq=parts%numeric%etq

!call cpu_time(t1)

!$acc update self(parts%x(1:2,ntotal))
!$acc update self(soil%x(1:2,ntotal2))
      
call create_link_list(parts)
call create_link_list(soil)

!$acc update device(parts%grid,parts%celldata(1:ntotal),parts%maxgridx,parts%mingridx) 
!$acc update device(parts%ngridx,parts%dgeomx,parts%ghsmlx,parts%dim)
!$acc update device(soil%grid, soil%celldata(1:ntotal2),soil%maxgridx,soil%mingridx)
!$acc update device(soil%ngridx,soil%dgeomx,soil%ghsmlx,soil%dim)

call cpu_time(t2)
!write(*,*) 'time:', t2-t1

ntotal  = parts%ntotal + parts%nvirt

!$omp parallel do private(i) shared(ntotal,H2O)

!!$acc parallel loop independent present(parts%p%r,parts%rho%r,parts%vof%r) &
!!$acc default(present) !num_gangs(ntotal)

!$acc kernels present(parts%p%r,parts%rho%r,parts%vof%r) &
!$acc default(present) num_gangs(ntotal) !vector_length(64) !num_workers(2)
!$acc loop independent vector
do i = 1, ntotal
  parts%p%r(i) = H2O%b*((parts%rho%r(i)/(H2O%rho0*parts%vof%r(i)))**H2O%gamma-1.d0)
enddo
!$omp end parallel do
!$omp barrier

!$acc end loop
!$acc end kernels

!$omp parallel do &
!$omp private(i,d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp) &
!$omp private(factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,dim,scale_k,hsml,gravity,H2O,alpha,beta,etq,parts)

!$acc loop independent &!present(parts%drho%r,parts%dvx%x%r,parts%dvx%y%r) &
!!$acc present(parts%x,parts%vx%x%r,parts%vx%y%r,parts%mass%r,parts%vof%r)  &
!!$acc present(parts%p%r,parts%rho%r,parts%dvof%r)  &
!!$acc present(parts%maxgridx,parts%mingridx,parts%ngridx,parts%dgeomx) &
!!$acc present(parts%ghsmlx,parts%grid,parts%celldata) &
!!$acc default(present) num_gangs(ntotal) & 
!!$acc num_gangs(ntotal) &
!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp) &
!$acc private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$acc vector


!!$acc kernels present(parts%drho%r,parts%dvx%x%r,parts%dvx%y%r) &
!!$acc present(parts%x,parts%vx%x%r,parts%vx%y%r,parts%mass%r,parts%vof%r)  &
!!$acc present(parts%p%r,parts%rho%r,parts%dvof%r)  &
!!$acc present(parts%maxgridx,parts%mingridx,parts%ngridx,parts%dgeomx) &
!!$acc present(parts%ghsmlx,parts%grid,parts%celldata) !&
!!$acc default(present) num_gangs(ntotal) !vector_length(64) !num_workers(2)
!!$acc num_gangs(100) !vector_length(512) !num_workers(2)

!!$acc loop independent private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp) &
!!$acc private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) !&
!!$acc worker

do i=1,ntotal

parts%drho%r(i) = 0.d0
parts%dvx%x%r(i) = 0.d0
parts%dvx%y%r(i) = 0.d0
parts%dvof%r(i) = 0.d0

call find_neibors(parts,i,scale_k,hsml,nneibor,neibor,max_neibor)

#if abc
      do d=1,3
        xgcell(d) = 1
      enddo
      do d=1,dim
        if ((parts%x(d,i).gt.parts%maxgridx(d)).or.(parts%x(d,i).lt.parts%mingridx(d))) then
          print *,' >>> ERROR <<< : Particle out of range'
          print *,'    Particle position: x(',i,d,') = ',parts%x(d,i)
          print *,'    Range: [xmin,xmax](',D,') =                    &
                 [',parts%mingridx(d),',',parts%maxgridx(d),']'
          stop
        else
          xgcell(d) = int(real(parts%ngridx(d))/parts%dgeomx(d)*      &
                      (parts%x(d,i)-parts%mingridx(d)) + 1.e0)
!          xgcell(d) = int((parts%x(d,i)-parts%mingridx(d))/parts%dx_grid(d)) + 1
        endif
      enddo
        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,dim
          dnxgcell(d) = xgcell(d) - parts%ghsmlx(d) - 1 
          dpxgcell(d) = xgcell(d) + parts%ghsmlx(d) + 1
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),parts%ngridx(d))
        enddo

!     Search grid:
        nneibor = 0; neibor = 0
!$acc loop seq
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = parts%grid(xcell,ycell,zcell)
              if(j==i)j=parts%celldata(j)
              do while(j.gt.0)
                dx(1) = parts%x(1,i) - parts%x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,dim
                  dx(d) = parts%x(d,i) - parts%x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                   nneibor = nneibor + 1 
                   if(nneibor>max_neibor) stop 'max_neibor reached!'
                   neibor(nneibor) = j
                endif
                j = parts%celldata(j)
                if(j==i) j=parts%celldata(j)
               enddo !(j.gt.0)
            enddo !xcell
          enddo !ycell
        enddo !zcell
!$acc end loop
#endif


do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = parts%x(d,i) - parts%x(d,j)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = parts%vx%x%r(i) - parts%vx%x%r(j)
   dvx(2) = parts%vx%y%r(i) - parts%vx%y%r(j)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)  
#if abc
          q=dr/hsml
          factor = 7.e0 / (4.e0*pi*hsml*hsml)
	if(q.ge.0.and.q.le.2) then
          tw = factor * ( (1-q/2)**4 *(1+2*q) )
          do d= 1, dim
            tdwdx(d) = factor*(-5+15*q/2-15*q**2/4+5*q**3/8)/    &
                      hsml**2*dx(d)
          enddo 
	else   
	  tw = 0.
          do d = 1, dim
            tdwdx(d) = 0.
          enddo  
        endif 
#endif
 
   parts%drho%r(i)=parts%drho%r(i) + (dvx(1)*tdwdx(1)+dvx(2)*tdwdx(2))*parts%mass%r(j)

!---Delta SPH

   tmp(1)=(dx(1)*tdwdx(1)+dx(2)*tdwdx(2))/(dr**2.d0+small)

!Without delta-sph, OK, but particles disordered. With delta-sph, particles well organized,
!e.g. line by line.
if(parts%vof%r(i)+parts%vof%r(j)>1.8.or.parts%vof%r(i)+parts%vof%r(j)<1.2)then
   parts%drho%r(i) = parts%drho%r(i) + 2.d0*delta*hsml*water_c*   &
                    (parts%rho%r(i)-parts%rho%r(j))*tmp(1)*parts%mass%r(j)/parts%rho%r(j)
   parts%dvof%r(i) = parts%dvof%r(i) + 2.d0*delta*hsml*water_c*   &
                    (parts%vof%r(i)-parts%vof%r(j))*tmp(1)*parts%mass%r(j)/parts%rho%r(j)
endif
   
   tmp(2) = -parts%vof%r(i)*parts%mass%r(j)*(parts%p%r(i)+parts%p%r(j))/   &
            (parts%rho%r(i)*parts%rho%r(j))
   parts%dvx%x%r(i) = parts%dvx%x%r(i) + tmp(2) * tdwdx(1)
   parts%dvx%y%r(i) = parts%dvx%y%r(i) + tmp(2) * tdwdx(2)

   tmp(2) = 2.0*H2O%viscosity*tmp(1)*parts%mass%r(j)/(parts%rho%r(i)*parts%rho%r(j))
   parts%dvx%x%r(i) = parts%dvx%x%r(i) + tmp(2) *                                &
         (parts%vof%r(i)*parts%vx%x%r(i)-parts%vof%r(j)*parts%vx%x%r(j))
   parts%dvx%y%r(i) = parts%dvx%y%r(i) + tmp(2) *                                &
         (parts%vof%r(i)*parts%vx%y%r(i)-parts%vof%r(j)*parts%vx%y%r(j))

!---artifical viscosity----
   tmp(3) = dx(1)*dvx(1) + dx(2)*dvx(2)
   if(tmp(3)<0.d0)then
      tmp(3) = hsml*tmp(3)/(dr**2.d0+hsml*hsml*etq*etq)
      tmp(3) = 2.d0*(beta*tmp(3)-alpha*water_c)*tmp(3)*parts%mass%r(j)/        &
              (parts%rho%r(i)+parts%rho%r(j))
      parts%dvx%x%r(i) = parts%dvx%x%r(i) - tmp(3)*tdwdx(1)
      parts%dvx%y%r(i) = parts%dvx%y%r(i) - tmp(3)*tdwdx(2)
   endif 
enddo ! neibor

   parts%dvx%y%r(i) = parts%dvx%y%r(i) + gravity

enddo !i
!$omp end parallel do
!$omp barrier

!$acc end loop
!!$acc end kernels

call cpu_time(t1)
write(*,*) 'time2: ', t1-t2

return
end subroutine


!----------------------------------------------------------------------      
      subroutine single_step_soil_GPU
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by using a sorting grid linked list  
!----------------------------------------------------------------------
implicit none

integer i, j, d, dim,scale_k, ntotal, niac, ntotal2    
real(dp) hsml,dr,r,dx(3),tdwdx(3),x(3),dvx(3),tmp(5)
real(dp) :: t1,t2,t3, tw, small = 1.0d-8, miu
integer, parameter :: max_neibor = 30
integer :: nneibor, neibor_j, neibor(max_neibor)
real(dp) :: dudx, dudy, dvdx, dvdy
real(dp) :: rrw,cf,sp,mprr,gravity
real(dp) :: delta, water_c, soil_c
real(dp) :: alpha,beta,etq
real(dp), parameter :: pi = 3.14159265358979323846
real(dp) :: factor,q
real(dp) :: phi_eq,J2

      integer xgcell(3),xcell, ycell,zcell,minxcell(3),maxxcell(3)
      integer dnxgcell(3),dpxgcell(3)

integer itmod
itmod = mod(parts%itimestep,10000)

dim = parts%dim
scale_k = parts%get_scale_k()
hsml = parts%hsml(1)             !!! hsml is same for all particles!
gravity = parts%gravity
ntotal = parts%ntotal + parts%nvirt
ntotal2= soil%ntotal+soil%nvirt
delta = parts%numeric%delta
water_c  = parts%material%c
soil_c = soil%material%c
alpha = parts%numeric%alpha; beta = parts%numeric%beta; etq=parts%numeric%etq

phi_eq = 0.58d0

!call cpu_time(t1)
!--------------------------internal force for soil----------------------------------

ntotal = soil%ntotal + soil%nvirt

!$omp parallel do &
!$omp private(d,i,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$omp private(dudx,dudy,dvdx,dvdy) &
!$omp private(factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,dim,scale_k,hsml,alpha,beta,etq,soil) &
!$omp private(J2)

!!$acc parallel loop independent present(soil%drho%r,soil%dvx%x%r) &
!!$acc present(soil%dvx%y%r,soil%x,soil%vx%x%r,soil%vx%y%r,soil%mass%r) &
!!$acc present(soil%vof%r,soil%p%r,soil%rho%r,soil%str%x%r,soil%str%xy%r) &
!!$acc present(soil%str%y%r,soil%tab%x%r,soil%tab%xy%r,soil%tab%y%r) &
!!$acc present(soil%vcc%r,soil%wxy%r,soil%dp%r,soil%dstr%x%r) &
!!$acc present(soil%dstr,soil%dstr%x,soil%dstr%xy,soil%dstr%y) &
!!$acc present(soil%dstr%xy%r,soil%dstr%y%r,soil%dvof%r,soil%psi%r) &
!!$acc default(present) &
!!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!!$acc private(dudx,dudy,dvdx,dvdy,J2)

!$acc parallel present(soil%drho%r,soil%dvx%x%r) &
!$acc present(soil%dvx%y%r,soil%x,soil%vx%x%r,soil%vx%y%r,soil%mass%r) &
!$acc present(soil%vof%r,soil%p%r,soil%rho%r,soil%str%x%r,soil%str%xy%r) &
!$acc present(soil%str%y%r,soil%tab%x%r,soil%tab%xy%r,soil%tab%y%r) &
!$acc present(soil%vcc%r,soil%wxy%r,soil%dp%r,soil%dstr%x%r) &
!$acc present(soil%dstr,soil%dstr%x,soil%dstr%xy,soil%dstr%y) &
!$acc present(soil%dstr%xy%r,soil%dstr%y%r,soil%dvof%r,soil%psi%r) &
!$acc default(present) num_gangs(ntotal)

!$acc loop independent private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$acc private(dudx,dudy,dvdx,dvdy,J2) vector
do i = 1, ntotal

   soil%drho%r(i) = 0.d0
   soil%dvx%x%r(i) = 0.d0
   soil%dvx%y%r(i) = 0.d0
   soil%dp%r(i) = 0.d0
   soil%dstr%x%r(i) = 0.d0
   soil%dstr%xy%r(i) = 0.d0
   soil%dstr%y%r(i) = 0.d0
   dudx = 0.d0; dudy = 0.d0; dvdx = 0.d0; dvdy = 0.d0

   call find_neibors(soil,i,scale_k,hsml,nneibor,neibor,max_neibor)


#if abc
      do d=1,3
        xgcell(d) = 1
      enddo
!      if(i==1) write(*,*)'adfad'
      do d=1,dim
        if ((soil%x(d,i).gt.soil%maxgridx(d)).or.(soil%x(d,i).lt.soil%mingridx(d))) then
          print *,' >>> ERROR <<< : Particle out of range'
          print *,'    Particle position: x(',i,d,') = ',soil%x(d,i)
          print *,'    Range: [xmin,xmax](',D,') =                    &
                 [',soil%mingridx(d),',',soil%maxgridx(d),']'
          stop
        else
          xgcell(d) = int(real(soil%ngridx(d))/soil%dgeomx(d)*      &
                      (soil%x(d,i)-soil%mingridx(d)) + 1.e0)
          !xgcell(d) = int((soil%x(d,i)-soil%mingridx(d))/soil%dx_grid(d)) + 1
        endif
      enddo
        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,dim
          dnxgcell(d) = xgcell(d)- soil%ghsmlx(d) - 1
          dpxgcell(d) = xgcell(d)+ soil%ghsmlx(d) + 1
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),soil%ngridx(d))
        enddo

!     Search grid:
        nneibor = 0; neibor = 0
!$acc loop seq
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = soil%grid(xcell,ycell,zcell)
              if(j==i)j=soil%celldata(j)
              do while(j.gt.0)
                dx(1) = soil%x(1,i) - soil%x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,dim
                  dx(d) = soil%x(d,i) - soil%x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                   nneibor = nneibor + 1 
                   if(nneibor>max_neibor) stop 'max_neibor reached!'
                   neibor(nneibor) = j
                endif
                j = soil%celldata(j)
                if(j==i) j=soil%celldata(j)
               enddo !(j.gt.0)
            enddo !xcell
          enddo !ycell
        enddo !zcell
!$acc end loop
#endif

!$acc loop seq   
do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = soil%x(d,i) - soil%x(d,j)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = soil%vx%x%r(i) - soil%vx%x%r(j)
   dvx(2) = soil%vx%y%r(i) - soil%vx%y%r(j)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)

!---continuous equation

    soil%drho%r(i)=soil%drho%r(i) + (dvx(1)*tdwdx(1)+dvx(2)*tdwdx(2))*soil%mass%r(j)

!---Delta SPH

   tmp(1)=(dx(1)*tdwdx(1)+dx(2)*tdwdx(2))/(dr**2.d0+small)

if(soil%vof%r(i)+soil%vof%r(j)>1.8d0.or.soil%vof%r(i)+soil%vof%r(j)<1.2d0)then

   soil%drho%r(i) = soil%drho%r(i) + 2.d0*delta*hsml*soil_c*   &
                    (soil%rho%r(i)-soil%rho%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    
   soil%dp%r(i) = soil%dp%r(i) + 2.d0*delta*hsml*soil_c*   &
                    (soil%p%r(i)-soil%p%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)   
   soil%dstr%x%r(i) = soil%dstr%x%r(i) + 2.d0*delta*hsml*soil_c*   &
              (soil%str%x%r(i)-soil%str%x%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    
   soil%dstr%xy%r(i) = soil%dstr%xy%r(i) + 2.d0*delta*hsml*soil_c*   &
              (soil%str%xy%r(i)-soil%str%xy%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    
   soil%dstr%y%r(i) = soil%dstr%y%r(i) + 2.d0*delta*hsml*soil_c*   &
              (soil%str%y%r(i)-soil%str%y%r(j))*tmp(1)*soil%mass%r(j)/soil%rho%r(j)    

endif            

!---momentum equation

    tmp(1) = soil%vof%r(i)*soil%p%r(i)/soil%rho%r(i)**2                  &
           + soil%vof%r(j)*soil%p%r(j)/soil%rho%r(j)**2

    tmp(2) = soil%vof%r(i)*soil%str%x%r(i)/soil%rho%r(i)**2              &
           + soil%vof%r(j)*soil%str%x%r(j)/soil%rho%r(j)**2   

    tmp(3) = soil%vof%r(i)*soil%str%xy%r(i)/soil%rho%r(i)**2             &
           + soil%vof%r(j)*soil%str%xy%r(j)/soil%rho%r(j)**2  
    
    tmp(4) = soil%vof%r(i)*soil%str%y%r(i)/soil%rho%r(i)**2              &
           + soil%vof%r(j)*soil%str%y%r(j)/soil%rho%r(j)**2

           
    soil%dvx%x%r(i) = soil%dvx%x%r(i) + soil%mass%r(j)*                  &
                    (-tmp(1)*tdwdx(1)+tmp(2)*tdwdx(1)+tmp(3)*tdwdx(2)) 

    soil%dvx%y%r(i) = soil%dvx%y%r(i) + soil%mass%r(j)*                  &
                    (-tmp(1)*tdwdx(2)+tmp(3)*tdwdx(1)+tmp(4)*tdwdx(2))

!---artifical viscosity----

   tmp(3) = dx(1)*dvx(1) + dx(2)*dvx(2)
   if(tmp(3)<0.d0)then
      tmp(3) = hsml*tmp(3)/(dr**2.d0+hsml*hsml*etq*etq)
      tmp(3) = 2.d0*(beta*tmp(3)-alpha*soil_c)*tmp(3)*soil%mass%r(j)/        &
              (soil%rho%r(i)+soil%rho%r(j))
      soil%dvx%x%r(i) = soil%dvx%x%r(i) - tmp(3)*tdwdx(1)
      soil%dvx%y%r(i) = soil%dvx%y%r(i) - tmp(3)*tdwdx(2)
   endif        

!---velocity gradient

    dudx = dudx - soil%mass%r(j)*dvx(1)*tdwdx(1)/soil%rho%r(j)
    dudy = dudy - soil%mass%r(j)*dvx(1)*tdwdx(2)/soil%rho%r(j)
    dvdx = dvdx - soil%mass%r(j)*dvx(2)*tdwdx(1)/soil%rho%r(j)
    dvdy = dvdy - soil%mass%r(j)*dvx(2)*tdwdx(2)/soil%rho%r(j)
enddo !neibor_j
!$acc end loop

    soil%dvx%y%r(i) = soil%dvx%y%r(i) + gravity
    
    soil%tab%x%r(i) = 2.d0/3.d0*(2.d0*dudx-dvdy)
    soil%tab%xy%r(i)= dudy + dvdx
    soil%tab%y%r(i) = 2.d0/3.d0*(2.d0*dvdy-dudx)

    soil%vcc%r(i)   = dudx + dvdy
    soil%wxy%r(i)   = 0.5d0*(dudy-dvdx)

! No critical state    
!    soil%vof%r(i) = soil%rho%r(i)/sio2%rho0

!Critical state
    J2 = soil%tab%x%r(i)**2.+2.*soil%tab%xy%r(i)**2.+soil%tab%y%r(i)**2.+  &
         (soil%tab%x%r(i)+soil%tab%y%r(i))**2.
    J2 = sqrt(J2/2)

         !if(parts%p%r(i)>1.d-3)then
         !   phi_eq = 0.58 - 25.d0*0.012*J2/parts%p%r(i)         !!! 0.60
         !else
            !phi_eq = 0.58
         !endif

!parts%psi%r(i) = atan(0.5109d0*(parts%vof%r(i)-phi_eq))  !!!K3=4.09,Pouliquen
soil%psi%r(i) = atan(4.09d0*(soil%vof%r(i)-phi_eq))  !!!K3=4.09,Pouliquen
soil%dvof%r(i) = -soil%vof%r(i)*tan(soil%psi%r(i))*J2

enddo  ! i
!$omp end parallel do
!$omp barrier

!$acc end loop
!$acc end parallel

!call cpu_time(t2)

!write(*,*) 'time3: ', t2-t1


call stress_rate_GPU(soil)


!-------------interaction for water---------------------
!call cpu_time(t1)
ntotal = parts%ntotal+parts%nvirt

!$omp parallel do &
!$omp private(d,i,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$omp private(rrw,cf,sp)  &
!$omp private(factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,scale_k,hsml,parts,soil)

!!$acc parallel loop independent present(parts%dvof%r,parts%x,soil%x) &
!!$acc present(parts%vx%x%r,parts%vx%y%r,soil%vx%x%r,soil%vx%y%r) &
!!$acc present(parts%rho%r,soil%rho%r,parts%vof%r,parts%dvx%x%r)  &
!!$acc present(parts%dvx%y%r,soil%mass%r,soil%vof%r)  &
!!$acc default(present) & !num_gangs(ntotal) &
!!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!!$acc private(rrw,cf,sp)

!$acc parallel present(parts%dvof%r,parts%x,soil%x) &
!$acc present(parts%vx%x%r,parts%vx%y%r,soil%vx%x%r,soil%vx%y%r) &
!$acc present(parts%rho%r,soil%rho%r,parts%vof%r,parts%dvx%x%r)  &
!$acc present(parts%dvx%y%r,soil%mass%r,soil%vof%r)  &
!$acc default(present) num_gangs(ntotal) 

!$acc loop independent private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$acc private(rrw,cf,sp) vector
do i = 1, ntotal 

call find_neibors2(parts,i,soil,scale_k,hsml,nneibor,neibor,max_neibor)

#if abc
      do d=1,3
        xgcell(d) = 1
      enddo

      do d=1,dim
        if ((parts%x(d,i).gt.soil%maxgridx(d)).or.(parts%x(d,i).lt.soil%mingridx(d))) then
          print *,' >>> ERROR <<< : Particle out of range'
          print *,'    Particle position: x(',i,d,') = ',parts%x(d,i)
          print *,'    Range: [xmin,xmax](',D,') =                    &
                 [',soil%mingridx(d),',',soil%maxgridx(d),']'
          stop
        else
          xgcell(d) = int(real(soil%ngridx(d))/soil%dgeomx(d)*      &
                      (parts%x(d,i)-soil%mingridx(d)) + 1.e0)
          !xgcell(d) = int((parts%x(d,i)-soil%mingridx(d))/soil%dx_grid(d)) + 1
        endif
      enddo

        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,dim
          dnxgcell(d) = xgcell(d) - soil%ghsmlx(d) - 1
          dpxgcell(d) = xgcell(d) + soil%ghsmlx(d) + 1
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),soil%ngridx(d))
        enddo
!     Search grid:
        nneibor = 0; neibor = 0
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = soil%grid(xcell,ycell,zcell)
              do while(j.gt.0)
                dx(1) = parts%x(1,i) - soil%x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,dim
                  dx(d) = parts%x(d,i) - soil%x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                   nneibor = nneibor + 1 
                   if(nneibor>max_neibor) stop 'max_neibor reached!'
                   neibor(nneibor) = j
                endif
                j = soil%celldata(j)
                !if(j==i) j=parts%celldata(j)
               enddo !(j.gt.0)
            enddo !xcell
          enddo !ycell
        enddo !zcell
#endif

do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = parts%x(d,i) - soil%x(d,j)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = parts%vx%x%r(i) - soil%vx%x%r(j)
   dvx(2) = parts%vx%y%r(i) - soil%vx%y%r(j)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)  !!!!!!!!
#if abc
          q=dr/hsml
          factor = 7.e0 / (4.e0*pi*hsml*hsml)
	if(q.ge.0.and.q.le.2) then
          tw = factor * ( (1-q/2)**4 *(1+2*q) )
          do d= 1, dim
            tdwdx(d) = factor*(-5+15*q/2-15*q**2/4+5*q**3/8)/    &
                      hsml**2*dx(d)
          enddo 
	else   
	  tw = 0.
          do d = 1, dim
            tdwdx(d) = 0.
          enddo  
        endif 
#endif

   rrw = tw/(parts%rho%r(i)*soil%rho%r(j))
   cf = 150.*1.e-3*(1-parts%vof%r(i))**2/(parts%vof%r(i)*0.225e-3**2)

   sp = cf*dvx(1)*rrw
   parts%dvx%x%r(i) = parts%dvx%x%r(i)-soil%mass%r(j)*sp
   sp = cf*dvx(2)*rrw
   parts%dvx%y%r(i) = parts%dvx%y%r(i)-soil%mass%r(j)*sp


   parts%dvof%r(i) = parts%dvof%r(i) -    &
   soil%mass%r(j)*soil%vof%r(j)*(dvx(1)*tdwdx(1)+dvx(2)*tdwdx(2))/soil%rho%r(j)

   enddo !neibor_j   
enddo ! i
!$omp end parallel do
!$omp barrier

!$acc end loop
!$acc end parallel

!call cpu_time(t2)
!write(*,*) 'time4: ', t2-t1

!--------------darcy_law for soil-------------------------------

ntotal = soil%ntotal + soil%nvirt
!$omp parallel do &
!$omp private(d,i,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$omp private(rrw,cf,sp,mprr,factor,q) &
!$omp private(xgcell,minxcell,maxxcell,dnxgcell,dpxgcell,zcell,ycell,xcell) &
!$omp shared(ntotal,dim,scale_k,parts,soil)

!!$acc parallel loop independent present(parts%x,soil%x,parts%vx%x%r)   &
!!$acc present(parts%vx%y%r,soil%vx%x%r,soil%vx%y%r,parts%rho%r,soil%rho%r) &
!!$acc present(parts%vof%r,soil%dvx%x%r,soil%dvx%y%r,parts%mass%r) &
!!$acc present(parts%p%r,soil%vof%r,soil%spp%r,soil%spp0%r) &
!!$acc default(present) &
!!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!!$acc private(rrw,cf,sp,mprr)


!$acc parallel present(parts%x,soil%x,parts%vx%x%r)   &
!$acc present(parts%vx%y%r,soil%vx%x%r,soil%vx%y%r,parts%rho%r,soil%rho%r) &
!$acc present(parts%vof%r,soil%dvx%x%r,soil%dvx%y%r,parts%mass%r) &
!$acc present(parts%p%r,soil%vof%r,soil%spp%r,soil%spp0%r) &
!$acc default(present) num_gangs(ntotal)

!$acc loop independent vector &
!$acc private(d,j,dx,dr,nneibor,neibor,neibor_j,tw,tdwdx,miu,dvx,tmp)  &
!$acc private(rrw,cf,sp,mprr)

do i = 1, ntotal

   call find_neibors2(soil,i,parts,scale_k,hsml,nneibor,neibor,max_neibor)

#if abc

      do d=1,3
        xgcell(d) = 1
      enddo

      do d=1,dim
        if ((soil%x(d,i).gt.parts%maxgridx(d)).or.(soil%x(d,i).lt.parts%mingridx(d))) then
          print *,' >>> ERROR <<< : Particle out of range'
          print *,'    Particle position: x(',i,d,') = ',soil%x(d,i)
          print *,'    Range: [xmin,xmax](',D,') =                    &
                 [',parts%mingridx(d),',',parts%maxgridx(d),']'
          stop
        else
          xgcell(d) = int(real(parts%ngridx(d))/parts%dgeomx(d)*      &
                      (soil%x(d,i)-parts%mingridx(d)) + 1.e0)
          !xgcell(d) = int((soil%x(d,i)-parts%mingridx(d))/parts%dx_grid(d)) + 1
        endif
      enddo

        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,dim
          dnxgcell(d) = xgcell(d) - parts%ghsmlx(d) - 1
          dpxgcell(d) = xgcell(d) + parts%ghsmlx(d) + 1
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),parts%ngridx(d))
        enddo
!     Search grid:
        nneibor = 0; neibor = 0
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = parts%grid(xcell,ycell,zcell)
              do while(j.gt.0)
                dx(1) = soil%x(1,i) - parts%x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,dim
                  dx(d) = soil%x(d,i) - parts%x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                   nneibor = nneibor + 1 
                   if(nneibor>max_neibor) stop 'max_neibor reached!'
                   neibor(nneibor) = j
                endif
                j = parts%celldata(j)
                !if(j==i) j=parts%celldata(j)
               enddo !(j.gt.0)
            enddo !xcell
          enddo !ycell
        enddo !zcell

#endif

soil%spp%r(i) = 0.d0
do neibor_j = 1, nneibor
   dr = 0.d0
   j = neibor(neibor_j)
   do d = 1, dim
      dx(d) = parts%x(d,j) - soil%x(d,i)
      dr = dr + dx(d)*dx(d)
   enddo
   dr = dsqrt(dr)   

   dvx(1) = parts%vx%x%r(j) - soil%vx%x%r(i)
   dvx(2) = parts%vx%y%r(j) - soil%vx%y%r(i)

!--- Kernel and derivations of kernel

   call kernel_GPU(dr,dx,hsml,tw,tdwdx,2,4)  !!!!!!!!
#if abc
          q=dr/hsml
          factor = 7.e0 / (4.e0*pi*hsml*hsml)
	if(q.ge.0.and.q.le.2) then
          tw = factor * ( (1-q/2)**4 *(1+2*q) )
          do d= 1, dim
            tdwdx(d) = factor*(-5+15*q/2-15*q**2/4+5*q**3/8)/    &
                      hsml**2*dx(d)
          enddo 
	else   
	  tw = 0.
          do d = 1, dim
            tdwdx(d) = 0.
          enddo  
        endif 
#endif   


   rrw = tw/(parts%rho%r(j)*soil%rho%r(i))
   cf = 150.*1.e-3*(1-parts%vof%r(j))**2/(parts%vof%r(j)*0.225e-3**2)

   sp = cf*dvx(1)*rrw
   soil%dvx%x%r(i) = soil%dvx%x%r(i)+parts%mass%r(j)*sp
   sp = cf*dvx(2)*rrw
   soil%dvx%y%r(i) = soil%dvx%y%r(i)+parts%mass%r(j)*sp

!--- Pore water pressure

   mprr = parts%mass%r(j)*(parts%p%r(j)-soil%spp0%r(i))/(parts%rho%r(j)*soil%rho%r(i))*soil%vof%r(i)
   soil%dvx%x%r(i) = soil%dvx%x%r(i) + mprr * tdwdx(1)
   soil%dvx%y%r(i) = soil%dvx%y%r(i) + mprr * tdwdx(2)

!--- Pore water pressure interpolation

    soil%spp%r(i) = soil%spp%r(i) +  &
                    parts%mass%r(j)*parts%p%r(j)*tw/parts%rho%r(j)

   enddo !neibor_j   
enddo ! i
!$omp end parallel do
!$omp barrier

!$acc end loop
!$acc end parallel

!$omp parallel do private(i) shared(ntotal)
!!$acc parallel loop independent present(soil%spp%r,soil%spp0%r) &
!!$acc default(present)

!$acc parallel present(soil%spp%r,soil%spp0%r) &
!$acc default(present) num_gangs(ntotal)

!$acc loop independent vector
do i = 1, ntotal
   soil%spp0%r(i) = soil%spp%r(i)
enddo
!$omp end parallel do
!$omp barrier

!$acc end loop
!$acc end parallel

      if(parts%itimestep>=parts%save_step_from.and.   &
         mod(parts%itimestep,parts%save_step).eq.0)then
         call output
      endif 

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

      call update_self(parts)
      call update_self(soil)

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
