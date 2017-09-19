    module m_particles
    !---------------------------
    use m_array
    implicit none

    type keyword
        character(len=72) name

    endtype

    type realarray
        real, pointer, dimension(:,:) :: r
    endtype

    type block
        real(dp) xl, yl
        integer  m,  n, np
        real(dp) dx, dy
        real(dp), pointer, dimension(:) :: x,y
        integer,  pointer, dimension(:) :: zone
        ! boundary point
        real(dp), pointer, dimension(:) :: bx, by
        integer,  pointer, dimension(:) :: bzone
        real(dp), pointer, dimension(:,:) :: bn
        real(dp), pointer, dimension(:) :: bs
        ! vertex point
        real(dp), pointer, dimension(:) :: vx,vy
        integer,  pointer, dimension(:) :: vzone

        type(particles), pointer :: p1=>null(), p2=>null()
    contains
    procedure :: set => set_block_sub
    procedure :: alloc => block_allocate_sub
    procedure :: cell_center => get_cell_center_sub
    procedure :: get_vertex => get_vertex_sub
    procedure :: get_boundary_nodes
    procedure :: get_segment_centers

    end type


    !-------------------------------------------------------------------------------------
    type material

    ! For water

    real(dp) rho0,b,gamma,c,viscosity

    ! For soil

    real(dp) k,porosity,permeability,G,E,niu,cohesion,phi,psi

    end type
    !------------------------------------------------------------------------------------

    type numerical

    ! Artificial viscosity
    ! alpha: shear viscosity; beta: bulk viscosity; etq: parameter to avoid sigularities   
    real(dp) :: alpha=0.1d0, beta=0.d0, etq=0.1d0

    ! Leonard_Johns repulsive force
    real(dp) :: dd = 0.1d0, p1 = 12, p2 = 4

    ! Delta-SPH
    real(dp) :: delta = 0.1d0

    ! Velocity average
    real(dp) :: epsilon = 0.001

    end type
    !-----------------------------------------------------------------------------------

    ! Particles in SPH method
    !---------------
    type, extends(dao) :: particles
        !---------------

        !Geometry parameters

        !dim : Dimension of the problem (1, 2 or 3)
        integer :: dim = 2        

        real(dp) dspp  ! initial particle interval  

        !Parameters for the computational geometry,  
        !x_maxgeom : Upper limit of allowed x-regime 
        !x_mingeom : Lower limit of allowed x-regime 
        !y_maxgeom : Upper limit of allowed y-regime 
        !y_mingeom : Lower limit of allowed y-regime 
        !z_maxgeom : Upper limit of allowed z-regime 
        !z_mingeom : Lower limit of allowed z-regime 

        real(dp) :: x_maxgeom, x_mingeom, y_maxgeom, y_mingeom, z_maxgeom, z_mingeom

        !Parameter used for sorting grid cells in the link list algorithm
        !maxngx  : Maximum number of sorting grid cells in x-direction
        !maxngy  : Maximum number of sorting grid cells in y-direction
        !maxngz  : Maximum number of sorting grid cells in z-direction
        !Determining maximum number of sorting grid cells:
        !(For an homogeneous particle distribution:)
        !1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
        !2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
        !3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)
        integer :: maxngx = 100,maxngy = 100, maxngz = 1

        integer, pointer, dimension(:,:,:) :: grid => null()
        integer, pointer, dimension(:,:) :: xgcell => null()
        integer, pointer, dimension(:) :: celldata => null() 
        integer  ngridx(3),ghsmlx(3)
        real(dp) mingridx(3),maxgridx(3),dgeomx(3)

        !Memory allocation

        !maxn: Maximum number of particles
        !max_interation : Maximum number of interaction pairs
        integer :: maxn = 200000, max_interaction = 20 * 200000

        !SPH algorithm

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

        logical :: single_phase = .true.
        integer :: integrate_scheme = 1  ! =1, LF; =2, Verlet
        logical :: summation_density  = .false.         
        logical :: average_velocity  = .true. 
        !logical :: visc  = .true.  
        logical :: ex_force  = .true.
        logical :: visc_artificial  = .true. 
        !logical :: nor_density  = .false.              
        !logical :: self_gravity  = .true.      

        ! water_pressure: 0->hydro-static, 1->dynamic
        integer :: water_pressure = 1
        ! Particle size: 1-> coarse; 2->fine
        integer :: particle_size = 1
        integer:: nsym = 0

        !integer :: soil_pressure = 2  ! =1, eos; =2, mean trace
        integer :: stress_integration = 1
        integer :: yield_criterion = 2
        integer :: plasticity = 3  ! =0 non; =1 Bui, =2 return mapping =3 Lopez
        integer :: flow_rule = 1  ! =1 associated; =2 non-associated
        logical :: artificial_density = .true.                  
        logical :: soil_artificial_stress = .true.

        logical :: volume_fraction = .true.
        logical :: water_artificial_volume = .true.
        logical :: volume_fraction_renorm = .true.

        integer :: critical_state = 0
        logical :: usaw = .false.
        logical :: particle_shift = .false.

        ! 0 ignor; 1 negative pressure to zero; 2 artficial stress
        !integer :: water_tension_instability = 0

        !Particle approximation (pa_sph)
        !pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
        !         2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
        integer :: pa_sph = 2 

        !Nearest neighbor particle searching (nnps) method
        !nnps = 1 : Simplest and direct searching
        !       2 : Sorting grid linked list
        !       3 : Tree algorithm
        integer :: nnps = 1 

        !Smoothing length evolution (sle) algorithm
        !sle = 0 : Keep unchanged,
        !      1 : h = fac * (m/rho)^(1/dim)
        !      2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
        !      3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) ) 
        integer :: sle = 0 

        !Smoothing kernel function 
        !skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
        !    = 2, Gauss kernel   (Gingold and Monaghan 1981) 
        !    = 3, Quintic kernel (Morris 1997)
        !    = 4, Wendland    
        integer :: skf = 4 

        !Material parameters
        character(len=32) :: imaterial

        !PGI FORTRAN doesn't support *
        class(material), pointer :: material => null()

        !class(*), pointer :: material => null()

        !Numerical parameters
        class(numerical), pointer :: numeric => null()

        !Physics parameter
        real(dp) :: gravity = -9.8

        !Time stepping
        real(dp) :: dt = 0.d0, time = 0.d0
        integer :: maxtimestep = 0 , itimestep = 0

        ! Control parameters for output 
        ! int_stat = .true. : Print statistics about SPH particle interactions.
        !                     including virtual particle information.
        ! print_step: Print Timestep (On Screen)
        ! save_step : Save Timestep    (To Disk File)
        ! moni_particle: The particle number for information monitoring.

        logical :: int_stat = .true.
        integer :: print_step, save_step, moni_particle = 264

        !Recorde time interval
        integer :: save_step_from = 0, save_step_to = 100

        !Paralle computing using OpenMP
        integer :: nthreads = 1
        integer, pointer, dimension(:) :: niac_start, niac_end

        integer :: niac  = 0
        integer :: ntotal = 0, nvirt = 0

        ! For symmetry boundary
        integer  :: nsymm = 0 

        ! Particle interaction pair
        integer :: maxp = 0, minp = 0
        integer :: maxiac = 0, miniac = 0
        integer :: sumiac = 0, noiac  = 0
        integer, pointer, dimension(:)  :: pair_i => null()
        integer, pointer, dimension(:)  :: pair_j => null()
        integer, pointer, dimension(:)  :: countiac=>null()

        ! Kernel and its derivative
        real(dp), pointer, dimension(:)   :: w    => null()
        real(dp), pointer, dimension(:,:) :: dwdx => null()

        ! Sherpard Filter
        type(array), pointer :: wi => null()



        ! Particle fundamental data
        integer,  pointer, dimension(:)   :: itype=> null()
        integer,  pointer, dimension(:)   :: iDeath => null()
        real(dp), pointer, dimension(:,:) :: x    => null()
        type(array), pointer :: vol  => null()
        type(array), pointer :: mass => null()
        real(dp), pointer, dimension(:)   :: hsml => null()
        integer,  pointer, dimension(:)   :: zone => null()

        ! Volume of Fraction
        type(array), pointer :: vof  => null()
        type(array), pointer :: vof2 =>null() ! phi_f= 1-phi_s
        type(array), pointer :: dvof  => null()
        type(array), pointer :: vof_min  => null()

        ! Dilatancy angle
        type(array), pointer :: psi => null()

        ! Field variables
        type(array), pointer :: rho  => null()  ! Density
        type(array), pointer :: rho0 => null()  ! True density, for critical state
        type(array), pointer :: vx   => null()  ! Velocity
        type(array), pointer :: p => null()     ! Pressure
        type(array), pointer :: c => null()     ! Sound speed

        type(array), pointer :: spp => null()   !!!! VOF pore pressure for soil particles

        ! Stress tensor 
        type(array), pointer :: str => null()

        ! artificial stress tensor
        type(array), pointer :: strp => null()   !stress tensor in the principal coordinate
        type(array), pointer :: thetai => null()
        type(array), pointer :: ast => null()    !artificial stress tensor
        type(array), pointer :: astp => null()   !artificial stress tensor in the principal coordinate

        ! Shear strain rate
        type(array), pointer :: tab => null()
        ! Bulk strain rate, i.e. divergency of velocity
        type(array), pointer :: vcc => null()

        ! Accumulative deviatoric strain
        type(array), pointer :: epsilon_p => null()
        ! Rotation tensor
        !   real(dp), pointer, dimension(:) :: wxy => null()
        type(array), pointer :: wxy => null()
        ! Old values
        !   real(dp), pointer, dimension(:)   :: rho_min => null()
        type(array), pointer :: rho_min => null()
        type(array), pointer :: v_min => null()
        type(array), pointer :: p_min   => null()
        type(array), pointer :: str_min => null()

        ! Acceleration
        type(array), pointer :: drho => null()
        type(array), pointer :: dvx  => null()
        type(array), pointer :: dp   => null()
        type(array), pointer :: av   => null()

        ! Particle shifting
        type(array), pointer :: r0 => null()
        type(array), pointer :: mt => null()
        type(array), pointer :: ps => null()

        type(array), pointer :: dstr => null()

        ! For plasticity corrector step (return mapping algorithm)
        type(array), pointer :: dstr2 => null()
        type(array), pointer :: dp2 => null()

        ! For soil failure
        integer :: nfail = 0
        integer, pointer, dimension(:) :: fail => null()

        ! Boundry particles defined as type particles
        type(particles), pointer :: bor
        real(dp), pointer, dimension(:,:) :: norm
        type(array), pointer :: bn => null()
        type(array), pointer :: bs => null()


    contains

    procedure :: write_particles
    procedure :: write_component 
    procedure :: append_virtual_particles
    procedure :: depend_virtual_particles
    procedure :: interaction_statistics
    procedure :: minimum_time_step
    procedure :: take_real => take_real_points1
    procedure :: take_virtual => take_virtual_points1
    procedure :: take_boundary_nodes
    procedure :: setup_itype
    procedure :: setup_ndim1
    procedure :: get_size => get_total_particles_number
    procedure :: get_maxn => get_maximum_particles_number
    procedure :: get_dim => get_dimension
    procedure :: get_scale_k
    procedure :: direct_find
    procedure :: init_grid
    procedure :: grid_geom
    procedure :: link_list
    procedure :: link_list2
    procedure :: find_pairs
    procedure :: calc_norm
    procedure :: kernel
    procedure :: initial_density
    procedure :: pressure
    procedure :: pressure_water
    procedure :: pressure_soil
    procedure :: tension_instability_water 
    procedure :: sum_density
    procedure :: con_density
    procedure :: newtonian_fluid
    procedure :: art_visc
    procedure :: art_visc_omp
    procedure :: delta_sph
    procedure :: delta_sph_omp
    procedure :: diffusive_sph_omp
    procedure :: av_vel
    procedure :: repulsive_force
    procedure :: repulsive_force_omp
    procedure :: df
    procedure :: df_omp
    procedure :: df_omp2
    procedure :: df2
    procedure :: df3
    procedure :: df3_omp
    procedure :: df4
    procedure :: df4_omp
    procedure :: div
    procedure :: div_omp
    procedure :: div2
    procedure :: div2_USAW
    procedure :: lap_omp
    procedure :: velocity_divergence
    procedure :: shear_strain_rate
    procedure :: artificial_stress
    procedure :: Jaumann_rate
    procedure :: mohr_coulomb_failure_criterion
    procedure :: drucker_prager_failure_criterion
    procedure :: plastic_or_not
    procedure :: plastic_flow_rule
    procedure :: plastic_flow_rule2
    procedure :: plastic_flow_rule3
    procedure :: stress_rate
    procedure :: scaling_back
    procedure :: get_num_threads
    procedure :: get_niac_start_end

    end type 

    !PGI fortran compiler doesn't support ::

    interface write_field
    !module procedure ::  write_scalar_field
    module procedure  write_scalar_field
    module procedure  write_vector2d_field
    end interface

    !=======
    contains
    !=======

#ifdef REMOVE
    !#if(.false.)
    !!DEC$IF(.FALSE.)
    !---------------------------------------------------------------------------------
    subroutine set_geo_sub(this,dl,dr,dh,ds,dj,D,L,H,immerse,       &
    mdl,mdr,mdh,mds,mdj,mD,mL,mH)
    !---------------------------------------------------------------------------------
    implicit none
    class(geo) this
    real(dp), optional :: dl,dr,dh,ds,dj,D,L,H,immerse
    integer,  optional :: mdl,mdr,mdh,mds,mdj,mD,mL,mH

    if(present(dl)) this%dl = dl  
    if(present(dr)) this%dr = dr  
    if(present(dh)) this%dh = dh  
    if(present(ds)) this%ds = ds  
    if(present(dj)) this%dj = dj  
    if(present(D)) this%D = D  
    if(present(L)) this%L = L  
    if(present(H)) this%H = H  
    if(present(immerse)) this%immerse = immerse

    if(present(mdl)) this%mdl = mdl  
    if(present(mdr)) this%mdr = mdr  
    if(present(mdh)) this%mdh = mdh  
    if(present(mds)) this%mds = mds  
    if(present(mdj)) this%mdj = mdj  
    if(present(mD)) this%mD = mD  
    if(present(mL)) this%mL = mL  
    if(present(mH)) this%mH = mH  

    return
    end subroutine

    !!DEC$ENDIF
#endif

    !--------------------------------------------
    subroutine set_block_sub(this,xl,yl,m,n)
    !--------------------------------------------
    implicit none
    class(block) this
    real(dp),optional :: xl,yl
    integer, optional :: m, n

    if(present(xl)) this%xl = xl
    if(present(yl)) this%yl = yl
    if(present(m))  this%m  = m
    if(present(n))  this%n  = n

    return
    end subroutine

    !------------------------------------------------
    subroutine block_allocate_sub(this)
    !------------------------------------------------
    implicit none
    class(block) this
    integer nbp

    nbp = (this%m+this%n)*2

    allocate(this%bx(nbp),this%by(nbp))
    allocate(this%bzone(nbp))

    allocate(this%bn(2,nbp))
    allocate(this%bs(nbp))

    return
    end subroutine

    !------------------------------------------------
    subroutine get_cell_center_sub(this)
    !------------------------------------------------
    implicit none
    class(block) this
    integer i,j,k

    this%dx = this%xl/this%m
    this%dy = this%yl/this%n

    k = 0
    do i = 1, this%m
        do j = 1, this%n
            k = k + 1
            this%x(k) = (i-1)*this%dx + this%dx/2.
            this%y(k) = (j-1)*this%dy + this%dy/2.
        enddo
    enddo

    return
    end subroutine


    !------------------------------------------------
    subroutine get_vertex_sub(this)
    !------------------------------------------------
    implicit none
    class(block) this
    integer i,j,k

    this%dx = this%xl/this%m
    this%dy = this%yl/this%n

    !! Boundary points!

    k = 0
    do j = 1, this%n
        k = k + 1
        this%x(k) = 0.0
        this%y(k) = this%yl-(j-1)*this%dy
    enddo
    do i = 1, this%m
        k = k + 1
        this%x(k) = (i-1)*this%dx
        this%y(k) = 0.0
    enddo
    do j = 1, this%n
        k = k + 1
        this%x(k) = this%xl
        this%y(k) = 0.d0 + (j-1)*this%dy
    enddo
    do i = 1, this%m
        k = k + 1
        this%x(k) = this%xl - (i-1)*this%dx
        this%y(k) = this%yl
    enddo

    !! Inner points

    do i = 1, this%m-1
        do j = 1, this%n-1
            k = k + 1
            this%x(k) = i*this%dx
            this%y(k) = j*this%dy
        enddo
    enddo

    return
    end subroutine

    !------------------------------------------------
    subroutine get_boundary_nodes(this)
    !------------------------------------------------
    implicit none
    class(block) this
    integer i,j,k

    this%dx = this%xl/this%m
    this%dy = this%yl/this%n

    !! Boundary points!

    k = 0
    do j = 1, this%n
        k = k + 1
        this%bx(k) = 0.0
        this%by(k) = this%yl-(j-1)*this%dy
    enddo
    do i = 1, this%m
        k = k + 1
        this%bx(k) = (i-1)*this%dx
        this%by(k) = 0.0
    enddo
    do j = 1, this%n
        k = k + 1
        this%bx(k) = this%xl
        this%by(k) = 0.d0 + (j-1)*this%dy
    enddo
    do i = 1, this%m
        k = k + 1
        this%bx(k) = this%xl - (i-1)*this%dx
        this%by(k) = this%yl
    enddo

    return
    end subroutine


    !------------------------------------------------
    subroutine get_segment_centers(this)
    !------------------------------------------------
    implicit none
    class(block) this
    integer i,j,k

    this%dx = this%xl/this%m
    this%dy = this%yl/this%n

    !! Boundary points!

    k = 0
    do j = 1, this%n
        k = k + 1
        this%bx(k) = 0.0
        this%by(k) = this%yl-(j-1)*this%dy-0.5*this%dy
        this%bn(1,k) = 1.d0
        this%bn(2,k) = 0.d0
        this%bs(k) = this%dy 
    enddo
    do i = 1, this%m
        k = k + 1
        this%bx(k) = (i-1)*this%dx+0.5*this%dx
        this%by(k) = 0.0
        this%bn(1,k) = 0.d0
        this%bn(2,k) = 1.d0
        this%bs(k) = this%dx       
    enddo
    do j = 1, this%n
        k = k + 1
        this%bx(k) = this%xl
        this%by(k) = 0.d0 + (j-1)*this%dy+0.5*this%dy
        this%bn(1,k) = -1.d0
        this%bn(2,k) = 0.d0
        this%bs(k) = this%dy       
    enddo
    do i = 1, this%m
        k = k + 1
        this%bx(k) = this%xl - (i-1)*this%dx-0.5*this%dx
        this%by(k) = this%yl
        this%bn(1,k) = 0.d0
        this%bn(2,k) = -1.d0
        this%bs(k) = this%dx       
    enddo

    return
    end subroutine

    !----------------------------------------
    subroutine write_particles(this) 
    !----------------------------------------
    implicit none
    class(particles) this
    integer i

    write(*,*) 'Real particles:'
    write(*,*) 'Num.  itype  x   y   vol  hsml'
    do i = 1, this%ntotal
        !write(*,*) i,  this%itype(i), this%x(1,i), this%x(2,i), this%vol(i), this%hsml(i) 
    end do

    write(*,*) 'Virtual particles:'
    write(*,*) 'Num.  itype  x   y   vol  hsml'
    do i = this%ntotal+1, this%ntotal+this%nvirt
        !write(*,*) i,  this%itype(i), this%x(1,i), this%x(2,i), this%vol(i), this%hsml(i) 
    end do

    return
    end subroutine

    !---------------------------------------------
    subroutine write_component(this,vars) 
    !---------------------------------------------
    implicit none
    class(particles) this
    character(len=*) vars
    integer ntotal, nvirt, i

    ntotal = this%ntotal
    nvirt  = this%nvirt

    select case (trim(vars))
    case ('x')
        write(*,*) 'Particle position:'
        do i = 1, ntotal+nvirt
            write(*,*) i,this%x(1,i),this%x(2,i)
        enddo
    case ('mass')
        write(*,*) 'Particle mass:'
        do i = 1, ntotal+nvirt
            write(*,*) i, this%mass%r(i)
        enddo
    case('w')
        write(*,*) 'Kernels:'
        write(*,*) (i,this%w(i), i = 1, this%niac)
    case('dwdx')
        write(*,*) 'Derivatives of Kernels:'
        write(*,*) (i,this%dwdx(1,i),this%dwdx(2,i),i = 1, this%niac)
    case ('pairs')
        write(*,*) 'Interaction pairs:'
        write(*,*) (i, this%pair_i(i), this%pair_j(i),'|',i=1,this%niac)
        case default
        write(*,*) 'Write this compoent is not implemented yet!'
    end select

    return
    end subroutine

    !---------------------------------------------------
    subroutine count_interaction(parts)
    !---------------------------------------------------
    implicit none
    class(particles) parts
    integer ntotal,i,j,k

    ntotal = parts%ntotal + parts%nvirt

    do i = 1, ntotal
        parts%countiac(i) = 0
    enddo
    do k = 1, parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        parts%countiac(i) = parts%countiac(i) + 1
        parts%countiac(j) = parts%countiac(j) + 1
    enddo

    return
    end subroutine

    ! Statistics for the interaction
    !-------------------------------------------------
    subroutine interaction_statistics(parts)
    !-------------------------------------------------
    implicit none
    class(particles) parts
    integer ntotal, i,j,d,k
    integer  sumiac, maxiac, miniac, noiac, maxp, minp
    logical :: dbg = .false.

    ntotal   =   parts%ntotal + parts%nvirt
    sumiac = 0
    maxiac = 0
    miniac = 1000
    noiac  = 0

    !For Particle Shifting, we need parts%countiac
    !do i = 1, ntotal
    !   parts%countiac(i) = 0
    !enddo
    !do k = 1, parts%niac
    !   i = parts%pair_i(k)
    !   j = parts%pair_j(k)
    !   parts%countiac(i) = parts%countiac(i) + 1
    !   parts%countiac(j) = parts%countiac(j) + 1
    !enddo
    !-----------------------------------------------

    do i=1,ntotal
        sumiac = sumiac + parts%countiac(i)
        if (parts%countiac(i).gt.maxiac) then
            maxiac = parts%countiac(i)
            maxp = i
        endif
        if (parts%countiac(i).lt.miniac) then 
            miniac = parts%countiac(i)
            minp = i
        endif
        if (parts%countiac(i).eq.0) noiac  = noiac + 1
    enddo

    parts%maxp = maxp;        parts%minp = minp  
    parts%maxiac = maxiac;    parts%miniac = miniac
    parts%sumiac = sumiac;    parts%noiac = noiac

    print *,' >> Statistics: interactions per particle:'
    print *,'**** Particle:',parts%maxp,                 &
    ' maximal interactions:', parts%maxiac
    print *,'**** Particle:',parts%minp,                 &
    ' minimal interactions:', parts%miniac
    print *,'**** Average :',real(parts%sumiac)/         &
    real(parts%ntotal+parts%nvirt)
    print *,'**** Total pairs : ',parts%niac
    print *,'**** Particles with no interactions:',parts%noiac

    return
    end subroutine

    !-------------------------------------------------
    subroutine minimum_time_step(parts)
    !-------------------------------------------------
    implicit none
    class(particles) parts

    type(material), pointer :: mat
    integer ntotal
    double precision c, miu 
    double precision mint1, mint21, mint22, mint3

    ntotal = parts%ntotal

    mat => parts%material
    c   = mat%c
    miu = mat%viscosity/mat%rho0

    ! CFL condition
    mint1 = 0.25*minval(parts%hsml(1:ntotal)/c)
    ! Acceleration 
    !mint21 = 0.25*minval(sqrt(parts%hsml(1:ntotal)/parts%dvx(1,1:ntotal)))
    !mint22 = 0.25*minval(sqrt(parts%hsml(1:ntotal)/parts%dvx(2,1:ntotal)))
    ! Viscous diffustion time
    mint3 = 0.125*minval(parts%hsml(1:ntotal)**2/miu)

    !write(*,*) 'Minimum time step is: '
    !write(*,*) ' mint1 = ', mint1
    !write(*,*) ' mint21 = ', mint21, ' mint22 = ', mint22
    !write(*,*) ' mint3 = ', mint3

    return
    end subroutine

    !-------------------------------------------------
    subroutine particle_monitor(this)
    !-------------------------------------------------
    implicit none
    class(particles) this

    !          write(*,*)
    !          write(*,*) '**** Information for particle ****', 
    !     &        		moni_particle         
    !          write(*,101)'internal a ','artifical a=',
    !     &         		'external a ','total a '   
    !          write(*,100)indvxdt%rr(1,moni_particle),
    !     &                ardvxdt%rr(1,moni_particle),
    !     &                exdvxdt%rr(1,moni_particle),
    !     &                parts%dvx(1,moni_particle)          
    !      endif
    !101   format(1x,4(2x,a12))      
    !100   format(1x,4(2x,e13.6))      

    return
    end subroutine

    !-------------------------------------------------
    subroutine write_scalar_field(field,parts)
    !-------------------------------------------------
    implicit none
    double precision, dimension(:) :: field
    type(particles) parts
    integer i,ntotal

    ntotal = parts%ntotal+parts%nvirt

    write(*,*) (i, field(i),i=1,ntotal)

    return
    end subroutine

    !-------------------------------------------------
    subroutine write_vector2d_field(field,parts)
    !-------------------------------------------------
    implicit none
    double precision, dimension(:,:) :: field
    type(particles) parts
    integer i,ntotal

    ntotal = parts%ntotal+parts%nvirt

    write(*,*) (i, field(1,i),field(2,i),i=1,ntotal)

    return
    end subroutine

    !---------------------------------------------------
    subroutine append_virtual_particles(this)
    !---------------------------------------------------
    implicit none
    class(particles) this

    type(particles), pointer :: bor
    integer ntotal, i, k

    bor    =>this%bor

    ! First, assembly all the real particles

    !ntotal = this%ntotal
    !k = 0 
    !do i = 1, bor%ntotal
    !   if(bor%itype(i)<0) cycle   ! real particles
    !   k = k + 1
    !   this%x(:,ntotal+k)  = bor%x(:,i)
    !   this%hsml(ntotal+k) = bor%hsml(i)
    !   this%vol(ntotal+k)  = bor%vol(i)
    !   this%itype(ntotal+k)= bor%itype(i)
    !   this%mass(ntotal+k) = bor%mass(i)
    !enddo
    !this%ntotal = this%ntotal + k

    ! Then assembly all the virtual particles

    ntotal = this%ntotal
    k = 0
    do i = 1, bor%ntotal
        if(bor%itype(i)>0) cycle  ! virtual particles
        k = k + 1
        this%x(:,ntotal+k)  = bor%x(:,i)
        this%hsml(ntotal+k) = bor%hsml(i)
        this%itype(ntotal+k)= bor%itype(i)
        this%vol%r(ntotal+k)  = bor%vol%r(i)
        this%mass%r(ntotal+k) = bor%mass%r(i)
    enddo
    this%nvirt = k
    !write(*,*) 'Appended virtual particle number: ', this%nvirt

    return
    end subroutine

    !---------------------------------------------------
    subroutine depend_virtual_particles(this)
    !---------------------------------------------------
    implicit none
    class(particles) this

    !type(particles), pointer :: bor
    !integer ntotal,nvirt, i,j,k

    !ntotal = this%ntotal
    !nvirt  = this%nvirt

    !bor => this%bor

    !k = 0
    !do i = ntotal+1, ntotal+nvirt

    !!   if(this%zone(i)==?) cycle      ! Symetery boundary virtural particles is skipped! 

    !   j = i - ntotal
    !   k = k + 1
    !   bor%x(:,j)  = this%x(:,i)
    !   bor%hsml(j) = this%hsml(i)
    !   bor%vol(j)  = this%vol(i)
    !   bor%zone(j) = this%zone(i)
    !enddo
    !this%nbor = k 

    return
    end subroutine


#ifdef REMOVE
    !!DEC$IF(.FALSE.)
    !-----------------------------------------------
    subroutine take_real_points(this,tank)
    !-----------------------------------------------
    implicit none

    class(particles) this
    type(block) tank

    integer i,j,k

    ! Take real particles in tank

    k = this%ntotal
    do j = 1, this%nreal_zone
        do i = 1, tank%m*tank%n
            if(tank%zone(i)==this%real_zone(j))then             
                k = k + 1
                this%x(1,k) = tank%x(i)
                this%x(2,k) = tank%y(i)
                this%zone(k) = tank%zone(i)
            endif
        enddo
    enddo
    this%ntotal = k

    return
    end subroutine

    !!DEC$ENDIF
#endif

    !--------------------------------------------------
    subroutine take_real_points1(this,tank,zone)
    !--------------------------------------------------
    implicit none

    class(particles) this
    type(block) tank
    integer zone

    integer i,j,k

    ! Take real particles in tank

    k = this%ntotal

    !   do i = 1, tank%m*tank%n
    do i = 1, tank%np
        if(tank%zone(i)== zone)then             
            k = k + 1
            this%x(1,k) = tank%x(i)
            this%x(2,k) = tank%y(i)
            this%zone(k) = tank%zone(i)
        endif
    enddo

    this%ntotal = k

    return
    end subroutine

#ifdef REMOVE
    !!DEC$IF(.FALSE.)
    !-----------------------------------------------------
    subroutine take_virtual_points(this,tank)
    !-----------------------------------------------------
    implicit none

    class(particles) this
    type(block) tank

    integer i,j,k

    ! Take virtual particles in tank

    k = this%ntotal
    do j = 1, this%nvirtual_zone
        do i = 1, tank%m*tank%n
            if(tank%zone(i)==this%virtual_zone(j))then             
                k = k + 1
                this%x(1,k) = tank%x(i)
                this%x(2,k) = tank%y(i)
                this%zone(k) = tank%zone(i)
            endif
        enddo
    enddo
    this%nvirt = k-this%ntotal

    return
    end subroutine

    !!DEC$ENDIF
#endif

    !-----------------------------------------------------
    subroutine take_virtual_points1(this,tank,zone)
    !-----------------------------------------------------
    implicit none

    class(particles) this
    type(block) tank
    integer zone

    integer i,j,k

    ! Take virtual particles in tank

    k = this%ntotal+this%nvirt
    !   do i = 1, tank%m*tank%n  
    do i = 1, tank%np  
        if(tank%zone(i)== zone)then             
            k = k + 1
            this%x(1,k) = tank%x(i)
            this%x(2,k) = tank%y(i)
            this%zone(k) = tank%zone(i)
        endif
    enddo
    this%nvirt = k-this%ntotal

    return
    end subroutine


    !-----------------------------------------------------
    subroutine take_boundary_nodes(this,tank,bzone)
    !-----------------------------------------------------
    implicit none

    class(particles) this
    type(block) tank
    integer bzone

    integer i,j,k

    ! Take virtual particles in tank

    k = this%ntotal+this%nvirt
    !   do i = 1, tank%m*tank%n  
    do i = 1, 2*(tank%m+tank%n)  
        if(tank%bzone(i)== bzone)then             
            k = k + 1
            this%x(1,k) = tank%bx(i)
            this%x(2,k) = tank%by(i)
            this%zone(k) = tank%bzone(i)
            this%bn%x%r(k) = tank%bn(1,i)
            this%bn%y%r(k) = tank%bn(2,i)
            this%bs%r(k) = tank%bs(i)
        endif
    enddo
    this%nvirt = k-this%ntotal

    return
    end subroutine


    !----------------------------------------
    subroutine setup_itype(parts)
    !----------------------------------------
    implicit none
    class(particles) parts
    integer ntotal, bntotal

    ntotal = parts%ntotal; bntotal = parts%nvirt

    if(trim(parts%imaterial)=='water')then
        parts%itype(1:ntotal) = 2
        parts%itype(ntotal+1:ntotal+bntotal) = -2
    elseif(trim(parts%imaterial)=='soil')then
        parts%itype(1:ntotal) = 3
        parts%itype(ntotal+1:ntotal+bntotal) = -3
    endif

    return
    end subroutine

    !----------------------------------------
    subroutine setup_ndim1(parts)
    !----------------------------------------
    implicit none
    class(particles) parts
    integer ntotal, bntotal

    ntotal = parts%ntotal + parts%nvirt

    if(associated(parts%rho))parts%rho%ndim1 = ntotal
    if(associated(parts%drho))parts%drho%ndim1 = ntotal
    if(associated(parts%mass))parts%mass%ndim1 = ntotal
    if(associated(parts%p))parts%p%ndim1 = ntotal
    !if(associated(parts%wxy))parts%wxy%ndim1 = ntotal
    if(associated(parts%vof))parts%vof%ndim1 = ntotal
    if(associated(parts%vof2))parts%vof2%ndim1 = ntotal
    if(associated(parts%dvof))parts%dvof%ndim1 = ntotal
    if(associated(parts%vof_min))parts%vof_min%ndim1 = ntotal
    if(associated(parts%str%x))parts%str%x%ndim1 = ntotal
    if(associated(parts%str%xy))parts%str%xy%ndim1 = ntotal
    if(associated(parts%str%y))parts%str%y%ndim1 = ntotal
    if(associated(parts%dstr))then
        if(associated(parts%dstr%x))parts%dstr%x%ndim1 = ntotal
        if(associated(parts%dstr%xy))parts%dstr%xy%ndim1 = ntotal
        if(associated(parts%dstr%y))parts%dstr%y%ndim1 = ntotal
    endif
    if(associated(parts%dstr2))then
        if(associated(parts%dstr2%x))parts%dstr2%x%ndim1 = ntotal
        if(associated(parts%dstr2%xy))parts%dstr2%xy%ndim1 = ntotal
        if(associated(parts%dstr2%y))parts%dstr2%y%ndim1 = ntotal
    endif
    !!!!!!!!
    if(associated(parts%tab%x))parts%tab%x%ndim1 = ntotal
    if(associated(parts%tab%xy))parts%tab%xy%ndim1 = ntotal
    if(associated(parts%tab%y))parts%tab%y%ndim1 = ntotal
    !!!!!!!!
    if(associated(parts%vcc))parts%vcc%ndim1 = ntotal
    if(associated(parts%v_min%x))parts%v_min%x%ndim1 = ntotal
    if(associated(parts%v_min%y))parts%v_min%y%ndim1 = ntotal
    if(associated(parts%epsilon_p))parts%epsilon_p%ndim1 = ntotal
    if(associated(parts%vx%x))parts%vx%x%ndim1 = ntotal
    if(associated(parts%vx%y))parts%vx%y%ndim1 = ntotal
    if(associated(parts%dvx%x))parts%dvx%x%ndim1 = ntotal
    if(associated(parts%dvx%y))parts%dvx%y%ndim1 = ntotal
    if(associated(parts%vol))parts%vol%ndim1 = ntotal
    if(associated(parts%rho_min))parts%rho_min%ndim1 = ntotal
    !if(associated(parts%u_min))parts%u_min%ndim1 = ntotal
    if(associated(parts%p_min))parts%p_min%ndim1 = ntotal
    if(associated(parts%str_min))then
        if(associated(parts%str_min%x))parts%str_min%x%ndim1 = ntotal
        if(associated(parts%str_min%y))parts%str_min%y%ndim1 = ntotal
        if(associated(parts%str_min%xy))parts%str_min%xy%ndim1 = ntotal
    endif
    if(associated(parts%dp))parts%dp%ndim1 = ntotal
    if(associated(parts%av%x))parts%av%x%ndim1 = ntotal
    if(associated(parts%av%y))parts%av%y%ndim1 = ntotal
    if(associated(parts%dp2))parts%dp2%ndim1 = ntotal
    !if(associated(parts%drhodx))parts%drhodx%ndim1 = ntotal
    if(associated(parts%c))parts%c%ndim1 = ntotal

    return
    end subroutine

    !-------------------------------------------------------------
    function get_total_particles_number(this) result(ntotal)
    !-------------------------------------------------------------
    implicit none
    class(particles) this
    integer ntotal
    ntotal = this%ntotal + this%nvirt
    end function

    !-----------------------------------------------------------        
    function get_maximum_particles_number(this) result(maxn)
    !-----------------------------------------------------------           
    implicit none
    class(particles) this
    integer maxn
    maxn = this%maxn
    end function
    !-----------------------------------------------------------        
    function get_dimension(this) result(dim)
    !-----------------------------------------------------------           
    implicit none
    class(particles) this
    integer dim
    dim = this%dim
    end function   

    !-------------------------------------------------
    function get_scale_k(this) result(scale_k)
    !-------------------------------------------------
    implicit none
    class(particles) this
    integer scale_k

    select case (this%skf)
    case(1)
        scale_k = 2
    case(2)
        scale_k = 3
    case(3)
        scale_k = 3
    case(4)
        scale_k = 2
        case default
        stop 'get_scale_k: No such kernel function!'
    end select 

    end function

    !----------------------------------------------------------------------
    subroutine direct_find(parts)
    !----------------------------------------------------------------------
    !   Subroutine to calculate the smoothing funciton for each particle and
    !   the interaction parameters used by the SPH algorithm. Interaction 
    !   pairs are determined by directly comparing the particle distance 
    !   with the corresponding smoothing length.
    !-----------------------------------------------------------------------
    implicit none

    class(particles) parts

    integer ntotal, niac,i,j,d,scale_k
    real(dp) dxiac(3), driac, r, mhsml, tdwdx(3)
    logical :: dbg = .false.

    if(dbg) write(*,*) 'In direct_find...'

    ntotal   =   parts%ntotal + parts%nvirt
    scale_k = parts%get_scale_k()

    do i=1,ntotal
        parts%countiac(i) = 0
    enddo

    niac = 0

    do i=1,ntotal-1     
        do j = i+1, ntotal
            dxiac(1) = parts%x(1,i) - parts%x(1,j)
            driac    = dxiac(1)*dxiac(1)
            do d=2,parts%dim
                dxiac(d) = parts%x(d,i) - parts%x(d,j)
                driac    = driac + dxiac(d)*dxiac(d)
            enddo
            mhsml = (parts%hsml(i)+parts%hsml(j))/2.
            if (sqrt(driac).lt.scale_k*mhsml) then
                if (niac.lt.parts%max_interaction) then    

                !     Neighboring pair list, and totalinteraction number and
                !     the interaction number for each particle 

                niac = niac + 1
                parts%pair_i(niac) = i
                parts%pair_j(niac) = j
                r = sqrt(driac)
                parts%countiac(i) = parts%countiac(i) + 1
                parts%countiac(j) = parts%countiac(j) + 1

                !     Kernel and derivations of kernel
                call parts%kernel(r,dxiac,mhsml,parts%w(niac),tdwdx)
                do d=1,parts%dim
                    parts%dwdx(d,niac) = tdwdx(d)
                enddo                                      
                else
                    print *,  ' >>> ERROR <<< : Too many interactions' 
                    stop
                endif
            endif
        enddo
    enddo  

    parts%niac = niac

    end subroutine

    !----------------------------------------------------------------------
    subroutine direct_find_2(parts,part2)
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts, part2
    integer ntotal, niac,i,j,d, scale_k
    real(dp) dxiac(3), driac, r, mhsml, tdwdx(3)
    logical :: dbg = .false.

    if(dbg) write(*,*) 'In direct_find2...'

    ntotal   =   parts%ntotal + parts%nvirt
    scale_k  =   parts%get_scale_k()

    do i=1,ntotal
        parts%countiac(i) = 0
    enddo
    do i = 1, part2%ntotal+part2%nvirt
        part2%countiac(i) = 0
    enddo

    niac = 0

    do i=1,ntotal     
        do j = 1, part2%ntotal+part2%nvirt
            dxiac(1) = parts%x(1,i) - part2%x(1,j)
            driac    = dxiac(1)*dxiac(1)
            do d=2,parts%dim
                dxiac(d) = parts%x(d,i) - part2%x(d,j)
                driac    = driac + dxiac(d)*dxiac(d)
            enddo
            mhsml = (parts%hsml(i)+part2%hsml(j))/2.
            if (sqrt(driac).lt.scale_k*mhsml) then
                if (niac.lt.parts%max_interaction) then    

                !     Neighboring pair list, and totalinteraction number and
                !     the interaction number for each particle 

                niac = niac + 1
                parts%pair_i(niac) = i
                parts%pair_j(niac) = j
                r = sqrt(driac)
                parts%countiac(i) = parts%countiac(i) + 1
                part2%countiac(j) = part2%countiac(j) + 1

                !     Kernel and derivations of kernel
                call parts%kernel(r,dxiac,mhsml,parts%w(niac),tdwdx)
                do d=1,parts%dim
                    parts%dwdx(d,niac) = tdwdx(d)
                enddo                                     
                else
                    print *, ' >>> ERROR <<< : Too many interactions' 
                    stop
                endif
            endif
        enddo
    enddo  

    parts%niac = niac

    end subroutine

    !   Subroutine to established a pair linked list by sorting grid cell.
    !   It is suitable for a homogeneous particle distribution with the 
    !   same smoothing length in an instant. A fixed number of particles
    !   lie in each cell.       
    !----------------------------------------------------------------------      
    subroutine init_grid(parts)
    !----------------------------------------------------------------------      
    implicit none
    class(particles) parts
    real(dp) hsml 
    integer i, j, k, d, dim, ntotal
    !     Averaged number of particles per grid cell
    real(dp), parameter :: nppg = 3.e0

    hsml = parts%hsml(1)   !!! hsml is same for all particles
    dim = parts%dim; ntotal = parts%ntotal + parts%nvirt

    !     Range of sorting grid

    parts%maxgridx(1) = parts%x_maxgeom;  parts%mingridx(1) = parts%x_mingeom
    parts%maxgridx(2) = parts%y_maxgeom;  parts%mingridx(2) = parts%y_mingeom
    parts%maxgridx(3) = parts%z_maxgeom;  parts%mingridx(3) = parts%z_mingeom
    parts%dgeomx = parts%maxgridx - parts%mingridx

    !     Number of grid cells in x-, y- and z-direction:

    if (dim.eq.1) then
        parts%ngridx(1) = min(int(ntotal/nppg) + 1,parts%maxngx)
    else if (dim.eq.2) then
        parts%ngridx(1) = &
        min(int(sqrt(ntotal*parts%dgeomx(1)/(parts%dgeomx(2)*nppg))) + 1,parts%maxngx)
        parts%ngridx(2) = &
        min(int(parts%ngridx(1)*parts%dgeomx(2)/parts%dgeomx(1)) + 1,parts%maxngy)
    else if (dim.eq.3) then
        parts%ngridx(1) = min(int((ntotal*parts%dgeomx(1)*parts%dgeomx(1)/   &
        (parts%dgeomx(2)*parts%dgeomx(3)*nppg))**(1.e0/3.e0)) + 1,parts%maxngx)
        parts%ngridx(2) =  &
        min(int(parts%ngridx(1)*parts%dgeomx(2)/parts%dgeomx(1)) + 1,parts%maxngy)
        parts%ngridx(3) =  &
        min(int(parts%ngridx(1)*parts%dgeomx(3)/parts%dgeomx(1)) + 1,parts%maxngz)
    endif

    !      parts%ngridx(1)=26; parts%ngridx(2)=13; parts%ngridx(3) = 1   !!!! VOF

    !     Smoothing Length measured in grid cells:
    do d=1,dim
        parts%ghsmlx(d) = int(real(parts%ngridx(d))*hsml/parts%dgeomx(d)) + 1
    enddo

    parts%grid = 0   ! Initialize grid  !! Better to move it out!

    end subroutine

    !   Subroutine to calculate the coordinates (xgcell) of the cell of 
    !   the sorting  grid, in which the particle with coordinates (x) lies.      
    !-----------------------------------------------------------------------      
    subroutine grid_geom(parts,i,x,xgcell)
    !-----------------------------------------------------------------------
    implicit none
    class(particles) parts
    integer i, xgcell(3)
    real(dp) x(3)
    integer d

    do d=1,3
        xgcell(d) = 1
    enddo

    do d=1,parts%dim
        if ((x(d).gt.parts%maxgridx(d)).or.(x(d).lt.parts%mingridx(d))) then
            print *,' >>> ERROR <<< : Particle out of range'
            print *,'    Particle position: x(',i,d,') = ',x(d)
            print *,'    Range: [xmin,xmax](',D,') =                    &
            [',parts%mingridx(d),',',parts%maxgridx(d),']'
            stop
        else
            xgcell(d) = int(real(parts%ngridx(d))/parts%dgeomx(d)*      &
            (x(d)-parts%mingridx(d)) + 1.e0)
        endif
    enddo

    end subroutine

    !----------------------------------------------------------------------      
    subroutine link_list(parts)
    !----------------------------------------------------------------------
    !   Subroutine to calculate the smoothing funciton for each particle and
    !   the interaction parameters used by the SPH algorithm. Interaction 
    !   pairs are determined by using a sorting grid linked list  
    !----------------------------------------------------------------------
    !      use ifport
    implicit none

    class(particles) parts
    integer i, j, d, scale_k, ntotal, niac    
    integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
    integer dnxgcell(3),dpxgcell(3)
    real(dp) hsml,dr,r,dx(3),tdwdx(3),x(3)

    INTEGER nthreads,n_per_threads, niac_per_threads, it
    integer,allocatable,dimension(:) :: n_start, n_end
    integer last
    real(dp) t1,t2,t3

    ntotal  = parts%ntotal + parts%nvirt
    scale_k = parts%get_scale_k()
    hsml = parts%hsml(1)             !!! hsml is same for all particles!

    do i=1,ntotal
        parts%countiac(i) = 0
    enddo

    call parts%init_grid

    !     Position particles on grid and create linked list:

    do i=1,ntotal
        do d = 1, parts%dim
            x(d) = parts%x(d,i)
        enddo
        call parts%grid_geom(i,x,gcell)
        do d=1,parts%dim
            parts%xgcell(d,i) = gcell(d)
        enddo
        parts%celldata(i) = parts%grid(gcell(1),gcell(2),gcell(3))
        parts%grid(gcell(1),gcell(2),gcell(3)) = i
    enddo

    !     Determine interaction parameters:
    !!!$omp parallel
    !!      nthreads = omp_get_num_threads()
    !!!$omp end parallel
    nthreads = parts%nthreads
    allocate(n_start(nthreads),n_end(nthreads))

    n_per_threads = ntotal/nthreads
    niac_per_threads = parts%max_interaction/nthreads
    do it = 1, nthreads
        n_start(it)=(it-1)*n_per_threads+1
        n_end(it) = it*n_per_threads
        parts%niac_start(it)=(it-1)*niac_per_threads
        parts%niac_end(it)=parts%niac_start(it)
    enddo   
    n_end(nthreads) = n_end(nthreads) + mod(ntotal,nthreads)

    do i =1, parts%max_interaction
        parts%pair_i(i) = 0; parts%pair_j(i)=0
    enddo
    !      t1 = rtc()
    !$omp parallel 
    !$omp do private(niac,i,d,minxcell,maxxcell,dnxgcell,dpxgcell,xcell,ycell,zcell,j,dx,dr,r,tdwdx)
    do it = 1, nthreads
        niac = parts%niac_start(it)
        do i = n_start(it),n_end(it)

        !      niac = 0
        !      do i=1,ntotal-1

        !     Determine range of grid to go through:

        do d=1,3
            minxcell(d) = 1
            maxxcell(d) = 1
        enddo
        do d=1,parts%dim
            dnxgcell(d) = parts%xgcell(d,i) - parts%ghsmlx(d)
            dpxgcell(d) = parts%xgcell(d,i) + parts%ghsmlx(d)
            minxcell(d) = max(dnxgcell(d),1)
            maxxcell(d) = min(dpxgcell(d),parts%ngridx(d))
        enddo

        !     Search grid:

        do zcell=minxcell(3),maxxcell(3)
            do ycell=minxcell(2),maxxcell(2)
                do xcell=minxcell(1),maxxcell(1)
                    j = parts%grid(xcell,ycell,zcell)
1                   if (j.gt.i) then
                        dx(1) = parts%x(1,i) - parts%x(1,j)
                        dr    = dx(1)*dx(1)
                        do d=2,parts%dim
                            dx(d) = parts%x(d,i) - parts%x(d,j)
                            dr    = dr + dx(d)*dx(d)
                        enddo
                        if (sqrt(dr).lt.scale_k*hsml) then
                            !if (niac.lt.parts%max_interaction) then
                            if (niac.lt.niac_per_threads*it) then


                            !   if(parts%x(2,i)>0.3.and.(parts%zone(i)==6.or.parts%zone(i)==100))then           
                            !   elseif(parts%x(2,j)>0.3.and.(parts%zone(j)==6.or.parts%zone(j)==100))then
                            !   else     

                            !     Neighboring pair list, and totalinteraction number and
                            !     the interaction number for each particle 

                            niac = niac + 1
                            parts%pair_i(niac) = i
                            parts%pair_j(niac) = j
                            r = sqrt(dr)
                            !parts%countiac(i) = parts%countiac(i) + 1
                            !parts%countiac(j) = parts%countiac(j) + 1

                            !--- Kernel and derivations of kernel

                            call parts%kernel(r,dx,hsml,parts%w(niac),tdwdx)  !!!!!!!!
                            do d = 1, parts%dim
                                parts%dwdx(d,niac)=tdwdx(d)
                            enddo                  

                            !   endif

                            else
                                print *, ' >>> Error <<< : too many interactions', &
                                it,niac,niac_per_threads
                                stop
                            endif
                        endif
                        j = parts%celldata(j)
                        goto 1
                    endif !(j.gt.i)
                enddo !xcell
            enddo !ycell
        enddo !zcell
        enddo !i

        parts%niac_start(it) = parts%niac_start(it)+1
        parts%niac_end(it) = niac
    enddo !it      
    !$omp end do
    !$omp barrier
    !$omp end parallel
    !      t2 = rtc()
    !      write(*,*) t2-t1
    parts%niac = 0
    do it = 1, nthreads
        parts%niac = parts%niac + parts%niac_end(it)-parts%niac_start(it)+1
    enddo

    last = parts%max_interaction
    do while(parts%pair_i(last)==0)
        last = last - 1
    enddo  
    i = 1 
    do while(i<last)
        if(parts%pair_i(i)/=0)then
            i = i + 1
            cycle
        endif        
        parts%pair_i(i) = parts%pair_i(last)
        parts%pair_j(i) = parts%pair_j(last)
        parts%w(i)      = parts%w(last)
        parts%dwdx(1,i) = parts%dwdx(1,last)
        parts%dwdx(2,i) = parts%dwdx(2,last)
        last = last - 1
        do while(parts%pair_i(last)==0)
            last = last - 1
        enddo
        i = i + 1
    enddo   
    !      t3 = rtc()
    !      write(*,*) t1,t2,t3
    !      write(*,*) t2-t1,t3-t2
    !      pause
    !      if(last/=parts%niac) stop 'adfasdf'
    !      do i = 1, last
    !         if(parts%pair_i(i)==0) stop 'adf'
    !      enddo   

    end subroutine

    !--------------------------------------------
    subroutine link_list2(parts,part2)
    !--------------------------------------------
    !      use ifport
    implicit none

    class(particles) parts,part2
    integer i, j, d, scale_k, ntotal, niac    
    integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
    integer dnxgcell(3),dpxgcell(3)
    real(dp) hsml,dr,r,dx(3),tdwdx(3),x(3)

    INTEGER nthreads,n_per_threads, niac_per_threads, it
    integer,allocatable,dimension(:) :: n_start, n_end
    integer last
    real(dp) t1,t2,t3

    ntotal  = parts%ntotal + parts%nvirt
    scale_k = parts%get_scale_k()
    hsml = parts%hsml(1)

    do i = 1,ntotal
        parts%countiac(i) = 0
    enddo
    do i = 1, part2%ntotal+part2%nvirt
        part2%countiac(i) = 0
    enddo

    call parts%init_grid
    parts%grid = 0; parts%xgcell = 0; parts%celldata = 0  !?? grid already
    ! initialized in
    ! init_grid

    !     Position particles on grid and create linked list:

    do i=1,ntotal
        do d = 1, parts%dim
            x(d) = parts%x(d,i)
        enddo
        call parts%grid_geom(i,x,gcell)
        do d=1,parts%dim
            parts%xgcell(d,i) = gcell(d)
        enddo
        parts%celldata(i) = parts%grid(gcell(1),gcell(2),gcell(3))
        parts%grid(gcell(1),gcell(2),gcell(3)) = i
    enddo

    !call part2%init_grid

    part2%maxgridx(1) = part2%x_maxgeom;  part2%mingridx(1) = part2%x_mingeom
    part2%maxgridx(2) = part2%y_maxgeom;  part2%mingridx(2) = part2%y_mingeom
    part2%maxgridx(3) = part2%z_maxgeom;  part2%mingridx(3) = part2%z_mingeom
    part2%dgeomx = part2%maxgridx - part2%mingridx      

    part2%ngridx = parts%ngridx; part2%ghsmlx = parts%ghsmlx   !Same grid
    part2%grid = 0; part2%xgcell = 0; part2%celldata = 0


    !     Position particles on grid and create linked list:

    do i=1,part2%ntotal+part2%nvirt
        do d = 1, part2%dim
            x(d) = part2%x(d,i)
        enddo      
        call part2%grid_geom(i,x,gcell)
        do d=1,part2%dim
            part2%xgcell(d,i) = gcell(d)
        enddo
        part2%celldata(i) = part2%grid(gcell(1),gcell(2),gcell(3))
        part2%grid(gcell(1),gcell(2),gcell(3)) = i
    enddo

    nthreads = parts%nthreads
    allocate(n_start(nthreads),n_end(nthreads))

    n_per_threads = ntotal/nthreads
    niac_per_threads = parts%max_interaction/nthreads
    do it = 1, nthreads
        n_start(it)=(it-1)*n_per_threads+1
        n_end(it) = it*n_per_threads
        parts%niac_start(it)=(it-1)*niac_per_threads
        parts%niac_end(it)=parts%niac_start(it)
    enddo   
    n_end(nthreads) = n_end(nthreads) + mod(ntotal,nthreads)

    do i =1, parts%max_interaction
        parts%pair_i(i) = 0; parts%pair_j(i) = 0
    enddo

    !      t1 = rtc()

    !!$omp parallel 
    !!$omp do private(niac,i,d,minxcell,maxxcell,dnxgcell,dpxgcell,xcell,ycell,zcell,j,dx,dr,r,tdwdx)
    !do it = 1, nthreads
    !       niac = parts%niac_start(it)
    !       do i = n_start(it),n_end(it)

    niac = 0                       ! No parallel
    do i=1,ntotal   !ntotal-1

    !     Determine range of grid to go through:

    do d=1,3
        minxcell(d) = 1
        maxxcell(d) = 1
    enddo
    do d=1,parts%dim
        dnxgcell(d) = parts%xgcell(d,i) - parts%ghsmlx(d)
        dpxgcell(d) = parts%xgcell(d,i) + parts%ghsmlx(d)
        minxcell(d) = max(dnxgcell(d),1)
        maxxcell(d) = min(dpxgcell(d),parts%ngridx(d))
    enddo

    !     Search grid:

    do zcell=minxcell(3),maxxcell(3)
        do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
                j = part2%grid(xcell,ycell,zcell)
                ! 1            if (j.gt.i) then
1               if (j/=0) then
                    !dx(1) = parts%x(1,i) - part2%x(1,j)
                    !dr    = dx(1)*dx(1)
                    dr = 0.d0
                    do d=1,parts%dim
                        dx(d) = parts%x(d,i) - part2%x(d,j)
                        dr    = dr + dx(d)*dx(d)
                    enddo
                    if (sqrt(dr).lt.scale_k*hsml) then
                        if (niac.lt.parts%max_interaction) then
                            !if (niac.lt.niac_per_threads) then

                            !     Neighboring pair list, and totalinteraction number and
                            !     the interaction number for each particle 

                            niac = niac + 1
                            parts%pair_i(niac) = i
                            parts%pair_j(niac) = j
                            r = sqrt(dr)
                            !parts%countiac(i) = parts%countiac(i) + 1
                            !parts%countiac(j) = parts%countiac(j) + 1

                            !--- Kernel and derivations of kernel     

                            call parts%kernel(r,dx,hsml,parts%w(niac),tdwdx)  !!!!!!!!
                            do d = 1, parts%dim
                                parts%dwdx(d,niac)=tdwdx(d)
                            enddo                  
                        else
                            print *, ' >>> Error <<< : too many interactions', &
                            it,niac,niac_per_threads
                            stop
                        endif
                    endif
                    j = part2%celldata(j)
                    goto 1
                endif
            enddo !xcell
        enddo !ycell
    enddo !zcell
    enddo !i

    parts%niac = niac

    parts%niac_start(it) = parts%niac_start(it)+1
    parts%niac_end(it) = niac
    !      enddo !it      
    !!$omp end do
    !!$omp end parallel
    return

    !      t2 = rtc()
    !      write(*,*) t2-t1
    parts%niac = 0
    do it = 1, nthreads
        parts%niac = parts%niac + parts%niac_end(it)-parts%niac_start(it)+1
    enddo

    last = parts%max_interaction
    do while(parts%pair_i(last)==0)
        last = last - 1
    enddo  
    i = 1 
    do while(i<last)
        if(parts%pair_i(i)/=0)then
            i = i + 1
            cycle
        endif        
        parts%pair_i(i) = parts%pair_i(last)
        parts%pair_j(i) = parts%pair_j(last)
        parts%w(i)      = parts%w(last)
        parts%dwdx(1,i) = parts%dwdx(1,last)
        parts%dwdx(2,i) = parts%dwdx(2,last)
        last = last - 1
        do while(parts%pair_i(last)==0)
            last = last - 1
        enddo
        i = i + 1
    enddo

    ! To check

    if(last/=parts%niac) stop 'adfasdf'
    do i = 1, last
        if(parts%pair_i(i)==0) stop 'adf'
    enddo   
    end subroutine

    !--------------------------------------------
    subroutine link_list2_omp(parts,part2)
    !--------------------------------------------
    !      use ifport
    implicit none

    class(particles) parts,part2
    integer i, j, d, scale_k, ntotal, niac    
    integer gcell(3),xcell,ycell,zcell,minxcell(3),maxxcell(3)
    integer dnxgcell(3),dpxgcell(3)
    real(dp) hsml,dr,r,dx(3),tdwdx(3), x(3)

    INTEGER nthreads,n_per_threads, niac_per_threads, it
    integer,allocatable,dimension(:) :: n_start, n_end
    integer last
    real(dp) t1,t2,t3

    ntotal  = parts%ntotal + parts%nvirt
    scale_k = parts%get_scale_k()
    hsml = parts%hsml(1)

    do i = 1,ntotal
        parts%countiac(i) = 0
    enddo
    do i = 1, part2%ntotal+part2%nvirt
        part2%countiac(i) = 0
    enddo

    call parts%init_grid
    parts%grid = 0; parts%xgcell = 0; parts%celldata = 0   ! grid already
    ! initialized in
    ! init_grid

    !     Position particles on grid and create linked list:

    do i=1,ntotal
        do d = 1, parts%dim
            x(d) = parts%x(d,i)
        enddo            
        call parts%grid_geom(i,x,gcell)
        do d=1,parts%dim
            parts%xgcell(d,i) = gcell(d)
        enddo
        parts%celldata(i) = parts%grid(gcell(1),gcell(2),gcell(3))
        parts%grid(gcell(1),gcell(2),gcell(3)) = i
    enddo

    !call part2%init_grid
    part2%maxgridx(1) = part2%x_maxgeom;  part2%mingridx(1) = part2%x_mingeom
    part2%maxgridx(2) = part2%y_maxgeom;  part2%mingridx(2) = part2%y_mingeom
    part2%maxgridx(3) = part2%z_maxgeom;  part2%mingridx(3) = part2%z_mingeom
    part2%dgeomx = part2%maxgridx - part2%mingridx

    part2%ngridx = parts%ngridx; part2%ghsmlx = parts%ghsmlx
    part2%grid = 0; part2%xgcell = 0; part2%celldata = 0


    !     Position particles on grid and create linked list:

    do i=1,part2%ntotal+part2%nvirt
        do d = 1, part2%dim
            x(d) = part2%x(d,i)
        enddo            
        call part2%grid_geom(i,x,gcell)
        do d=1,part2%dim
            part2%xgcell(d,i) = gcell(d)
        enddo
        part2%celldata(i) = part2%grid(gcell(1),gcell(2),gcell(3))
        part2%grid(gcell(1),gcell(2),gcell(3)) = i
    enddo

    nthreads = parts%nthreads
    allocate(n_start(nthreads),n_end(nthreads))

    n_per_threads = ntotal/nthreads
    niac_per_threads = parts%max_interaction/nthreads
    do it = 1, nthreads
        n_start(it)=(it-1)*n_per_threads+1
        n_end(it) = it*n_per_threads
        parts%niac_start(it)=(it-1)*niac_per_threads
        parts%niac_end(it)=parts%niac_start(it)
    enddo   
    n_end(nthreads) = n_end(nthreads) + mod(ntotal,nthreads)

    do i = 1, parts%max_interaction
        parts%pair_i(i) = 0; parts%pair_j(i) = 0
    enddo

    !      t1 = rtc()

    !$omp parallel 
    !$omp do private(niac,i,d,minxcell,maxxcell,dnxgcell,dpxgcell,xcell,ycell,zcell,j,dx,dr,r,tdwdx)
    do it = 1, nthreads
        niac = parts%niac_start(it)
        do i = n_start(it),n_end(it)

        !      niac = 0
        !      do i=1,ntotal   !ntotal-1

        !     Determine range of grid to go through:

        do d=1,3
            minxcell(d) = 1
            maxxcell(d) = 1
        enddo
        do d=1,parts%dim
            dnxgcell(d) = parts%xgcell(d,i) - parts%ghsmlx(d)
            dpxgcell(d) = parts%xgcell(d,i) + parts%ghsmlx(d)
            minxcell(d) = max(dnxgcell(d),1)
            maxxcell(d) = min(dpxgcell(d),parts%ngridx(d))
        enddo

        !     Search grid:

        do zcell=minxcell(3),maxxcell(3)
            do ycell=minxcell(2),maxxcell(2)
                do xcell=minxcell(1),maxxcell(1)
                    j = part2%grid(xcell,ycell,zcell)
                    ! 1            if (j.gt.i) then
1                   if (j/=0) then
                        !dx(1) = parts%x(1,i) - part2%x(1,j)
                        !dr    = dx(1)*dx(1)
                        dr = 0.d0
                        do d=1,parts%dim
                            dx(d) = parts%x(d,i) - part2%x(d,j)
                            dr    = dr + dx(d)*dx(d)
                        enddo
                        if (sqrt(dr).lt.scale_k*hsml) then
                            !if (niac.lt.parts%max_interaction) then
                            if (niac.lt.niac_per_threads*it) then

                            !     Neighboring pair list, and totalinteraction number and
                            !     the interaction number for each particle 

                            niac = niac + 1
                            parts%pair_i(niac) = i
                            parts%pair_j(niac) = j
                            r = sqrt(dr)
                            !parts%countiac(i) = parts%countiac(i) + 1
                            !parts%countiac(j) = parts%countiac(j) + 1

                            !--- Kernel and derivations of kernel     

                            call parts%kernel(r,dx,hsml,parts%w(niac),tdwdx)  !!!!!!!!
                            do d = 1, parts%dim
                                parts%dwdx(d,niac)=tdwdx(d)
                            enddo                  
                            else
                                print *, ' >>> Error <<< : too many interactions', &
                                it,niac,niac_per_threads
                                stop
                            endif
                        endif
                        j = part2%celldata(j)
                        goto 1
                    endif
                enddo !xcell
            enddo !ycell
        enddo !zcell
        enddo !i

        !      parts%niac = niac

        parts%niac_start(it) = parts%niac_start(it)+1
        parts%niac_end(it) = niac
    enddo !it      
    !$omp end do
    !$omp barrier
    !$omp end parallel

    !      t2 = rtc()
    !      write(*,*) t2-t1
    parts%niac = 0
    do it = 1, nthreads
        parts%niac = parts%niac + parts%niac_end(it)-parts%niac_start(it)+1
    enddo

    last = parts%max_interaction
    do while(parts%pair_i(last)==0)
        last = last - 1
    enddo  
    i = 1 
    do while(i<last)
        if(parts%pair_i(i)/=0)then
            i = i + 1
            cycle
        endif        
        parts%pair_i(i) = parts%pair_i(last)
        parts%pair_j(i) = parts%pair_j(last)
        parts%w(i)      = parts%w(last)
        parts%dwdx(1,i) = parts%dwdx(1,last)
        parts%dwdx(2,i) = parts%dwdx(2,last)
        last = last - 1
        do while(parts%pair_i(last)==0)
            last = last - 1
        enddo
        i = i + 1
    enddo

    ! To check

    if(last/=parts%niac) stop 'adfasdf'
    do i = 1, last
        if(parts%pair_i(i)==0) stop 'adf'
    enddo   
    end subroutine

    !--------------------------------------------
    subroutine parts2cell(parts,fp,f)
    !--------------------------------------------
    implicit none
    class(particles) parts
    type(array) fp
    type(realarray) f
    integer i,j,k,zcell,ycell,xcell
    !write(*,*) 'c', f%r(8,2)
    !        do zcell=1,1   !parts%ngridx(3)
    do ycell=1,parts%ngridx(2)
        do xcell=1,parts%ngridx(1)
            j = parts%grid(xcell,ycell,1)
            f%r(xcell,ycell) = 0.d0
            k = 0
1           if (j.gt.0) then
                k = k + 1  
                f%r(xcell,ycell) = f%r(xcell,ycell) +  fp%r(j)
                j = parts%celldata(j)
                goto 1
            endif
            if(k>0)then
                f%r(xcell,ycell) = f%r(xcell,ycell)/k
            endif
        enddo !xcell
    enddo !ycell
    !        enddo !zcell
    !write(*,*) 'c', f%r(8,2)
    return
    end subroutine

    !--------------------------------------------
    subroutine cell2parts(parts,f,fp)
    !--------------------------------------------
    implicit none
    class(particles) parts
    type(realarray) f
    type(array) fp
    integer i,j,k,zcell,ycell,xcell
    !write(*,*) fp%r(1:10)
    do i = 1, parts%ntotal + parts%nvirt
        fp%r(i) = f%r(parts%xgcell(1,i),parts%xgcell(2,i))
    enddo
    !write(*,*) parts%xgcell(1:2,9),f%r(2,6)
    !write(*,*) fp%r(1:10)
    !write(*,*) ((f%r(i,j),i = 1,2),j=1,6)
    !write(*,*) 'ddd'
    return
    end subroutine      

    !--------------------------------------------
    subroutine find_pairs(parts)
    !--------------------------------------------
    implicit none
    class(particles) parts
    integer nnps

    nnps = parts%nnps
    if(nnps == 1)then
        call direct_find(parts)
    elseif(nnps == 2)then
        call link_list(parts)
        !   call link_list(parts%itimestep, parts%ntotal+parts%nvirt,parts%hsml(1),  &
        !   parts%x,parts%niac,parts%pair_i,parts%pair_j,parts%w,parts%dwdx,parts%countiac)
        !   stop 'find_pairs: link_list method not implemented yet!'
    elseif(nnps == 3)then
        !    call tree_search(parts)
        stop 'find_pairs: tree_search method not implemented yet!'
    endif

    return
    end subroutine

    !----------------------------------------------------------------------
    subroutine kernel(this,r,dx,hsml,w,dwdx)   
    !----------------------------------------------------------------------
    !   Subroutine to calculate the smoothing kernel wij and its 
    !   derivatives dwdxij.
    !     if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
    !            = 2, Gauss kernel   (Gingold and Monaghan 1981) 
    !            = 3, Quintic kernel (Morris 1997)
    !            = 4, Wendland kernel
    !----------------------------------------------------------------------
    implicit none
    class(particles) this
    real(dp) r, dx(3), hsml, w, dwdx(3)
    integer i, j, d, dim, skf      
    real(dp) q, dw, factor
    real(dp), parameter :: pi = 3.14159265358979323846 
    logical :: dbg = .false.

    if(dbg) write(*,*) 'In kernel...'

    dim = this%dim; skf = this%skf

    if(dabs(hsml)<1.0d-10) stop 'kernel: Too small hsml!'
    q = r/hsml 
    w = 0.e0
    do d=1,dim         
        dwdx(d) = 0.e0
    enddo   

    if (skf.eq.1) then     
        if (dim.eq.1) then
            factor = 1.e0/hsml
        elseif (dim.eq.2) then
            factor = 15.e0/(7.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
            factor = 3.e0/(2.e0*pi*hsml*hsml*hsml)
        else
            print *,' >>> Error <<< : Wrong dimension: Dim =',dim
            stop
        endif                                           
        if (q.ge.0.and.q.le.1.e0) then          
            w = factor * (2./3. - q*q + q**3 / 2.)
            do d = 1, dim
                dwdx(d) = factor * (-2.+3./2.*q)/hsml**2 * dx(d)       
            enddo   
        else if (q.gt.1.e0.and.q.le.2) then          
            w = factor * 1.e0/6.e0 * (2.-q)**3 
            do d = 1, dim
                dwdx(d) =-factor * 1.e0/6.e0 * 3.*(2.-q)**2/hsml * (dx(d)/r)        
            enddo              
        else
            w=0.
            do d= 1, dim
                dwdx(d) = 0.
            enddo             
        endif     

    else if (skf.eq.2) then

    factor = 1.e0 / (hsml**dim * pi**(dim/2.))      
    if(q.ge.0.and.q.le.3) then
        w = factor * exp(-q*q)
        do d = 1, dim
            dwdx(d) = w * ( -2.* dx(d)/hsml/hsml)
        enddo 
    else
        w = 0.
        do d = 1, dim
            dwdx(d) = 0.
        enddo 	   
    endif	       

    else if (skf.eq.3) then	

    if (dim.eq.1) then
        factor = 1.e0 / (120.e0*hsml)
    elseif (dim.eq.2) then
        factor = 7.e0 / (478.e0*pi*hsml*hsml)
    elseif (dim.eq.3) then
        factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml)
    else
        print *,' >>> Error <<< : Wrong dimension: Dim =',dim
        stop
    endif              
    if(q.ge.0.and.q.le.1) then
        w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
        do d= 1, dim
            dwdx(d) = factor * ( (-120 + 120*q - 50*q**2)      &
            / hsml**2 * dx(d) )
        enddo 
    else if(q.gt.1.and.q.le.2) then
        w = factor * ( (3-q)**5 - 6*(2-q)**5 )
        do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4 + 30*(2-q)**4)     &
            / hsml * (dx(d)/r) 
        enddo 
    else if(q.gt.2.and.q.le.3) then
        w = factor * (3-q)**5 
        do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4) / hsml * (dx(d)/r) 
        enddo 
    else   
        w = 0.
        do d = 1, dim
            dwdx(d) = 0.
        enddo  
    endif                      

    elseif(skf.eq.4)then    ! Wendland. refer to SPHysic manual

    if (dim.eq.1) then
        factor = 0.
    elseif (dim.eq.2) then
        factor = 7.e0 / (4.e0*pi*hsml*hsml)
    elseif (dim.eq.3) then
        factor = 0.
    else
        print *,' >>> Error <<< : Wrong dimension: Dim =',dim
        stop
    endif              
    if(q.ge.0.and.q.le.2) then
        w = factor * ( (1-q/2)**4 *(1+2*q) )
        do d= 1, dim
            dwdx(d) = factor*(-5+15*q/2-15*q**2/4+5*q**3/8)/    &
            hsml**2*dx(d)
        enddo 
    else   
        w = 0.
        do d = 1, dim
            dwdx(d) = 0.
        enddo  
    endif 

    endif 

    end subroutine

    !--------------------------------------------------------------------------
    subroutine calc_norm(parts)
    !--------------------------------------------------------------------------
    implicit none
    class(particles) parts
    real(dp) ds
    integer i,j,k,d

    do i = parts%ntotal+1, parts%ntotal + parts%nvirt-1 
        k = i-parts%ntotal
        ds = sqrt((parts%x(1,i)-parts%x(1,i+1))**2.+(parts%x(2,i)-parts%x(2,i+1))**2.)
        do d = 1, parts%dim
            parts%norm(1,k) = (parts%x(2,i+1)-parts%x(2,i))/ds
            parts%norm(2,k) = (parts%x(1,i)-parts%x(1,i+1))/ds
        enddo
    enddo

    return
    end subroutine

    !--------------------------------------------------------------------------
    subroutine repulsive_force(parts)
    !--------------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision dx(3), rr, f, rr0, dd, p1, p2,small     
    integer i, j, k, d, ii

    rr0 = parts%dspp; dd = parts%numeric%dd
    p1 = parts%numeric%p1; p2 = parts%numeric%p2
    small = 1.d-6

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)*parts%itype(j)>0)cycle    ! inlet boundary condition
        !   if(parts%itype(i).gt.0.and.parts%itype(j).lt.0)then  
        rr = 0.      
        do d=1,parts%dim
            dx(d) =  parts%x(d,i) -  parts%x(d,j)
            rr = rr + dx(d)*dx(d)
        enddo  
        rr = sqrt(rr) + small  !!!!!
        !if(rr.lt.rr0)then
        if(rr.gt.rr0)cycle
        f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2

        ii = i
        if(parts%itype(i)<0)then
            ii=j; f = -f       ! direction
        endif        

        !do d = 1, parts%dim
        !   parts%dvx(d, ii) = parts%dvx(d, ii) + dd*dx(d)*f
        !enddo

        parts%dvx%x%r(ii) = parts%dvx%x%r(ii) + dd * dx(1)*f    ! direction?
        parts%dvx%y%r(ii) = parts%dvx%y%r(ii) + dd * dx(2)*f

        !endif
        !endif        
    enddo   

    return
    end subroutine

    !--------------------------------------------------------------------------
    subroutine repulsive_force_omp(parts)
    !--------------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision dx(3), rr, f, rr0, dd, p1, p2,small
    real(dp),allocatable,dimension(:,:,:) :: local
    integer i, j, k, d, ii,ntotal,nthreads,it,dim

    rr0 = parts%dspp; dd = parts%numeric%dd
    p1 = parts%numeric%p1; p2 = parts%numeric%p2
    ntotal = parts%ntotal+parts%nvirt; dim = parts%dim
    nthreads = parts%nthreads
    small = 1.d-6

    allocate(local(dim,ntotal,nthreads))
    call parts%get_niac_start_end

    !$omp parallel
    !$omp do private(i,j,rr,d,dx,f,ii,k) 
    do it = 1,nthreads
        do i = 1,ntotal
            local(1:dim,i,it) = 0.d0
        enddo
        do k=parts%niac_start(it),parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            if(parts%itype(i)*parts%itype(j)>0)cycle    ! inlet boundary condition
            !   if(parts%itype(i).gt.0.and.parts%itype(j).lt.0)then  
            rr = 0.      
            do d=1,dim
                dx(d) =  parts%x(d,i) -  parts%x(d,j)
                rr = rr + dx(d)*dx(d)
            enddo  
            rr = sqrt(rr) + small !!!!!!
            !if(rr.lt.rr0)then
            if(rr.gt.rr0)cycle
            f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2

            ii = i
            if(parts%itype(i)<0)then
                ii=j; f = -f         ! direction
            endif        

            do d = 1,dim
                !parts%dvx(d, ii) = parts%dvx(d, ii) + dd*dx(d)*f
                local(d,ii,it) = local(d,ii,it) + dd * dx(d)*f     ! direction?
            enddo
            !endif
            !endif        
        enddo  !k
    enddo  !it
    !$omp end do

    !$omp barrier
    !$omp do private(it)
    do i = 1,ntotal
        do it = 1,nthreads
            parts%dvx%x%r(i) = parts%dvx%x%r(i) + local(1,i,it)
            parts%dvx%y%r(i) = parts%dvx%y%r(i) + local(2,i,it)
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
    end subroutine


    ! Calculate partial derivatives of a field
    !--------------------------------------------
    function df(parts,f,x) result(res)
    !--------------------------------------------
    implicit none

    type(array) :: f
    character(len=1) x
    class(particles),target :: parts
    type(array) :: res
    real(dp), pointer, dimension(:) :: dwdx
    real(dp) fwx
    integer i, j, k, ntotal
    ntotal = parts%ntotal +parts%nvirt

    allocate(res%r(ntotal))
    res%parts => parts
    res%ndim1 = ntotal; res = 0.d0

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        fwx = (f%r(i)+f%r(j))*dwdx(k)
        res%r(i) = res%r(i) + parts%mass%r(j)/parts%rho%r(j)*fwx
        res%r(j) = res%r(j) - parts%mass%r(i)/parts%rho%r(i)*fwx
    enddo
    end function


    ! Calculate partial derivatives of a field
    !-----------------------------------------------------
    function df_omp(parts,f,x)  result(fun)
    !-----------------------------------------------------
    implicit none

    type(array) f
    character(len=1) x
    class(particles),target :: parts
    type(array) :: fun
    real(dp), allocatable, dimension(:,:) :: local
    real(dp), pointer, dimension(:) :: dwdx
    real(dp) fwx
    integer i, j, k, ntotal, it, nthreads

    !write(*,*) 'In df_omp...'

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads

    fun%parts => parts
    fun%ndim1 = ntotal
    allocate(fun%r(ntotal))
    if(nthreads>=1)then
        allocate(local(ntotal,nthreads))
        call parts%get_niac_start_end
    endif   

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)

    !$omp parallel
    !$omp do private(i,j,k,fwx) 
    do it = 1, nthreads
        do i = 1, ntotal
            local(i,it) = 0.d0
        enddo 
        do k = parts%niac_start(it), parts%niac_end(it)
            !do k=1,parts%niac
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            fwx = (f%r(i)+f%r(j))*dwdx(k)
            local(i,it) = local(i,it) + parts%mass%r(j)/parts%rho%r(j)*fwx
            local(j,it) = local(j,it) - parts%mass%r(i)/parts%rho%r(i)*fwx
        enddo !k
    enddo !it
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        fun%r(i) = 0.d0
        do it = 1, nthreads
            fun%r(i) = fun%r(i) + local(i,it)
        enddo
    enddo   
    !$omp end do
    !$omp end parallel
    end function

    ! Calculate partial derivatives of a field
    !----------------------------------------------
    function df_omp2(parts,f,x) result(val)
    !----------------------------------------------
    implicit none

    type(array) :: f
    character(len=1) x
    class(particles),target :: parts
    type(array) :: val
    !real(dp), allocatable, dimension(:,:) :: df_local
    real(dp), pointer, dimension(:) :: dwdx
    real(dp) fwx
    integer i, j, k, ntotal, it, nthreads, ii, jj

    !write(*,*) 'In df_omp...'

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads
    !write(*,*) 'sadf', nthreads

    val%parts => parts
    allocate(val%r(ntotal*nthreads))
    !if(nthreads>=1)then
    !   allocate(df_local(ntotal,nthreads))
    call parts%get_niac_start_end
    !endif   

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)

    !$omp parallel
    !$omp do private(i,j,k,fwx,ii,jj)
    do it = 1, parts%nthreads
        do i = 1, ntotal
            ii = i + (it-1)*ntotal
            val%r(ii) = 0.d0
        enddo
        do k = parts%niac_start(it), parts%niac_end(it)


        !do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        fwx = (f%r(i)+f%r(j))*dwdx(k)
        ii = i + (it-1)*ntotal; jj = j + (it-1)*ntotal
        val%r(ii) = val%r(ii) + parts%mass%r(j)/parts%rho%r(j)*fwx
        val%r(jj) = val%r(jj) - parts%mass%r(i)/parts%rho%r(i)*fwx   
        !   df_local(i,it) = df_local(i,it) + parts%mass(j)/parts%rho(j)*fwx
        !   df_local(j,it) = df_local(j,it) - parts%mass(i)/parts%rho(i)*fwx
        enddo !k
    enddo !it
    !$omp end do
    !$omp barrier

    !$omp do private(it,ii)
    do i = 1, ntotal
        !val(i) = 0.d0
        do it = 2, nthreads
            ii = i + (it-1)*ntotal
            val%r(i) = val%r(i)+val%r(ii)
        enddo
    enddo   
    !$omp end do
    !$omp end parallel

    val%ndim1 = ntotal !!!
    end function

    ! Calculate partial derivatives of a field
    !-------------------------------------------
    function df2(parts,f,x)
    !-------------------------------------------
    implicit none

    type(array) :: f
    character(len=1) x
    class(particles), target :: parts
    type(array) :: df2
    real(dp), pointer, dimension(:) :: dwdx
    real(dp) fwx
    integer i, j, k, ntotal
    ntotal = parts%ntotal +parts%nvirt

    df2%parts => parts
    allocate(df2%r(ntotal))
    df2%ndim1 = ntotal; df2 = 0.d0

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        fwx = (f%r(j)-f%r(i))*dwdx(k)
        df2%r(i) = df2%r(i) + parts%mass%r(j)*fwx
        df2%r(j) = df2%r(j) + parts%mass%r(i)*fwx
    enddo

    do i = 1, parts%ntotal + parts%nvirt 
        df2%r(i) = df2%r(i)/parts%rho%r(i)
    enddo   
    end function


    ! Calculate partial derivatives of a field
    !-------------------------------------------
    function df3(parts,f,x)
    !------------------------------------------- 
    implicit none

    type(array) :: f
    character(len=1) x
    class(particles),target :: parts
    type(array) :: df3
    real(dp), pointer, dimension(:) :: dwdx
    real(dp) fwx
    integer i,j,k,ntotal
    ntotal = parts%ntotal + parts%nvirt

    df3%parts => parts
    allocate(df3%r(ntotal))
    df3%ndim1 = ntotal; df3 = 0.d0

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        fwx = ((f%r(i)/parts%rho%r(i)**2)+(f%r(j)/parts%rho%r(j)**2))*dwdx(k)
        df3%r(i) = df3%r(i) + parts%mass%r(j)*fwx
        df3%r(j) = df3%r(j) - parts%mass%r(i)*fwx
    enddo

    do i = 1, parts%ntotal + parts%nvirt
        df3%r(i) = df3%r(i)*parts%rho%r(i)
    enddo
    end function

    !-------------------------------------------
    function df3_omp(parts,f,x)
    !------------------------------------------- 
    implicit none

    type(array) f
    character(len=1) x
    class(particles),target :: parts
    type(array) :: df3_omp
    real(dp), allocatable, dimension(:,:) :: local
    !real(dp) df3_local(parts%maxn,8)
    real(dp), pointer, dimension(:) :: dwdx, bn
    real(dp) fwx
    integer i,j,k,it,ntotal,nthreads

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads

    df3_omp%parts => parts
    df3_omp%ndim1 = ntotal
    allocate(df3_omp%r(ntotal))

    if(nthreads>=1)then
        allocate(local(ntotal,nthreads))
        call parts%get_niac_start_end
    endif 

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)
    if(parts%usaw)then
        if(x=='x')bn=>parts%bn%x%r
        if(x=='y')bn=>parts%bn%y%r
    endif

    !$omp parallel
    !$omp do private(i,j,k,fwx)
    do it = 1, parts%nthreads
        do i = 1, ntotal
            local(i,it) = 0.d0
        enddo

        do k = parts%niac_start(it), parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)

            ! USAW________________________________________________________________

            if(parts%usaw.and.parts%itype(i)>0.and.parts%itype(j)<0)then   
                fwx = (f%r(i)/parts%rho%r(i)**2)+(f%r(j)/parts%rho%r(j)**2)
                local(i,it) = local(i,it) - fwx*parts%w(k)*parts%bs%r(j)*bn(j)*parts%rho%r(j)
                cycle
            endif        
            if(parts%usaw.and.parts%itype(i)<0.and.parts%itype(j)<0)cycle   
            !---------------------------------------------------------------------

            fwx = ((f%r(i)/parts%rho%r(i)**2)+(f%r(j)/parts%rho%r(j)**2))*dwdx(k)
            local(i,it)= local(i,it) + parts%mass%r(j)*fwx
            local(j,it)= local(j,it) - parts%mass%r(i)*fwx
            !   df3_omp(i) = df3_omp(i) + parts%mass(j)*fwx
            !   df3_omp(j) = df3_omp(j) - parts%mass(i)*fwx
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        df3_omp%r(i) = 0.d0
        do it = 1, nthreads
            df3_omp%r(i) = df3_omp%r(i)+local(i,it)
        enddo
        df3_omp%r(i) = df3_omp%r(i)*parts%rho%r(i)
        !USAW
        if(parts%usaw.and.parts%itype(i)>0) df3_omp%r(i) = df3_omp%r(i)/parts%wi%r(i)
    enddo   
    !$omp end do
    !$omp end parallel


    !do i = 1, ntotal    !moved up
    !   df3_omp%r(i) = df3_omp%r(i)*parts%rho%r(i)
    !enddo


    end function

    ! Calculate partial derivatives of a field
    !-------------------------------------------
    function df4(parts,f,x)
    !-------------------------------------------
    implicit none

    type(array) :: f
    character(len=1) x
    class(particles),target :: parts
    type(array),allocatable :: df4
    real(dp), pointer, dimension(:) :: dwdx
    real(dp) fwx
    integer i, j, k, ntotal
    ntotal = parts%ntotal +parts%nvirt

    allocate(df4); allocate(df4%r(ntotal))
    df4%parts => parts
    df4%ndim1 = ntotal; df4 = 0.d0

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        fwx = (f%r(j)-f%r(i))*dwdx(k)
        df4%r(i) = df4%r(i) + parts%mass%r(j)/parts%rho%r(j)*fwx
        df4%r(j) = df4%r(j) + parts%mass%r(i)/parts%rho%r(i)*fwx
    enddo
    end function

    !-------------------------------------------------
    function df4_omp(parts,f,x) result(fun)
    !-------------------------------------------------
    implicit none

    type(array) f
    character(len=1) x
    class(particles),target :: parts
    type(array) :: fun
    real(dp), allocatable, dimension(:,:) :: local
    real(dp), pointer, dimension(:) :: dwdx, bn
    real(dp) fwx
    integer i, j, k, ntotal, it, nthreads

    !write(*,*) 'In df_omp...'

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads

    fun%parts => parts
    fun%ndim1 = ntotal
    allocate(fun%r(ntotal))
    if(nthreads>=1)then
        allocate(local(ntotal,nthreads))
        call parts%get_niac_start_end
    endif   

    if(x=='x')dwdx=>parts%dwdx(1,:)
    if(x=='y')dwdx=>parts%dwdx(2,:)
    if(parts%usaw)then
        if(x=='x')bn=>parts%bn%x%r
        if(x=='y')bn=>parts%bn%y%r
    endif
    !$omp parallel
    !$omp do private(i,j,k,fwx) 
    do it = 1, parts%nthreads
        do i = 1, ntotal
            local(i,it) = 0.d0
        enddo  
        do k = parts%niac_start(it), parts%niac_end(it)
            !do k=1,parts%niac
            i = parts%pair_i(k)
            j = parts%pair_j(k)

            ! USAW________________________________________________________________

            if(parts%usaw.and.parts%itype(i)>0.and.parts%itype(j)<0)then   
                fwx = f%r(j)-f%r(i)
                local(i,it) = local(i,it) - fwx*parts%w(k)*parts%bs%r(j)*bn(j)
                cycle
            endif        
            if(parts%usaw.and.parts%itype(i)<0.and.parts%itype(j)<0)cycle   
            !---------------------------------------------------------------------

            fwx = (f%r(j)-f%r(i))*dwdx(k)
            local(i,it) = local(i,it) + parts%mass%r(j)/parts%rho%r(j)*fwx
            local(j,it) = local(j,it) + parts%mass%r(i)/parts%rho%r(i)*fwx
        enddo !k
    enddo !it
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        fun%r(i) = 0.d0
        do it = 1, nthreads
            fun%r(i) = fun%r(i) + local(i,it)
        enddo
        if(parts%usaw.and.parts%itype(i)>0) fun%r(i)=fun%r(i)/parts%wi%r(i)   ! USAW

    enddo   
    !$omp end do
    !$omp end parallel
    end function

    !-------------------------------------------
    function div(parts,f)
    !-------------------------------------------
    implicit none

    type(array) :: f
    class(particles), target :: parts
    type(array),allocatable :: div
    real(dp) df(3)
    real(dp) :: hdiv
    real(dp),pointer,dimension(:,:) :: dwdx
    integer, pointer, dimension(:) :: pair_i, pair_j
    integer i,j,k,ntotal,nthreads,niac,dim,d
    type(p2r) f_i(3),f_j(3)

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads
    niac = parts%niac; dim = parts%dim
    dwdx =>parts%dwdx

    div%parts => parts
    allocate(div);allocate(div%r(ntotal))
    div%ndim1 = ntotal;div = 0.d0

    do k = 1,niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        f_i = f%cmpt(i); f_j = f%cmpt(j)
        do d=1,dim
            df(d) = f_j(d)%p - f_i(d)%p
        enddo
        hdiv = df(1)*dwdx(1,k)
        do d = 2,dim
            hdiv = hdiv + df(d)*dwdx(d,k)
        enddo
        div%r(i) = div%r(i) + parts%mass%r(j)*hdiv/parts%rho%r(j)
        div%r(j) = div%r(j) + parts%mass%r(i)*hdiv/parts%rho%r(i)
    enddo
    end function

    !--------------------------------------------------
    function div_omp(parts,f) result(fun)
    !--------------------------------------------------
    implicit none

    type(array) :: f
    class(particles),target :: parts
    type(array),allocatable :: fun
    real(dp) df(3)
    real(dp) :: hdiv
    real(dp),pointer,dimension(:,:) :: dwdx
    integer, pointer, dimension(:) :: pair_i, pair_j
    real(dp), allocatable, dimension(:,:) :: local
    integer i,j,k,ntotal,nthreads,niac,dim,d,it
    type(p2r) f_i(3),f_j(3)

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads
    niac = parts%niac; dim = parts%dim
    dwdx =>parts%dwdx

    allocate(fun);allocate(fun%r(ntotal))
    fun%ndim1 = ntotal
    fun%parts => parts

    if(nthreads>=1)then
        allocate(local(ntotal,nthreads))
        call parts%get_niac_start_end
    endif   

    !$omp parallel
    !$omp do private(i,j,d,f_i,f_j,df,k,hdiv)
    do it = 1, nthreads
        do i = 1, ntotal
            local(i,it) = 0.d0
        enddo   
        do k = parts%niac_start(it), parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            f_i = f%cmpt(i); f_j = f%cmpt(j)
            do d=1,dim
                df(d) = f_j(d)%p - f_i(d)%p
            enddo

            ! USAW________________________________________________________________

            if(parts%usaw.and.parts%itype(i)>0.and.parts%itype(j)<0)then   
                hdiv = df(1)*parts%bn%x%r(j)+df(2)*parts%bn%y%r(j)
                local(i,it) = local(i,it) - hdiv*parts%w(k)*parts%bs%r(j)
                cycle
            endif        
            if(parts%usaw.and.parts%itype(i)<0.and.parts%itype(j)<0)cycle   
            !---------------------------------------------------------------------


            hdiv = df(1)*dwdx(1,k)
            do d = 2,dim
                hdiv = hdiv + df(d)*dwdx(d,k)
            enddo
            local(i,it) = local(i,it) + parts%mass%r(j)*hdiv/parts%rho%r(j)
            local(j,it) = local(j,it) + parts%mass%r(i)*hdiv/parts%rho%r(i)
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        fun%r(i) = 0.d0
        do it = 1,nthreads
            fun%r(i) =  fun%r(i) + local(i,it) 
        enddo  
        if(parts%usaw.and.parts%itype(i)>0) fun%r(i)=fun%r(i)/parts%wi%r(i)   ! USAW

    enddo   
    !$omp end do
    !$omp end parallel

    end function

    !--------------------------------------------------
    function lap_omp(parts,miu,f) result(fun)
    !--------------------------------------------------
    implicit none

    type(array) :: f
    class(particles),target :: parts
    type(array),allocatable :: fun
    real(dp) miu
    real(dp) df,dx(3)
    real(dp) :: hdiv,rr,small = 1.0d-8
    real(dp),pointer,dimension(:,:) :: dwdx
    integer, pointer, dimension(:) :: pair_i, pair_j
    real(dp), allocatable, dimension(:,:) :: local
    integer i,j,k,ntotal,nthreads,niac,dim,d,it
    type(p2r) f_i(3),f_j(3)

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads
    niac = parts%niac; dim = parts%dim
    dwdx =>parts%dwdx

    allocate(fun);allocate(fun%r(ntotal))
    fun%ndim1 = ntotal
    fun%parts => parts

    if(nthreads>=1)then
        allocate(local(ntotal,nthreads))
        call parts%get_niac_start_end
    endif   

    !$omp parallel
    !$omp do private(i,j,d,f_i,f_j,df,k,hdiv,rr,dx)
    do it = 1, nthreads
        do i = 1, ntotal
            local(i,it) = 0.d0
        enddo   
        do k = parts%niac_start(it), parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            df = f%r(i) - f%r(j)
            hdiv = 0.d0
            rr = 0.d0
            do d = 1, dim
                dx(d) = parts%x(d,i) - parts%x(d,j)
                rr = rr + dx(d)*dx(d)
                hdiv = hdiv + dx(d)*dwdx(d,k)
            enddo
            hdiv = 2.*miu*hdiv*df/(rr+small)/(parts%rho%r(i)*parts%rho%r(j))
            local(i,it) = local(i,it) + parts%mass%r(j)*hdiv
            local(j,it) = local(j,it) - parts%mass%r(i)*hdiv
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        fun%r(i) = 0.d0
        do it = 1,nthreads
            fun%r(i) =  fun%r(i) + local(i,it) 
        enddo
        fun%r(i) = fun%r(i)*parts%rho%r(i)  
    enddo   
    !$omp end do
    !$omp end parallel

    end function



    !-----------------------------------------------------------------------
    subroutine velocity_divergence(parts)
    !-----------------------------------------------------------------------
    implicit none

    class(particles) parts      
    integer, pointer, dimension(:) :: pair_i, pair_j
    double precision, pointer, dimension(:) :: mass, rho, vcc
    double precision, pointer, dimension(:,:) :: vx, dwdx
    integer i,j,k,d, dim, ntotal, niac
    double precision dvx(3), hvcc
    type(p2r) vx_i(3), vx_j(3)

    pair_i => parts%pair_i
    pair_j => parts%pair_j
    mass   => parts%mass%r
    rho    => parts%rho%r
    vcc    => parts%vcc%r
    !vx     => parts%vx
    dwdx   => parts%dwdx

    ntotal = parts%ntotal + parts%nvirt
    niac   = parts%niac; dim = parts%dim 

    do i=1,ntotal
        vcc(i) = 0.e0         
    enddo

    do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
            !dvx(d) = vx(d,j) - vx(d,i) 
            dvx(d) = vx_j(d)%p - vx_i(d)%p 
        enddo        
        hvcc = dvx(1)*dwdx(1,k)
        do d=2,dim
            hvcc = hvcc + dvx(d)*dwdx(d,k)
        enddo    
        vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
        vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
    enddo  

    return
    end subroutine

    !-------------------------------------------
    function div2(parts,f) result(res)
    !-------------------------------------------
    implicit none

    type(array) :: f
    class(particles),target :: parts
    type(array),allocatable :: res
    real(dp) df(3)
    real(dp) :: temp
    real(dp),pointer,dimension(:,:) :: dwdx
    integer, pointer, dimension(:) :: pair_i, pair_j
    real(dp), allocatable, dimension(:,:) :: local
    integer i,j,k,ntotal,nthreads,niac,dim,d, it
    type(p2r) f_i(3),f_j(3)

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads
    niac = parts%niac; dim = parts%dim
    dwdx =>parts%dwdx

    allocate(res);allocate(res%r(ntotal))
    res%ndim1 = ntotal
    res%parts=>parts

    allocate(local(ntotal,nthreads))
    call parts%get_niac_start_end

    !$omp parallel
    !$omp do private(k,i,j,f_i,f_j,d,df,temp) 
    do it = 1, nthreads
        do i = 1, ntotal
            local(i,it) = 0.d0
        enddo    
        do k = parts%niac_start(it),parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            f_i = f%cmpt(i); f_j = f%cmpt(j)
            temp = 0.d0
            do d=1,dim
                df(d) = f_j(d)%p - f_i(d)%p
                temp = temp + df(d)*dwdx(d,k)
            enddo
            local(i,it) = local(i,it) + parts%mass%r(j)*temp
            local(j,it) = local(j,it) + parts%mass%r(i)*temp
        enddo !k
    enddo !it
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        res%r(i) = 0
        do it = 1, nthreads
            res%r(i) = res%r(i) + local(i,it)
        enddo
        res%r(i) = res%r(i)/parts%rho%r(i)
    enddo
    !$omp end do
    !$omp end parallel


    end function


    !----------------------------------------------------
    function div2_USAW(parts,f) result(res)
    !----------------------------------------------------
    implicit none

    type(array) :: f
    class(particles),target :: parts
    type(array),allocatable :: res
    real(dp) df(3)
    real(dp) :: temp
    real(dp),pointer,dimension(:,:) :: dwdx
    integer, pointer, dimension(:) :: pair_i, pair_j
    real(dp), allocatable, dimension(:,:) :: local
    integer i,j,k,ntotal,nthreads,niac,dim,d, it
    type(p2r) f_i(3),f_j(3),bn(3)

    ntotal = parts%ntotal + parts%nvirt
    nthreads = parts%nthreads
    niac = parts%niac; dim = parts%dim
    dwdx =>parts%dwdx

    allocate(res);allocate(res%r(ntotal))
    res%ndim1 = ntotal
    res%parts=>parts

    allocate(local(ntotal,nthreads))
    call parts%get_niac_start_end

    !$omp parallel
    !$omp do private(k,i,j,f_i,f_j,d,df,temp) 
    do it = 1, nthreads
        do i = 1, ntotal
            local(i,it) = 0.d0
        enddo    
        do k = parts%niac_start(it),parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            f_i = f%cmpt(i); f_j = f%cmpt(j)
            temp = 0.d0
            do d=1,dim
                df(d) = f_j(d)%p - f_i(d)%p
                temp = temp + df(d)*dwdx(d,k)
            enddo
            local(i,it) = local(i,it) + parts%mass%r(j)*temp
            local(j,it) = local(j,it) + parts%mass%r(i)*temp

            ! USAW
            if(parts%usaw)then
                if(parts%itype(i)<0.and.parts%itype(j)>0)then
                    temp = 0.d0
                    bn = parts%bn%cmpt(i)
                    do d = 1, dim
                        temp = temp - df(d)*parts%w(k)*parts%bs%r(i)*bn(d)%p*parts%rho%r(i)
                    enddo
                    local(j,it) = local(j,it) - temp
                endif
                if(parts%itype(i)>0.and.parts%itype(j)<0)then
                    temp = 0.d0
                    bn = parts%bn%cmpt(j)
                    do d = 1, dim
                        temp = temp + df(d)*parts%w(k)*parts%bs%r(j)*bn(d)%p*parts%rho%r(j)
                    enddo
                    local(i,it) = local(i,it) - temp
                endif
            endif 


        enddo !k
    enddo !it
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        res%r(i) = 0
        do it = 1, nthreads
            res%r(i) = res%r(i) + local(i,it)
        enddo
        if(.not.parts%usaw) res%r(i) = res%r(i)/parts%rho%r(i)
        if(parts%usaw.and.parts%itype(i)>0) res%r(i)=res%r(i)/parts%rho%r(i)/parts%wi%r(i)  !USAW
    enddo
    !$omp end do
    !$omp end parallel


    end function

    !---------------------------------------------------------------------
    subroutine initial_density(parts)
    !---------------------------------------------------------------------
    implicit none

    class(particles) parts

    type(material), pointer :: water, soil
    integer ntotal

    ntotal = parts%ntotal+parts%nvirt

    if(parts%imaterial=='water')then

    water => parts%material
    parts%rho = water%rho0*(parts%p/water%b+1.d0)**(1/water%gamma)

    elseif(parts%imaterial=='soil')then

    soil => parts%material
    !parts%rho(1:ntotal) = soil%rho0*(parts%p(1:ntotal)/soil%k+1)

    endif

    return
    end subroutine

    !----------------------------------------------------------------------
    subroutine pressure(parts)
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts
    type(material), pointer :: water, soil
    integer ntotal, i

    ntotal = parts%ntotal+parts%nvirt

    if(parts%imaterial=='water')then

    water => parts%material
    parts%p = water%b*((parts%rho/(water%rho0*parts%vof))**water%gamma-1.d0) !!! False density  
    !if(parts%water_pressure==0)then   !!!Hydrostatic pressure!
    !   do i = 1, ntotal
    !      parts%p%r(i) = water%rho0*(-9.81)*(parts%x(2,i)-0.1)
    !   enddo     
    !endif
    ! Tension instability
    !                              if(water_tension_instability==1)then
    !         do i = 1, ntotal
    !            if(parts%p(i)<0)then
    !               parts%p(i) = 0.d0
    !               parts%rho(i) = water%rho0*parts%vof(i)
    !            endif
    !         enddo
    !                              endif

    !parts%c(1:ntotal) = water%c         
    parts%c = water%c         
    !parts%c = water%c*(parts%rho/(water%rho0*parts%vof))**3.0         

    elseif(parts%imaterial=='soil')then

    soil => parts%material
    !         if(parts%soil_pressure==1)          &
    !         parts%p(1:ntotal) = soil%k*(parts%rho(1:ntotal)/soil%rho0-1)
    !!         parts%p(1:ntotal) = parts%p(1:ntotal)
    !!     &                      -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
    parts%c = soil%c

    endif

    return      
    end subroutine

    !--------------------------------------------------------------
    subroutine pressure_water(parts)  !not used
    !--------------------------------------------------------------
    implicit none

    class(particles) parts
    type(material), pointer :: water, soil
    integer ntotal, i

    ntotal = parts%ntotal+parts%nvirt

    water => parts%material
    parts%p%r(1:ntotal) = water%b*((parts%rho%r(1:ntotal)/(water%rho0  &
    *parts%vof%r(1:ntotal))) &   !!! False density
    **water%gamma-1)

    parts%c%r(1:ntotal) = water%c*(parts%rho%r(1:ntotal)/(water%rho0*parts%vof%r(1:ntotal)))**3.0         

    return
    end subroutine

    !--------------------------------------------------------------
    subroutine tension_instability_water(parts)   !not used
    !--------------------------------------------------------------
    implicit none

    class(particles) parts
    type(material), pointer :: water, soil
    integer ntotal, i

    ntotal = parts%ntotal+parts%nvirt
    water => parts%material

    do i = 1, ntotal
        if(parts%p%r(i)<0)then
            parts%p%r(i) = 0.d0
            !parts%rho(i) = water%rho0*parts%vof(i)
        endif
    enddo

    return
    end subroutine


    !------------------------------------------------------------
    subroutine pressure_soil(parts) !not used
    !------------------------------------------------------------
    implicit none

    class(particles) parts
    type(material), pointer :: water, soil
    integer ntotal, i

    ntotal = parts%ntotal+parts%nvirt

    soil => parts%material
    parts%p%r(1:ntotal) = soil%k*(parts%rho%r(1:ntotal)/soil%rho0-1)
    ! parts%p(1:ntotal) = parts%p(1:ntotal)  &
    !                    -soil%k*parts%vcc(1:ntotal)*0.000005   !*dt
    parts%c%r(1:ntotal) = soil%c

    return      
    end subroutine

    !Subroutine to calculate the density with SPH summation algorithm.
    !----------------------------------------------------------------------
    subroutine sum_density(parts) 
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts
    integer ntotal, i, j, k, d      
    real(dp) selfdens, hv(3), r
    real(dp), allocatable, dimension(:) :: wi     

    ntotal = parts%ntotal + parts%nvirt

    !     wi(maxn)---integration of the kernel itself
    allocate(wi(parts%maxn))

    hv = 0.d0

    !     Self density of each particle: Wii (Kernel for distance 0)
    !     and take contribution of particle itself:

    r=0.d0

    !     Firstly calculate the integration of the kernel over the space

    do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        wi(i)=selfdens*parts%mass%r(i)/parts%rho%r(i)
    enddo

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        wi(i) = wi(i) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k)
        wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
    enddo

    !     Secondly calculate the rho integration over the space

    do i=1,ntotal
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        parts%rho%r(i) = selfdens*parts%mass%r(i)
    enddo

    !     Calculate SPH sum for rho:
    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        parts%rho%r(i) = parts%rho%r(i) + parts%mass%r(j)*parts%w(k)
        parts%rho%r(j) = parts%rho%r(j) + parts%mass%r(i)*parts%w(k)
    enddo

    !     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
    !      if (nor_density) then 
    do i=1, ntotal
        parts%rho%r(i)=parts%rho%r(i)/wi(i)
    enddo
    !      endif 

    end subroutine

    !Subroutine to interpolate fields between two particle systems
    !----------------------------------------------------------------------
    subroutine interpolation(water,fwater,soil,fsoil) 
    !----------------------------------------------------------------------
    implicit none

    type(particles) water,soil
    type(array) fwater, fsoil
    real(dp), allocatable, dimension(:,:) :: local
    integer ntotal, i, j, k, d, it, nthreads      

    nthreads = water%nthreads
    ntotal = soil%ntotal + soil%nvirt
    allocate(local(ntotal,nthreads))
    call water%get_niac_start_end

    !$omp parallel
    !$omp do private(i,j,k) 
    do it =1,nthreads
        do j = 1,ntotal
            local(j,it) = 0.d0
        enddo 
        do k = water%niac_start(it),water%niac_end(it)
            i = water%pair_i(k)   ! water
            j = water%pair_j(k)   ! soil
            local(j,it) = local(j,it) + water%mass%r(i)*fwater%r(i)  &  
            *water%w(k)/water%rho%r(i)
        enddo  
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do j = 1,soil%ntotal+soil%nvirt
        fsoil%r(j) = 0.d0
        do it = 1,water%nthreads
            fsoil%r(j) = fsoil%r(j) + local(j,it)
        enddo  
    enddo   
    !$omp end do
    !$omp end parallel
    return
    end subroutine 



    !!DEC$IF(.FALSE.)
    !----------------------------------------------------------------------
    subroutine Sherpard_filter(parts) 
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts
    integer ntotal, i, j, k, d      
    real(dp) selfdens, hv(3), r
    !      real(dp), allocatable, dimension(:) :: wi     
    real(dp), pointer, dimension(:) :: wi     

    ntotal = parts%ntotal + parts%nvirt

    !     wi(maxn)---integration of the kernel itself
    !      allocate(wi(parts%maxn))
    wi => parts%wi%r

    hv = 0.d0

    !     Self density of each particle: Wii (Kernel for distance 0)
    !     and take contribution of particle itself:

    r=0.d0

    !     Firstly calculate the integration of the kernel over the space

    do i=1,ntotal
        wi(i) = 0.d0  ! initialized to 0
        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
        if(parts%itype(i)>0)then
            wi(i)=selfdens*parts%mass%r(i)/parts%rho%r(i)
        endif
    enddo

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(j)>0)wi(i) = wi(i) + parts%mass%r(j)/parts%rho%r(j)*parts%w(k)
        if(parts%itype(i)>0)wi(j) = wi(j) + parts%mass%r(i)/parts%rho%r(i)*parts%w(k)
    enddo

    !      do i=1,ntotal
    !      if(wi(i) < 1.d-8) wi(i)=1.d10
    !      enddo

    !     Secondly calculate the rho integration over the space

    !      do i=1,ntotal
    !        call parts%kernel(r,hv,parts%hsml(i),selfdens,hv)
    !        parts%rho%r(i) = selfdens*parts%mass%r(i)
    !      enddo

    !     Calculate SPH sum for rho:
    !      do k=1,parts%niac
    !        i = parts%pair_i(k)
    !        j = parts%pair_j(k)
    !        parts%rho%r(i) = parts%rho%r(i) + parts%mass%r(j)*parts%w(k)
    !        parts%rho%r(j) = parts%rho%r(j) + parts%mass%r(i)*parts%w(k)
    !      enddo

    !     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
    !!      if (nor_density) then 
    !        do i=1, ntotal
    !          parts%rho%r(i)=parts%rho%r(i)/wi(i)
    !        enddo
    !!      endif 

    end subroutine


    !!DEC$ENDIF


    !----------------------------------------------------------------------
    subroutine get_boundary_field(parts,f) 
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts
    type(array) f
    integer ntotal, i, j, k, d      

    ntotal = parts%ntotal + parts%nvirt

    do i = 1, ntotal
        if(parts%itype(i)>0)cycle
        f%r(i)=0.d0
    enddo

    do k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        if(parts%itype(i)<0.and.parts%itype(j)>0)then
            f%r(i) = f%r(i) +  parts%mass%r(j)*f%r(j)*parts%w(k)/parts%rho%r(j)
        endif

        if(parts%itype(i)>0.and.parts%itype(j)<0)then
            f%r(j) = f%r(j) +  parts%mass%r(i)*f%r(i)*parts%w(k)/parts%rho%r(i)
        endif
    enddo

    !     Calculate the normalized rho, rho=sum(rho)/sum(w)
    !      if (nor_density) then 
    do i=1, ntotal
        if(parts%wi%r(i)<1.d-8)cycle
        if(parts%itype(i)>0)cycle
        f%r(i)=f%r(i)/parts%wi%r(i)
    enddo
    !      endif 

    end subroutine


    ! Subroutine to calculate the density with SPH continuiity approach.
    !----------------------------------------------------------------------      
    subroutine con_density(parts)
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts      

    integer ntotal,i,j,k,d    
    double precision vcc, dvx(3), dwdx(3) 
    type(p2r) vx_i(3), vx_j(3) 

    ntotal = parts%ntotal + parts%nvirt

    do i = 1, ntotal
        parts%drho%r(i) = 0.
    enddo

    do k=1,parts%niac      
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1, parts%dim
            dwdx(d) = parts%dwdx(d,k) 
            dvx(d) = vx_i(d)%p - vx_j(d)%p 
        enddo        
        !vcc = dvx(1)* parts%dwdx(1,k)        
        vcc = dvx(1)*dwdx(1)        
        do d=2, parts%dim
            !vcc = vcc + dvx(d)*parts%dwdx(d,k)
            vcc = vcc + dvx(d)*dwdx(d)
        enddo    
        parts%drho%r(i) = parts%drho%r(i) + parts%mass%r(j)*vcc
        parts%drho%r(j) = parts%drho%r(j) + parts%mass%r(i)*vcc       
    enddo    

    end subroutine

    !----------------------------------------------------------------------
    subroutine newtonian_fluid(parts)
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts
    type(material), pointer :: liquid
    integer ntotal

    ntotal = parts%ntotal + parts%nvirt
    liquid => parts%material
    !parts%tab%x%ndim1 = ntotal; parts%tab%y%ndim1 = ntotal; parts%tab%xy%ndim1 = ntotal

    parts%str%x = liquid%viscosity*parts%tab%x
    parts%str%y = liquid%viscosity*parts%tab%y
    parts%str%xy = liquid%viscosity*parts%tab%xy

    return
    end subroutine

    ! Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
    !----------------------------------------------------------------------      
    subroutine art_visc(parts)
    !----------------------------------------------------------------------
    implicit none

    class(particles), target :: parts

    type(numerical), pointer :: numeric
    real(dp) dx, dvx(3), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml
    type(p2r) vx_i(3), vx_j(3)
    integer i,j,k,d,dim,ntotal,niac

    ntotal   =  parts%ntotal + parts%nvirt
    niac     =  parts%niac; dim = parts%dim

    numeric  => parts%numeric      
    alpha = numeric%alpha; beta = numeric%beta; etq = numeric%etq

    !     Calculate SPH sum for artificial viscosity

    do k=1,niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        mhsml= (parts%hsml(i)+parts%hsml(j))/2.
        vr = 0.e0
        rr = 0.e0
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
            !dvx(d) = parts%vx(d,i) - parts%vx(d,j)
            dvx(d) = vx_i(d)%p - vx_j(d)%p
            dx     = parts%x(d,i)  - parts%x(d,j)
            vr     = vr + dvx(d)*dx
            rr     = rr + dx*dx
        enddo

        !     Artificial viscous force only if v_ij * r_ij < 0

        if (vr.lt.0.e0) then

        !     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )

        muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)

        !     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

        mc   = 0.5e0*(parts%c%r(i) + parts%c%r(j))
        mrho = 0.5e0*(parts%rho%r(i) + parts%rho%r(j))
        piv  = (beta*muv - alpha*mc)*muv/mrho              

        !     Calculate SPH sum for artificial viscous force

        !do d=1,dim
        !  h = -piv*parts%dwdx(d,k)
        !  parts%dvx(d,i) = parts%dvx(d,i) + parts%mass(j)*h
        !  parts%dvx(d,j) = parts%dvx(d,j) - parts%mass(i)*h
        !enddo

        h = -piv*parts%dwdx(1,k)
        parts%dvx%x%r(i) = parts%dvx%x%r(i) + parts%mass%r(j)*h
        parts%dvx%x%r(j) = parts%dvx%x%r(j) - parts%mass%r(i)*h
        if(dim.ge.2)then
            h = -piv*parts%dwdx(2,k)
            parts%dvx%y%r(i) = parts%dvx%y%r(i) + parts%mass%r(j)*h
            parts%dvx%y%r(j) = parts%dvx%y%r(j) - parts%mass%r(i)*h
        endif


        endif
    enddo

    return
    end subroutine

    ! Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
    !----------------------------------------------------------------------      
    subroutine art_visc_omp(parts)
    !----------------------------------------------------------------------
    implicit none

    class(particles),target :: parts

    type(numerical), pointer :: numeric
    real(dp) dx, dvx(3), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml
    real(dp), allocatable, dimension(:,:,:) :: local
    type(p2r) vx_i(3), vx_j(3)
    integer i,j,k,d,dim,ntotal,niac,nthreads,it

    ntotal   =  parts%ntotal + parts%nvirt
    niac     =  parts%niac; dim = parts%dim
    nthreads =  parts%nthreads

    allocate(local(dim,ntotal,nthreads))
    call parts%get_niac_start_end

    numeric  => parts%numeric      
    alpha = numeric%alpha; beta = numeric%beta; etq = numeric%etq

    !     Calculate SPH sum for artificial viscosity

    !$omp parallel 
    !$omp do private(i,j,k,mhsml,vr,rr,vx_i,vx_j,d,dvx,dx,muv,mc,mrho,piv,h) 
    do it = 1, nthreads 
        do i = 1, ntotal
            local(1:dim,i,it) = 0.d0
        enddo
        do k=parts%niac_start(it),parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            mhsml= (parts%hsml(i)+parts%hsml(j))/2.
            vr = 0.e0
            rr = 0.e0
            vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
            do d=1,dim
                !dvx(d) = parts%vx(d,i) - parts%vx(d,j)
                dvx(d) = vx_i(d)%p - vx_j(d)%p
                dx     = parts%x(d,i)  - parts%x(d,j)
                vr     = vr + dvx(d)*dx
                rr     = rr + dx*dx
            enddo

            !     Artificial viscous force only if v_ij * r_ij < 0

            if (vr.lt.0.e0) then

            !     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )

            muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)

            !     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

            mc   = 0.5e0*(parts%c%r(i) + parts%c%r(j))
            mrho = 0.5e0*(parts%rho%r(i) + parts%rho%r(j))
            piv  = (beta*muv - alpha*mc)*muv/mrho              

            !     Calculate SPH sum for artificial viscous force

            do d=1,dim
                h = -piv*parts%dwdx(d,k)
                !  parts%dvx(d,i) = parts%dvx(d,i) + parts%mass(j)*h
                !  parts%dvx(d,j) = parts%dvx(d,j) - parts%mass(i)*h
                local(d,i,it) = local(d,i,it) + parts%mass%r(j)*h
                local(d,j,it) = local(d,j,it) - parts%mass%r(i)*h
            enddo
            endif
        enddo !k
    enddo !it
    !$omp end do
    !$omp barrier

    if(.not.parts%usaw)then
        !$omp do private(it)
        do i = 1, ntotal
            do it = 1 ,nthreads
                parts%dvx%x%r(i) = parts%dvx%x%r(i) + local(1,i,it)
                if (dim.ge.2) parts%dvx%y%r(i) = parts%dvx%y%r(i) +local(2,i,it)
            enddo  
        enddo
        !$omp end do
    elseif(parts%usaw)then

    !$omp do private(it)
    do i = 1, ntotal
        do it = 1 ,nthreads
            if(parts%itype(i)>0)then                               
                parts%dvx%x%r(i) = parts%dvx%x%r(i) + local(1,i,it)/parts%wi%r(i)
                if (dim.ge.2) parts%dvx%y%r(i) = parts%dvx%y%r(i)     &
                +local(2,i,it)/parts%wi%r(i)
            endif
        enddo  
    enddo                                     
    !$omp end do
    endif
    !$omp end parallel

    return
    end subroutine

    !-------------------------------------------------------------------
    subroutine delta_sph(parts,f,df)
    !-------------------------------------------------------------------
    implicit none

    class(particles) parts
    real(dp), dimension(:) :: f,df

    real(dp) dx(3),delta, muv, rr, h, small
    integer i,j,k,d,ntotal,niac,dim

    !      write(*,*) 'In art_density...'

    ntotal   =  parts%ntotal + parts%nvirt
    niac     =  parts%niac; dim = parts%dim            
    delta    = parts%numeric%delta
    small    = 1.d-8

    do k=1,niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        rr = 0.e0
        do d=1,dim
            dx(d)  =  parts%x(d,i) -  parts%x(d,j)
            rr     = rr + dx(d)*dx(d)
        enddo
        rr = rr + small   !!!!  
        muv = 2.0*(f(i)-f(j))/rr

        h = 0.d0
        do d=1,dim
            !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
            h = h + dx(d)*muv*parts%dwdx(d,k)
        enddo
        df(i) = df(i) + delta*parts%hsml(i)*parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
        df(j) = df(j) - delta*parts%hsml(j)*parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)
    enddo

    return
    end subroutine

    !     Subroutine to calculate the average velocity to correct velocity
    !     for preventing.penetration (monaghan, 1992)      

    !-------------------------------------------------------------------
    subroutine delta_sph_omp(parts,f,df)
    !-------------------------------------------------------------------
    implicit none

    class(particles) parts
    type(array) f
    type(array) df
    real(dp), allocatable, dimension(:,:) :: local
    real(dp) dx(3),delta, muv, rr, h, small,sum_vof
    integer i,j,k,d,ntotal,niac,dim,it,nthreads

    !      write(*,*) 'In art_density...'

    ntotal   =  parts%ntotal + parts%nvirt
    niac     =  parts%niac; dim = parts%dim; nthreads = parts%nthreads            
    delta    = parts%numeric%delta
    small    = 1.d-8

    if(nthreads>=1)then
        allocate(local(ntotal,nthreads))
        call parts%get_niac_start_end
    endif   

    !$omp parallel 
    !$omp do private(i,j,d,k,rr,muv,dx,h,sum_vof)
    do it = 1, nthreads 
        do i =1,ntotal
            local(i,it) = 0.d0
        enddo     
        do k = parts%niac_start(it),parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)

            sum_vof = parts%vof%r(i) + parts%vof%r(j)  !!! Jump condition
            if(sum_vof<1.8.and.sum_vof>1.2) cycle      !!! remedied except dstr

            !if(parts%wi%r(i)<1.017.or.parts%wi%r(j)<1.017)cycle !!! boundary parts

            rr = 0.e0
            do d=1,dim
                dx(d)  =  parts%x(d,i) -  parts%x(d,j)
                rr     = rr + dx(d)*dx(d)
            enddo

            rr = rr + small    !!!!
            muv = 2.0*(f%r(i)-f%r(j))/rr

            h = 0.d0
            do d=1,dim
                !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
                h = h + dx(d)*muv*parts%dwdx(d,k)
            enddo
            !          df(i) = df(i) + delta*parts%hsml(i)*parts%c(i)*parts%mass(j)*h/parts%rho(j)
            !          df(j) = df(j) - delta*parts%hsml(j)*parts%c(j)*parts%mass(i)*h/parts%rho(i)
            if(parts%usaw)then
                if(parts%itype(j)>0)local(i,it) = local(i,it) + delta*parts%hsml(i)*parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
                if(parts%itype(i)>0)local(j,it) = local(j,it) - delta*parts%hsml(j)*parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)


            else
                local(i,it) = local(i,it) + delta*parts%hsml(i)*parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
                local(j,it) = local(j,it) - delta*parts%hsml(j)*parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)
            endif
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1,ntotal

    !if(parts%wi%r(i)<1.017)cycle  !!! remedy dstr of soil

    do it = 1,nthreads
        !df%r(i) =  df%r(i) + local(i,it)   ! Original
        if(parts%usaw.and.parts%itype(i)>0)df%r(i)=df%r(i)+local(i,it)/parts%wi%r(i) !!??
        if(.not.parts%usaw) df%r(i) =  df%r(i) + local(i,it)
    enddo
    !if(parts%usaw.and.parts%itype(i)>0)df%r(i)=df%r(i)/parts%wi%r(i) !!??
    enddo
    !$omp end do
    !$omp end parallel
    return
    end subroutine


    !-------------------------------------------------------------------
    subroutine diffusive_sph_omp(parts,f,df)
    !-------------------------------------------------------------------
    implicit none

    class(particles) parts
    type(array) f
    type(array) df
    real(dp), allocatable, dimension(:,:) :: local
    real(dp) dx(3),delta, muv, rr, h, small
    integer i,j,k,d,ntotal,niac,dim,it,nthreads

    !      write(*,*) 'In art_density...'

    ntotal   =  parts%ntotal + parts%nvirt
    niac     =  parts%niac; dim = parts%dim; nthreads = parts%nthreads            
    delta    = parts%numeric%delta
    small    = 1.d-8

    if(nthreads>=1)then
        allocate(local(ntotal,nthreads))
        call parts%get_niac_start_end
    endif   

    !$omp parallel 
    !$omp do private(i,j,d,k,rr,muv,dx,h)
    do it = 1, nthreads 
        do i =1,ntotal
            local(i,it) = 0.d0
        enddo     
        do k = parts%niac_start(it),parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)
            rr = 0.e0
            do d=1,dim
                dx(d)  =  parts%x(d,i) -  parts%x(d,j)
                rr     = rr + dx(d)*dx(d)
            enddo

            rr = rr + small    !!!!
            muv = (f%r(i)-f%r(j))/sqrt(rr)

            h = 0.d0
            do d=1,dim
                !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
                h = h + dx(d)*muv*parts%dwdx(d,k)
            enddo
            !          df(i) = df(i) + delta*parts%hsml(i)*parts%c(i)*parts%mass(j)*h/parts%rho(j)
            !          df(j) = df(j) - delta*parts%hsml(j)*parts%c(j)*parts%mass(i)*h/parts%rho(i)
            if(parts%usaw)then
                if(parts%itype(j)>0)local(i,it) = local(i,it) + parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
                if(parts%itype(i)>0)local(j,it) = local(j,it) - parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)


            else
                local(i,it) = local(i,it) + parts%c%r(i)*parts%mass%r(j)*h/parts%rho%r(j)
                local(j,it) = local(j,it) - parts%c%r(j)*parts%mass%r(i)*h/parts%rho%r(i)
            endif
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1,ntotal
        do it = 1,nthreads
            !df%r(i) =  df%r(i) + local(i,it)   ! Original
            if(parts%usaw.and.parts%itype(i)>0)df%r(i)=df%r(i)+local(i,it)/parts%wi%r(i) !!??
            if(.not.parts%usaw) df%r(i) =  df%r(i) + local(i,it)
        enddo
        !if(parts%usaw.and.parts%itype(i)>0)df%r(i)=df%r(i)/parts%wi%r(i) !!??
    enddo
    !$omp end do
    !$omp end parallel
    return
    end subroutine


    !----------------------------------------------------------------------      
    subroutine av_vel(parts)
    !----------------------------------------------------------------------   
    implicit none

    class(particles),target :: parts
    real(dp) vcc, dvx(3), epsilon, mrho
    type(p2r) vx_i(3), vx_j(3), av_i(3), av_j(3)
    integer i,j,k,d,ntotal,niac    

    ntotal = parts%ntotal
    niac   = parts%niac
    epsilon = parts%numeric%epsilon

    parts%av%x = 0.d0; parts%av%y = 0.d0

    do k=1,niac       
        i = parts%pair_i(k)
        j = parts%pair_j(k)       
        mrho = (parts%rho%r(i)+parts%rho%r(j))/2.0
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        av_i = parts%av%cmpt(i); av_j = parts%av%cmpt(j)
        do d=1,parts%dim
            !dvx(d) = parts%vx(d,i) - parts%vx(d,j)            
            dvx(d) = vx_i(d)%p - vx_j(d)%p            
            av_i(d)%p = av_i(d)%p - parts%mass%r(j)*dvx(d)/mrho*parts%w(k)
            av_j(d)%p = av_j(d)%p + parts%mass%r(i)*dvx(d)/mrho*parts%w(k)       
        enddo                    
    enddo  

    !      do i = 1, ntotal
    !         do d = 1, parts%dim
    !            av_i(d)%p = epsilon * av_i(d)%p
    !         enddo 
    !      enddo

    parts%av%x = epsilon * parts%av%x; parts%av%y = epsilon * parts%av%y

    return
    end subroutine

    !----------------------------------------------------------------------      
    subroutine av_vel_omp(parts)
    !----------------------------------------------------------------------   
    implicit none

    class(particles),target :: parts
    real(dp) vcc, dvx(3), epsilon, mrho
    type(p2r) vx_i(3), vx_j(3), av_i(3), av_j(3)
    real(dp),allocatable,dimension(:,:,:) :: local
    integer i,j,k,d,ntotal,niac,it,nthreads  

    ntotal = parts%ntotal + parts%nvirt   !!! originally parts%ntotal
    niac   = parts%niac
    epsilon = parts%numeric%epsilon
    nthreads = parts%nthreads
    allocate(local(parts%dim,parts%ntotal+parts%nvirt,nthreads))
    call parts%get_niac_start_end

    parts%av%x = 0.d0; parts%av%y = 0.d0

    !$omp parallel
    !$omp do private(k,i,j,mrho,vx_i,vx_j,av_i,av_j,d,dvx) 
    do it = 1,nthreads    
        do i = 1, ntotal
            do d= 1, parts%dim
                local(d,i,it)=0.d0
            enddo
        enddo
        do k=parts%niac_start(it),parts%niac_end(it)
            i = parts%pair_i(k)
            j = parts%pair_j(k)       
            mrho = (parts%rho%r(i)+parts%rho%r(j))/2.0

            if(parts%usaw.and.mrho<1.e-8)cycle   !!! USAW

            vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
            !        av_i = parts%av%cmpt(i); av_j = parts%av%cmpt(j)
            do d=1,parts%dim
                !dvx(d) = parts%vx(d,i) - parts%vx(d,j)            
                dvx(d) = vx_i(d)%p - vx_j(d)%p        
                local(d,i,it) = local(d,i,it) - parts%mass%r(j)*dvx(d)/mrho*parts%w(k)
                local(d,j,it) = local(d,j,it) + parts%mass%r(i)*dvx(d)/mrho*parts%w(k)
                !            av_i(d)%p = av_i(d)%p - parts%mass%r(j)*dvx(d)/mrho*parts%w(k)
                !            av_j(d)%p = av_j(d)%p + parts%mass%r(i)*dvx(d)/mrho*parts%w(k)       
            enddo                    
        enddo  
    enddo
    !$omp end do
    !$omp barrier
    !$omp do private(it)
    do i = 1, ntotal
        do it = 1,nthreads
            parts%av%x%r(i) = parts%av%x%r(i) + local(1,i,it)
            parts%av%y%r(i) = parts%av%y%r(i) + local(2,i,it)
        enddo   
    enddo
    !$omp end do
    !$omp end parallel
    !      do i = 1, ntotal
    !         do d = 1, parts%dim
    !            av_i(d)%p = epsilon * av_i(d)%p
    !         enddo 
    !      enddo

    parts%av%x = epsilon * parts%av%x; parts%av%y = epsilon * parts%av%y

    return
    end subroutine


    !----------------------------------------------------------------------
    subroutine shear_strain_rate(parts)  ! not used
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts
    real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz
    type(p2r) vx_i(3), vx_j(3)
    integer i, j, k, d, dim, ntotal, niac

    ntotal = parts%ntotal + parts%nvirt
    niac = parts%niac; dim = parts%dim

    parts%tab%x%r(1:ntotal) = 0.e0
    if(dim>=2)then
        parts%tab%y%r(1:ntotal) = 0.e0
        parts%tab%xy%r(1:ntotal) = 0.e0
    endif
    if(dim==3)then
        parts%tab%z%r(1:ntotal) = 0.e0
        parts%tab%xz%r(1:ntotal) = 0.e0
        parts%tab%yz%r(1:ntotal) = 0.e0
    endif

    !     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

    do k=1,niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)
        vx_i = parts%vx%cmpt(i); vx_j = parts%vx%cmpt(j)
        do d=1,dim
            !dvx(d) = parts%vx(d,j) - parts%vx(d,i)
            dvx(d) = vx_j(d)%p - vx_i(d)%p
        enddo
        if (dim.eq.1) then 
            hxx = 2.e0*dvx(1)*parts%dwdx(1,k)        
        else if (dim.eq.2) then           
            hxx = 2.e0*dvx(1)*parts%dwdx(1,k) - dvx(2)*parts%dwdx(2,k) 
            hxy = dvx(1)*parts%dwdx(2,k) + dvx(2)*parts%dwdx(1,k)
            hyy = 2.e0*dvx(2)*parts%dwdx(2,k) - dvx(1)*parts%dwdx(1,k)
        else if (dim.eq.3) then
            hxx = 2.e0*dvx(1)*parts%dwdx(1,k) - dvx(2)*parts%dwdx(2,k)   &
            - dvx(3)*parts%dwdx(3,k) 
            hxy = dvx(1)*parts%dwdx(2,k) + dvx(2)*parts%dwdx(1,k)
            hxz = dvx(1)*parts%dwdx(3,k) + dvx(3)*parts%dwdx(1,k)          
            hyy = 2.e0*dvx(2)*parts%dwdx(2,k) - dvx(1)*parts%dwdx(1,k)   &    
            - dvx(3)*parts%dwdx(3,k)
            hyz = dvx(2)*parts%dwdx(3,k) + dvx(3)*parts%dwdx(2,k)
            hzz = 2.e0*dvx(3)*parts%dwdx(3,k) - dvx(1)*parts%dwdx(1,k)   &
            - dvx(2)*parts%dwdx(2,k)
        endif                              
        hxx = 2.e0/3.e0*hxx
        hyy = 2.e0/3.e0*hyy
        hzz = 2.e0/3.e0*hzz
        if (dim.eq.1) then 
            parts%tab%x%r(i) = parts%tab%x%r(i) + parts%mass%r(j)*hxx/parts%rho%r(j)
            parts%tab%x%r(j) = parts%tab%x%r(j) + parts%mass%r(i)*hxx/parts%rho%r(i)             
        else if (dim.eq.2) then           
            parts%tab%x%r(i) = parts%tab%x%r(i) + parts%mass%r(j)*hxx/parts%rho%r(j)
            parts%tab%x%r(j) = parts%tab%x%r(j) + parts%mass%r(i)*hxx/parts%rho%r(i)   
            parts%tab%xy%r(i) = parts%tab%xy%r(i) + parts%mass%r(j)*hxy/parts%rho%r(j)
            parts%tab%xy%r(j) = parts%tab%xy%r(j) + parts%mass%r(i)*hxy/parts%rho%r(i)            
            parts%tab%y%r(i) = parts%tab%y%r(i) + parts%mass%r(j)*hyy/parts%rho%r(j)
            parts%tab%y%r(j) = parts%tab%y%r(j) + parts%mass%r(i)*hyy/parts%rho%r(i)          
        else if (dim.eq.3) then
            !             txx(i) = txx(i) + mass(j)*hxx/rho(j)
            !             txx(j) = txx(j) + mass(i)*hxx/rho(i)   
            !             txy(i) = txy(i) + mass(j)*hxy/rho(j)
            !             txy(j) = txy(j) + mass(i)*hxy/rho(i) 
            !             txz(i) = txz(i) + mass(j)*hxz/rho(j)
            !             txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
            !             tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
            !             tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
            !             tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
            !             tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
            !             tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
            !             tzz(j) = tzz(j) + mass(i)*hzz/rho(i)                 
        endif                              
    enddo

    return
    end subroutine      

    !---------------------------------------------------------------------
    subroutine artificial_stress(soil)
    !----------------------------------------------------------------------
    implicit none

    class(particles) soil
    integer :: i,j,ntotal,dim,niac,d,k
    real(dp) :: epsilon = 0.7
    real(dp) :: n = 2.55 !exponent
    real(dp), parameter :: pi = 3.14159265358979323846
    real(dp) :: f, piv, hsml, q, w, cos2,sin2,sincos
    type(p2r) :: strp_i(3), astp_i(3)

    niac   = soil%niac
    dim    = soil%dim
    ntotal = soil%ntotal+soil%nvirt

    !    write(*,*) 'aaaaaaaaaaaaaaaa1',ntotal

    do i = 1, ntotal
        !    write(*,*) 'aaaaaaaaaaaaaaaa1'
        soil%thetai%r(i) = 0.d0
        if(abs(soil%str%x%r(i)-soil%str%y%r(i))<0.000001)cycle
        !        soil%thetai%r(i) = 0.5*atan(2*soil%str%xy%r(i)/(soil%str%x%r(i)-soil%str%y%r(i)+0.00000001))
        soil%thetai%r(i) = 0.5*atan(2.*soil%str%xy%r(i)/(soil%str%x%r(i)-soil%str%y%r(i)))
    enddo

    !  write(*,*) 'aaaaaaaaaaaaaaaa2'

    do i = 1, ntotal
        cos2 = cos(soil%thetai%r(i))**2.0
        sin2 = sin(soil%thetai%r(i))**2.0
        sincos = cos(soil%thetai%r(i))*sin(soil%thetai%r(i))
        !        soil%strp%x%r(i) = cos(soil%thetai%r(i))**2*soil%str%x%r(i)+sin(soil%thetai%r(i))**2*soil%str%y%r(i) &
        !                             +2.0*cos(soil%thetai%r(i))*sin(soil%thetai%r(i))*soil%str%xy%r(i)+soil%p%r(i)
        !        soil%strp%y%r(i) = sin(soil%thetai%r(i))**2*soil%str%x%r(i)+cos(soil%thetai%r(i))**2*soil%str%y%r(i) &
        !                             -2.0*cos(soil%thetai%r(i))*sin(soil%thetai%r(i))*soil%str%xy%r(i)+soil%p%r(i)

        soil%strp%x%r(i) = cos2*soil%str%x%r(i)+sin2*soil%str%y%r(i) &
        + 2.0*sincos*soil%str%xy%r(i)-soil%p%r(i)       ! +p LU
        soil%strp%y%r(i) = sin2*soil%str%x%r(i)+cos2*soil%str%y%r(i) &
        - 2.0*sincos*soil%str%xy%r(i)-soil%p%r(i)       ! +p
        !      enddo 

        !    write(*,*) 'aaaaaaaaaaaaaaaa3' 
        !      do i = 1, ntotal
        strp_i = soil%strp%cmpt(i)
        astp_i = soil%astp%cmpt(i)
        do d = 1, dim
            if (strp_i(d)%p > 0) then                                  ! >0 LU
                astp_i(d)%p = -epsilon*strp_i(d)%p/soil%rho%r(i)**2
            else
                astp_i(d)%p = 0.d0
            endif
        enddo
        soil%ast%x%r(i) = soil%astp%x%r(i)*cos2 + soil%astp%y%r(i)*sin2
        soil%ast%y%r(i) = soil%astp%x%r(i)*sin2 + soil%astp%y%r(i)*cos2
        soil%ast%xy%r(i) = (soil%astp%x%r(i)-soil%astp%y%r(i))*sincos

    enddo

    !   write(*,*) (soil%ast%x%r(i),i=1,10)
    !   write(*,*) (soil%ast%y%r(i),i=1,10)
    !   write(*,*) (soil%ast%xy%r(i),i=1,10)
    !    write(*,*) 'aaaaaaaaaaaaaaaa4'

    q = 0.8333333
    hsml = 1.01/808*1.2
    w = 7.0/(4.0*pi*hsml**2)*((1-q/2)**4*(1+2*q))  !kernel(skf = 4)  w = w(delta d, h)

    !    write(*,*) 'aaaaaaaaaaaaaaaa5'

    do k = 1, niac
        i = soil%pair_i(k)
        j = soil%pair_j(k)
        f = soil%w(k)/w

        soil%dvx%x%r(i) = soil%dvx%x%r(i) + soil%mass%r(j)*((soil%ast%x%r(i)+soil%ast%x%r(j))*soil%dwdx(1,k) &
        +(soil%ast%xy%r(i)+soil%ast%xy%r(j))*soil%dwdx(2,k))*f**n  
        soil%dvx%x%r(j) = soil%dvx%x%r(j) - soil%mass%r(i)*((soil%ast%x%r(i)+soil%ast%x%r(j))*soil%dwdx(1,k) &
        +(soil%ast%xy%r(i)+soil%ast%xy%r(j))*soil%dwdx(2,k))*f**n

        soil%dvx%y%r(i) = soil%dvx%y%r(i) + soil%mass%r(j)*((soil%ast%xy%r(i)+soil%ast%xy%r(j))*soil%dwdx(1,k) &
        +(soil%ast%y%r(i)+soil%ast%y%r(j))*soil%dwdx(2,k))*f**n
        soil%dvx%y%r(j) = soil%dvx%y%r(j) - soil%mass%r(i)*((soil%ast%xy%r(i)+soil%ast%xy%r(j))*soil%dwdx(1,k) &
        +(soil%ast%y%r(i)+soil%ast%y%r(j))*soil%dwdx(2,k))*f**n

    end do

    !   write(*,*) 'aaaaaaaaaaaaaaaa6'

    return
    end subroutine

    !---------------------------------------------------------------------
    subroutine Jaumann_rate(parts)
    !----------------------------------------------------------------------
    implicit none

    class(particles) parts
    type(material), pointer :: soil
    real(dp) dvx(3), hxx, hyy, hzz, hxy, hxz, hyz, G
    type(p2r) vx_i(3), vx_j(3) 
    integer i, j, k, d, dim, ntotal, niac

    ntotal = parts%ntotal + parts%nvirt
    niac = parts%niac; dim = parts%dim
    soil => parts%material
    G = soil%e/(2.0*(1+soil%niu))

    !     Hook's law

    !$omp parallel do
    do i = 1, ntotal
        parts%dstr%x%r(i) = parts%dstr%x%r(i)+G*parts%tab%x%r(i)   ! No accumulation origionaly
        parts%dstr%xy%r(i) = parts%dstr%xy%r(i)+G*parts%tab%xy%r(i)
        parts%dstr%y%r(i) = parts%dstr%y%r(i)+G*parts%tab%y%r(i)
    enddo

    !         if(parts%soil_pressure==2)then
    !$omp parallel do
    do i = 1, ntotal
        parts%dp%r(i)   = parts%dp%r(i) - soil%k*parts%vcc%r(i)  !!! simultaneous pressure  
    enddo
    !         endif         

    !     spin tensor
    parts%wxy = 0.5*(parts%df4_omp(parts%vx%x,'y')-parts%df4_omp(parts%vx%y,'x'))

    ! USAW      
    if(parts%usaw) call get_boundary_field(parts,parts%wxy)

    !   Jaumann rate

    !$omp parallel do
    do i = 1, ntotal
        parts%dstr%x%r(i) = parts%dstr%x%r(i)+2.0*parts%str%xy%r(i)*parts%wxy%r(i)
        parts%dstr%xy%r(i) = parts%dstr%xy%r(i)-(parts%str%x%r(i)-parts%str%y%r(i))*parts%wxy%r(i)
        parts%dstr%y%r(i) = parts%dstr%y%r(i)-2.0*parts%str%xy%r(i)*parts%wxy%r(i)
    enddo         

    return
    end subroutine

    !------------------------------------------------------------
    subroutine mohr_coulomb_failure_criterion(soil)
    !------------------------------------------------------------
    implicit none

    class(particles) soil
    type(material), pointer :: property
    double precision yield, phi, skale, cohesion, tmax
    integer i, k, ntotal

    ntotal = soil%ntotal+soil%nvirt
    property => soil%material
    cohesion = property%cohesion
    phi      = property%phi

    k = 0
    soil%fail = 0
    do i = 1, ntotal
        tmax  = sqrt(((soil%str%x%r(i)-soil%str%y%r(i))/2)**2+soil%str%xy%r(i)**2)
        if(tmax<1.e-6)cycle
        yield = cohesion*cos(phi)+soil%p%r(i)*sin(phi)

        if(yield<=0.)then    ! <
            yield=0.; soil%p%r(i)=-cohesion*tan(phi)**(-1.0)
            soil%rho%r(i) = property%rho0
        endif

        !         if(yield<=0.)then   ! Collapse
        !            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            p(i)   = 0.
        !         endif

        if(tmax>yield)then
            !         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = yield/tmax
            soil%str%x%r(i) = skale * soil%str%x%r(i)
            soil%str%xy%r(i) = skale * soil%str%xy%r(i)
            soil%str%y%r(i) = skale * soil%str%y%r(i)
        endif
        !sxx(i) = 0.; sxy(i) = 0.; syy(i) = 0.
    enddo
    soil%nfail = k

    return
    end subroutine      

    !------------------------------------------------------------
    subroutine scaling_back(soil)
    !------------------------------------------------------------
    implicit none

    class(particles) soil
    type(material), pointer :: property
    double precision yield, phi, skale, cohesion, tmax, alpha1,I1,J2,kc
    integer i, k, ntotal

    !if(soil%itimestep==82)write(*,*) 'in drucker...'

    ntotal = soil%ntotal+soil%nvirt
    property => soil%material
    cohesion = property%cohesion
    phi      = property%phi
    alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    kc       = 3.*cohesion/sqrt(9.+12.*tan(phi)**2.)

    k = 0
    !soil%fail = 0
    !$omp parallel do private(I1,J2,skale) reduction(+:k)
    do i = 1, ntotal
        !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
        !if(tmax<1.e-6)cycle
        !yield = cohesion*cos(phi)+p(i)*sin(phi)

        if(soil%p%r(i)<-kc/(3.*alpha1))then    ! <
            !k = k + 1
            !soil%fail(i) = 1
            soil%p%r(i)=-kc/(3.*alpha1)
            !soil%rho(i) = property%rho0
        endif

        !         if(yield<=0.)then   ! Collapse
        !            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            p(i)   = 0.
        !         endif

        I1 = 3.*soil%p%r(i) 

        J2 = soil%str%x%r(i)**2.+2.*soil%str%xy%r(i)**2.+soil%str%y%r(i)**2.+(soil%str%x%r(i)+soil%str%y%r(i))**2.
        J2 = sqrt(J2/2.)+1.d-6
        !if(J2<1.e-6)cycle

        if(J2>alpha1*I1+kc)then
            !         if(yield>0.and.tmax>yield)then
            k = k + 1
            !soil%fail(i) = 1 
            skale = (alpha1*I1+kc)/J2
            soil%str%x%r(i) = skale * soil%str%x%r(i)
            soil%str%xy%r(i) = skale * soil%str%xy%r(i)
            soil%str%y%r(i) = skale * soil%str%y%r(i)
        endif

    enddo
    !$omp end parallel do
    !soil%nfail = k

    return
    end subroutine

    !------------------------------------------------------------
    subroutine drucker_prager_failure_criterion(soil)
    !------------------------------------------------------------
    implicit none

    class(particles) soil
    type(material), pointer :: property
    double precision yield, phi, skale, cohesion, tmax, alpha1,I1,J2,kc
    integer i, k, ntotal

    !if(soil%itimestep==82)write(*,*) 'in drucker...'

    ntotal = soil%ntotal+soil%nvirt
    property => soil%material
    cohesion = property%cohesion
    phi      = property%phi
    alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    kc       = 3.*cohesion/sqrt(9.+12.*tan(phi)**2.)

    k = 0
    soil%fail = 0
    !$omp parallel do private(I1,J2,skale) reduction(+:k)
    do i = 1, ntotal
        !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
        !if(tmax<1.e-6)cycle
        !yield = cohesion*cos(phi)+p(i)*sin(phi)

        if(soil%p%r(i)<-kc/(3.*alpha1))then    ! <
            !k = k + 1
            !soil%fail(i) = 1
            soil%p%r(i)=-kc/(3.*alpha1)
            !soil%rho(i) = property%rho0
        endif

        !         if(yield<=0.)then   ! Collapse
        !            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            p(i)   = 0.
        !         endif

        I1 = 3.*soil%p%r(i) 

        J2 = soil%str%x%r(i)**2.+2.*soil%str%xy%r(i)**2.+soil%str%y%r(i)**2.+(soil%str%x%r(i)+soil%str%y%r(i))**2.
        J2 = sqrt(J2/2.)+1.d-6
        !if(J2<1.e-6)cycle

        if(J2>alpha1*I1+kc)then
            !         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = (alpha1*I1+kc)/J2
            soil%str%x%r(i) = skale * soil%str%x%r(i)
            soil%str%xy%r(i) = skale * soil%str%xy%r(i)
            soil%str%y%r(i) = skale * soil%str%y%r(i)
        endif

    enddo
    !$omp end parallel do
    soil%nfail = k

    return
    end subroutine


    !------------------------------------------------------------
    subroutine drucker_prager_failure_criterion_critical_state(soil)
    !------------------------------------------------------------
    implicit none

    class(particles) soil
    type(material), pointer :: property
    double precision yield, phi, skale, cohesion, tmax, alpha1,I1,J2
    integer i, k, ntotal

    !if(soil%itimestep==82)write(*,*) 'in drucker...'

    ntotal = soil%ntotal+soil%nvirt
    property => soil%material
    cohesion = property%cohesion
    phi      = property%phi
    !alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)

    k = 0
    soil%fail = 0
    !$omp parallel do private(I1,J2,skale) reduction(+:k)
    do i = 1, ntotal
        !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
        !if(tmax<1.e-6)cycle
        !yield = cohesion*cos(phi)+p(i)*sin(phi)


        if(soil%p%r(i)<0.)then    ! <
            !k = k + 1
            !soil%fail(i) = 1
            soil%p%r(i)=0.d0
            !soil%rho(i) = property%rho0
        endif

        !         if(yield<=0.)then   ! Collapse
        !            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            p(i)   = 0.
        !         endif

        alpha1   = tan(phi+soil%psi%r(i))/sqrt(9.+12.*tan(phi+soil%psi%r(i))**2.)

        I1 = 3.*soil%p%r(i) 

        J2 = soil%str%x%r(i)**2.+2.*soil%str%xy%r(i)**2.+soil%str%y%r(i)**2.+(soil%str%x%r(i)+soil%str%y%r(i))**2.
        J2 = sqrt(J2/2.)+1.d-6
        !if(J2<1.e-6)cycle

        if(J2>alpha1*I1)then
            !         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = alpha1*I1/J2
            soil%str%x%r(i) = skale * soil%str%x%r(i)
            soil%str%xy%r(i) = skale * soil%str%xy%r(i)
            soil%str%y%r(i) = skale * soil%str%y%r(i)
        endif

    enddo
    !$omp end parallel do
    soil%nfail = k

    return
    end subroutine
    !------------------------------------------------------------
    subroutine plastic_or_not(soil)
    !------------------------------------------------------------
    implicit none

    class(particles) soil
    type(material), pointer :: property
    real(dp) yield, phi, skale, cohesion, tmax, alpha1,I1,J2,kc
    integer i, k, ntotal

    !if(soil%itimestep==82)write(*,*) 'in drucker...'

    ntotal = soil%ntotal+soil%nvirt
    property => soil%material
    cohesion = property%cohesion
    phi      = property%phi
    alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    kc       = 3.*cohesion/sqrt(9.+12.*tan(phi)**2.)

    k = 0
    soil%fail = 0

    !$omp parallel do private(I1,J2) reduction(+:k)
    do i = 1, ntotal
        !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
        !if(tmax<1.e-6)cycle
        !yield = cohesion*cos(phi)+p(i)*sin(phi)

        !if(p(i)<0.)then    ! <
        !   !k = k + 1
        !   !soil%fail(i) = 1
        !   p(i)=0.d0
        !   !soil%rho(i) = property%rho0
        !endif

        !         if(yield<=0.)then   ! Collapse
        !            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            p(i)   = 0.
        !         endif

        I1 = 3.*soil%p%r(i) 

        J2 = soil%str%x%r(i)**2.+2.*soil%str%xy%r(i)**2.+soil%str%y%r(i)**2.+(soil%str%x%r(i)+soil%str%y%r(i))**2.
        J2 = sqrt(J2/2.)+1.d-6
        !if(J2<1.e-6)cycle


        !         if(J2>alpha1*I1)then
        if(abs(J2-alpha1*I1-kc)<1.0e-5)then
            k = k + 1
            soil%fail(i) = 1
            !skale = alpha1*I1/J2
            !sxx(i) = skale * sxx(i)
            !sxy(i) = skale * sxy(i)
            !syy(i) = skale * syy(i)
        endif

    enddo
    !$omp end parallel do
    soil%nfail = k

    return
    end subroutine

    !------------------------------------------------------------
    subroutine plastic_or_not_critical_state(soil)
    !------------------------------------------------------------
    implicit none

    class(particles) soil
    type(material), pointer :: property
    real(dp) yield, phi, skale, cohesion, tmax, alpha1,I1,J2
    integer i, k, ntotal

    !if(soil%itimestep==82)write(*,*) 'in drucker...'

    ntotal = soil%ntotal+soil%nvirt
    property => soil%material
    cohesion = property%cohesion
    phi      = property%phi
    !alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)

    k = 0
    soil%fail = 0

    !$omp parallel do private(I1,J2) reduction(+:k)
    do i = 1, ntotal
        !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
        !if(tmax<1.e-6)cycle
        !yield = cohesion*cos(phi)+p(i)*sin(phi)

        !if(p(i)<0.)then    ! <
        !   !k = k + 1
        !   !soil%fail(i) = 1
        !   p(i)=0.d0
        !   !soil%rho(i) = property%rho0
        !endif

        !         if(yield<=0.)then   ! Collapse
        !            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
        !            p(i)   = 0.
        !         endif

        alpha1   = tan(phi+soil%psi%r(i))/sqrt(9.+12.*tan(phi+soil%psi%r(i))**2.)

        I1 = 3.*soil%p%r(i) 

        J2 = soil%str%x%r(i)**2.+2.*soil%str%xy%r(i)**2.+soil%str%y%r(i)**2.+(soil%str%x%r(i)+soil%str%y%r(i))**2.
        J2 = sqrt(J2/2.)+1.d-6
        !if(J2<1.e-6)cycle

        !         if(J2>alpha1*I1)then
        if(abs(J2-alpha1*I1)<1.0e-5)then
            k = k + 1
            soil%fail(i) = 1 
            !skale = alpha1*I1/J2
            !sxx(i) = skale * sxx(i)
            !sxy(i) = skale * sxy(i)
            !syy(i) = skale * syy(i)
        endif

    enddo
    !$omp end parallel do
    soil%nfail = k

    return
    end subroutine


    !---------------------------------------------------------------------
    subroutine plastic_flow_rule(parts)
    !---------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
    double precision, pointer, dimension(:) :: sxx, sxy, syy
    type(material), pointer :: property
    double precision alpha, phi, K, G,e,niu, J2, sde, dlambda
    double precision exx, exy, eyy                ! total strain rate
    double precision :: small_value = 1.d-10
    integer i, ntotal

    ntotal = parts%ntotal + parts%nvirt

    dsxx => parts%dstr%x%r
    dsxy => parts%dstr%xy%r
    dsyy => parts%dstr%y%r
    vcc  => parts%vcc%r
    sxx => parts%str%x%r
    sxy => parts%str%xy%r
    syy => parts%str%y%r

    property => parts%material
    phi = property%phi
    k   = property%k; e = property%E; niu = property%niu
    alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    G = e/(2.0*(1+niu))
    !$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
    do i = 1, ntotal
        !if(parts%fail(i)==1)then

        exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
        exy = parts%tab%xy%r(i)/2.
        eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

        sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
        J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
        J2 = J2 + small_value

        !G = parts%eta(i)
        dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
        !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
        if(dlambda<0)cycle
        !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
        dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
        dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
        !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
        dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

        !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
        parts%dp%r(i) = parts%dp%r(i) + 3.*k*alpha*dlambda

        ! Accumulative deviatoric strain

        parts%epsilon_p%r(i) = parts%epsilon_p%r(i)    &
        + dlambda*sxy(i)/(2*sqrt(J2))*parts%dt

        !endif ! Fail
    enddo
    !$omp end parallel do
    return
    end subroutine

    !---------------------------------------------------------------------
    subroutine plastic_flow_rule2(parts)
    !---------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision, pointer, dimension(:) :: dsxx2,dsxy2,dsyy2,vcc
    double precision, pointer, dimension(:) :: sxx, sxy, syy
    type(material), pointer :: property
    double precision alpha, phi, K, G,e,niu, J2, sde, dlambda
    double precision exx, exy, eyy                ! total strain rate
    double precision :: small_value = 1.d-10
    integer i, ntotal

    ntotal = parts%ntotal + parts%nvirt

    dsxx2 => parts%dstr2%x%r
    dsxy2 => parts%dstr2%xy%r
    dsyy2 => parts%dstr2%y%r
    vcc  => parts%vcc%r
    sxx => parts%str%x%r
    sxy => parts%str%xy%r
    syy => parts%str%y%r

    property => parts%material
    phi = property%phi
    k   = property%k; e = property%E; niu = property%niu
    alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    G = e/(2.0*(1+niu))
    !$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
    do i = 1, ntotal
        !if(parts%fail(i)==1)then

        exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
        exy = parts%tab%xy%r(i)/2.
        eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

        sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
        J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
        J2 = J2 + small_value

        !G = parts%eta(i)
        dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)

        !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
        dsxx2(i) = -dlambda*(G/sqrt(J2)*sxx(i))
        dsxy2(i) = -dlambda*(G/sqrt(J2)*sxy(i))
        !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
        dsyy2(i) = -dlambda*(G/sqrt(J2)*syy(i))

        !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
        parts%dp2%r(i) =  3.*k*alpha*dlambda

        ! Accumulative deviatoric strain

        !      parts%epsilon_p(i) = parts%epsilon_p(i) 
        !     &                   + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

        !endif ! Fail
    enddo
    !$omp end parallel do

    return
    end subroutine

    !---------------------------------------------------------------------
    subroutine plastic_flow_rule3(parts)
    !---------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
    double precision, pointer, dimension(:) :: sxx, sxy, syy
    type(material), pointer :: property
    double precision alpha, phi, K, G,e,niu, J2, sde, dlambda
    double precision exx, exy, eyy                ! total strain rate
    double precision :: small_value = 1.d-10
    integer i, ntotal

    ntotal = parts%ntotal + parts%nvirt

    dsxx => parts%dstr%x%r
    dsxy => parts%dstr%xy%r
    dsyy => parts%dstr%y%r
    vcc  => parts%vcc%r
    sxx => parts%str%x%r
    sxy => parts%str%xy%r
    syy => parts%str%y%r

    property => parts%material
    phi = property%phi
    k   = property%k; e = property%E; niu = property%niu
    alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    G = e/(2.0*(1+niu))
    !$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
    do i = 1, ntotal
        if(parts%fail(i)==1)then

        exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
        exy = parts%tab%xy%r(i)/2.
        eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

        sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
        J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
        J2 = J2 + small_value

        !G = parts%eta(i)
        dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
        !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
        if(dlambda<0)cycle
        !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
        dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
        dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
        !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
        dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

        !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
        parts%dp%r(i) = parts%dp%r(i) + 3.*k*alpha*dlambda

        ! Accumulative deviatoric strain

        parts%epsilon_p%r(i) = parts%epsilon_p%r(i)       & 
        + dlambda*sxy(i)/(2*sqrt(J2))*parts%dt

        endif ! Fail
    enddo
    !$omp end parallel do

    return
    end subroutine


    !!DEC$IF(.FALSE.)
    !---------------------------------------------------------------------
    subroutine particle_shifting(parts)
    !---------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision :: small_value = 1.d-10
    double precision :: rij,xij1,xij2,vmax = 2.d0
    double precision :: beta = 0.01    
    integer i,j,k, ntotal

    ntotal = parts%ntotal + parts%nvirt

    parts%r0 = 0.d0
    parts%mt = 0.d0

    do  k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)  
        rij = sqrt((parts%x(1,i)-parts%x(1,j))**2+(parts%x(2,i)-parts%x(2,j))**2)
        parts%r0%r(i) = parts%r0%r(i) + rij
        parts%r0%r(j) = parts%r0%r(j) + rij
        parts%mt%r(i) = parts%mt%r(i) + parts%mass%r(j)
        parts%mt%r(j) = parts%mt%r(j) + parts%mass%r(i)
    enddo

    do i = 1, ntotal
        if(parts%countiac(i)==0)cycle
        parts%r0%r(i)=parts%r0%r(i)/parts%countiac(i)
    enddo

    do i = 1, ntotal
        parts%mt%r(i)=parts%mt%r(i)+parts%mass%r(i)
    enddo

    parts%ps%x = 0.d0; parts%ps%y = 0.d0

    do  k=1,parts%niac
        i = parts%pair_i(k)
        j = parts%pair_j(k)  
        rij = sqrt((parts%x(1,i)-parts%x(1,j))**2+(parts%x(2,i)-parts%x(2,j))**2)
        xij1 = parts%x(1,i) - parts%x(1,j)
        xij2 = parts%x(2,i) - parts%x(2,j)
        !if(rij<1.0d-10)cycle
        parts%ps%x%r(i) = parts%ps%x%r(i) + parts%mass%r(j)*xij1/rij**3.0
        parts%ps%y%r(i) = parts%ps%y%r(i) + parts%mass%r(j)*xij2/rij**3.0
        parts%ps%x%r(j) = parts%ps%x%r(j) - parts%mass%r(i)*xij1/rij**3.0
        parts%ps%y%r(j) = parts%ps%y%r(j) - parts%mass%r(i)*xij2/rij**3.0
    enddo

    do i = 1, ntotal
        if(parts%countiac(i)==0)cycle

        if(parts%zone(i)/=100)then    !!! Only particles near nozzle shift
            parts%ps%x%r(i)=0.d0
            parts%ps%y%r(i)=0.d0
            cycle
        endif

        parts%ps%x%r(i) = beta*parts%ps%x%r(i)*parts%r0%r(i)**2.0*vmax/parts%mt%r(i)
        parts%ps%y%r(i) = beta*parts%ps%y%r(i)*parts%r0%r(i)**2.0*vmax/parts%mt%r(i)
    enddo

    return
    end subroutine


    !!DEC$ENDIF
    !---------------------------------------------------------------------
    subroutine non_associated_plastic_flow_rule3(parts)
    !---------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
    double precision, pointer, dimension(:) :: sxx, sxy, syy
    type(material), pointer :: property
    double precision alpha, phi, K, G,e,niu, J2, sde, dlambda, psi
    double precision exx, exy, eyy                ! total strain rate
    double precision :: small_value = 1.d-10
    integer i, ntotal

    ntotal = parts%ntotal + parts%nvirt

    dsxx => parts%dstr%x%r
    dsxy => parts%dstr%xy%r
    dsyy => parts%dstr%y%r
    vcc  => parts%vcc%r
    sxx => parts%str%x%r
    sxy => parts%str%xy%r
    syy => parts%str%y%r

    property => parts%material
    phi = property%phi
    !psi = 0.   !0.01744
    psi = property%psi
    k   = property%k; e = property%E; niu = property%niu
    alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    G = e/(2.0*(1+niu))
    !$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
    do i = 1, ntotal
        if(parts%fail(i)==1)then

        exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
        exy = parts%tab%xy%r(i)/2.
        eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

        sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
        J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
        J2 = J2 + small_value

        !G = parts%eta(i)
        !!      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
        dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(27.*alpha*K*sin(psi)+G)
        !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
        if(dlambda<0)cycle
        !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
        dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
        dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
        !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
        dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

        !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
        !!      parts%dp%r(i) = parts%dp%r(i) + 3.*k*alpha*dlambda
        parts%dp%r(i) = parts%dp%r(i) + 9.*K*sin(psi)*dlambda

        ! Accumulative deviatoric strain

        parts%epsilon_p%r(i) = parts%epsilon_p%r(i)       & 
        + dlambda*sxy(i)/(2*sqrt(J2))*parts%dt

        endif ! Fail
    enddo
    !$omp end parallel do

    return
    end subroutine


    !---------------------------------------------------------------------
    subroutine non_associated_flow_rule_critical_state(parts)
    !---------------------------------------------------------------------
    implicit none

    class(particles) parts
    double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
    double precision, pointer, dimension(:) :: sxx, sxy, syy
    type(material), pointer :: property
    double precision alpha, phi, K, G,e,niu, J2, sde, dlambda, psi
    double precision exx, exy, eyy                ! total strain rate
    double precision :: small_value = 1.d-10
    integer i, ntotal

    ntotal = parts%ntotal + parts%nvirt

    dsxx => parts%dstr%x%r
    dsxy => parts%dstr%xy%r
    dsyy => parts%dstr%y%r
    vcc  => parts%vcc%r
    sxx => parts%str%x%r
    sxy => parts%str%xy%r
    syy => parts%str%y%r

    property => parts%material
    phi = property%phi
    !psi = 0.   !0.01744
    psi = property%psi
    k   = property%k; e = property%E; niu = property%niu
    !alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
    G = e/(2.0*(1+niu))
    !$omp parallel do private(exx,exy,eyy,sde,J2,dlambda)
    do i = 1, ntotal

    if(parts%fail(i)==1)then
        alpha = tan(phi+parts%psi%r(i))/sqrt(9.+12.*tan(phi+parts%psi%r(i))**2.)
        exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.   ! Due to this, this should before Jaumman
        exy = parts%tab%xy%r(i)/2.
        eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3.

        sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
        J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
        J2 = J2 + small_value

        !G = parts%eta(i)
        !!      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)

        !      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(27.*alpha*K*sin(psi)+G)
        dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(27.*alpha*K*sin(parts%psi%r(i))+G)


        !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
        if(dlambda<0)cycle
        !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
        dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
        dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
        !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
        dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

        !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
        !!      parts%dp%r(i) = parts%dp%r(i) + 3.*k*alpha*dlambda

        !      parts%dp%r(i) = parts%dp%r(i) + 9.*K*sin(psi)*dlambda
        parts%dp%r(i) = parts%dp%r(i) + 9.*K*sin(parts%psi%r(i))*dlambda

        ! Accumulative deviatoric strain

        parts%epsilon_p%r(i) = parts%epsilon_p%r(i)       & 
        + dlambda*sxy(i)/(2*sqrt(J2))*parts%dt

    endif ! Fail
    enddo
    !$omp end parallel do

    return
    end subroutine      


    !---------------------------------------------------------------
    subroutine darcy_law(water, soil)
    !---------------------------------------------------------------
    !      use param, only: volume_fraction
    !      use m_particles
    implicit none

    type(particles),target :: water, soil
    double precision dx(3), ks, ns, gw, cf, sp, rrw
    type(p2r) vx_i(3), vx_j(3), dvx_i(3), dvx_j(3)
    type(material), pointer :: h2o,sio2  
    type(numerical), pointer :: numeric
    double precision gravity   
    integer i, j, k, d, dim

    h2o => water%material
    sio2=>  soil%material
    numeric => water%numeric
    gravity = water%gravity
    dim = water%dim

    gw = h2o%rho0*(-gravity); ns = sio2%porosity
    ks = sio2%permeability
    cf = gw*ns/ks
    !cf = 6.e6

    do  k=1,water%niac
        i = water%pair_i(k)
        j = water%pair_j(k)  
        rrw = water%w(k)/(water%rho%r(i)*soil%rho%r(j))

        ! For staturated soil
        !        if(volume_fraction) cf = water%vof(i)*water%rho(i)*(-gravity)/ks
        if(water%volume_fraction) cf = water%vof%r(i)*soil%vof%r(j)*water%rho%r(i)*(-gravity)/ks
        !        cf = water%vof%r(i)*water%rho%r(i)*(-gravity)/ks
        !


        vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
        dvx_i = water%dvx%cmpt(i); dvx_j = soil%dvx%cmpt(j)
        do d=1,dim
            !sp = cf*(water%vx(d,i)-soil%vx(d,j))*rrw
            !water%dvx(d,i) = water%dvx(d,i) - soil%mass(j)*sp
            !soil%dvx(d,j)  = soil%dvx(d,j) + water%mass(i)*sp   

            !Ergun

            cf = 150.*1.e-3*(1-water%vof%r(i))**2/(water%vof%r(i)*0.225e-3**2)+1.75*(1-water%vof%r(i))*water%rho%r(i)*abs(vx_i(d)%p-vx_j(d)%p)/0.225e-3

            sp = cf*(vx_i(d)%p-vx_j(d)%p)*rrw
            dvx_i(d)%p = dvx_i(d)%p - soil%mass%r(j)*sp
            dvx_j(d)%p = dvx_j(d)%p + water%mass%r(i)*sp   
        enddo

        !sp = cf*(water%vx%x%r(i)-soil%vx%x%r(j))*rrw
        !water%dvx%x%r(i) = water%dvx%x%r(i) - soil%mass(j)*sp
        !soil%dvx%x%r(j)  =  soil%dvx%x%r(j) + water%mass(i)*sp   
        !sp = cf*(water%vx%y%r(i)-soil%vx%y%r(j))*rrw
        !water%dvx%y%r(i) = water%dvx%y%r(i) - soil%mass(j)*sp
        !soil%dvx%y%r(j)  =  soil%dvx%y%r(j) + water%mass(i)*sp   

    enddo

    return
    end subroutine 

    !---------------------------------------------------------------
    subroutine darcy_law_omp(water, soil)
    !---------------------------------------------------------------
    !      use param, only: volume_fraction
    !      use m_particles
    implicit none

    type(particles), target :: water, soil
    real(dp), allocatable, dimension(:,:,:) :: wlocal, slocal
    real(dp) dx(3), ks, ns, gw, cf, sp, rrw, gravity
    type(p2r) vx_i(3), vx_j(3), dvx_i(3), dvx_j(3)
    type(material), pointer :: h2o,sio2  
    type(numerical), pointer :: numeric
    integer i, j, k, d, dim, nthreads,it, wntotal, sntotal

    nthreads = water%nthreads
    wntotal = water%ntotal + water%nvirt
    sntotal = soil%ntotal + soil%nvirt

    h2o => water%material
    sio2=>  soil%material
    numeric => water%numeric
    gravity = water%gravity
    dim = water%dim

    gw = h2o%rho0*(-gravity); ns = sio2%porosity
    ks = sio2%permeability
    cf = gw*ns/ks
    !cf = 6.e6

    allocate(wlocal(dim,wntotal,nthreads))
    allocate(slocal(dim,sntotal,nthreads))
    call water%get_niac_start_end

    !$omp parallel
    !$omp do private(i,j,rrw,cf,vx_i,vx_j,d,sp,k) 
    do it = 1,nthreads
        do i = 1, wntotal
            do d = 1,dim
                wlocal(d,i,it) = 0.d0
            enddo  
        enddo
        do j = 1, sntotal
            do d = 1,dim
                slocal(d,j,it) = 0.d0
            enddo  
        enddo          
        do k=water%niac_start(it),water%niac_end(it)
            i = water%pair_i(k)
            j = water%pair_j(k)  
            rrw = water%w(k)/(water%rho%r(i)*soil%rho%r(j))

            ! For staturated soil
            !        if(volume_fraction) cf = water%vof(i)*water%rho(i)*(-gravity)/ks
            if(water%volume_fraction) cf = water%vof%r(i)*soil%vof%r(j)*water%rho%r(i)*(-gravity)/ks
            !        cf = water%vof%r(i)*water%rho%r(i)*(-gravity)/ks

            vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
            !          dvx_i = water%dvx%cmpt(i); dvx_j = soil%dvx%cmpt(j)
            do d=1,dim
                !sp = cf*(water%vx(d,i)-soil%vx(d,j))*rrw
                !water%dvx(d,i) = water%dvx(d,i) - soil%mass(j)*sp
                !soil%dvx(d,j)  = soil%dvx(d,j) + water%mass(i)*sp  

                !Ergun

                !        cf = 150.*1.e-3*(1-water%vof%r(i))**2/(water%vof%r(i)*0.225e-3**2)+1.75*(1-water%vof%r(i))*water%rho%r(i)*abs(vx_i(d)%p-vx_j(d)%p)/0.225e-3

                sp = cf*(vx_i(d)%p-vx_j(d)%p)*rrw
                !dvx_i(d)%p = dvx_i(d)%p - soil%mass%r(j)*sp
                !dvx_j(d)%p = dvx_j(d)%p + water%mass%r(i)*sp 
                wlocal(d,i,it) = wlocal(d,i,it) - soil%mass%r(j)*sp
                slocal(d,j,it) = slocal(d,j,it) + water%mass%r(i)*sp
            enddo
        enddo
    enddo     
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1,wntotal
        !   dvx_i = water%dvx%cmpt(i)
        !      do d=1,dim
        !      dvx_i(d)%p = local_water(d,i)+dvx_i(d)%p
        !      enddo
        do it =1,nthreads
            water%dvx%x%r(i) = water%dvx%x%r(i) + wlocal(1,i,it)
            water%dvx%y%r(i) = water%dvx%y%r(i) + wlocal(2,i,it)
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do j = 1,sntotal
        !   dvx_j = soil%dvx%cmpt(j)
        !      do d=1,dim
        !      dvx_j(d)%p = local_soil(d,j)+dvx_j(d)%p
        !      enddo
        do it = 1, nthreads
            soil%dvx%x%r(j) = soil%dvx%x%r(j) + slocal(1,j,it)
            soil%dvx%y%r(j) = soil%dvx%y%r(j) + slocal(2,j,it)      
        enddo 
    enddo
    !$omp end do
    !$omp end parallel

    return
    end subroutine 


    !-----------------------------------------------------
    subroutine pore_water_pressure(water,soil)
    !-----------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) water, soil
    double precision mprr
    integer i, j, k, d

    do k = 1, water%niac
        i = water%pair_i(k)   ! water
        j = water%pair_j(k)   ! soil
        mprr = water%mass%r(i)*water%p%r(i)/(water%rho%r(i)*soil%rho%r(j))
        !         mprr = water%mass(i)*(water%p(i)+soil%p(j))/       &      Bui2014
        !                (water%rho(i)*soil%rho(j))
        !         do d = 1, water%dim
        !            soil%dvx(d,j) = soil%dvx(d,j) + mprr*water%dwdx(d,k)  &  !+
        !                            *soil%vof%r(j)  
        !         enddo

        soil%dvx%x%r(j) = soil%dvx%x%r(j) + mprr*water%dwdx(1,k)  &  !+
        *soil%vof%r(j)  
        soil%dvx%y%r(j) = soil%dvx%y%r(j) + mprr*water%dwdx(2,k)  &  !+
        *soil%vof%r(j)  

        ! saturated soil
        !if(water%volume_fraction)then
        !do d = 1, dim
        !water%dvx(d,i) = water%dvx(d,i) -                     & ! Must be -
        !soil%mass(j)*water%p(i)*soil%vof(j)*water%dwdx(d,k)/  &
        !(water%rho(i)*soil%rho(j))
        !enddo
        !endif
    enddo

    return
    end subroutine       

    !-----------------------------------------------------
    subroutine pore_water_pressure_omp(water,soil)
    !-----------------------------------------------------
    !      use param
    !      use m_particles
    implicit none
    real(dp), allocatable, dimension(:,:,:) :: local
    type(particles) water, soil
    real(dp) mprr
    integer i, j, k, d, ntotal, it, nthreads, dim

    nthreads = water%nthreads
    ntotal = soil%ntotal+soil%nvirt; dim = water%dim

    allocate(local(dim,ntotal,nthreads))
    call water%get_niac_start_end

    !$omp parallel
    !$omp do private(i,j,mprr,k) 
    do it =1,nthreads
        do j = 1,ntotal
            local(1:dim,j,it) = 0.d0
        enddo 
        do k = water%niac_start(it),water%niac_end(it)
            i = water%pair_i(k)   ! water
            j = water%pair_j(k)   ! soil

            !mprr = water%mass%r(i)*water%p%r(i)/(water%rho%r(i)*soil%rho%r(j)) !!! Original

            ! Symmetric form, Added on 25 July 2016.
            mprr = water%mass%r(i)*(water%p%r(i)-soil%spp%r(j))/(water%rho%r(i)*soil%rho%r(j)) ! Bui2011

            !         mprr = water%mass(i)*(water%p(i)+soil%p(j))/       &      Bui2014
            !                (water%rho(i)*soil%rho(j))
            do d = 1, dim
                !            soil%dvx(d,j) = soil%dvx(d,j) + mprr*water%dwdx(d,k)  &  !+
                !                            *soil%vof%r(j)  
                local(d,j,it) = local(d,j,it) + mprr*water%dwdx(d,k)  &  !+
                *soil%vof%r(j)  
            enddo

            ! saturated soil
            !if(water%volume_fraction)then
            !do d = 1, dim
            !water%dvx(d,i) = water%dvx(d,i) -                     & ! Must be -
            !soil%mass(j)*water%p(i)*soil%vof(j)*water%dwdx(d,k)/  &
            !(water%rho(i)*soil%rho(j))
            !enddo
            !endif
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do j = 1,soil%ntotal+soil%nvirt
        do it = 1,water%nthreads
            soil%dvx%x%r(j) = soil%dvx%x%r(j) + local(1,j,it)
            soil%dvx%y%r(j) = soil%dvx%y%r(j) + local(2,j,it)
        enddo  
    enddo   
    !$omp end do
    !$omp end parallel
    return
    end subroutine 

    !-------------------------------------------------------------
    subroutine volume_fraction_water(water, soil)
    !-------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) water, soil
    integer i,j,k,d, ntotal
    type(material), pointer :: sio2

    sio2 => soil%material
    ntotal = water%ntotal+water%nvirt

    water%vof2 = 0.d0
    do k = 1, water%niac
        i = water%pair_i(k)
        j = water%pair_j(k)
        water%vof2%r(i) = water%vof2%r(i)+soil%mass%r(j)*water%w(k)
    enddo

    do k = 1, water%ntotal+water%nvirt
        water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0
    enddo

    return
    end subroutine 

    !-------------------------------------------------------------
    subroutine volume_fraction_water_omp(water, soil)
    !-------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) water, soil
    integer i,j,k,d, ntotal,nthreads,it
    real(dp), allocatable, dimension(:,:) :: local
    type(material), pointer :: sio2

    sio2 => soil%material
    ntotal = water%ntotal+water%nvirt
    nthreads =water%nthreads

    allocate(local(ntotal,nthreads))
    call water%get_niac_start_end

    water%vof2 = 0.d0

    !$omp parallel
    !$omp do private(i,j,k)
    do it = 1, nthreads 
        do i = 1,ntotal
            local(i,it) = 0.d0
        enddo  
        do k = water%niac_start(it),water%niac_end(it)
            i = water%pair_i(k)
            j = water%pair_j(k)
            !         water%vof2%r(i) = water%vof2%r(i)+soil%mass%r(j)*water%w(k)
            !         local(i,it) = local(i,it)+soil%mass%r(j)*water%w(k) ! due to critical
            !         state
            local(i,it) = local(i,it)+soil%mass%r(j)*soil%vof%r(j)*water%w(k)/soil%rho%r(j)
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do k = 1, water%ntotal+water%nvirt
        !         water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0
        do it= 1, nthreads 
            water%vof2%r(k) = water%vof2%r(k) + local(k,it)
        enddo 
        !water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0 ! due to critical
        !state
        water%vof2%r(k) = 1.d0 - water%vof2%r(k)

    enddo
    !$omp end do
    !$omp end parallel
    return
    end subroutine 


    !-------------------------------------------------------------
    subroutine volume_fraction_water_omp_with_filter(water, soil)
    !-------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) water, soil
    integer i,j,k,d, ntotal,nthreads,it
    real(dp), allocatable, dimension(:,:) :: local, local_wi
    type(material), pointer :: sio2

    sio2 => soil%material
    ntotal = water%ntotal+water%nvirt
    nthreads =water%nthreads

    allocate(local(ntotal,nthreads))
    allocate(local_wi(ntotal,nthreads))
    call water%get_niac_start_end

    water%vof2 = 0.d0
    water%wi%r = 0.d0
    !write(*,*) 'adfdf',size(water%wi%r)

    !$omp parallel
    !$omp do private(i,j,k)
    do it = 1, nthreads 
        do i = 1,ntotal
            local(i,it) = 0.d0
            local_wi(i,it) = 0.d0
        enddo  
        do k = water%niac_start(it),water%niac_end(it)
            i = water%pair_i(k)
            j = water%pair_j(k)
            !         water%vof2%r(i) = water%vof2%r(i)+soil%mass%r(j)*water%w(k)
            local(i,it) = local(i,it)+soil%mass%r(j)*water%w(k)
            local_wi(i,it) = local_wi(i,it) +                     &
            soil%mass%r(j)/soil%rho%r(j)*water%w(k)
        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do k = 1, water%ntotal+water%nvirt
        !         water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0
        do it= 1, nthreads 
            water%vof2%r(k) = water%vof2%r(k) + local(k,it)
            water%wi%r(k) = water%wi%r(k) + local_wi(k,it)
        enddo 
        !water%vof2%r(k) = 1.d0 - water%vof2%r(k)/sio2%rho0
        !write(*,*) 'adfaf'
        if(abs(water%wi%r(k))>1.0d-5)then
            !        write(*,*) 'water%wi%r(k)=', water%wi%r(k), k
            water%vof2%r(k) = 1.d0 - water%vof2%r(k)/water%wi%r(k)/sio2%rho0
        else
            water%vof2%r(k) = 1.d0
        endif
    enddo
    !$omp end do
    !$omp end parallel
    return
    end subroutine 


    !-------------------------------------------------------------
    subroutine volume_fraction_water2(water, soil)
    !-------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) water, soil
    integer i,j,k,d, ntotal
    type(material), pointer :: sio2
    double precision dvx(3),tmp
    type(p2r) vx_i(3), vx_j(3)

    sio2 => soil%material
    ntotal = water%ntotal+water%nvirt

    !      water%dvof = 0.d0

    do k = 1, water%niac
        i = water%pair_i(k)
        j = water%pair_j(k)
        vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
        do d = 1, water%dim
            !dvx(d) = water%vx(d,i)-soil%vx(d,j)
            dvx(d) = vx_i(d)%p - vx_j(d)%p
        enddo 
        tmp = dvx(1)*water%dwdx(1,k)+dvx(2)*water%dwdx(2,k)
        water%dvof%r(i) = water%dvof%r(i)-soil%mass%r(j)*tmp/sio2%rho0
    enddo

    return
    end subroutine 

    !-------------------------------------------------------------
    subroutine volume_fraction_water2_omp(water, soil)
    !-------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) water, soil
    integer i,j,k,d, ntotal,nthreads,it
    type(material), pointer :: sio2
    real(dp), allocatable, dimension(:,:) :: local
    double precision dvx(3),tmp
    type(p2r) vx_i(3), vx_j(3)

    sio2 => soil%material
    ntotal = water%ntotal+water%nvirt
    nthreads = water%nthreads

    allocate(local(ntotal,nthreads))
    call water%get_niac_start_end

    !      water%dvof = 0.d0

    !$omp parallel
    !$omp do private(i,j,k,vx_i,vx_j,d,dvx,tmp)
    do it = 1, nthreads
        do i = 1,ntotal
            local(i,it) = 0.d0
        enddo  
        do k = water%niac_start(it),water%niac_end(it)
            i = water%pair_i(k)
            j = water%pair_j(k)
            vx_i = water%vx%cmpt(i); vx_j = soil%vx%cmpt(j)
            do d = 1, water%dim
                !dvx(d) = water%vx(d,i)-soil%vx(d,j)
                dvx(d) = vx_i(d)%p - vx_j(d)%p
            enddo 
            tmp = dvx(1)*water%dwdx(1,k)+dvx(2)*water%dwdx(2,k)
            !         water%dvof%r(i) = water%dvof%r(i)-soil%mass%r(j)*tmp/sio2%rho0
            !         local(i,it) = local(i,it)-soil%mass%r(j)*tmp/sio2%rho0  ! due to
            !         critical state
            local(i,it) = local(i,it)-soil%mass%r(j)*soil%vof%r(i)*tmp/soil%rho%r(j)

        enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do private(it)
    do i = 1, ntotal
        do it = 1,nthreads
            water%dvof%r(i) = water%dvof%r(i) + local(i,it)
        enddo   
    enddo
    !$omp end do
    !$omp end parallel

    return
    end subroutine 

    !--------------------------------------------------------------------
    subroutine volume_fraction_soil(parts)
    !--------------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) parts
    integer i,j,k, ntotal
    type(material), pointer :: sio2

    sio2 => parts%material
    ntotal = parts%ntotal + parts%nvirt

    do i = 1, ntotal
        parts%vof%r(i) = parts%rho%r(i)/sio2%rho0
    enddo

    return
    end subroutine

    !--------------------------------------------------------------------
    subroutine volume_fraction_soil_omp(parts)
    !--------------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) parts
    integer i,j,k, ntotal
    type(material), pointer :: sio2

    sio2 => parts%material
    ntotal = parts%ntotal + parts%nvirt

    !$omp parallel do
    do i = 1, ntotal
        parts%vof%r(i) = parts%rho%r(i)/sio2%rho0
    enddo
    !$omp end parallel do 

    return
    end subroutine


    !--------------------------------------------------------------------
    subroutine volume_fraction_soil_critical_state(parts)
    !--------------------------------------------------------------------
    !      use param
    !      use m_particles
    implicit none

    type(particles) parts
    integer i,j,k, ntotal
    real(dp) J2, phi_eq
    type(array), pointer :: tab
    type(material), pointer :: sio2

    sio2 => parts%material      

    ntotal = parts%ntotal + parts%nvirt
    tab => parts%tab

    !$omp parallel do private(J2,phi_eq)
    do i = 1, ntotal

    !parts%psi%r(i) = atan(15.09d0*(parts%vof%r(i)-0.60))    !!! 4.09

    J2 = tab%x%r(i)**2.+2.*tab%xy%r(i)**2.+tab%y%r(i)**2.+(tab%x%r(i)+tab%y%r(i))**2.
    J2 = sqrt(J2/2)
    !if(parts%p%r(i)>1.d-3)then
    !   phi_eq = 0.58 - 25.d0*0.012*J2/parts%p%r(i)         !!! 0.60
    !else
    phi_eq = 0.58
    !endif

    parts%psi%r(i) = atan(4.09d0*(parts%vof%r(i)-phi_eq))  !!!K3=4.09,Pouliquen
    parts%dvof%r(i) = -parts%vof%r(i)*tan(parts%psi%r(i))*J2


    !if only theta replaced by theta+psi, i.e. original idea, not critical_state
    !theory
    !         parts%psi%r(i) = -0.35    
    !         parts%dvof%r(i) = -parts%vof%r(i)*tan(parts%psi%r(i))*J2

    enddo
    !$omp end parallel do 

    return
    end subroutine      



    !--------------------------------------------------
    subroutine get_num_threads(this)
    !--------------------------------------------------
    implicit none
    class(particles) this
    integer i,j,k
    integer, external :: OMP_GET_NUM_THREADS

    !$omp parallel
    this%nthreads = OMP_GET_NUM_THREADS()
    !$omp end parallel

    !allocate(this%niac_start(this%nthreads),this%niac_end(this%nthreads))

    end subroutine

    !--------------------------------------------------
    subroutine get_niac_start_end(this)
    !--------------------------------------------------
    implicit none
    class(particles) this
    integer i,j,k,it

    i = this%niac/this%nthreads
    j = mod(this%niac,this%nthreads)

    do it = 1, this%nthreads
        this%niac_start(it) = (it-1)*i + 1
        this%niac_end(it) = it*i
    enddo

    this%niac_end(this%nthreads) = this%niac_end(this%nthreads)+j

    end subroutine

    !--------------------------------------------------

    Subroutine stress_rate(parts)
    !--------------------------------------------------
    implicit none

    class(particles) parts
    type(material),pointer:: property
    double precision,pointer,dimension(:) :: vcc
    double precision:: sxx_trial,sxy_trial,syy_trial,sp_trial
    double precision:: dt,alpha,phi,K,G,E,niu,I1,J2,sde,dlambda, &
    Kc,cohesion,f,I1_trial,J2_trial,f_trial,exx,exy,eyy
    integer ntotal,i

    ntotal = parts%ntotal
    dt     = parts%dt
    vcc    => parts%vcc%r
    property => parts%material
    K = property%k ; E = property%E ; niu = property%niu
    phi = property%phi ; cohesion = property%cohesion
    G = E/(2.d0*(1.d0+niu))
    alpha = tan(phi)/dsqrt(9.d0+12.d0*tan(phi)**2.d0)
    Kc    = 3.d0*cohesion/dsqrt(9.d0+12.d0*tan(phi)**2.d0)

    do i=1, ntotal
        !Elastic stress rate
        parts%dstr%x%r(i)  = G*parts%tab%x%r(i)
        parts%dstr%xy%r(i) = G*parts%tab%xy%r(i)
        parts%dstr%y%r(i)  = G*parts%tab%y%r(i)
        parts%dp%r(i)      = - K*parts%vcc%r(i)

        !Jaumann stress rate
        parts%wxy = 0.5*(parts%df4_omp(parts%vx%x,'y')-parts%df4_omp(parts%vx%y,'x'))
        parts%dstr%x%r(i) = parts%dstr%x%r(i)+2.0*parts%str%xy%r(i)*parts%wxy%r(i)
        parts%dstr%xy%r(i) = parts%dstr%xy%r(i)-(parts%str%x%r(i)-parts%str%y%r(i))*parts%wxy%r(i)
        parts%dstr%y%r(i) = parts%dstr%y%r(i)-2.0*parts%str%xy%r(i)*parts%wxy%r(i)

        !Checking yielding condition of present stress tensor.
        I1 = 3.*parts%p%r(i)
        J2 = 0.5d0*(parts%str%x%r(i)**2.+2.*parts%str%xy%r(i)**2.+parts%str%y%r(i)**2.+(parts%str%x%r(i)+parts%str%y%r(i))**2.)
        f = DSQRT(J2) + alpha*I1-Kc

        !Update stress tensor using current trial stress rate.
        sxx_trial = parts%str%x%r(i)  + parts%dstr%x%r(i)*dt
        sxy_trial = parts%str%xy%r(i) + parts%dstr%xy%r(i)*dt
        syy_trial = parts%str%y%r(i)  + parts%dstr%y%r(i)*dt
        sp_trial = parts%p%r(i) + parts%dp%r(i)*dt
        I1_trial = 3.d0*sp_trial
        J2_trial = 0.5d0*(sxx_trial**2.+2.*sxy_trial**2.+syy_trial**2.+(sxx_trial+syy_trial)**2.)
        f_trial  = DSQRT(J2_trial) + alpha*I1_trial - Kc

        !Check if yield and loading criterion to determine the plastic part of stress rate.
        if ((f>1.d-20).and.(f_trial>1.d-20)) then
            exx = parts%tab%x%r(i)/2.+parts%vcc%r(i)/3.
            exy = parts%tab%xy%r(i)/2.
            eyy = parts%tab%y%r(i)/2.+parts%vcc%r(i)/3. 
            sde = parts%str%x%r(i)*exx+2.*parts%str%xy%r(i)*exy+parts%str%y%r(i)*eyy
            dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
            parts%dstr%x%r(i) = parts%dstr%x%r(i)-dlambda*(G/sqrt(J2)*parts%str%x%r(i))
            parts%dstr%xy%r(i) = parts%dstr%xy%r(i)-dlambda*(G/sqrt(J2)*parts%str%xy%r(i))
            parts%dstr%y%r(i) = parts%dstr%y%r(i)-dlambda*(G/sqrt(J2)*parts%str%y%r(i))
            parts%dp%r(i) = parts%dp%r(i) + 3.*k*alpha*dlambda
        endif
    enddo
    end subroutine    


    end module