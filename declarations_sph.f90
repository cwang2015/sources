!-----------------------------
    module declarations_sph
!-----------------------------
use m_particles
implicit none

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
type(geo) wasserjet
logical with_nozzle

! Particle size: large, small?
character(len=32) particle_size
character(len=32) operating_mode

logical single_phase
logical dry_soil

! Particle type : water, soil or mixture?
character(len=32) particle_type

type(material), target :: H2O, SiO2
type(numerical), target :: numeric

! Particles 
type(particles), target :: parts, soil

! Working arrays
double precision, pointer, dimension(:) :: txx=>null(), tyy=>null(), txy=>null()

! Physical variables carried by the particles
!double precision, pointer, dimension(:)   ::  rho => null()
!double precision, pointer, dimension(:,:) ::   vx => null()
!double precision, pointer, dimension(:)   ::    p => null()
!double precision, pointer, dimension(:)   ::    c => null()
!double precision, pointer, dimension(:)   ::    t => null()
!double precision, pointer, dimension(:)   ::    u => null()
!double precision, pointer, dimension(:)   ::    s => null()
!double precision, pointer, dimension(:)   ::    e => null()
!double precision, pointer, dimension(:)   ::  eta => null()

! Temporal variables
!double precision, pointer, dimension(:)   ::  drho => null()
!double precision, pointer, dimension(:,:) ::   dvx => null()
!double precision, pointer, dimension(:)   ::    du => null()
!double precision, pointer, dimension(:,:) ::    av => null()

!double precision, pointer, dimension(:,:) :: indvxdt=>null()
!double precision, pointer, dimension(:,:) :: ardvxdt=>null()
!double precision, pointer, dimension(:,:) :: exdvxdt=>null()
!double precision, pointer, dimension(:)   ::  avdudt=>null()
!double precision, pointer, dimension(:)   ::  ahdudt=>null()

!type(array) indvxdt, ardvxdt, exdvxdt, indudt, avdudt, ahdudt

!Numerical paramters
double precision :: dt, time = 0.d0
integer :: maxtimestep = 0 , itimestep = 0

end module
