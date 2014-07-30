module particle_header

  use ace_header,      only: MaterialMacroXS, NuclideMicroXS, Nuclide, &
                             SAlphaBeta 
  use bank_header,     only: Bank
  use constants,       only: NEUTRON, ONE, NONE, ZERO
  use geometry_header, only: BASE_UNIVERSE
  use material_header, only: Material

  implicit none

!===============================================================================
! LOCALCOORD describes the location of a particle local to a single
! universe. When the geometry consists of nested universes, a particle will have
! a list of coordinates in each level
!===============================================================================

  type LocalCoord
    ! Indices in various arrays for this level
    integer :: cell      = NONE
    integer :: universe  = NONE
    integer :: lattice   = NONE
    integer :: lattice_x = NONE
    integer :: lattice_y = NONE
    integer :: lattice_z = NONE

    ! Particle position and direction for this level
    real(8) :: xyz(3)
    real(8) :: uvw(3)

    ! Is this level rotated?
    logical :: rotated = .false.

    ! Pointer to next (more local) set of coordinates
    type(LocalCoord), pointer :: next => null()
  end type LocalCoord

!===============================================================================
! PARTICLE describes the state of a particle being transported through the
! geometry
!===============================================================================

  type Particle
    ! Basic data
    integer(8) :: id            ! Unique ID
    integer    :: type          ! Particle type (n, p, e, etc)

    ! Particle coordinates
    type(LocalCoord), pointer :: coord0 => null() ! coordinates on universe 0
    type(LocalCoord), pointer :: coord  => null() ! coordinates on lowest universe

    ! Other physical data
    real(8)    :: wgt           ! particle weight
    real(8)    :: E             ! energy
    real(8)    :: mu            ! angle of scatter
    logical    :: alive         ! is particle alive?

    ! Pre-collision physical data
    real(8)    :: last_xyz(3)   ! previous coordinates
    real(8)    :: last_uvw(3)   ! previous direction coordinates
    real(8)    :: last_wgt      ! pre-collision particle weight
    real(8)    :: last_E        ! pre-collision energy
    real(8)    :: absorb_wgt    ! weight absorbed for survival biasing

    ! What event last took place
    logical    :: fission       ! did the particle cause implicit fission
    integer    :: event         ! scatter, absorption
    integer    :: event_nuclide ! index in nuclides array
    integer    :: event_MT      ! reaction MT

    ! How far did the neutron travel
    real(8), pointer :: dist

    ! Pointer to global memory
    type(Bank), pointer :: fission_bank(:) => null()
    type(MaterialMacroXS), pointer :: material_xs => null()
    type(NuclideMicroXS), pointer :: micro_xs(:) => null()
    type(Nuclide), pointer :: nuclides(:) => null()
    type(SAlphaBeta), pointer :: sab_tables(:) => null()
    logical, pointer :: survival_biasing => null()
    real(8), pointer :: weight_survive => null()
    real(8), pointer :: weight_cutoff => null()

    ! Post-collision physical data
    integer :: nu ! number of fission sites banked
    integer(8), pointer :: n_bank ! current number of fission sites banked
    real(8) :: wgt_bank ! weight of fission sites banked
    real(8), pointer :: keff ! keff used in physics module

    ! Indices for various arrays
    integer    :: surface       ! index for surface particle is on
    integer    :: cell_born     ! index for cell particle was born in
    type(Material), pointer :: material => null() ! current material
    type(Material), pointer :: last_material => null() ! current material

    ! Statistical data
    integer    :: n_collision   ! # of collisions

    ! Track output
    logical    :: write_track = .false.

  contains
    procedure :: initialize => initialize_particle
    procedure :: clear => clear_particle
    procedure :: copy => copy_particle 
    generic   :: assignment(=) => copy
  end type Particle

contains

!===============================================================================
! DEALLOCATE_COORD removes all levels of coordinates below a given level. This
! is used in distance_to_boundary when the particle moves from a lower universe
! to a higher universe since the data for the lower one is not needed anymore.
!===============================================================================

  recursive subroutine deallocate_coord(coord)

    type(LocalCoord), pointer :: coord

    if (associated(coord)) then
      ! recursively deallocate lower coordinates
      if (associated(coord % next)) call deallocate_coord(coord%next)

      ! deallocate original coordinate
      deallocate(coord)
    end if

  end subroutine deallocate_coord

!===============================================================================
! INITIALIZE_PARTICLE sets default attributes for a particle from the source
! bank
!===============================================================================

  subroutine initialize_particle(this)

    class(Particle) :: this

    ! Clear coordinate lists
    call this % clear()

    ! Set particle to neutron that's alive
    this % type  = NEUTRON
    this % alive = .true.

    ! clear attributes
    this % surface       = NONE
    this % cell_born     = NONE
    this % material      => null()
    this % last_material => null()
    this % wgt           = ONE
    this % last_wgt      = ONE
    this % absorb_wgt    = ZERO
    this % nu            = 0
    this % wgt_bank      = ZERO
    this % n_collision   = 0
    this % fission       = .false.

    ! Set up base level coordinates
    allocate(this % coord0)
    this % coord0 % universe = BASE_UNIVERSE
    this % coord             => this % coord0

  end subroutine initialize_particle

!===============================================================================
! COPY_PARTICLE copies all particle attributes, however uses pointers
! that associate to pointers in Particle. Therefore, do not use clear_particle
! when finished. This is currently used in tallies when performing fake
! collisions. 
!===============================================================================

  subroutine copy_particle(this, p)

    class(Particle), intent(inout) :: this
    class(Particle), intent(in) :: p

    ! Copy attributes 
    this % id = p % id
    this % type = p % type
    this % coord0 => p % coord0
    this % coord => p % coord
    this % wgt = p % wgt
    this % E = p % E
    this % mu = p % mu
    this % alive = p % alive
    this % last_xyz = p % last_xyz
    this % last_uvw = p % last_uvw
    this % last_wgt = p % last_wgt
    this % last_E = p % last_E
    this % dist => p % dist
    this % material_xs => p % material_xs
    this % micro_xs => p % micro_xs
    this % nuclides => p % nuclides
    this % sab_tables => p % sab_tables
    this % surface = p % surface
    this % cell_born = p % cell_born
    this % material => p % material
    this % last_material => p % last_material

  end subroutine copy_particle

!===============================================================================
! CLEAR_PARTICLE frees all memory associated with a particle
!===============================================================================

  subroutine clear_particle(this)

    class(Particle) :: this

    ! remove any coordinate levels
    call deallocate_coord(this % coord0)

    ! Make sure coord pointer is nullified
    nullify(this % coord)

  end subroutine clear_particle

end module particle_header
