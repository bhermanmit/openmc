module tally_filter_class

  use constants
  use mesh_header,     only: StructuredMesh
  use mesh,            only: get_mesh_bin
  use particle_header, only: Particle
  use search,          only: binary_search

  implicit none
  private

  ! General tally filter type
  type, abstract, public :: TallyFilterClass
    private
    character(len=MAX_WORD_LEN) :: type    ! Type of filter
    integer :: n_bins
    integer, allocatable :: int_bins(:)
    real(8), allocatable :: real_bins(:)
    contains
      procedure, public :: get_n_bins
      procedure, public :: destroy => tally_filter_destroy
      procedure, public :: write => write_filter 
      procedure(filter_index_interface), deferred :: get_filter_index
  end type TallyFilterClass

  ! Interface for deferred procedures
  abstract interface
    function filter_index_interface(self, p) result(filter_index)
      import TallyFilterClass
      import Particle
      class(TallyFilterClass) :: self
      type(Particle) :: p
      integer :: filter_index
    end function filter_index_interface
  end interface

  ! Tally filter pointer
  type, public :: TallyFilter_p
    class(TallyFilterClass), pointer :: p => null()
  end type TallyFilter_p

  ! Energy filter
  type, extends(TallyFilterClass), public :: EnergyFilterClass
    private
    contains
      procedure, public :: get_filter_index => energy_filter_get_index
      procedure, public :: set_bins => energy_filter_set_bins
  end type EnergyFilterClass
  interface EnergyFilterClass
    module procedure energy_filter_init
  end interface

  ! Energy-out filter
  type, extends(TallyFilterClass), public :: EnergyOutFilterClass
    private
    contains
      procedure, public :: get_filter_index => energyout_filter_get_index
      procedure, public :: set_bins => energyout_filter_set_bins
  end type EnergyOutFilterClass
  interface EnergyOutFilterClass
    module procedure energyout_filter_init
  end interface

  ! Mesh filter
  type, extends(TallyFilterClass), public :: MeshFilterClass
    private
    type(StructuredMesh), pointer :: mesh
    contains
      procedure, public :: get_filter_index => mesh_filter_get_index
      procedure, public :: set_bins => mesh_filter_set_bins
  end type MeshFilterClass
  interface MeshFilterClass
    module procedure mesh_filter_init
  end interface

  contains

!*******************************************************************************
!*******************************************************************************
! General abstract tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! GET_N_BINS returns the number of bins for TallyFilterClass instance
!===============================================================================

  function get_n_bins(self) result(n_bins)

    class(TallyFilterClass) :: self
    integer :: n_bins

    n_bins = self % n_bins

  end function get_n_bins

!===============================================================================
! SET_INT_BINS allocates and sets filter bins that are integers
!===============================================================================

  subroutine set_int_bins(self, n_bins, bins)

    class(TallyFilterClass), intent(inout) :: self
    integer, intent(in) :: n_bins
    integer, intent(in) :: bins(:)

    self % n_bins = n_bins
    allocate(self % int_bins(n_bins))
    self % int_bins = bins

  end subroutine set_int_bins

!===============================================================================
! SET_REAL_BINS allocates and sets filter bins that are reals
!===============================================================================

  subroutine set_real_bins(self, n_bins, bins)

    class(TallyFilterClass), intent(inout) :: self
    integer, intent(in) :: n_bins
    real(8), intent(in) :: bins(:)

    self % n_bins = n_bins
    allocate(self % real_bins(n_bins))
    self % real_bins = bins

  end subroutine set_real_bins

!===============================================================================
! TALLY_FILTER_DESTROY deallocates all members of TallyFilterClass
!===============================================================================

  subroutine tally_filter_destroy(self)

    class(TallyFilterClass), intent(inout) :: self

    if (allocated(self % int_bins)) deallocate(self % int_bins)
    if (allocated(self % real_bins)) deallocate(self % real_bins)

  end subroutine tally_filter_destroy

!===============================================================================
! WRITE_FILTER
!===============================================================================

  subroutine write_filter(self, unit)

    class(TallyFilterClass), intent(inout) :: self
    integer :: unit

    ! Write filter information
    write(unit, *) "    Type:", trim(self % type)
    write(unit, *) "    Number of bins:", self % n_bins
    if (allocated(self % int_bins)) write(unit, *) "    BINS:", self % int_bins
    if (allocated(self % real_bins)) write(unit, *) "    BINS:", self % real_bins

  end subroutine write_filter

!*******************************************************************************
!*******************************************************************************
! Energy filter methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! ENERGY_FILTER_INIT allocates and sets up an EnergyFilterClass instance
!===============================================================================

  function energy_filter_init() result(self)

    class(EnergyFilterClass), pointer :: self

    ! Create object
    allocate(self)

    ! Set type of filter
    self % type = 'energyin'

  end function energy_filter_init

!===============================================================================
! ENERGY_FILTER_GET_INDEX returns the index for an energy filter
!===============================================================================

  function energy_filter_get_index(self, p) result(filter_index)

    class(EnergyFilterClass) :: self
    type(Particle) :: p
    integer :: filter_index

    ! perform binary search
    filter_index = binary_search(self % real_bins, self % n_bins + 1, &
                   p % last_E)

  end function energy_filter_get_index

!===============================================================================
! ENERGY_FILTER_SET_BINS is a special routine to set energy filter bins
!===============================================================================

  subroutine energy_filter_set_bins(self, n_bins, bins)

    class(EnergyFilterClass), intent(inout) :: self
    integer, intent(in) :: n_bins
    real(8), intent(in) :: bins(:)

    ! Set and allocate bins
    self % n_bins = n_bins
    allocate(self % real_bins(n_bins + 1))
    self % real_bins = bins

  end subroutine energy_filter_set_bins

!*******************************************************************************
!*******************************************************************************
! Energy-out filter methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! ENERGYOUT_FILTER_INIT allocates and sets up an EnergyOutFilterClass instance
!===============================================================================

  function energyout_filter_init() result(self)

    class(EnergyOutFilterClass), pointer :: self

    ! Create object
    allocate(self)

    ! Set type of filter
    self % type = 'energy-out'

  end function energyout_filter_init

!===============================================================================
! ENERGYOUT_FILTER_GET_INDEX returns the index for an energy-out filter
!===============================================================================

  function energyout_filter_get_index(self, p) result(filter_index)

    class(EnergyOutFilterClass) :: self
    type(Particle) :: p
    integer :: filter_index

    ! perform binary search
    filter_index = binary_search(self % real_bins, self % n_bins + 1, &
                   p % E)

  end function energyout_filter_get_index

!===============================================================================
! ENERGYOUT_FILTER_SET_BINS is a special routine to set energy-out filter bins
!===============================================================================

  subroutine energyout_filter_set_bins(self, n_bins, bins)

    class(EnergyOutFilterClass), intent(inout) :: self
    integer, intent(in) :: n_bins
    real(8), intent(in) :: bins(:)

    ! Set and allocate bins
    self % n_bins = n_bins
    allocate(self % real_bins(n_bins + 1))
    self % real_bins = bins

  end subroutine energyout_filter_set_bins

!*******************************************************************************
!*******************************************************************************
! Mesh filter methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! MESH_FILTER_INIT allocates and sets up a MeshFilterClass instance
!===============================================================================

  function mesh_filter_init() result(self)

    class(MeshFilterClass), pointer :: self

    ! Create object
    allocate(self)

    ! Set type of filter
    self % type = 'mesh'

  end function mesh_filter_init

!===============================================================================
! MESH_FILTER_GET_INDEX returns the index for a mesh filter
!===============================================================================

  function mesh_filter_get_index(self, p) result(filter_index)

    class(MeshFilterClass) :: self
    type(Particle) :: p
    integer :: filter_index

    integer :: bin_offset

    ! Get mesh filter bin
    call get_mesh_bin(self % mesh, p % coord0 % xyz, filter_index)

  end function mesh_filter_get_index

!===============================================================================
! MESH_FILTER_SET_BINS is a special routine to set mesh filter bins
!===============================================================================

  subroutine mesh_filter_set_bins(self, mesh)

    class(MeshFilterClass), intent(inout) :: self
    type(StructuredMesh), target, intent(in) :: mesh

    integer :: n_bins

    ! Calculate total number of bins in mesh
    self % n_bins = product(mesh % dimension)

    ! Associate mesh
    self % mesh => mesh

  end subroutine mesh_filter_set_bins

end module tally_filter_class
