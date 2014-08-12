module tally_filter_class

  use constants
  use mesh_header,     only: StructuredMesh
  use mesh,            only: get_mesh_bin, mesh_intersects_2d, &
                             mesh_intersects_3d, get_mesh_indices, &
                             mesh_indices_to_bin
  use particle_header, only: Particle
  use search,          only: binary_search

  implicit none
  private

  ! General tally filter type
  type, abstract, public :: TallyFilterClass
    private
    integer :: type
    integer :: n_bins
    integer, allocatable :: int_bins(:)
    real(8), allocatable :: real_bins(:)
    contains
      procedure, public :: get_n_bins
      procedure, public :: get_int_bin
      procedure, public :: get_real_bin
      procedure, public :: set_int_bin
      procedure, public :: get_real_bins_pointer
      procedure, public :: get_int_bins_pointer
      procedure, public :: get_type => get_tally_filter_type
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
    integer :: n_cross ! number of mesh crossings
    integer, allocatable :: bins(:)
    real(8), allocatable :: dists(:)
    type(StructuredMesh), pointer :: mesh
    contains
      procedure, public :: get_filter_index => mesh_filter_get_index
      procedure, public :: get_mesh_pointer => mesh_filter_get_pointer
      procedure, public :: calculate_distances => mesh_filter_calc_distances
      procedure, public :: get_distance => mesh_filter_get_distance
      procedure, public :: get_mesh_bin => mesh_filter_get_bin
      procedure, public :: get_n_cross => mesh_filter_get_n_cross
      procedure, public :: set_bins => mesh_filter_set_bins
  end type MeshFilterClass
  interface MeshFilterClass
    module procedure mesh_filter_init
  end interface

  ! Surface filter
  type, extends(TallyFilterClass), public :: SurfaceFilterClass
    private
    contains
      procedure, public :: get_filter_index => surface_filter_get_index
      procedure, public :: set_bins => surface_filter_set_bins
  end type SurfaceFilterClass
  interface SurfaceFilterClass
    module procedure surface_filter_init
  end interface

  contains

!*******************************************************************************
!*******************************************************************************
! General abstract tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! GET_TALLY_FILTER_TYPE returns the type of tally filter
!===============================================================================

  function get_tally_filter_type(self) result(filter_type)

    class(TallyFilterClass) :: self
    integer :: filter_type 

    filter_type = self % type

  end function get_tally_filter_type

!===============================================================================
! GET_N_BINS returns the number of bins for TallyFilterClass instance
!===============================================================================

  function get_n_bins(self) result(n_bins)

    class(TallyFilterClass) :: self
    integer :: n_bins

    n_bins = self % n_bins

  end function get_n_bins

!===============================================================================
! SET_INT_BIN
!===============================================================================

  subroutine set_int_bin(self, idx, val)

    class(TallyFilterClass), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(in) :: val

    self % int_bins(idx) = val

  end subroutine set_int_bin

!===============================================================================
! GET_INT_BIN
!===============================================================================

  function get_int_bin(self, bin) result(int_bin)

    class(TallyFilterClass) :: self
    integer :: bin
    integer :: int_bin

    int_bin = self % int_bins(bin)

  end function get_int_bin

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
! GET_REAL_BIN
!===============================================================================

  function get_real_bin(self, bin) result(real_bin)

    class(TallyFilterClass) :: self
    integer :: bin
    real(8) :: real_bin

    real_bin = self % real_bins(bin)

  end function get_real_bin

!===============================================================================
! GET_INT_BINS_POINTER
!===============================================================================

  function get_int_bins_pointer(self) result(int_bins)

    class(TallyFilterClass), target :: self
    integer, pointer :: int_bins(:)

    int_bins => self % int_bins

  end function get_int_bins_pointer

!===============================================================================
! GET_REAL_BINS_POINTER
!===============================================================================

  function get_real_bins_pointer(self) result(real_bins)

    class(TallyFilterClass), target :: self
    real(8), pointer :: real_bins(:)

    real_bins => self % real_bins

  end function get_real_bins_pointer

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
    write(unit, *) "    Type:", self % type
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
    self % type = FILTER_ENERGYIN

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
    self % type = FILTER_ENERGYOUT

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
    self % type = FILTER_MESH

  end function mesh_filter_init

!===============================================================================
! MESH_FILTER_CALC_DISTANCES
!===============================================================================

  subroutine mesh_filter_calc_distances(self, p, dists, bins, n_cross)

    class(MeshFilterClass), intent(inout) :: self
    type(Particle), intent(in) :: p
    real(8), intent(inout), pointer :: dists(:)
    integer, intent(inout), pointer :: bins(:)
    integer, intent(inout) :: n_cross

    integer :: j                    ! loop index for dimension
    integer :: k                    ! loop index for mesh cell crossings
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: ijk_cross(3)         ! indices of mesh cell crossed
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross(3)         ! coordinates of next boundary
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: distance             ! distance traveled in mesh cell
    logical :: found_bin            ! was a scoring bin found?
    logical :: start_in_mesh        ! starting coordinates inside mesh?
    logical :: end_in_mesh          ! ending coordinates inside mesh?
    type(StructuredMesh), pointer, save :: m => null()

    ! Reset n_cross
    n_cross = 0

    ! Copy starting and ending location of particle
    xyz0 = p % coord0 % xyz - (p % dist - TINY_BIT) * p % coord0 % uvw
    xyz1 = p % coord0 % xyz  - TINY_BIT * p % coord0 % uvw

    ! Determine indices for starting and ending location
    m => self % mesh 
    call get_mesh_indices(m, xyz0, ijk0(:m % n_dimension), start_in_mesh)
    call get_mesh_indices(m, xyz1, ijk1(:m % n_dimension), end_in_mesh)

    ! Check if start or end is in mesh -- if not, check if track still
    ! intersects with mesh
    if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
      if (m % n_dimension == 2) then
        if (.not. mesh_intersects_2d(m, xyz0, xyz1)) return
      else
        if (.not. mesh_intersects_3d(m, xyz0, xyz1)) return
      end if
    end if

    ! Reset starting and ending location
    xyz0 = p % coord0 % xyz - p % dist * p % coord0 % uvw
    xyz1 = p % coord0 % xyz

    ! Calculate number of surface crossings
    n_cross = sum(abs(ijk1(:m % n_dimension) - ijk0(:m % n_dimension))) + 1

    ! Allocate dists and bins and set to 0
    allocate(dists(n_cross))
    allocate(bins(n_cross))
    dists = ZERO
    bins = 0

    ! Copy particle's direction
    uvw = p % coord0 % uvw

    ! Bounding coordinates
    do j = 1, m % n_dimension
      if (uvw(j) > 0) then
        xyz_cross(j) = m % lower_left(j) + ijk0(j) * m % width(j)
      else
        xyz_cross(j) = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
      end if
    end do

    MESH_LOOP: do k = 1, n_cross
      found_bin = .false.

      ! Calculate distance to each bounding surface. We need to treat special
      ! case where the cosine of the angle is zero since this would result in a
      ! divide-by-zero.

      if (k == n_cross) xyz_cross = xyz1

      do j = 1, m % n_dimension
        if (uvw(j) == 0) then
          d(j) = INFINITY
        else
          d(j) = (xyz_cross(j) - xyz0(j))/uvw(j)
        end if
      end do

      ! Determine the closest bounding surface of the mesh cell by calculating
      ! the minimum distance

      j = minloc(d(:m % n_dimension), 1)
      distance = d(j)

      ! Now use the minimum distance and direction of the particle to determine
      ! which surface was crossed

      if (all(ijk0(:m % n_dimension) >= 1) .and. all(ijk0(:m % n_dimension) <= m % dimension)) then
        ijk_cross = ijk0
        found_bin = .true.
      end if

      ! Increment indices and determine new crossing point
      if (uvw(j) > 0) then
        ijk0(j) = ijk0(j) + 1
        xyz_cross(j) = xyz_cross(j) + m % width(j)
      else
        ijk0(j) = ijk0(j) - 1
        xyz_cross(j) = xyz_cross(j) - m % width(j)
      end if

      ! Save values
      if (found_bin) then
        dists(k) = distance
        bins(k) = mesh_indices_to_bin(m, ijk_cross)
      end if

      ! Calculate new coordinates
      xyz0 = xyz0 + distance * uvw

    end do MESH_LOOP

  end subroutine mesh_filter_calc_distances

!===============================================================================
! MESH_FILTER_GET_DISTANCE
!===============================================================================

  function mesh_filter_get_distance(self, k) result(distance)

    class(MeshFilterClass) :: self
    integer :: k
    real(8) :: distance

    distance = self % dists(k)

  end function mesh_filter_get_distance

!===============================================================================
! MESH_FILTER_GET_BIN
!===============================================================================

  function mesh_filter_get_bin(self, k) result(bin)

    class(MeshFilterClass) :: self
    integer :: k
    integer :: bin

    bin = self % bins(k)

  end function mesh_filter_get_bin

!===============================================================================
! MESH_FILTER_GET_N_CROSS
!===============================================================================

  function mesh_filter_get_n_cross(self) result(n_cross)

    class(MeshFilterClass) :: self
    integer :: n_cross

    n_cross = self % n_cross

  end function mesh_filter_get_n_cross

!===============================================================================
! MESH_FILTER_GET_INDEX returns the index for a mesh filter
!===============================================================================

  function mesh_filter_get_index(self, p) result(filter_index)

    class(MeshFilterClass) :: self
    type(Particle) :: p
    integer :: filter_index

    ! Get mesh filter bin
    call get_mesh_bin(self % mesh, p % coord0 % xyz, filter_index)

  end function mesh_filter_get_index

!===============================================================================
! MESH_FILTER_SET_BINS is a special routine to set mesh filter bins
!===============================================================================

  subroutine mesh_filter_set_bins(self, mesh, surface)

    class(MeshFilterClass), intent(inout) :: self
    type(StructuredMesh), target, intent(in) :: mesh
    logical, intent(in) :: surface

    ! Calculate total number of bins in mesh
    if (surface) then
      self % n_bins = product(mesh % dimension + 1)
    else
      self % n_bins = product(mesh % dimension)
    end if

    ! Associate mesh
    self % mesh => mesh

    ! Set mesh as bins
    allocate(self % int_bins(1))
    self % int_bins(1) = self % mesh % id

  end subroutine mesh_filter_set_bins

!===============================================================================
! MESH_FILTER_GET_POINTER
!===============================================================================

  function mesh_filter_get_pointer(self) result(m)

    class(MeshFilterClass) :: self
    type(StructuredMesh), pointer :: m

    m => self % mesh

  end function mesh_filter_get_pointer

!*******************************************************************************
!*******************************************************************************
! Surface filter methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! SURFACE_FILTER_INIT allocates and sets up a SurfaceFilterClass instance
!===============================================================================

  function surface_filter_init() result(self)

    class(SurfaceFilterClass), pointer :: self

    ! Create object
    allocate(self)

    ! Set type of filter
    self % type = FILTER_SURFACE

  end function surface_filter_init

!===============================================================================
! SURFACE_FILTER_SET_BINS
!===============================================================================

  subroutine surface_filter_set_bins(self, bins, n_bins)

    class(SurfaceFilterClass), intent(inout) :: self
    integer :: bins(:)
    integer :: n_bins

    self % n_bins = n_bins
    allocate(self % int_bins(n_bins))
    self % int_bins = bins

  end subroutine surface_filter_set_bins

!===============================================================================
! SURFACE_FILTER_GET_INDEX returns the index for a surface filter
!===============================================================================

  function surface_filter_get_index(self, p) result(filter_index)

    class(SurfaceFilterClass) :: self
    type(Particle) :: p
    integer :: filter_index

  end function surface_filter_get_index

end module tally_filter_class
