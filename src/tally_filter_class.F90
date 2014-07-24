module tally_filter_class

  use constants
  use search,          only: binary_search
  use particle_header

  implicit none
  private

  ! General tally filter type
  type, abstract, public :: TallyFilterClass
    private
    integer :: type    ! Type of filter from constants
    integer :: n_bins
    integer, allocatable :: int_bins(:)
    real(8), allocatable :: real_bins(:)
    contains
      procedure, public :: get_n_bins
      procedure, public :: get_type
      procedure, public :: set_type
      procedure :: set_int_bins
      procedure :: set_real_bins 
      generic, public :: set_bins => set_int_bins, set_real_bins
      procedure, public :: destroy => tally_filter_destroy
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
      procedure, public :: set_real_bins => energy_filter_set_bins
  end type EnergyFilterClass
  interface EnergyFilterClass
    module procedure energy_filter_init
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
! SET_TYPE sets the member type in TallyClass instance
!===============================================================================

  subroutine set_type(self, type)

    class(TallyFilterClass), intent(inout) :: self
    integer, intent(in) :: type

    ! Set the type to instance
    self % type = type

  end subroutine set_type

!===============================================================================
! GET_TYPE returns the member type from TallyClass instance
!===============================================================================

  function get_type(self) result(type)

    class(TallyFilterClass) :: self
    integer :: type

    ! Get the type from instance
    type = self % type

  end function get_type

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
    call self % set_type(FILTER_ENERGYIN)

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

    class(EnergyFilterClass) :: self
    integer :: n_bins
    real(8) :: bins(:)

    ! Set and allocate bins
    self % n_bins = n_bins
    allocate(self % real_bins(n_bins + 1))
    self % real_bins = bins

  end subroutine energy_filter_set_bins

end module tally_filter_class
