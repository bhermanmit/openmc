module tally_filter_class

  use constants

  implicit none
  private
  public :: TallyFilter_p, TallyFilterClass, EnergyFilterClass

  ! General tally filter type
  type, abstract :: TallyFilterClass
    private
    integer :: type    ! Type of filter from constants
    integer :: n_bins
    integer, allocatable :: int_bins(:)
    real(8), allocatable :: real_bins(:)
    contains
      procedure, public :: set_type
      procedure, public :: get_type
      generic, public :: set_bins => set_int_bins, set_real_bins
      procedure :: set_int_bins
      procedure :: set_real_bins 
      procedure, public :: destroy => tally_filter_destroy
  end type TallyFilterClass

  ! Tally filter pointer
  type :: TallyFilter_p
    class(TallyFilterClass), pointer :: p => null()
  end type TallyFilter_p

  ! Energy filter
  type, extends(TallyFilterClass) :: EnergyFilterClass
    private
!     procedure, public :: destroy => energy_filter_destroy
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
! SET_TYPE sets the member type in TallyClass instance
!===============================================================================

  subroutine set_type(self, type)

    class(TallyFilterClass) :: self
    integer :: type

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

    class(TallyFilterClass) :: self
    integer :: n_bins
    integer :: bins(:)

    self % n_bins = n_bins
    allocate(self % int_bins(n_bins))
    self % int_bins = bins

  end subroutine set_int_bins

!===============================================================================
! SET_REAL_BINS allocates and sets filter bins that are reals
!===============================================================================

  subroutine set_real_bins(self, n_bins, bins)

    class(TallyFilterClass) :: self
    integer :: n_bins
    real(8) :: bins(:)

    self % n_bins = n_bins
    allocate(self % real_bins(n_bins))
    self % real_bins = bins
    print *, "In real bins"
  end subroutine set_real_bins

!===============================================================================
! TALLY_FILTER_DESTROY deallocates all members of TallyFilterClass
!===============================================================================

  subroutine tally_filter_destroy(self)

    class(TallyFilterClass) :: self

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

end module tally_filter_class
