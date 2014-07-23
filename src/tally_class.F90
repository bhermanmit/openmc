module tally_class

  use constants
  use tally_filter_class 
  use tally_result_class

  implicit none
  private
  public :: Tally_p, TallyClass, AnalogTallyClass

  ! General tally
  type, abstract :: TallyClass
    private
    integer :: i_filter = 1  ! Current filter
    integer :: n_filters ! Number of filters
    integer :: n_scores ! Number of scores
    integer :: total_score_bins ! Total number of score bins
    integer :: total_filter_bins ! Total number of filter bins
    integer :: type ! Type of tally from constants
    type(TallyFilter_p), allocatable :: filters(:) ! Polymorphic array of filter objects
    type(TallyResultClass), allocatable :: results(:,:) ! Array of result objects
    contains
      procedure, public :: add_filter
      procedure, public :: allocate_filters
      procedure, public :: destroy => tally_destroy
      procedure :: allocate_results
      procedure, public :: set_type
  end type TallyClass

  ! Tally pointer type
  type :: Tally_p
    class(TallyClass), pointer :: p => null()
  end type Tally_p

  ! Analog tally 
  type, extends(TallyClass) :: AnalogTallyClass
    private
  end type AnalogTallyClass

  ! Constructor call
  interface AnalogTallyClass
    module procedure analog_tally_init
  end interface

  contains

!*******************************************************************************
!*******************************************************************************
! General abstract tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! ADD_FILTER adds a filter to tally array
!===============================================================================

  subroutine add_filter(self, filter)

    class(TallyClass), intent(inout) :: self
    class(TallyFilterClass), pointer, intent(in) :: filter

    ! Set filter to array
    self % filters(self % i_filter) % p => filter
    self % i_filter = self % i_filter + 1

  end subroutine add_filter

!===============================================================================
! ALLOCATE_FILTERS allocates the filters array in TallyClass instance
!===============================================================================

  subroutine allocate_filters(self, n_filters)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: n_filters

    self % n_filters = n_filters
    allocate(self % filters(n_filters))
    
  end subroutine allocate_filters

!===============================================================================
! ALLOCATE_RESULTS allocates the results array
!===============================================================================

  subroutine allocate_results(self)

    class(TallyClass), intent(inout) :: self

    allocate(self % results(self % total_score_bins, self % total_filter_bins))
    
  end subroutine allocate_results

!===============================================================================
! TALLY_DESTROY frees all memory used by TallyClass
!===============================================================================

  subroutine tally_destroy(self)

    class(TallyClass), intent(inout) :: self

    integer :: i

    ! Loop through all filters and destroy
    do i = 1, self % n_filters
      call self % filters(i) % p % destroy()
      if (associated(self % filters(i) % p)) nullify(self % filters(i) % p)
    end do

    ! Destroy tally results
    if (allocated(self % results)) deallocate(self % results)

  end subroutine tally_destroy

!===============================================================================
! SET_TYPE sets the member type in TallyClass instance
!===============================================================================

  subroutine set_type(self, type)

    class(TallyClass), intent(inout) :: self
    integer :: type
    self % type = type

  end subroutine set_type

!*******************************************************************************
!*******************************************************************************
! Analog tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! ANALOG_TALLY_INIT initializes an AnalogTallyClass
!===============================================================================

  function analog_tally_init() result(self)

    class(AnalogTallyClass), pointer :: self

    ! Allocate
    allocate(self)

    ! Set the tally type
    call self % set_type(ESTIMATOR_ANALOG)

  end function analog_tally_init

end module tally_class
