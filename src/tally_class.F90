module tally_class

  use constants
  use tally_filter_class 

  implicit none
  private
  public :: Tally_p, TallyClass, AnalogTallyClass

  ! General tally
  type, abstract :: TallyClass
    private
    integer :: i_filter = 1  ! Current filter
    integer :: n_filters ! Number of filters
    integer :: type ! Type of tally from constants
    type(TallyFilter_p), allocatable :: filters(:) ! Array of filter objects
    contains
      procedure, public :: add_filter
      procedure, public :: allocate_filters
      procedure, public :: set_type
      procedure (tally_destroy), deferred :: destroy 
  end type TallyClass

  ! Interface for deferred procedures
  abstract interface
    subroutine tally_destroy(self)
      import TallyClass
      class(TallyClass) :: self
    end subroutine tally_destroy
  end interface

  ! Tally pointer type
  type :: Tally_p
    class(TallyClass), pointer :: p => null()
  end type Tally_p

  ! Analog tally 
  type, extends(TallyClass) :: AnalogTallyClass
    private
    contains
      procedure, public :: destroy => analog_tally_destroy
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
! ALLOCATE_FILTERS allocates the filters array in TallyClass instance
!===============================================================================

  subroutine allocate_filters(self, n_filters)

    class(TallyClass) :: self
    integer :: n_filters
    self % n_filters = n_filters
    allocate(self % filters(n_filters))
    
  end subroutine allocate_filters

!===============================================================================
! ADD_FILTER adds a filter to tally array
!===============================================================================

  subroutine add_filter(self, filter)

    class(TallyClass) :: self
    class(TallyFilterClass), pointer :: filter
    self % filters(self % i_filter) % p => filter
    self % i_filter = self % i_filter + 1

  end subroutine add_filter

!===============================================================================
! SET_TYPE sets the member type in TallyClass instance
!===============================================================================

  subroutine set_type(self, type)

    class(TallyClass) :: self
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

!===============================================================================
! ANALOG_TALLY_DESTROY frees all memory used by AnalogTallyClass
!===============================================================================

  subroutine analog_tally_destroy(self)
    class(AnalogTallyClass) :: self
  end subroutine analog_tally_destroy

end module tally_class
