module tally_class

  use constants
  use tally_filter_class 

  implicit none
  private
  public :: Tally_p, TallyClass, AnalogTallyClass

  ! General tally
  type, abstract :: TallyClass
    private
    integer :: type    ! Type of tally from constants
    contains
      procedure, public :: set_type
      procedure (tally_destroy), deferred :: destroy 
  end type TallyClass

  abstract interface
    subroutine tally_destroy(self)
      import TallyClass
      class(TallyClass) :: self
    end subroutine tally_destroy
  end interface

  ! Tally pointer
  type :: Tally_p
    class(TallyClass), pointer :: p => null()
  end type Tally_p

  ! Analog tally 
  type, extends(TallyClass) :: AnalogTallyClass
    private
    contains
      procedure, public :: destroy => analog_tally_destroy
  end type AnalogTallyClass
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
! SET_TYPE sets the member type in TallyClass instance
!===============================================================================

  subroutine set_type(self, type)
    class(TallyClass) :: self
    integer :: type
    self % type = type
  end subroutine set_type

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
