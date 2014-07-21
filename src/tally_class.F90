module tally_class

  use constants
  use tally_filter_class 

  implicit none
  private
  public :: TallyClass

  ! General tally
  type, abstract :: TallyClass
    private
    integer :: type    ! Type of tally from constants
    contains
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

  contains

!===============================================================================
! ANALOG_TALLY_DESTROY frees all memory used by AnalogTallyClass
!===============================================================================

  subroutine analog_tally_destroy(self)
    class(AnalogTallyClass) :: self
  end subroutine analog_tally_destroy

end module tally_class
