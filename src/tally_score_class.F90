module tally_score_class

  use constants

  implicit none
  private
  public :: TallyScore_p, TallyScoreClass, TotalScoreClass

  ! General tally score type
  type, abstract :: TallyScoreClass
    private
    integer :: type    ! Type of score from constants
    contains
      procedure, public :: set_type
      procedure, public :: get_type
  end type TallyScoreClass

  ! Tally score pointer
  type :: TallyScore_p
    class(TallyScoreClass), pointer :: p => null()
  end type TallyScore_p

  ! Total score type
  type, extends(TallyScoreClass) :: TotalScoreClass
    private
  end type TotalScoreClass
  interface TotalScoreClass
    module procedure total_score_init
  end interface

  contains

!*******************************************************************************
!*******************************************************************************
! General abstract tally score methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! SET_TYPE sets the member type in TallyScoreClass instance
!===============================================================================

  subroutine set_type(self, type)

    class(TallyScoreClass), intent(inout) :: self
    integer, intent(in) :: type

    ! Set the type to instance
    self % type = type

  end subroutine set_type

!===============================================================================
! GET_TYPE returns the member type from TallyScoreClass instance
!===============================================================================

  function get_type(self) result(type)

    class(TallyScoreClass) :: self
    integer :: type

    ! Get the type from instance
    type = self % type

  end function get_type

!===============================================================================
! TALLY_SCORE_DESTROY deallocates all members of TallyScoreClass
!===============================================================================

  subroutine tally_score_destroy(self)

    class(TallyScoreClass), intent(inout) :: self

  end subroutine tally_score_destroy

!*******************************************************************************
!*******************************************************************************
! Total score methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! TOTAL_SCORE_INIT allocates and sets up a TotalScoreClass instance
!===============================================================================

  function total_score_init() result(self)

    class(TotalScoreClass), pointer :: self

    ! Create object
    allocate(self)

    ! Set type of filter
    call self % set_type(SCORE_TOTAL)

  end function total_score_init

end module tally_score_class
