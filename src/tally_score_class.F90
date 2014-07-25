module tally_score_class

  use constants
  use particle_header, only: Particle

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
      procedure(score_match_interface), deferred :: score_match
      procedure(get_score_interface), deferred :: get_score
  end type TallyScoreClass

  ! Abstract interface for deferred procedures
  abstract interface
    function score_match_interface(self, event) result(match)
      import TallyScoreClass
      class(TallyScoreClass) :: self
      integer :: event
      logical :: match
    end function score_match_interface
    function get_score_interface(self, p) result(score)
      import TallyScoreClass
      import Particle
      class(TallyScoreClass) :: self
      type(Particle) :: p
      real(8) :: score
    end function get_score_interface
  end interface

  ! Tally score pointer
  type :: TallyScore_p
    class(TallyScoreClass), pointer :: p => null()
  end type TallyScore_p

  ! Total score type
  type, extends(TallyScoreClass) :: TotalScoreClass
    private
    contains
      procedure, public :: score_match => total_score_match
      procedure, public :: get_score => total_get_score
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

!===============================================================================
! TOTAL_SCORE_MATCH results a true such that this tally is always scored
!===============================================================================

  function total_score_match(self, event) result(match)

    class(TotalScoreClass) :: self
    integer :: event
    logical :: match

    ! Should evaluate to true always
    match = (event == EVENT_SCATTER) .or. (event /= EVENT_SCATTER)

  end function total_score_match

!===============================================================================
! TOTAL_GET_SCORE returns the score for a TotalScoreClass instance
!===============================================================================

  function total_get_score(self, p) result(score)

    class(TotalScoreClass) :: self
    type(Particle) :: p
    real(8) :: score

    score = p % last_wgt

  end function total_get_score

end module tally_score_class
