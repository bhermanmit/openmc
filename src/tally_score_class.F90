module tally_score_class

  use ace_header,      only: material_xs, micro_xs
  use constants
  use particle_header, only: Particle

  implicit none
  private

  ! General tally score type
  type, abstract, public :: TallyScoreClass
    private
    integer :: type
    integer :: moment_order = 0
    real(8), pointer :: flux
    real(8), pointer :: response
    contains
      procedure, public :: get_type => get_score_type
      procedure, public :: get_moment_order
      procedure, public :: set_moment_order
      procedure, public :: write => write_score
      procedure(score_match_interface), deferred, nopass :: score_match
      procedure(get_weight_interface), deferred, nopass :: get_weight
      procedure(get_response_interface), deferred, nopass :: get_response
  end type TallyScoreClass

  ! Abstract interface for deferred procedures
  abstract interface
    function score_match_interface(p) result(match)
      import Particle
      type(Particle) :: p
      logical :: match
    end function score_match_interface
    function get_response_interface(i_nuclide) result(response)
      integer :: i_nuclide
      real(8) :: response
    end function get_response_interface
    function get_weight_interface(p) result(weight)
      import Particle
      type(Particle) :: p
      real(8) :: weight
    end function get_weight_interface
  end interface

  ! Tally score pointer
  type, public :: TallyScore_p
    class(TallyScoreClass), pointer :: p => null()
  end type TallyScore_p

  ! Total score type
  type, extends(TallyScoreClass), public :: TotalScoreClass
    private
    contains
      procedure, public, nopass :: score_match => total_score_match
      procedure, public, nopass :: get_response => total_get_response
      procedure, public, nopass :: get_weight => total_get_weight
  end type TotalScoreClass
  interface TotalScoreClass
    module procedure total_score_init
  end interface

  ! Nu-fission score type
  type, extends(TallyScoreClass), public :: NuFissionScoreClass
    private
    contains
      procedure, public, nopass :: score_match => nufission_score_match
      procedure, public, nopass :: get_response => nufission_get_response
      procedure, public, nopass :: get_weight => nufission_get_weight
  end type NuFissionScoreClass
  interface NuFissionScoreClass
    module procedure nufission_score_init
  end interface

  contains

!*******************************************************************************
!*******************************************************************************
! General abstract tally score methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! GET_SCORE_TYPE
!===============================================================================

  function get_score_type(self) result(type)

    class(TallyScoreClass) :: self
    integer :: type

    type = self % type

  end function get_score_type

!===============================================================================
! GET_MOMENT_ORDER
!===============================================================================

  function get_moment_order(self) result(moment_order)

    class(TallyScoreClass) :: self
    integer :: moment_order

    moment_order = self % moment_order

  end function get_moment_order

!===============================================================================
! SET_MOMENT_ORDER
!===============================================================================

  subroutine set_moment_order(self, moment_order)

    class(TallyScoreClass), intent(inout) :: self
    integer, intent(in) :: moment_order

    self % moment_order = moment_order

  end subroutine set_moment_order

!===============================================================================
! TALLY_SCORE_DESTROY deallocates all members of TallyScoreClass
!===============================================================================

  subroutine tally_score_destroy(self)

    class(TallyScoreClass), intent(inout) :: self

  end subroutine tally_score_destroy

!===============================================================================
! WRITE_SCORE writes out score information to a unit
!===============================================================================

  subroutine write_score(self, unit)

    class(TallyScoreClass), intent(inout) :: self
    integer, intent(in) :: unit

    write(unit, *) "    Type:", self % type

  end subroutine write_score

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
    self % type = SCORE_TOTAL

  end function total_score_init

!===============================================================================
! TOTAL_SCORE_MATCH results a true such that this tally is always scored
!===============================================================================

  function total_score_match(p) result(match)

    type(Particle) :: p
    logical :: match

    ! Should evaluate to true always
    match = (p % event == EVENT_SCATTER) .or. (p % event /= EVENT_SCATTER)

  end function total_score_match

!===============================================================================
! TOTAL_GET_RESPONSE gets the total macro xs
!===============================================================================

  function total_get_response(i_nuclide) result(response)

    integer :: i_nuclide
    real(8) :: response

    if (i_nuclide == MATERIAL_TOTAL) then
      response = material_xs % total
    else
      response = micro_xs(i_nuclide) % total
    end if

  end function total_get_response

!===============================================================================
! TOTAL_GET_WEIGHT returns the weight for a TotalScoreClass instance
!===============================================================================

  function total_get_weight(p) result(weight)

    type(Particle) :: p
    real(8) :: weight

    weight = p % last_wgt

  end function total_get_weight

!*******************************************************************************
!*******************************************************************************
! Nu-fission score methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! NUFISSION_SCORE_INIT allocates and sets up a NuFissionScoreClass instance
!===============================================================================

  function nufission_score_init() result(self)

    class(NuFissionScoreClass), pointer :: self

    ! Create object
    allocate(self)

    ! Set type of filter
    self % type = SCORE_NU_FISSION

  end function nufission_score_init

!===============================================================================
! NUFISSION_SCORE_MATCH checks for an implicit fission event
!===============================================================================

  function nufission_score_match(p) result(match)

    type(Particle) :: p
    logical :: match

    ! Get fission logical from particle
    match = p % fission

  end function nufission_score_match

!===============================================================================
! NUFISSION_GET_RESPONSE gets the nu-fission macro xs
!===============================================================================

  function nufission_get_response(i_nuclide) result(response)

    integer :: i_nuclide
    real(8) :: response

    if (i_nuclide == MATERIAL_TOTAL) then
      response = material_xs % nu_fission
    else
      response = micro_xs(i_nuclide) % nu_fission
    end if

  end function nufission_get_response

!===============================================================================
! NUFISSION_GET_WEIGHT returns the weight for a NuFissionScoreClass instance
!===============================================================================

  function nufission_get_weight(p) result(weight)

    type(Particle) :: p
    real(8) :: weight

    weight = p % wgt

  end function nufission_get_weight

end module tally_score_class
