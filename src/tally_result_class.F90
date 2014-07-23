module tally_result_class

  use constants

  implicit none

  type TallyResultClass
    private
    real(8) :: value    = ZERO
    real(8) :: sum      = ZERO
    real(8) :: sum_sq   = ZERO
    contains
      procedure, public :: accumulate => accumulate_result
  end type TallyResultClass

  contains

!===============================================================================
! ACCUMULATE_RESULT accumulates results from many histories (or many generations)
! into a single realization of a random variable.
!===============================================================================

  elemental subroutine accumulate_result(self, total_weight)

    class(TallyResultClass), intent(inout) :: self
    real(8), intent(in) :: total_weight

    real(8) :: val

    ! Add the sum and square of the sum of contributions from a tally result to
    ! the variables sum and sum_sq. This will later allow us to calculate a
    ! variance on the tallies.

    val = self % value/total_weight
    self % sum    = self % sum    + val
    self % sum_sq = self % sum_sq + val*val

    ! Reset the single batch estimate
    self % value = ZERO

  end subroutine accumulate_result

end module tally_result_class
