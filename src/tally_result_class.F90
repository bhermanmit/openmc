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
      procedure, public :: add => add_result
      procedure, public :: get_sum
      procedure, public :: get_sum_sq
      procedure, public :: get_value
      procedure, public :: set_value
      procedure, public :: statistics => statistics_result
      procedure, public :: reset => reset_result
      procedure, public :: write => write_result
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

!===============================================================================
! GET_SUM
!===============================================================================

  function get_sum(self) result(sum_)

    class(TallyResultClass) :: self
    real(8) :: sum_

    sum_ = self % sum

  end function get_sum

!===============================================================================
! GET_SUM_SQ
!===============================================================================

  function get_sum_sq(self) result(sum_sq)

    class(TallyResultClass) :: self
    real(8) :: sum_sq

    sum_sq = self % sum_sq

  end function get_sum_sq

!===============================================================================
! GET_VALUE
!===============================================================================

  elemental function get_value(self) result(val)

    class(TallyResultClass), intent(in) :: self
    real(8) :: val

    val = self % value

 end function get_value

!===============================================================================
! GET_VALUE
!===============================================================================

  elemental subroutine set_value(self, val)

    class(TallyResultClass), intent(inout) :: self
    real(8), intent(in) :: val

    self % value = val

 end subroutine set_value

!===============================================================================
! STATISTICS_RESULT determines the sample mean and the standard deviation of the
! mean for a TallyResult.
!===================`===========================================================

  elemental subroutine statistics_result(self, n)

    class(TallyResultClass), intent(inout) :: self
    integer, intent(in)    :: n

    ! Calculate sample mean and standard deviation of the mean -- note that we
    ! have used Bessel's correction so that the estimator of the variance of the
    ! sample mean is unbiased.

    self % sum    = self % sum/real(n, 8)
    self % sum_sq = sqrt((self % sum_sq/real(n, 8) - self % sum * &
         self % sum) / real(n - 1, 8))

  end subroutine statistics_result

!===============================================================================
! RESET_RESULT zeroes out the value and accumulated sum and sum-squared for a
! single TallyResult.
!===============================================================================

  elemental subroutine reset_result(self)

    class(TallyResultClass), intent(inout) :: self

    self % value    = ZERO
    self % sum      = ZERO
    self % sum_sq   = ZERO

  end subroutine reset_result

!===============================================================================
! ADD_RESULT adds a score to the results array
!===============================================================================

  subroutine add_result(self, score)

    class(TallyResultClass), intent(inout) :: self
    real(8), intent(in) :: score

!$omp critical
    self % value = self % value + score
!$omp end critical

  end subroutine add_result

!===============================================================================
! WRITE_RESULT
!===============================================================================

  subroutine write_result(self, unit)

    class(TallyResultClass), intent(inout) :: self
    integer, intent(in) :: unit

    ! Write result information
    write(unit, *) self % sum, self % sum_sq

  end subroutine write_result 

end module tally_result_class
