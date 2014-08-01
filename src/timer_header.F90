module timer_header

  use constants, only: ZERO

  implicit none
  private

!===============================================================================
! TIMER represents a timer that can be started and stopped to measure how long
! different routines run. The intrinsic routine system_clock is used to measure
! time rather than cpu_time.
!===============================================================================

  type, public :: Timer
    private
    logical :: running      = .false. ! is timer running?
    integer :: start_counts = 0       ! counts when started
    real(8), public :: elapsed = ZERO ! total time elapsed in seconds
   contains
     procedure, public :: start     => timer_start
     procedure, public :: get_value => timer_get_value
     procedure, public :: stop      => timer_stop
     procedure, public :: reset     => timer_reset
  end type Timer

  ! Timing global variables
  type(Timer), save, public :: time_total         ! timer for total run
  type(Timer), save, public :: time_initialize    ! timer for initialization
  type(Timer), save, public :: time_read_xs       ! timer for reading cross sections
  type(Timer), save, public :: time_unionize      ! timer for unionizing energy grid
  type(Timer), save, public :: time_bank          ! timer for fission bank synchronization
  type(Timer), save, public :: time_bank_sample   ! timer for fission bank sampling
  type(Timer), save, public :: time_bank_sendrecv ! timer for fission bank SEND/RECV
  type(Timer), save, public :: time_tallies       ! timer for accumulate tallies
  type(Timer), save, public :: time_inactive      ! timer for inactive batches
  type(Timer), save, public :: time_active        ! timer for active batches
  type(Timer), save, public :: time_transport     ! timer for transport only
  type(Timer), save, public :: time_finalize      ! timer for finalization
  type(Timer), save, public :: time_cmfd          ! timer for whole cmfd calculation
  type(Timer), save, public :: time_cmfdbuild     ! timer for matrix build
  type(Timer), save, public :: time_cmfdsolve     ! timer for solver 

contains

!===============================================================================
! TIMER_START starts running a timer and measures the current time
!===============================================================================

  subroutine timer_start(self)

    class(Timer), intent(inout) :: self

    ! Turn timer on and measure starting time
    self % running = .true.
    call system_clock(self % start_counts)

  end subroutine timer_start

!===============================================================================
! TIMER_GET_VALUE returns the current value of the timer
!===============================================================================

  function timer_get_value(self) result(elapsed)

    class(Timer), intent(in) :: self   ! the timer
    real(8)                 :: elapsed ! total elapsed time

    integer :: end_counts   ! current number of counts
    integer :: count_rate   ! system-dependent counting rate
    real    :: elapsed_time ! elapsed time since last start

    if (self % running) then
      call system_clock(end_counts, count_rate)
      elapsed_time = real(end_counts - self % start_counts)/real(count_rate)
      elapsed = self % elapsed + elapsed_time
    else
      elapsed = self % elapsed
    end if

  end function timer_get_value

!===============================================================================
! TIMER_STOP stops the timer and sets the elapsed time
!===============================================================================

  subroutine timer_stop(self)

    class(Timer), intent(inout) :: self

    ! Check to make sure timer was running
    if (.not. self % running) return

    ! Stop timer and add time
    self % elapsed = timer_get_value(self)
    self % running = .false.

  end subroutine timer_stop

!===============================================================================
! TIMER_RESET resets a timer to have a zero value
!===============================================================================

  subroutine timer_reset(self)
    
    class(Timer), intent(inout) :: self

    self % running      = .false.
    self % start_counts = 0
    self % elapsed      = ZERO

  end subroutine timer_reset

end module timer_header
