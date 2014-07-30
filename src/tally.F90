module tally

  use error,              only: fatal_error, warning, write_message
  use global
  use particle_header,    only: Particle
  use string,             only: lower_case
  use tally_class
  use tally_filter_class
  use tally_score_class
  use xml_interface

  implicit none

  character(2*MAX_LINE_LEN) :: message

  contains

!===============================================================================
! SETUP_ACTIVE_USERTALLIES_NEW
!===============================================================================

  subroutine setup_active_usertallies()

    integer :: i ! loop counter

    do i = 1, n_user_tallies
      ! Add tally to active tallies
      call active_tallies % add(i)

      ! Check what type of tally this is and add it to the appropriate list
      select type(t => tallies(i) % p)

      type is (AnalogTallyClass)
        call active_analog_tallies % add(i)

      type is (TracklengthTallyClass)
        call active_tracklength_tallies % add(i)

      type is (CollisionTallyClass)
        call active_collision_tallies % add(i)

      end select
    end do

  end subroutine setup_active_usertallies

!===============================================================================
! SCORE_ANALOG_TALLIES_NEW
!===============================================================================

  subroutine score_analog_tallies(p)

    type(Particle) :: p

    integer :: i
    integer :: i_tally

    ! Loop around tallies and score
    do i = 1, active_analog_tallies % size()
      i_tally = active_analog_tallies % get_item(i)
      call tallies(i_tally) % p % score(p)
    end do

  end subroutine score_analog_tallies

!===============================================================================
! SCORE_TRACKLENGTH_TALLIES
!===============================================================================

  subroutine score_tracklength_tallies(p)

    type(Particle) :: p

    integer :: i
    integer :: i_tally

    ! Loop around tallies and score
    do i = 1, active_tracklength_tallies % size()
      i_tally = active_tracklength_tallies % get_item(i)
      call tallies(i_tally) % p % score(p)
    end do

  end subroutine score_tracklength_tallies

!===============================================================================
! SCORE_COLLISION_TALLIES
!===============================================================================

  subroutine score_collision_tallies(p)

    type(Particle) :: p

    integer :: i
    integer :: i_tally

    ! Loop around tallies and score
    do i = 1, active_collision_tallies % size()
      i_tally = active_collision_tallies % get_item(i)
      call tallies(i_tally) % p % score(p)
    end do

  end subroutine score_collision_tallies

!===============================================================================
! DESTROY_TALLIES frees all memory associated with tallies
!===============================================================================

  subroutine destroy_tallies()

    integer :: i

    do i = 1, n_user_tallies
      call tallies(i) % p % destroy()
    end do

  end subroutine destroy_tallies

!===============================================================================
! TALLY_STATISTICS
!===============================================================================

  subroutine tally_statistics()

    integer :: i

    do i = 1, n_user_tallies
      call tallies(i) % p % statistics
    end do

  end subroutine tally_statistics

!===============================================================================
! RESET_TALLIES
!===============================================================================

  subroutine reset_tallies()

    integer :: i

    do i = 1, n_user_tallies
      call tallies(i) % p % reset()
    end do

  end subroutine reset_tallies

!===============================================================================
! SYNCHRONIZE_TALLIES accumulates the sum of the contributions from each history
! within the batch to a new random variable
!===============================================================================

  subroutine synchronize_tallies()

    integer :: i
    real(8) :: k_col ! Copy of batch collision estimate of keff
    real(8) :: k_abs ! Copy of batch absorption estimate of keff
    real(8) :: k_tra ! Copy of batch tracklength estimate of keff

    ! Accumulate active tallies
    do i = 1, active_tallies % size()
      call tallies(i) % p % accumulate(total_weight)
    end do

    ! Accumulate global tallies
    if (run_mode == MODE_EIGENVALUE) then
      if (active_batches) then
        ! Accumulate products of different estimators of k
        k_col = global_tallies(K_COLLISION) % get_value() / total_weight
        k_abs = global_tallies(K_ABSORPTION) % get_value() / total_weight
        k_tra = global_tallies(K_TRACKLENGTH) % get_value() / total_weight
        k_col_abs = k_col_abs + k_col * k_abs
        k_col_tra = k_col_tra + k_col * k_tra
        k_abs_tra = k_abs_tra + k_abs * k_tra
      end if
    end if

    ! Accumulate global tallies
    call global_tallies % accumulate(total_weight)
    n_realizations = n_realizations + 1

  end subroutine synchronize_tallies

!===============================================================================
! WRITE_TALLIES
!===============================================================================

  subroutine write_tallies()

    integer :: i ! index in tallies array
    character(MAX_FILE_LEN) :: filename ! name of output file

    ! Create filename for tally output
    filename = trim(path_output) // "tallies.out"

    ! Open tally file for writing
    open(FILE=filename, UNIT=UNIT_TALLY, STATUS='replace', ACTION='write')

    ! Loop around tallies and write
    do i = 1, n_user_tallies
      call tallies(i) % p % write_output(nuclides, xs_listings)
    end do

    ! Close tally file
    close(UNIT=UNIT_TALLY)

  end subroutine write_tallies

end module tally
