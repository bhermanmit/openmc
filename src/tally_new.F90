module tally_new

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
! READ_TALLIES_NEW reads in tallies from XML to new tally variable
!===============================================================================

  subroutine read_tallies_new()

    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    class(TallyClass), pointer :: t => null()
    class(TallyFilterClass), pointer :: f => null()
    class(TallyScoreClass), pointer :: s => null()
    integer :: i
    integer :: int_scalar    ! temporary scalar integer
    integer :: j
    integer :: n_filters     ! number of filters
    integer :: n_words       ! number of words read
    integer :: n_scores      ! number of scores
    integer, allocatable :: int_bins(:)
    logical :: file_exists
    real(8), allocatable :: real_bins(:)
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_filt => null()
    type(Node), pointer :: node_tal => null()
    type(NodeList), pointer :: node_filt_list => null()
    type(NodeList), pointer :: node_tal_list => null()

    ! Check if tallies.xml exists
    filename = trim(path_input) // "tallies_new.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! Since a tallies.xml file is optional, no error is issued here
      return
    end if

    ! Display output message
    message = "Reading tallies XML file again..."
    call write_message(message, 5)

    ! Parse tallies.xml file
    call open_xmldoc(doc, filename)

    ! Get pointer list to XML <tally>
    call get_node_list(doc, "tally", node_tal_list)

    ! Check for user tallies
    n_user_tallies = get_list_size(node_tal_list)
    if (n_user_tallies == 0) then
      message = "No tallies present in tallies.xml file!"
      if (master) call warning(message)
    end if

    ! Allocate tally array
    allocate(tallies_new(n_user_tallies))

    ! Process tally input
    READ_TALLIES: do i = 1, n_user_tallies

      ! Get pointer to tally xml node
      call get_list_item(node_tal_list, i, node_tal)

      ! Check if user specified estimator
      if (check_for_node(node_tal, "estimator")) then
        temp_str = ''
        call get_node_value(node_tal, "estimator", temp_str)
      else
        temp_str = 'tracklength'
      end if

      ! Allocate tally pointer
      select case(trim(temp_str))
      case ('analog')
        tallies_new(i) % p => AnalogTallyClass()
      case ('tracklength', 'track-length', 'pathlength', 'path-length')
        tallies_new(i) % p => TracklengthTallyClass()
      case ('collision')
        tallies_new(i) % p => CollisionTallyClass()
      case default
        message = "Invalid estimator '" // trim(temp_str) &
             // "' on tally "
        call fatal_error(message)
      end select
      t => tallies_new(i) % p

      ! Copy material id
      if (check_for_node(node_tal, "id")) then
        call get_node_value(node_tal, "id", int_scalar)
      else
        message = "Must specify id for tally in tally XML file."
        call fatal_error(message)
      end if
      call t % set_id(int_scalar)

      ! Get pointer list to XML <filter> and get number of filters
      call get_node_list(node_tal, "filter", node_filt_list)
      n_filters = get_list_size(node_filt_list)

      ! Process filters
      if (n_filters /= 0) then

        ! Allocate filters in tally instance
        call t % allocate_filters(n_filters)

        READ_FILTERS: do j = 1, n_filters

          ! Get pointer to filter xml node
          call get_list_item(node_filt_list, j, node_filt)

          ! Convert filter type to lower case
          temp_str = ''
          if (check_for_node(node_filt, "type")) &
            call get_node_value(node_filt, "type", temp_str)
          call lower_case(temp_str)

          ! Determine number of bins
          if (check_for_node(node_filt, "bins")) then
            if (trim(temp_str) == 'energy' .or. &
                trim(temp_str) == 'energyout') then
              n_words = get_arraysize_double(node_filt, "bins")
            else
              n_words = get_arraysize_integer(node_filt, "bins")
            end if
          else
            message = "Bins not set in filter on tally "
            call fatal_error(message)
          end if

          ! Allocate type of filter
          select case (temp_str)

          case ('energy')
            f => EnergyFilterClass()

          case ('energyout')
            f => EnergyOutFilterClass()

          case ('mesh')
            f => MeshFilterClass()

          case default

            ! Specified tally filter is invalid, raise error
            message = "Unknown filter type '" // &
                 trim(temp_str) // "' on tally "
            call fatal_error(message)

          end select

          ! Set up filters
          select type (f)

          type is (EnergyFilterClass)

            ! Read in bins and set to filter
            allocate(real_bins(n_words))
            call get_node_array(node_filt, "bins", real_bins)
            call f % set_bins(n_words - 1, real_bins)
            deallocate(real_bins)

          type is (EnergyOutFilterClass)

            ! Read in bins and set to filter
            allocate(real_bins(n_words))
            call get_node_array(node_filt, "bins", real_bins)
            call f % set_bins(n_words - 1, real_bins)
            deallocate(real_bins)

          type is (MeshFilterClass)

            ! Check to make sure only one mesh index is specified
            if (n_words /= 1) then
              message = "Only one mesh can be specified in a filter."
              call fatal_error(message)
            end if

            ! Read in mesh index
            allocate(int_bins(1))
            call get_node_array(node_filt, "bins", int_bins)

            ! Get mesh index from id
            int_bins(1) = mesh_dict % get_key(int_bins(1))

            ! Set up mesh bins
            call f % set_bins(meshes(int_bins(1)))
            deallocate(int_bins)

          end select

          ! Add filter to tally instance
          call t % add_filter(f)

        end do READ_FILTERS
      end if

      ! Process tally scores
      if (check_for_node(node_tal, "scores")) then

        ! Read string array from input
        n_words = get_arraysize_string(node_tal, "scores")
        allocate(sarray(n_words))
        call get_node_array(node_tal, "scores", sarray)

        ! Handle moment stuff here in future

        ! Allocate scores
        n_scores = n_words
        call t % allocate_scores(n_scores)

        ! Loop arounds scores and initialize
        READ_SCORES: do j = 1, n_scores

          call lower_case(sarray(j))
          select case(sarray(j))
          case('total')

            ! Allocate a total score
            s => TotalScoreClass() 

          case('nu-fission')

            ! Allocate a nu-fission score
            s => NuFissionScoreClass()

          case default

            ! Specified tally score is invalid, raise error
            message = "Unknown score type '" // &
                 trim(sarray(j)) // "' on tally "
            call fatal_error(message)
          
          end select

          ! Add score to tally instance
          call t % add_score(s)

        end do READ_SCORES

        ! Deallocate temporary string array of scores
        deallocate(sarray)

      else

        ! Error if there are no scores specified for a tally
        message = "No <scores> specified on tally "
        call fatal_error(message)

      end if

      ! Finish tally setup
      call t % finish_setup()

    end do READ_TALLIES

    ! Close xml file
    call close_xmldoc(doc)

  end subroutine read_tallies_new

!===============================================================================
! SETUP_ACTIVE_USERTALLIES_NEW
!===============================================================================

  subroutine setup_active_usertallies_new()

    integer :: i ! loop counter

    do i = 1, n_user_tallies
      ! Add tally to active tallies
      call active_tallies_new % add(i)

      ! Check what type of tally this is and add it to the appropriate list
      select type(t => tallies_new(i) % p)

      type is (AnalogTallyClass)
        call active_analog_tallies_new % add(i)

      type is (TracklengthTallyClass)
        call active_tracklength_tallies_new % add(i)

      type is (CollisionTallyClass)
        call active_collision_tallies_new % add(i)

      end select
    end do

  end subroutine setup_active_usertallies_new

!===============================================================================
! SCORE_ANALOG_TALLIES_NEW
!===============================================================================

  subroutine score_analog_tallies_new(p)

    type(Particle) :: p

    integer :: i
    integer :: i_tally

    ! Loop around tallies and score
    do i = 1, active_analog_tallies_new % size()
      i_tally = active_analog_tallies_new % get_item(i)
      call tallies_new(i_tally) % p % score(p)
    end do

  end subroutine score_analog_tallies_new

!===============================================================================
! SCORE_TRACKLENGTH_TALLIES_NEW
!===============================================================================

  subroutine score_tracklength_tallies_new(p)

    type(Particle) :: p

    integer :: i
    integer :: i_tally

    ! Loop around tallies and score
    do i = 1, active_tracklength_tallies_new % size()
      i_tally = active_tracklength_tallies_new % get_item(i)
      call tallies_new(i_tally) % p % score(p)
    end do

  end subroutine score_tracklength_tallies_new

!===============================================================================
! SCORE_COLLISION_TALLIES_NEW
!===============================================================================

  subroutine score_collision_tallies_new(p)

    type(Particle) :: p

    integer :: i
    integer :: i_tally

    ! Loop around tallies and score
    do i = 1, active_collision_tallies_new % size()
      i_tally = active_collision_tallies_new % get_item(i)
      call tallies_new(i_tally) % p % score(p)
    end do

  end subroutine score_collision_tallies_new

!===============================================================================
! DESTROY_TALLIES_NEW frees all new tallies memory
!===============================================================================

  subroutine destroy_tallies_new()

    integer :: i

    do i = 1, n_user_tallies
      call tallies_new(i) % p % destroy()
    end do

  end subroutine destroy_tallies_new

!===============================================================================
! TALLY_NEW_STATISTICS
!===============================================================================

  subroutine tally_new_statistics()

    integer :: i

    do i = 1, n_user_tallies
      call tallies_new(i) % p % statistics
    end do

  end subroutine tally_new_statistics

!===============================================================================
! RESET_NEW_TALLIES
!===============================================================================

  subroutine reset_new_tallies()

    integer :: i

    do i = 1, n_user_tallies
      call tallies_new(i) % p % reset()
    end do

  end subroutine reset_new_tallies

!===============================================================================
! WRITE_NEW_TALLIES
!===============================================================================

  subroutine write_new_tallies()

    integer :: i            ! index in tallies array
    character(MAX_FILE_LEN) :: filename                    ! name of output file

    ! Create filename for tally output
    filename = trim(path_output) // "tallies_new.out"

    ! Open tally file for writing
    open(FILE=filename, UNIT=UNIT_TALLY, STATUS='replace', ACTION='write')

    ! Loop around tallies and write
    do i = 1, n_user_tallies
      call tallies_new(i) % p % write(UNIT_TALLY)
    end do

    ! Close tally file
    close(UNIT=UNIT_TALLY)

  end subroutine write_new_tallies

end module tally_new
