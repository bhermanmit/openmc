module tally_new

  use error,            only: fatal_error, warning
  use global
  use output,           only: write_message
  use string
  use tally_class
  use tally_filter_class
  use xml_interface

  implicit none

  contains

!===============================================================================
! READ_TALLIES_NEW reads in tallies from XML to new tally variable
!===============================================================================

  subroutine read_tallies_new()

    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: temp_str
    integer :: i
    integer :: j
    integer :: n_filters     ! number of filters
    integer :: n_words       ! number of words read
    logical :: file_exists
    real(8), allocatable :: real_bins(:)
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_mesh => null()
    type(Node), pointer :: node_tal => null()
    type(Node), pointer :: node_filt => null()
    type(NodeList), pointer :: node_mesh_list => null()
    type(NodeList), pointer :: node_tal_list => null()
    type(NodeList), pointer :: node_filt_list => null()
    type(Tally_p), pointer :: t => null()
    class(TallyFilterClass), pointer :: f => null()

    ! Check if tallies.xml exists
    filename = trim(path_input) // "tallies.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! Since a tallies.xml file is optional, no error is issued here
      return
    end if

    ! Display output message
    message = "Reading tallies XML file again..."
    call write_message(5)

    ! Parse tallies.xml file
    call open_xmldoc(doc, filename)

    ! Get pointer list to XML <tally>
    call get_node_list(doc, "tally", node_tal_list)

    ! Check for user tallies
    n_user_tallies = get_list_size(node_tal_list)
    if (n_user_tallies == 0) then
      message = "No tallies present in tallies.xml file!"
      call warning()
    end if

    ! Allocate tally array
    allocate(tallies_new(n_user_tallies))

    ! Process tally input
    READ_TALLIES: do i = 1, n_user_tallies

      ! Get pointer to tally xml node
      call get_list_item(node_tal_list, i, node_tal)

      ! Set pointer to tallies_new position
      t => tallies_new(i)

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
        t % p  => AnalogTallyClass()
      case ('tracklength', 'track-length', 'pathlength', 'path-length')
        message = "Invalid estimator '" // trim(temp_str) &
             // "' on tally "
        call fatal_error()
      case default
        message = "Invalid estimator '" // trim(temp_str) &
             // "' on tally "
        call fatal_error()
      end select

      ! Get pointer list to XML <filter> and get number of filters
      call get_node_list(node_tal, "filter", node_filt_list)
      n_filters = get_list_size(node_filt_list)

      ! Process filters
      if (n_filters /= 0) then

        ! Allocate filters in tally instance
        call t % p % allocate_filters(n_filters)

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
            call fatal_error()
          end if

          ! Determine type of filter
          select case (temp_str)
          case ('energy')

            ! Allocate an energy filter
            f => EnergyFilterClass()

            ! Read in bins and set to filter
            allocate(real_bins(n_words))
            call get_node_array(node_filt, "bins", real_bins)
            call f % set_bins(n_words, real_bins)
            deallocate(real_bins)
 
          case default
            ! Specified tally filter is invalid, raise error
            message = "Unknown filter type '" // &
                 trim(temp_str) // "' on tally "
            call fatal_error()

          end select

          ! Add filter to tally instance
          call t % p % add_filter(f)

        end do READ_FILTERS
      end if
 
    end do READ_TALLIES

    ! Close xml file
    call close_xmldoc(doc)

  end subroutine read_tallies_new

end module tally_new
