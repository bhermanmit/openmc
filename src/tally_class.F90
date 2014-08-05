module tally_class

  use ace_header,         only: Nuclide, Reaction, XSListing, nuclides, &
                                xs_listings, material_xs, n_nuclides_total
  use bank_header,        only: fission_bank, n_bank
  use endf,               only: reaction_name
  use constants
  use dict_header,        only: DictIntInt
  use material_header,    only: Material, materials
  use mesh_header,        only: StructuredMesh
  use mesh,               only: bin_to_mesh_indices
  use particle_header
  use physics,            only: sample_nuclide, sample_fission, &
                                sample_fission_energy, keff
  use random_lcg,         only: prn_set_stream, STREAM_TRACKING, &
                                STREAM_TALLIES
  use set_header,         only: SetInt
  use string,             only: to_str, upper_case
  use tally_filter_class 
  use tally_result_class
  use tally_score_class 

  implicit none
  private

  ! General tally
  type, abstract, public :: TallyClass
    private
    character(len=MAX_LINE_LEN) :: label = ""
    integer :: estimator ! type of tally estimator
    integer :: id ! ID of tally
    integer :: i_filter = 1  ! Current filter
    integer :: i_score = 1  ! Current score
    integer :: n_filters = 0 ! Number of filters
    integer :: n_realizations = 0 ! Number of tally realizations
    integer :: n_scores = 0 ! Number of scores
    integer :: total_nuclide_bins = 1 ! Total number of nuclide bins
    integer :: total_score_bins = 0 ! Total number of score bins
    integer :: total_filter_bins = 1 ! Total number of filter bins
    integer :: find_filter(N_FILTER_TYPES)
    integer, allocatable :: stride(:)
    integer, allocatable :: matching_bins(:)
    integer, allocatable :: nuclides(:)
    logical :: has_eout_filter = .false.
    logical :: has_mesh_filter = .false.
    logical :: score_all_nuclides = .false.
    type(MeshFilterClass), pointer :: mesh_filter => null()
    type(TallyFilter_p), allocatable :: filters(:) ! Polymorphic array of filter objects
    type(TallyResultClass), allocatable :: results(:,:,:) ! Array of result objects
    type(TallyScore_p), allocatable :: scores(:) ! Polymorphic array of filter objects
    contains
      procedure, public :: accumulate => accumulate_results
      procedure, public :: add_filter
      procedure, public :: add_score
      procedure, public :: allocate_filters
      procedure         :: allocate_results
      procedure, public :: allocate_scores
      procedure, public :: allocate_nuclides
      procedure, public :: destroy => tally_destroy
      procedure, public :: finish_setup
      procedure         :: get_filter
      procedure         :: get_filter_index
      procedure         :: get_label
      procedure         :: get_n_filters
      procedure         :: get_total_score_bins
      procedure         :: get_total_filter_bins
      procedure         :: get_total_nuclide_bins
      procedure, public :: get_results_pointer
      procedure, public :: get_n_realizations
      procedure         :: set_filter_index
      procedure         :: setup_filter_indices
      procedure, public :: get_nuclide_bin
      procedure, public :: set_nuclide_bin
      procedure, public :: reset => tally_reset
      procedure         :: setup_stride
      procedure, public :: set_id
      procedure, public :: get_id
      procedure, public :: set_score_all_nuclides
      procedure, public :: statistics => tally_statistics
      procedure, public :: confidence => tally_confidence
      procedure, public :: write => write_tally
      procedure, public :: write_output => write_tally_output
      procedure(score_interface), deferred :: score
  end type TallyClass

  ! Interface for deferred procedure
  abstract interface
    subroutine score_interface(self, p)
      import TallyClass
      import Particle
      class(TallyClass), intent(inout) :: self
      type(Particle), target, intent(in) :: p
    end subroutine score_interface
  end interface

  ! Tally pointer type
  type, public :: Tally_p
    class(TallyClass), pointer :: p => null()
  end type Tally_p

  ! Analog tally 
  type, extends(TallyClass), public :: AnalogTallyClass
    private
    contains
      procedure :: analog_nufission_eout
      procedure, public :: score => analog_tally_score
  end type AnalogTallyClass

  ! Constructor call for analog tally
  interface AnalogTallyClass
    module procedure analog_tally_init
  end interface

  ! Tracklength tally 
  type, extends(TallyClass), public :: TracklengthTallyClass
    private
    contains
      procedure, nopass :: get_flux => tracklength_get_flux
      procedure, public :: score => tracklength_tally_score
  end type TracklengthTallyClass

  ! Constructor call for tracklength tally
  interface TracklengthTallyClass
    module procedure tracklength_tally_init
  end interface

  ! Collision tally 
  type, extends(TallyClass), public :: CollisionTallyClass
    private
    contains
      procedure, nopass :: get_flux => collision_get_flux
      procedure, public :: score => collision_tally_score
  end type CollisionTallyClass

  ! Constructor call for collision tally
  interface CollisionTallyClass
    module procedure collision_tally_init
  end interface

  !====================================
  ! Global variables

  ! Main tally arrays
  integer, save, public :: n_tallies      = 0 ! # of tallies
  integer, save, public :: n_cmfd_tallies = 3 ! # of cmfd tallies
  integer, save, public :: n_user_tallies = 0 ! # of user tallies
  real(8), save, public :: total_weight       ! total starting particle weight in realization
  type(DictIntInt), save, public :: tally_dict
  type(Tally_p), save, public, allocatable, target :: tallies(:)
  type(Tally_p), save, public, pointer :: user_tallies(:) => null()
  type(Tally_p), save, public, pointer :: cmfd_tallies(:) => null()

  ! Location of first tally in tallies
  integer, save, public :: i_user_tallies = -1
  integer, save, public :: i_cmfd_tallies = -1

  ! Active tally lists
  type(SetInt), save, public :: active_analog_tallies
  type(SetInt), save, public :: active_tracklength_tallies
  type(SetInt), save, public :: active_collision_tallies
  type(SetInt), save, public :: active_current_tallies
  type(SetInt), save, public :: active_tallies
!$omp threadprivate(active_analog_tallies, active_tracklength_tallies, &
!$omp&              active_collision_tallies, active_current_tallies, &
!$omp&              active_tallies)

  integer, save, public :: n_realizations = 0 ! # of independent realizations

  ! Tally map structure
! type(TallyMap), save, public, allocatable :: tally_maps(:)

  ! Tally options
  logical, save, public :: assume_separate = .false. ! tallies overlap?
  logical, save, public :: confidence_intervals = .false. ! calculate 95% conf. intervals?
  logical, save, public :: tallies_on = .false. ! tallies on?
  logical, save, public :: reduce_tallies = .true.

  contains

!*******************************************************************************
!*******************************************************************************
! General abstract tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! ACCUMULATE_RESULTS
!===============================================================================

  subroutine accumulate_results(self, total_weight)

    class(TallyClass), intent(inout) :: self
    real(8), intent(in) :: total_weight

    self % n_realizations = self % n_realizations + 1
    call self % results % accumulate(total_weight)

  end subroutine accumulate_results

!===============================================================================
! ADD_FILTER adds a filter to tally filter array
!===============================================================================

  subroutine add_filter(self, filter)

    class(TallyClass), intent(inout) :: self
    class(TallyFilterClass), pointer, intent(in) :: filter

    ! Set filter to array
    self % filters(self % i_filter) % p => filter

    ! Set lookup id in find filter (for manual editing of filter indices)
    self % find_filter(filter % get_type()) = self % i_filter

    ! Increment counter for next filter
    self % i_filter = self % i_filter + 1

    ! Check if its a mesh filter and associate shortcut
    select type(filter)
    type is (EnergyOutFilterClass)
      self % has_eout_filter = .true.
    type is (MeshFilterClass)
      self % has_mesh_filter = .true.
      self % mesh_filter => filter
    end select

  end subroutine add_filter

!===============================================================================
! ADD_SCORE adds a score to tally score array
!===============================================================================

  subroutine add_score(self, score)

    class(TallyClass), intent(inout) :: self
    class(TallyScoreClass), pointer, intent(in) :: score

    ! Set filter to array
    self % scores(self % i_score) % p => score
    self % i_score = self % i_score + 1

  end subroutine add_score

!===============================================================================
! ALLOCATE_FILTERS allocates the filters array in TallyClass instance
!===============================================================================

  subroutine allocate_filters(self, n_filters)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: n_filters

    self % n_filters = n_filters
    allocate(self % filters(n_filters))
    
  end subroutine allocate_filters

!===============================================================================
! ALLOCATE_RESULTS allocates the results array
!===============================================================================

  subroutine allocate_results(self)

    class(TallyClass), intent(inout) :: self

    allocate(self % results(self % total_nuclide_bins, &
       self % total_score_bins, self % total_filter_bins))
    
  end subroutine allocate_results

!===============================================================================
! ALLOCATE_SCORES allocates the scores array in TallyClass instance
!===============================================================================

  subroutine allocate_scores(self, n_scores)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: n_scores

    self % n_scores = n_scores
    allocate(self % scores(n_scores))
 
  end subroutine allocate_scores

!===============================================================================
! ALLOCATE_NUCLIDES allocates the nuclides array
!===============================================================================

  subroutine allocate_nuclides(self, n_nuclides)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: n_nuclides

    self % total_nuclide_bins = n_nuclides
    allocate(self % nuclides(n_nuclides))

  end subroutine allocate_nuclides

!===============================================================================
! TALLY_DESTROY frees all memory used by TallyClass
!===============================================================================

  subroutine tally_destroy(self)

    class(TallyClass), intent(inout) :: self

    integer :: i

    ! Loop through all filters and destroy
    do i = 1, self % n_filters
      call self % filters(i) % p % destroy()
      if (associated(self % filters(i) % p)) deallocate(self % filters(i) % p)
    end do

    ! Destroy tally results
    if (allocated(self % results)) deallocate(self % results)

    ! Destroy stride array
    if (allocated(self % stride)) deallocate(self % stride)

    ! Destroy matching bins array
    if (allocated(self % matching_bins)) deallocate(self % matching_bins)

  end subroutine tally_destroy

!===============================================================================
! FINISH_SETUP finishes the rest of initialization for TallyClass
!===============================================================================

  subroutine finish_setup(self)

    class(TallyClass) :: self

    integer :: i

    ! Set total number of filter bins and scores
    do i = 1, self % n_filters 
      self % total_filter_bins = self % total_filter_bins * &
                                 self % filters(i) % p % get_n_bins()
    end do 
    self % total_score_bins = self % n_scores

    ! Allocate results array
    call self % allocate_results()

    ! Set up stride array
    call self % setup_stride()

  end subroutine finish_setup

!===============================================================================
! SETUP_FILTER_INDICES sets up matching_bins array
!===============================================================================

  subroutine setup_filter_indices(self, p)

    class(TallyClass), intent(inout) :: self
    type(Particle), intent(in) :: p

    integer :: i

    ! Loop through tally filters
    do i = 1, self % n_filters
      self % matching_bins(i) = self % filters(i) % p % get_filter_index(p)
    end do

  end subroutine setup_filter_indices

!===============================================================================
! SET_FILTER_INDEX sets an index for one filter manually
!===============================================================================

  subroutine set_filter_index(self, filter_type, filter_index)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: filter_type
    integer, intent(in) :: filter_index

    integer :: idx ! temporary index

    ! Look up array index from find_filter array
    idx = self % find_filter(filter_type)

    ! Set value in matching_bins array
    self % matching_bins(idx) = filter_index

  end subroutine set_filter_index

!===============================================================================
! GET_FILTER_INDEX calculates filter index from matching_bins array
!===============================================================================

  function get_filter_index(self) result(filter_index)

    class(TallyClass) :: self
    integer :: filter_index

    ! Calculate overall filter index
    filter_index = sum((self % matching_bins - 1) * self % stride) + 1

  end function get_filter_index

!===============================================================================
! GET_N_FILTERS
!===============================================================================

  function get_n_filters(self) result(n_filters)

    class(TallyClass) :: self
    integer :: n_filters

    n_filters = self % n_filters

  end function get_n_filters

!===============================================================================
! GET_TOTAL_SCORE_BINS
!===============================================================================

  function get_total_score_bins(self) result(total_score_bins)

    class(TallyClass) :: self
    integer :: total_score_bins

    total_score_bins = self % total_score_bins

  end function get_total_score_bins

!===============================================================================
! GET_TOTAL_FILTER_BINS
!===============================================================================

  function get_total_filter_bins(self) result(total_filter_bins)

    class(TallyClass) :: self
    integer :: total_filter_bins

    total_filter_bins = self % total_filter_bins

  end function get_total_filter_bins

!===============================================================================
! GET_TOTAL_NUCLIDE_BINS
!===============================================================================

  function get_total_nuclide_bins(self) result(total_nuclide_bins)

    class(TallyClass) :: self
    integer :: total_nuclide_bins

    total_nuclide_bins = self % total_nuclide_bins

  end function get_total_nuclide_bins

!===============================================================================
! GET_FILTERS
!===============================================================================

  function get_filter(self, i) result(filter)

    class(TallyClass) :: self
    integer :: i
    class(TallyFilterClass), pointer :: filter

    filter => self % filters(i) % p

  end function get_filter

!===============================================================================
! GET_RESULTS_POINTER
!===============================================================================

  subroutine get_results_pointer(self, results)

    class(TallyClass), intent(inout), target:: self
    class(TallyResultClass), pointer, intent(out) :: results(:,:,:)

    results => self % results

  end subroutine get_results_pointer

!===============================================================================
! GET_N_REALIZATIONS
!===============================================================================

  function get_n_realizations(self) result(n_realizations)

    class(TallyClass) :: self
    integer :: n_realizations

    n_realizations = self % n_realizations

  end function get_n_realizations

!===============================================================================
! GET_NUCLIDE_BIN
!===============================================================================

  function get_nuclide_bin(self, i) result(bin)

    class(TallyClass) :: self
    integer :: i
    integer :: bin

    bin = self % nuclides(i)

  end function get_nuclide_bin

!===============================================================================
! SET_NUCLIDE_BIN
!===============================================================================

  subroutine set_nuclide_bin(self, i, bin)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: i
    integer, intent(in) :: bin

    self % nuclides(i) = bin

  end subroutine set_nuclide_bin

!===============================================================================
! SETUP_STRIDE
!===============================================================================

  subroutine setup_stride(self)

    class(TallyClass) :: self

    integer :: j
    integer :: n

    ! Allocate stride array and matching bins
    allocate(self % stride(self % n_filters))
    allocate(self % matching_bins(self % n_filters))

    ! The filters are traversed in opposite order so that the last filter has
    ! the shortest stride in memory and the first filter has the largest
    ! stride

    n = 1
    do j = self % n_filters, 1, -1
      self % stride(j) = n
      n = n * self % filters(j) % p % get_n_bins()
    end do
    
  end subroutine setup_stride

!===============================================================================
! SET_ID sets the member id in TallyClass instance
!===============================================================================

  subroutine set_id(self, id)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: id

    self % id = id

  end subroutine set_id

!===============================================================================
! GET_ID
!===============================================================================

  function get_id(self) result(id)

    class(TallyClass) :: self
    integer :: id

    id = self % id

  end function get_id

!===============================================================================
! SET_SCORE_ALL_NUCLIDES
!===============================================================================

  subroutine set_score_all_nuclides(self, score_all)

    class(TallyClass), intent(inout) :: self
    logical :: score_all

    self % score_all_nuclides = score_all

  end subroutine set_score_all_nuclides

!===============================================================================
! TALLY_RESET
!===============================================================================

  subroutine tally_reset(self)

    class(TallyClass), intent(inout) :: self

    self % n_realizations = 0
    call self % results % reset()

  end subroutine tally_reset

!===============================================================================
! TALLY_STATISTICS
!===============================================================================

  subroutine tally_statistics(self)

    class(TallyClass), intent(inout) :: self

    call self % results % statistics(self % n_realizations) 

  end subroutine tally_statistics

!===============================================================================
! TALLY_CONFIDENCE
!===============================================================================

  subroutine tally_confidence(self, t_value)

    class(TallyClass), intent(inout) :: self
    real(8), intent(in) :: t_value

    call self % results % confidence(t_value) 

  end subroutine tally_confidence

!===============================================================================
! WRITE_TALLY
!===============================================================================

  subroutine write_tally(self, unit)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: unit

    integer :: i
    integer :: j
    integer :: k

    ! Write tally information
    write(unit, *) "Tally ID:", self % id
    write(unit, *) "  Estimator:", self % estimator
    write(unit, *) "  Number of realizations:", self % n_realizations
    write(unit, *) "  Filter Information:"

    ! Write filter information
    do i = 1, self % n_filters
      call self % filters(i) % p % write(unit)
    end do

    ! Write score information
    do i = 1, self % n_scores
      call self % scores(i) % p % write(unit)
    end do

    ! Write out results
    do i = 1, self % total_filter_bins
      do j = 1, self % total_score_bins
        do k = 1, self % total_filter_bins
          call self % results(k,j,i) % write(unit)
        end do
      end do
    end do

  end subroutine write_tally

!===============================================================================
! WRITE_TALLY_OUTPUT
!===============================================================================

  subroutine write_tally_output(self)

    class(TallyClass), intent(inout) :: self

    integer :: j            ! level in tally hierarchy
    integer :: k            ! loop index for scoring bins
    integer :: n            ! loop index for nuclides
    integer :: l            ! loop index for user scores
    integer :: type         ! type of tally filter
    integer :: indent       ! number of spaces to preceed output
    integer :: filter_index ! index in results array for filters
    integer :: score_index  ! scoring bin index
    integer :: i_nuclide    ! index in nuclides array
    integer :: i_listing    ! index in xs_listings array
    integer :: n_order      ! loop index for moment orders
    integer :: nm_order     ! loop index for Ynm moment orders
    character(15)           :: filter_name(N_FILTER_TYPES) ! names of tally filters
    character(36)           :: score_names(N_SCORE_TYPES)  ! names of scoring function
    character(36)           :: score_name                  ! names of scoring function

    integer :: nuclide_index

    nuclide_index = 1

    ! Initialize names for tally filter types
    filter_name(FILTER_UNIVERSE)  = "Universe"
    filter_name(FILTER_MATERIAL)  = "Material"
    filter_name(FILTER_CELL)      = "Cell"
    filter_name(FILTER_CELLBORN)  = "Birth Cell"
    filter_name(FILTER_SURFACE)   = "Surface"
    filter_name(FILTER_MESH)      = "Mesh"
    filter_name(FILTER_ENERGYIN)  = "Incoming Energy"
    filter_name(FILTER_ENERGYOUT) = "Outgoing Energy"

    ! Initialize names for scores
    score_names(abs(SCORE_FLUX))          = "Flux"
    score_names(abs(SCORE_TOTAL))         = "Total Reaction Rate"
    score_names(abs(SCORE_SCATTER))       = "Scattering Rate"
    score_names(abs(SCORE_NU_SCATTER))    = "Scattering Production Rate"
    score_names(abs(SCORE_TRANSPORT))     = "Transport Rate"
    score_names(abs(SCORE_N_1N))          = "(n,1n) Rate"
    score_names(abs(SCORE_ABSORPTION))    = "Absorption Rate"
    score_names(abs(SCORE_FISSION))       = "Fission Rate"
    score_names(abs(SCORE_NU_FISSION))    = "Nu-Fission Rate"
    score_names(abs(SCORE_KAPPA_FISSION)) = "Kappa-Fission Rate"
    score_names(abs(SCORE_EVENTS))        = "Events"
    score_names(abs(SCORE_FLUX_YN))       = "Flux Moment"
    score_names(abs(SCORE_TOTAL_YN))      = "Total Reaction Rate Moment"
    score_names(abs(SCORE_SCATTER_N))     = "Scattering Rate Moment"
    score_names(abs(SCORE_SCATTER_PN))    = "Scattering Rate Moment"
    score_names(abs(SCORE_SCATTER_YN))    = "Scattering Rate Moment"
    score_names(abs(SCORE_NU_SCATTER_N))  = "Scattering Prod. Rate Moment"
    score_names(abs(SCORE_NU_SCATTER_PN)) = "Scattering Prod. Rate Moment"
    score_names(abs(SCORE_NU_SCATTER_YN)) = "Scattering Prod. Rate Moment"

    ! Write header block (put labels back in later)
    if (self % label == "") then
      call header("TALLY " // trim(to_str(self % id)), unit=UNIT_TALLY, &
           level=3)
    else
      call header("TALLY " // trim(to_str(self % id)) // ": " &
           // trim(self % label), unit=UNIT_TALLY, level=3)
    endif

    ! WARNING: Admittedly, the logic for moving for printing results is
    ! extremely confusing and took quite a bit of time to get correct. The
    ! logic is structured this way since it is not practical to have a do
    ! loop for each filter variable (given that only a few filters are likely
    ! to be used for a given tally.

    ! Initialize bins, filter level, and indentation
    self % matching_bins(1:self % n_filters) = 0
    j = 1
    indent = 0

    print_bin: do
      find_bin: do
        ! Check for no filters
        if (self % n_filters == 0) exit find_bin

        ! Increment bin combination
        self % matching_bins(j) = self % matching_bins(j) + 1

        ! =================================================================
        ! REACHED END OF BINS FOR THIS FILTER, MOVE TO NEXT FILTER

        if (self % matching_bins(j) > self % filters(j) % p % get_n_bins()) then
          ! If this is the first filter, then exit
          if (j == 1) exit print_bin

          self % matching_bins(j) = 0
          j = j - 1
          indent = indent - 2

          ! =================================================================
          ! VALID BIN -- WRITE FILTER INFORMATION OR EXIT TO WRITE RESULTS

        else
          ! Check if this is last filter
          if (j == self % n_filters) exit find_bin

          ! Print current filter information
          type = self % filters(j) % p % get_type()
          write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
               trim(filter_name(type)), trim(self % get_label(j))
          indent = indent + 2
          j = j + 1
        end if

      end do find_bin

      ! Print filter information
      if (self % n_filters > 0) then
        type = self % filters(j) % p % get_type()
        write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
             trim(filter_name(type)), trim(self % get_label(j))
      end if

      ! Determine scoring index for this bin combination -- note that unlike
      ! in the score_tally subroutine, we have to use max(bins,1) since all
      ! bins below the lowest filter level will be zeros

      if (self % n_filters > 0) then
        filter_index = sum((max(self % matching_bins(1:self%n_filters),1) - 1) * self % stride) + 1
      else
        filter_index = 1
      end if

      ! Write results for this filter bin combination
      if (self % n_filters > 0) indent = indent + 2
      do n = 1, self % total_nuclide_bins
        ! Write label for nuclide
        i_nuclide = self % nuclides(n)
        if (i_nuclide == -1) then
          write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
               "Total Material"
        else
          i_listing = nuclides(i_nuclide) % listing
          write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
               trim(xs_listings(i_listing) % alias)
        end if

        indent = indent + 2
        k = 0
        do l = 1, self % total_score_bins
          k = k + 1
          select case(self % scores(k) % p % get_type())
          case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
            score_name = 'P' // trim(to_str(self % scores(k) % p % &
                get_moment_order())) // " " // &
              score_names(abs(self % scores(k) % p % get_type()))
            write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
              repeat(" ", indent), score_name, &
              to_str(self % results(n, k, filter_index) % get_sum()), &
              trim(to_str(self % results(n, k, filter_index) % get_sum_sq()))
          case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
            score_index = score_index - 1
            do n_order = 0, self % scores(k) % p % get_moment_order()
              score_index = score_index + 1
              score_name = 'P' // trim(to_str(n_order)) //  " " //&
                score_names(abs(self % scores(k) % p % get_type()))
              write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
                repeat(" ", indent), score_name, &
                to_str(self % results(n, k, filter_index) % get_sum()), &
                trim(to_str(self % results(n, k, filter_index) % get_sum_sq()))
            end do
            k = k + self % scores(k) % p % get_moment_order()
          case (SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN, SCORE_FLUX_YN, &
                SCORE_TOTAL_YN)
            score_index = score_index - 1
            do n_order = 0, self % scores(k) % p % get_moment_order()
              do nm_order = -n_order, n_order
                score_index = score_index + 1
                score_name = 'Y' // trim(to_str(n_order)) // ',' // &
                  trim(to_str(nm_order)) // " " // score_names(abs(self % scores(k) % p % get_type()))
                write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
                  repeat(" ", indent), score_name, &
                  to_str(self % results(n, k, filter_index) % get_sum()), &
                  trim(to_str(self % results(n, k, filter_index) % get_sum_sq()))
              end do
            end do
            k = k + (self % scores(k) % p % get_moment_order() + 1)**2 - 1
          case default
            if (self % scores(k) % p % get_type() > 0) then
              score_name = reaction_name(self % scores(k) % p % get_type())
            else
              score_name = score_names(abs(self % scores(k) % p % get_type()))
            end if
            write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
              repeat(" ", indent), score_name, &
              to_str(self % results(n, k, filter_index) % get_sum()), &
              trim(to_str(self % results(n, k, filter_index) % get_sum_sq()))
          end select
        end do
        indent = indent - 2

      end do
      indent = indent - 2

      if (self % n_filters == 0) exit print_bin

    end do print_bin

  end subroutine write_tally_output

!===============================================================================
! GET_LABEL returns a label for a cell/surface/etc given a tally, filter type,
! and corresponding bin
!===============================================================================

  function get_label(self, i_filter) result(label)

    class(TallyClass) :: self ! tally object
    integer :: i_filter ! index in filters array
    character(30) :: label ! user-specified identifier

    integer :: i      ! index in cells/surfaces/etc array
    integer :: bin
    integer, allocatable :: ijk(:) ! indices in mesh
    real(8)              :: E0     ! lower bound for energy bin
    real(8)              :: E1     ! upper bound for energy bin
    type(StructuredMesh), pointer :: m => null()

    bin = self % matching_bins(i_filter)

    select case(self % filters(i_filter) % p % get_type())
    case (FILTER_UNIVERSE)
      i = self % filters(i_filter) % p % get_int_bin(bin)
      label = "" 
    case (FILTER_MATERIAL)
      i = self % filters(i_filter) % p % get_int_bin(bin)
      label = "" 
    case (FILTER_CELL, FILTER_CELLBORN)
      i = self % filters(i_filter) % p % get_int_bin(bin)
      label = "" 
    case (FILTER_SURFACE)
      i = self % filters(i_filter) % p % get_int_bin(bin)
      label = "" 
    case (FILTER_MESH)
      m => self % mesh_filter % get_mesh_pointer()
      allocate(ijk(m % n_dimension))
      call bin_to_mesh_indices(m, bin, ijk)
      if (m % n_dimension == 2) then
        label = "Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ")"
      elseif (m % n_dimension == 3) then
        label = "Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
      end if
    case (FILTER_ENERGYIN, FILTER_ENERGYOUT)
      E0 = self % filters(i_filter) % p % get_real_bin(bin)
      E1 = self % filters(i_filter) % p % get_real_bin(bin + 1)
      label = "[" // trim(to_str(E0)) // ", " // trim(to_str(E1)) // ")"
    end select

  end function get_label

!*******************************************************************************
!*******************************************************************************
! Analog tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! ANALOG_TALLY_INIT initializes an AnalogTallyClass
!===============================================================================

  function analog_tally_init() result(self)

    class(AnalogTallyClass), pointer :: self

    ! Allocate
    allocate(self)

    ! Set the tally type
    self % estimator = ESTIMATOR_ANALOG

  end function analog_tally_init

!===============================================================================
! ANALOG_NUFISSION_EOUT scores nu-fission energy-out filter
!===============================================================================

  subroutine analog_nufission_eout(self, p, p_fiss, score_index, nuclide_index)

    class(AnalogTallyClass), intent(inout) :: self
    integer :: score_index
    type(Particle), intent(in) :: p
    type(Particle), pointer, intent(inout) :: p_fiss

    integer(8) :: k
    integer :: filter_index
    integer :: nuclide_index
    real(8) :: score

    ! Loop around particles in bank
    do k = n_bank - int(p % nu, 8) + 1, n_bank

      ! Check to create particle
      if (.not. associated(p_fiss)) allocate(p_fiss)
        call p_fiss % initialize()

      ! Copy bank attributes to fake fission particle
      p_fiss % last_E = p % last_E
      p_fiss % E = fission_bank(k) % E
      p_fiss % wgt = fission_bank(k) % wgt
      p_fiss % coord0 % xyz = fission_bank(k) % xyz

      ! Get filter index
      call self % setup_filter_indices(p_fiss)
      filter_index = self % get_filter_index()

      ! Get standard analog score
      score = keff * self % scores(score_index) % p % get_weight(p_fiss)

      ! Add score to results array
      call self % results(nuclide_index, score_index, filter_index) % add(score)

    end do

  end subroutine analog_nufission_eout

!===============================================================================
! ANALOG_TALLY_SCORE performs a score to an analog tally
!===============================================================================

  subroutine analog_tally_score(self, p)

    class(AnalogTallyClass), intent(inout) :: self
    type(Particle), target, intent(in) :: p

    integer :: filter_index
    integer :: i_nuclide
    integer :: nuclide_index
    integer :: j
    integer :: k
    real(8) :: score ! the score to record
    type(Particle), pointer :: p_fiss => null()

    ! Get filter index - this will be the filter index for every score
    ! but where there is nu-fission with an energyout filter. This
    ! filter index will be recalculated for each procedure neutron.
    call self % setup_filter_indices(p)
    filter_index = self % get_filter_index()

    ! Need to get nuclide index
    k = 0
    NUCLIDE_LOOP: do while (k < n_nuclides_total)

      ! Increment the index in the list of nuclide bins
      k = k + 1

      if (self % score_all_nuclides) then
        ! In the case that the user has requested to tally all nuclides, we
        ! can take advantage of the fact that we know exactly how nuclide
        ! bins correspond to nuclide indices.
        if (k == 1) then
          ! If we just entered, set the nuclide bin index to the index in
          ! the nuclides array because this will match the index in the
          ! nuclide bin array
          k = p % event_nuclide
          nuclide_index = k
        else if (k == p % event_nuclide + 1) then
          ! After we've tallied the individual nuclide bin, we also need
          ! to contribute to the total material bin which is the last bin
          k = n_nuclides_total + 1
          nuclide_index = k
        else
          ! After we've tallied in both the individual nuclide bin and the
          ! total material bin, we are done
          exit
        end if
      else
        ! If the user has explicitly specified nuclides (or specified
        ! none), we need to search through the nulide bin list one by
        ! one. First we need to get the value of the nuclide bin
        i_nuclide = self % nuclides(k)
        nuclide_index = k

        ! Now compare the value against that of the colliding nuclide
        if (i_nuclide /= p % event_nuclide .and. &
            i_nuclide /= MATERIAL_TOTAL) cycle
      end if

      ! Loop around score bins now that we now the nuclide or material
      SCORE_LOOP: do j = 1, self % n_scores

        ! Check if event matches score type
        if (.not. self % scores(j) % p % score_match(p)) cycle

        ! Get standard analog score
        score = self % scores(j) % p % get_weight(p)

        ! Special cases
        select type(s => self % scores(j) % p)

        ! Nu-fission score needs how many neutrons were produced
        type is (NuFissionScoreClass)

          ! For energy-out need to perform for each  
          if (self % has_eout_filter) then
            call self % analog_nufission_eout(p, p_fiss, j, nuclide_index)
            cycle
          else
            score = keff * p % wgt_bank
          end if

        end select

        ! Add score to results array
        call self % results(nuclide_index, j, filter_index) % add(score)

      end do SCORE_LOOP

    end do NUCLIDE_LOOP

    ! Deallocate particle pointers if associated
    if (associated(p_fiss)) then
      call p_fiss % clear()
      deallocate(p_fiss)
    end if

  end subroutine analog_tally_score

!*******************************************************************************
!*******************************************************************************
! Tracklength tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! TRACKLENGTH_TALLY_INIT initializes an TracklengthTallyClass
!===============================================================================

  function tracklength_tally_init() result(self)

    class(TracklengthTallyClass), pointer :: self

    ! Allocate
    allocate(self)

    ! Set the tally type
    self % estimator = ESTIMATOR_TRACKLENGTH

  end function tracklength_tally_init

!===============================================================================
! TRACKLENGTH_GET_FLUX gets the flux estimator for a TracklengthTallyClass
!===============================================================================

  function tracklength_get_flux(p) result(flux)

    type(Particle) :: p
    real(8) :: flux

    flux = p % dist

  end function tracklength_get_flux

!===============================================================================
! TRACKLENGTH_TALLY_SCORE performs a score to a tracklength tally
!===============================================================================

  subroutine tracklength_tally_score(self, p)

    class(TracklengthTallyClass), intent(inout) :: self
    type(Particle), target, intent(in) :: p

    integer :: bin ! mesh bin
    integer :: filter_index
    integer :: i
    integer :: ii
    integer :: j
    integer :: k
    integer :: i_nuclide
    integer :: n_cross ! number of surface crossings
    integer :: nuclide_index
    integer :: bins_left
    logical :: bins_done
    logical :: found_bin
    real(8) :: score ! the score to record
    real(8) :: flux ! estimate of the flux
    real(8) :: response ! tally response
    real(8) :: weight ! some form of a neutron statistical weight
    real(8) :: number_density
    type(Particle), pointer :: p_fiss => null()
    type(Particle), pointer :: p_score => null()
    type(Material), pointer :: mat => null()

    ! Point to material
    mat => materials(p % material)

    ! Set up bins left to score to
    bins_left = self % total_nuclide_bins
    bins_done = .false.

    ! Loop around nuclides in material
    NUCLIDE_MAT_LOOP: do i = 1, mat % n_nuclides + 1

      ! Just reset i_nuclide to 1 so that we know for a fact that it gets
      ! set to MATERIAL_TOTAL
      i_nuclide = 1

      ! Check for any bins left
      if (bins_left == 0) exit NUCLIDE_MAT_LOOP

      ! Check for material total score
      if (i == mat % n_nuclides + 1) then
        if (self % nuclides(self % total_nuclide_bins) == MATERIAL_TOTAL) then
          i_nuclide = MATERIAL_TOTAL
          nuclide_index = self % total_nuclide_bins
          number_density = ONE
        else
          exit NUCLIDE_MAT_LOOP
        end if
      end if

      ! Check if scoring to all nuclide bins (i_nuclide matches nuclide_index)
      if (i_nuclide /= MATERIAL_TOTAL) then
        if (self % score_all_nuclides) then
          i_nuclide = mat % nuclide(i)
          nuclide_index = mat % nuclide(i)
          number_density = mat % atom_density(i)
        else
          NUCLIDE_BIN_LOOP: do ii = 1, self % total_nuclide_bins

            ! Check to see if there is a match
            if (mat % nuclide(i) == self % nuclides(ii)) then
              i_nuclide = self % nuclides(ii)
              nuclide_index = ii
              number_density = mat % atom_density(i)
              exit NUCLIDE_BIN_LOOP
            end if

            ! Went through all nuclide bins
            if (ii == self % total_nuclide_bins) then
              cycle NUCLIDE_MAT_LOOP
           end if

          end do NUCLIDE_BIN_LOOP
        end if
      end if

      ! Going to score to a nuclide bin
      bins_left = bins_left - 1

      ! Loop around score bins
      SCORE_LOOP: do j = 1, self % n_scores

        ! Set scoring particle to be the actual particle
        p_score => p

        ! Special cases for enery out filter
        if (self % has_eout_filter) then
          select type (s => self % scores(j) % p)

          type is (NuFissionScoreClass)

            ! If no fission macro, don't score
            if (material_xs % fission == ZERO) cycle

            ! Check for nuclide fission
            if (i_nuclide /= MATERIAL_TOTAL) then
              if (.not. nuclides(i_nuclide) % fissionable) cycle
            end if

            ! Check to see if we need to sample a fission reaction
            if (.not. associated(p_fiss)) p_fiss => sample_fake_fission(p, i_nuclide)
            p_score => p_fiss

          end select 

        end if

        ! Check if there is a mesh filter
        if (self % has_mesh_filter) then

          ! Get filter index
          call self % setup_filter_indices(p_score)

          ! Get tally score response and weight
          response = self % scores(j) % p % get_response(i_nuclide)
          weight = self % scores(j) % p % get_weight(p_score)

          ! Get number of surface crossings
          call self % mesh_filter % get_crossings(p_score, n_cross)

          ! Loop around number of crossings and record
          ! Maybe this shouldn't be on the inner most loop because
          ! different score types will have the same flux contributions
          do k = 1, n_cross

            ! Get the distance traveled through a bin
            call self % mesh_filter % get_next_distance(p_score, flux, bin, found_bin)

            ! Don't continue if no bin was found
            if (.not. found_bin) cycle

            ! Alter filter indices
            call self % set_filter_index(FILTER_MESH, bin)

            ! Get overall filter index
            filter_index = self % get_filter_index()

            ! Calculate score
            score = weight * response * flux

            ! Add score to results array
            call self % results(nuclide_index, j, filter_index) % add(score)

          end do 

        else

          ! Get filter index
          call self % setup_filter_indices(p_score)
          filter_index = self % get_filter_index()

          ! Calculate standard tracklength score
          flux = self % get_flux(p_score)
          response = self % scores(j) % p % get_response(i_nuclide)
          weight = self % scores(j) % p % get_weight(p_score)
          score = weight * number_density * response * flux

          ! Add score to results array
          call self % results(nuclide_index, j, filter_index) % add(score)

        end if

      end do SCORE_LOOP

    end do NUCLIDE_MAT_LOOP

    ! Deallocate particle pointers if associated
    if (associated(p_fiss)) then
      deallocate(p_fiss)
    end if

  end subroutine tracklength_tally_score

!*******************************************************************************
!*******************************************************************************
! Collision tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! COLLISION_TALLY_INIT initializes an CollisionTallyClass
!===============================================================================

  function collision_tally_init() result(self)

    class(CollisionTallyClass), pointer :: self

    ! Allocate
    allocate(self)

    ! Set the tally type
    self % estimator = ESTIMATOR_COLLISION

  end function collision_tally_init

!===============================================================================
! COLLISION_GET_FLUX gets the flux estimator for a CollisionTallyClass
!===============================================================================

  function collision_get_flux() result(flux)

    real(8) :: flux

    flux = ONE/material_xs % total

  end function collision_get_flux

!===============================================================================
! COLLISION_TALLY_SCORE performs a score to a collision tally
!===============================================================================

  subroutine collision_tally_score(self, p)

    class(CollisionTallyClass), intent(inout) :: self
    type(Particle), target, intent(in) :: p

    integer :: filter_index
    integer :: i_nuclide
    integer :: nuclide_index
    integer :: j
    real(8) :: score ! the score to record
    real(8) :: flux ! estimate of the flux
    real(8) :: response ! tally response
    real(8) :: weight ! some form of a neutron statistical weight
    type(Particle), pointer :: p_fiss => null()
    type(Particle), pointer :: p_score => null()

    i_nuclide = MATERIAL_TOTAL
    nuclide_index = 1

    ! Loop around score bins
    SCORE_LOOP: do j = 1, self % n_scores

      ! Set scoring particle to be the actual particle
      p_score => p

      ! Special cases for enery out filter
      if (self % has_eout_filter) then
        select type (s => self % scores(j) % p)

        type is (NuFissionScoreClass)

          ! If no fission macro, don't score
          if (material_xs % fission == ZERO) cycle

          ! Check for nuclide fission
          if (i_nuclide /= MATERIAL_TOTAL) then
            if (.not. nuclides(i_nuclide) % fissionable) cycle
          end if

          ! Check to see if we need to sample a fission reaction
          if (.not. associated(p_fiss)) p_fiss => sample_fake_fission(p, i_nuclide)
          p_score => p_fiss

        end select 

      end if

      ! Calculate score
      flux = self % get_flux()
      response = self % scores(j) % p % get_response(i_nuclide)
      weight = self % scores(j) % p % get_weight(p_score)
      score = weight * response * flux

      ! Get filter index
      call self % setup_filter_indices(p_score)
      filter_index = self % get_filter_index()

      ! Add score to results array
      call self % results(nuclide_index, j, filter_index) % add(score)

    end do SCORE_LOOP

    ! Deallocate particle pointers if associated
    if (associated(p_fiss)) then
      deallocate(p_fiss)
    end if

  end subroutine collision_tally_score

!*******************************************************************************
!*******************************************************************************
! Private routines
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! SAMPLE_FAKE_FISSION samples a fake fission reaction
!===============================================================================

  function sample_fake_fission(p, i_nuclide) result(p_fiss)

    integer :: i_nuclide
    type(Particle) :: p
    type(Particle), pointer :: p_fiss

    integer :: i_nuclide_rxn
    integer :: i_reaction
    type(Nuclide),  pointer :: nuc
    type(Reaction), pointer :: rxn

    ! Switch over random number stream to tallies
    call prn_set_stream(STREAM_TALLIES)

    ! Allocate fake fission particle and initialize
    allocate(p_fiss)
    call p_fiss % initialize()

    ! Copy particle attributes over
    p_fiss = p

    ! Sample nuclide for fission reaction
    if (i_nuclide == MATERIAL_TOTAL) then
      i_nuclide_rxn = sample_nuclide(p_fiss, 'fission')
    else
      i_nuclide_rxn = i_nuclide
    end if

    ! Sample a fake fission
    call sample_fission(i_nuclide_rxn, i_reaction)

    ! Set up pointers and get fission energy
    nuc => nuclides(i_nuclide_rxn)
    rxn => nuc % reactions(i_reaction)
    p_fiss % E = sample_fission_energy(nuc, rxn, p_fiss % last_E)

    ! Switch over random number stream to tallies
    call prn_set_stream(STREAM_TRACKING)

  end function sample_fake_fission

!===============================================================================
! HEADER displays a header block according to a specified level. If no level is
! specified, it is assumed to be a minor header block (H3).
!===============================================================================

  subroutine header(msg, unit, level)

    use, intrinsic :: ISO_FORTRAN_ENV 

    character(*), intent(in) :: msg ! header message
    integer, optional :: unit       ! unit to write to
    integer, optional :: level      ! specified header level

    integer :: n            ! number of = signs on left
    integer :: m            ! number of = signs on right
    integer :: unit_        ! unit to write to
    integer :: header_level ! actual header level
    character(MAX_LINE_LEN) :: line

    ! set default level
    if (present(level)) then
      header_level = level
    else
      header_level = 3
    end if

    ! set default unit
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! determine how many times to repeat '=' character
    n = (63 - len_trim(msg))/2
    m = n
    if (mod(len_trim(msg),2) == 0) m = m + 1

    ! convert line to upper case
    line = msg
    call upper_case(line)

    ! print header based on level
    select case (header_level)
    case (1)
      write(UNIT=unit_, FMT='(/3(1X,A/))') repeat('=', 75), &
           repeat('=', n) // '>     ' // trim(line) // '     <' // &
           repeat('=', m), repeat('=', 75)
    case (2)
      write(UNIT=unit_, FMT='(/2(1X,A/))') trim(line), repeat('-', 75)
    case (3)
      write(UNIT=unit_, FMT='(/1X,A/)') repeat('=', n) // '>     ' // &
           trim(line) // '     <' // repeat('=', m)
    end select

  end subroutine header

end module tally_class
