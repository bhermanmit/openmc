module tally_class

  use constants
  use particle_header
  use tally_filter_class 
  use tally_result_class
  use tally_score_class 

  implicit none
  private

  ! General tally
  type, abstract, public :: TallyClass
    private
    character(len=MAX_LINE_LEN) :: label = ""
    character(len=MAX_WORD_LEN) :: estimator = ""
    integer :: id ! ID of tally
    integer :: i_filter = 1  ! Current filter
    integer :: i_score = 1  ! Current score
    integer :: n_filters = ZERO ! Number of filters
    integer :: n_realizations = ZERO ! Number of tally realizations
    integer :: n_scores = ZERO ! Number of scores
    integer :: total_score_bins = ZERO ! Total number of score bins
    integer :: total_filter_bins = ZERO ! Total number of filter bins
    integer, allocatable :: stride(:)
    integer, allocatable :: matching_bins(:)
    type(TallyFilter_p), allocatable :: filters(:) ! Polymorphic array of filter objects
    type(TallyResultClass), allocatable :: results(:,:) ! Array of result objects
    type(TallyScore_p), allocatable :: scores(:) ! Polymorphic array of filter objects
    contains
      procedure, public :: accumulate => accumulate_results
      procedure, public :: add_filter
      procedure, public :: add_score
      procedure, public :: allocate_filters
      procedure :: allocate_results
      procedure, public :: allocate_scores
      procedure, public :: destroy => tally_destroy
      procedure, public :: finish_setup
      procedure :: get_filter_index
      procedure, public :: reset => tally_reset
      procedure :: setup_stride
      procedure, public :: set_id
      procedure, public :: statistics => tally_statistics
      procedure, public :: write => write_tally
      procedure, public :: score => tally_score
  end type TallyClass

  ! Tally pointer type
  type, public :: Tally_p
    class(TallyClass), pointer :: p => null()
  end type Tally_p

  ! Analog tally 
  type, extends(TallyClass), public :: AnalogTallyClass
    private
  end type AnalogTallyClass

  ! Constructor call for analog tally
  interface AnalogTallyClass
    module procedure analog_tally_init
  end interface

  ! Tracklength tally 
  type, extends(TallyClass), public :: TracklengthTallyClass
    private
    contains
      procedure :: get_flux => tracklength_get_flux
  end type TracklengthTallyClass

  ! Constructor call for tracklength tally
  interface TracklengthTallyClass
    module procedure tracklength_tally_init
  end interface

  ! Collision tally 
  type, extends(TallyClass), public :: CollisionTallyClass
    private
    contains
      procedure :: get_flux => collision_get_flux
  end type CollisionTallyClass

  ! Constructor call for collision tally
  interface CollisionTallyClass
    module procedure collision_tally_init
  end interface

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
    self % i_filter = self % i_filter + 1

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

    allocate(self % results(self % total_score_bins, self % total_filter_bins))
    
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
! TALLY_DESTROY frees all memory used by TallyClass
!===============================================================================

  subroutine tally_destroy(self)

    class(TallyClass), intent(inout) :: self

    integer :: i

    ! Loop through all filters and destroy
    do i = 1, self % n_filters
      call self % filters(i) % p % destroy()
      if (associated(self % filters(i) % p)) nullify(self % filters(i) % p)
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
      self % total_filter_bins = self % total_filter_bins + &
                                 self % filters(i) % p % get_n_bins()
    end do 
    self % total_score_bins = self % n_scores

    ! Allocate results array
    call self % allocate_results()

    ! Set up stride array
    call self % setup_stride()

  end subroutine finish_setup

!===============================================================================
! GET_FILTER_INDEX
!===============================================================================

 function get_filter_index(self, p) result(filter_index)

    class(TallyClass) :: self
    type(Particle) :: p
    integer :: filter_index

    integer :: i

    ! Loop through tally filters
    do i = 1, self % n_filters
      self % matching_bins(i) = self % filters(i) % p % get_filter_index(p)
    end do

    ! Calculate overall filter index
    filter_index = sum((self % matching_bins - 1) * self % stride) + 1

 end function get_filter_index

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
! TALLY_RESET
!===============================================================================

  subroutine tally_reset(self)

    class(TallyClass), intent(inout) :: self

    self % n_realizations = 0
    call self % results % reset()

  end subroutine tally_reset

!===============================================================================
! TALLY_SCORE performs a score to a tally
!===============================================================================

  subroutine tally_score(self, p)

    class(TallyClass) :: self
    type(Particle) :: p

    integer :: filter_index
    integer :: j
    real(8) :: score
    real(8) :: flux
    real(8) :: response
    real(8) :: weight

    ! Get filter index
    filter_index = self % get_filter_index(p)

    ! Loop around score bins
    do j = 1, self % n_scores

      ! Get appropriate particle weight
      weight = self % scores(j) % p % get_weight(p)

      ! Calculate appropriate score 
      select type(self)

      type is (AnalogTallyClass)
        score = weight

      type is (TracklengthTallyClass)
        flux = self % get_flux(p)
        response = self % scores(j) % p % get_response(p)
        score = weight * response * flux

      type is (CollisionTallyClass)
        flux = self % get_flux(p)
        response = self % scores(j) % p % get_response(p)
        score = weight * response * flux

      end select

      ! Add score to results array
      call self % results(j, filter_index) % add(score)

    end do

  end subroutine tally_score

!===============================================================================
! TALLY_STATISTICS
!===============================================================================

  subroutine tally_statistics(self)

    class(TallyClass), intent(inout) :: self

    call self % results % statistics(self % n_realizations) 

  end subroutine tally_statistics

!===============================================================================
! WRITE_TALLY
!===============================================================================

  subroutine write_tally(self, unit)

    class(TallyClass), intent(inout) :: self
    integer, intent(in) :: unit

    integer :: i
    integer :: j

    ! Write tally information
    write(unit, *) "Tally ID:", self % id
    write(unit, *) "  Estimator:", trim(self % estimator)
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
        call self % results(j,i) % write(unit)
      end do
    end do

  end subroutine write_tally

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
    self % estimator = 'analog'

  end function analog_tally_init

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
    self % estimator = 'tracklength'

  end function tracklength_tally_init

!===============================================================================
! TRACKLENGTH_GET_FLUX gets the flux estimator for a TracklengthTallyClass
!===============================================================================

  function tracklength_get_flux(self, p) result(flux)

    class(TracklengthTallyClass) :: self
    type(Particle) :: p
    real(8) :: flux

    flux = p % dist

  end function tracklength_get_flux

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
    self % estimator = 'collision'

  end function collision_tally_init

!===============================================================================
! COLLISION_GET_FLUX gets the flux estimator for a CollisionTallyClass
!===============================================================================

  function collision_get_flux(self, p) result(flux)

    class(CollisionTallyClass) :: self
    type(Particle) :: p
    real(8) :: flux

    flux = ONE/p % material_xs % total

  end function collision_get_flux

end module tally_class
