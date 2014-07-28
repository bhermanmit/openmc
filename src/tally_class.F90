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
    integer :: find_filter(N_FILTER_TYPES)
    integer, allocatable :: stride(:)
    integer, allocatable :: matching_bins(:)
    logical :: has_eout_filter = .false.
    logical :: has_mesh_filter = .false.
    type(MeshFilterClass), pointer :: mesh_filter => null()
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
      procedure :: set_filter_index
      procedure :: setup_filter_indices
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
    contains
      procedure :: analog_nufission_eout
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
    class(TallyFilterClass), pointer :: f => null()

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

    class(TallyScoreClass), pointer :: s => null()
    integer :: bin ! mesh bin
    integer :: filter_index
    integer :: j
    integer :: k
    integer :: n_cross ! number of surface crossings
    logical :: found_bin
    real(8) :: score ! the score to record
    real(8) :: flux ! estimate of the flux
    real(8) :: response ! tally response
    real(8) :: weight ! some form of a neutron statistical weight
    type(Particle), pointer :: p_fiss => null()

    ! Loop around score bins
    SCORE_LOOP: do j = 1, self % n_scores

      ! Calculate appropriate score depending on TallyClass type
      select type(self)

      type is (AnalogTallyClass)

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

            call self % analog_nufission_eout(p, p_fiss, j)

          else

            ! Get filter index
            call self % setup_filter_indices(p)
            filter_index = self % get_filter_index()

            ! Get standard nu-fission score
            score = p % keff * p % wgt_bank

          end if

        end select

        ! Add score to results array
        call self % results(j, filter_index) % add(score)

      type is (TracklengthTallyClass)
        ! Special cases
        select type (s => self % scores(j) % p)

        type is (NuFissionScoreClass)

          ! Check to see if we need to sample a fission reaction
          if (.not. associated(p_fiss)) call sample_fake_fission(p_fiss)

        end select 

        ! Check if there is a mesh filter
        if (self % has_mesh_filter) then

          ! Get filter index
          call self % setup_filter_indices(p)

          ! Get tally score response and weight
          response = self % scores(j) % p % get_response(p)
          weight = self % scores(j) % p % get_weight(p)

          ! Get number of surface crossings
          call self % mesh_filter % get_crossings(p, n_cross)

          ! Loop around number of crossings and record
          ! Maybe this shouldn't be on the inner most loop because
          ! different score types will have the same flux contributions
          do k = 1, n_cross

            ! Get the distance traveled through a bin
            call self % mesh_filter % get_next_distance(p, flux, bin, found_bin)

            ! Don't continue if no bin was found
            if (.not. found_bin) cycle

            ! Alter filter indices
            call self % set_filter_index(FILTER_MESH, bin)

            ! Get overall filter index
            filter_index = self % get_filter_index()

            ! Calculate score
            score = weight * response * flux

            ! Add score to results array
            call self % results(j, filter_index) % add(score)

          end do 

        else

          ! Get filter index
          call self % setup_filter_indices(p)
          filter_index = self % get_filter_index()

          ! Calculate standard tracklength score
          flux = self % get_flux(p)
          response = self % scores(j) % p % get_response(p)
          weight = self % scores(j) % p % get_weight(p)
          score = weight * response * flux

          ! Add score to results array
          call self % results(j, filter_index) % add(score)

        end if

      type is (CollisionTallyClass)

        ! Calculate score
        flux = self % get_flux(p)
        response = self % scores(j) % p % get_response(p)
        weight = self % scores(j) % p % get_weight(p)
        score = weight * response * flux

        ! Get filter index
        call self % setup_filter_indices(p)
        filter_index = self % get_filter_index()

        ! Add score to results array
        call self % results(j, filter_index) % add(score)

      end select

    end do SCORE_LOOP

    ! Deallocate particle pointers if associated
    if (associated(p_fiss)) then
      call p_fiss % clear()
      deallocate(p_fiss)
    end if

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

!===============================================================================
! ANALOG_NUFISSION_EOUT scores nu-fission energy-out filter
!===============================================================================

  subroutine analog_nufission_eout(self, p, p_fiss, score_index)

    class(AnalogTallyClass), intent(inout) :: self
    integer :: score_index
    type(Particle), intent(in) :: p
    type(Particle), pointer, intent(inout) :: p_fiss

    integer :: k
    integer :: filter_index
    real(8) :: score

    ! Loop around particles in bank
    do k = 1, p % n_bank

      ! Check to create particle
      if (.not. associated(p_fiss)) allocate(p_fiss)
        call p_fiss % initialize()

      ! Copy bank attributes to fake fission particle
      p_fiss % last_E = p % last_E
      p_fiss % E = p % fission_bank(k) % E
      p_fiss % wgt = p % fission_bank(k) % wgt
      p_fiss % coord0 % xyz = p % fission_bank(k) % xyz

      ! Get filter index
      call self % setup_filter_indices(p_fiss)
      filter_index = self % get_filter_index()

      ! Get standard analog score
      score = p % keff * self % scores(score_index) % p % get_weight(p_fiss)

      ! Add score to results array
      call self % results(score_index, filter_index) % add(score)

    end do

  end subroutine analog_nufission_eout

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

!===============================================================================
! SAMPLE_FAKE_FISSION samples a fake fission reaction
!===============================================================================

  subroutine sample_fake_fission(p)

    type(Particle) :: p

    ! Sample nuclide for fission reaction
!   i_nuclide_rxn = sample_nuclide(p_fake, 'fission')

    ! Sample a fake fission
!   call sample_fission(i_nuclide_rxn, i_reaction)

    ! Set up pointers and get fission energy
!   nuc => nuclides(i_nuclide_rxn)
!   rxn => nuc % reactions(i_reaction)
!   p_fake % E = sample_fission_energy(nuc, rxn, p_fake % last_E)

  end subroutine sample_fake_fission

end module tally_class
