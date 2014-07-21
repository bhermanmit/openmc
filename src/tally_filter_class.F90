module tally_filter_class

  use constants

  implicit none
  private
  public :: TallyFilter_p, TallyFilterClass, EnergyFilterClass

  ! General tally filter type
  type, abstract :: TallyFilterClass
    private
    integer :: type    ! Type of filter from constants
    contains
      procedure, public :: set_type
      procedure, public :: get_type
      procedure (filter_destroy), deferred :: destroy 
  end type TallyFilterClass

  abstract interface
    subroutine filter_destroy(self)
      import TallyFilterClass
      class(TallyFilterClass) :: self
    end subroutine filter_destroy
  end interface

  ! Tally filter pointer
  type :: TallyFilter_p
    class(TallyFilterClass), pointer :: p => null()
  end type TallyFilter_p

  ! Energy filter
  type, extends(TallyFilterClass) :: EnergyFilterClass
    private
    integer :: n_bins
    real(8), allocatable :: bins(:)
    contains
      procedure, public :: set_bins => set_energy_bins
      procedure, public :: destroy => energy_filter_destroy
  end type EnergyFilterClass
  interface EnergyFilterClass
    module procedure energy_filter_init
  end interface

  contains

!*******************************************************************************
!*******************************************************************************
! General abstract tally methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! SET_TYPE sets the member type in TallyClass instance
!===============================================================================

  subroutine set_type(self, type)

    class(TallyFilterClass) :: self
    integer :: type

    ! Set the type to instance
    self % type = type

  end subroutine set_type

!===============================================================================
! GET_TYPE returns the member type from TallyClass instance
!===============================================================================

  function get_type(self) result(type)

    class(TallyFilterClass) :: self
    integer :: type

    ! Get the type from instance
    type = self % type

  end function get_type

!*******************************************************************************
!*******************************************************************************
! Energy filter methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! ENERGY_FILTER_INIT
!===============================================================================

  function energy_filter_init() result(self)

    class(EnergyFilterClass), pointer :: self

    ! Create object
    allocate(self)

    ! Set type of filter
    call self % set_type(FILTER_ENERGYIN)

  end function energy_filter_init

!===============================================================================
! SET_ENERGY_BINS
!===============================================================================

  subroutine set_energy_bins(self, n_bins, bins)

    class(EnergyFilterClass) :: self
    integer :: n_bins
    real(8) :: bins(:)

    ! Allocate number of bins and set information to instance
    allocate(self % bins(n_bins))
    self % n_bins = n_bins
    self % bins = bins

  end subroutine set_energy_bins

!===============================================================================
! ENERGY_FILTER_DESTROY
!===============================================================================

  subroutine energy_filter_destroy(self)

    class(EnergyFilterClass) :: self

    ! Free memory associated with energy filter
    deallocate(self % bins)

  end subroutine energy_filter_destroy

end module tally_filter_class
