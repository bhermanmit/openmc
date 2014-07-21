module tally_filter_class

  use constants

  implicit none

  ! General tally filter type
  type, abstract :: TallyFilter
    private
    integer :: type    ! Type of filter from constants
    integer :: stride  ! Stride in tally results
    contains
      procedure, public :: set_type
      procedure, public :: get_type
      procedure, public :: set_stride
      procedure, public :: get_stride
      procedure (tally_filter_destroy), deferred :: destroy 
  end type TallyFilter

  abstract interface
    subroutine tally_filter_destroy(self)
      import TallyFilter
      class(TallyFilter) :: self
    end subroutine tally_filter_destroy
  end interface

  ! Energy filter
  type, extends(TallyFilter) :: EnergyFilter
    private
    integer :: n_bins
    real(8), allocatable :: bins(:)
    contains
      procedure, public :: set_bins => set_energy_bins
      procedure, public :: destroy => energy_filter_destroy
  end type EnergyFilter

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
    class(TallyFilter) :: self
    integer :: type
    self % type = type
  end subroutine set_type

!===============================================================================
! GET_TYPE returns the member type from TallyClass instance
!===============================================================================

  integer function get_type(self) result(type)
    class(TallyFilter) :: self
    type = self % type
  end function get_type

!===============================================================================
! SET_STRIDE sets the member stride in TallyClass instance
!===============================================================================

  subroutine set_stride(self, stride)
    class(TallyFilter) :: self
    integer :: stride
    self % stride = stride
  end subroutine set_stride

!===============================================================================
! GET_STRIDE returns the member stride from TallyClass instance
!===============================================================================

  integer function get_stride(self) result(stride)
    class(TallyFilter) :: self
    stride = self % stride
  end function get_stride

!*******************************************************************************
!*******************************************************************************
! Energy filter methods
!*******************************************************************************
!*******************************************************************************

!===============================================================================
! INIT_ENERGY_FILTER
!===============================================================================

  type(EnergyFilter) function init_energy_filter() result(self)

    ! Set up the type
    call self % set_type(FILTER_ENERGYIN)
    
  end function init_energy_filter

!===============================================================================
! SET_ENERGY_BINS
!===============================================================================

  subroutine set_energy_bins(self, n_bins, bins)

    class(EnergyFilter) :: self
    integer :: n_bins
    real(8) :: bins(:)

    ! Set up energy bins
    allocate(self % bins(n_bins))
    self % bins = bins

  end subroutine set_energy_bins

!===============================================================================
! ENERGY_FILTER_DESTROY
!===============================================================================

  subroutine energy_filter_destroy(self)

    class(EnergyFilter) :: self

    ! Deallocate bin information
    deallocate(self % bins)

  end subroutine energy_filter_destroy

end module tally_filter_class
