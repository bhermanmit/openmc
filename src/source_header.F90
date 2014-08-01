module source_header

  implicit none

!===============================================================================
! EXTSOURCE describes an external source of neutrons for a fixed-source problem
! or for the starting source in a k eigenvalue problem
!===============================================================================

  type ExtSource
    integer :: type_space              ! spacial distribution, e.g. 'box' or 'point'
    integer :: type_angle              ! angle distribution, e.g. 'isotropic'
    integer :: type_energy             ! energy distribution, e.g. 'Watt'
    real(8), allocatable :: params_space(:)  ! parameters for spatial distribution
    real(8), allocatable :: params_angle(:)  ! parameters for angle distribution
    real(8), allocatable :: params_energy(:) ! parameters for energy distribution
  end type ExtSource

  ! External source
  type(ExtSource), save, target :: external_source

  ! Trace
  logical    :: trace
  integer    :: trace_batch
  integer    :: trace_gen
  integer(8) :: trace_particle
!$omp threadprivate(trace)


end module source_header
