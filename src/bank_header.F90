module bank_header

  implicit none

!===============================================================================
! BANK is used for storing fission sites in eigenvalue calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type Bank
    ! The 'sequence' attribute is used here to ensure that the data listed
    ! appears in the given order. This is important for MPI purposes when bank
    ! sites are sent from one processor to another.
    sequence

    real(8)    :: wgt    ! weight of bank site
    real(8)    :: xyz(3) ! location of bank particle
    real(8)    :: uvw(3) ! diretional cosines
    real(8)    :: E      ! energy
  end type Bank

  ! Source and fission bank
  type(Bank), save, allocatable, target :: source_bank(:)
  type(Bank), save, allocatable, target :: fission_bank(:)
#ifdef _OPENMP
  type(Bank), save, allocatable, target :: master_fission_bank(:)
#endif
  integer(8), save, target :: n_bank       ! # of sites in fission bank
  integer(8), save :: work         ! number of particles per processor
  integer(8), save, allocatable :: work_index(:) ! starting index in source bank for each process

!$omp threadprivate(fission_bank, n_bank)

end module bank_header
