module mpi_interface

#ifdef MPI
  use mpi
#endif

  implicit none

  ! The defaults set here for the number of processors, rank, and master and
  ! mpi_enabled flag are for when MPI is not being used at all, i.e. a serial
  ! run. In this case, these variables are still used at times.

  integer, save :: n_procs     = 1 ! number of processes
  integer, save :: rank        = 0 ! rank of process
  integer, save :: MPI_BANK ! MPI datatype for fission bank
  integer, save :: MPI_TALLYRESULT ! MPI datatype for TallyResult
  integer, save :: mpi_err ! error code from MPI calls
  logical, save :: master = .true. ! master process?
  logical, save :: mpi_enabled = .false. ! is MPI in use and initialized?

end module mpi_interface
