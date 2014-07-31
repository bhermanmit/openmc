module mpi_interface

#ifdef MPI
  use mpi
#endif

  implicit none

  integer, save :: mpi_err ! error code from MPI calls
  logical :: master = .true. ! master process?

end module mpi_interface
