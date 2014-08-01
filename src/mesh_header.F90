module mesh_header

  use dict_header, only: DictIntInt

  implicit none

!===============================================================================
! STRUCTUREDMESH represents a tessellation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type StructuredMesh
    integer :: id                          ! user-specified id
    integer :: type                        ! rectangular, hexagonal
    integer :: n_dimension                 ! rank of mesh
    real(8) :: volume_frac                 ! volume fraction of each cell
    integer, allocatable :: dimension(:)   ! number of cells in each direction
    real(8), allocatable :: lower_left(:)  ! lower-left corner of mesh
    real(8), allocatable :: upper_right(:) ! upper-right corner of mesh
    real(8), allocatable :: width(:)       ! width of each mesh cell
  end type StructuredMesh

  ! Meshes
  integer, save :: n_meshes       = 0 ! # of structured meshes
  integer, save :: n_cmfd_meshes  = 1 ! # of structured meshes
  integer, save :: n_user_meshes  = 0 ! # of structured user meshes
  type(DictIntInt), save :: mesh_dict
  type(StructuredMesh), save, allocatable, target :: meshes(:)

  ! Entropy mesh
  type(StructuredMesh), save, pointer :: entropy_mesh => null()

  ! Uniform fission source weighting
  logical, save :: ufs = .false.
  type(StructuredMesh), save, pointer :: ufs_mesh => null()
  real(8), save, allocatable :: source_frac(:,:,:,:)

end module mesh_header
