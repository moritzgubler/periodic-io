program fingerprint_distance
  !! reads two strunctures and calculates their fingerprint distance based on the ovelap matrix
  implicit none
  integer :: nat1, nat2
  real*8, allocatable, dimension(:,:) :: r1, r2
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: lat1, lat2
  !! lattice vectors (bohr)
  character(len=250) :: file1, file2
  !! input/ output file
  character(len=2), allocatable, dimension(:) :: symb1, symb2
  !! chemical name of all the atoms
  character(len=80) :: comment
  !! contents of first comment line
  integer :: stat

  logical :: is_periodic = .TRUE., f1periodic = .TRUE., f2periodic = .TRUE.
  character(len=10) :: filetype
  integer, parameter :: natx_sphere = 120, ns = 1, np = 1
  real(8) :: fp_dist
  real(8), allocatable, dimension(:,:) :: fp1, fp2
  real(8), parameter :: width_cutoff = 3.d0

  call get_command_argument(1,file1, status=stat)
  if ( stat/= 0 ) call help
  call get_command_argument(2,file2, status=stat)
  if ( stat /= 0 ) call help

  call get_file_type(file1, filetype)
  if ( trim(filetype) == 'xyz' ) f1periodic = .FALSE.
  call get_file_type(file2, filetype)
  if ( trim(filetype) == 'xyz' ) f2periodic = .FALSE.

  if ( f1periodic .neqv. f2periodic ) then
    print*, 'You are trying to compare a periodic file format with a non periodic one. This makes no sense.'
    print*, ''
    call help
  end if
  is_periodic = f1periodic

  if ( is_periodic ) then
    call get_nat_periodic(file1, nat1)
    call get_nat_periodic(file2, nat2)
  else
    call get_nat_xyz(file1, nat1)
    call get_nat_xyz(file2, nat2)
  end if

  if (nat1 /= nat2) then
    print*, 'input files contain different number of atoms'
    print*, ''
    call help
  end if
  allocate(r1(3, nat1), r2(3,nat2), symb1(nat1), symb2(nat2))
  if ( is_periodic ) then
    call read_periodic(trim(file1), nat1, r1, lat1, symb1, comment)
    call read_periodic(trim(file2), nat2, r2, lat2, symb2, comment)
  else
    call read_xyz(file1, nat1, r1, symb1, comment)
    call read_xyz(file2, nat2, r2, symb2, comment)
    lat1 = 0.d0
    lat2 = 0.d0
  end if


  allocate(fp1((ns + 3*np)*natx_sphere, nat1), fp2((ns + 3*np)*natx_sphere, nat2))

  if ( is_periodic ) then
    call back2cell(nat1, r1, lat1)
  end if
  call fingerprint(nat1, natx_sphere, ns, np, width_cutoff, lat1, r1, symb1, fp1)

  if ( is_periodic ) then
    call back2cell(nat2, r2, lat2)
  end if
  call fingerprint(nat2, natx_sphere, ns, np, width_cutoff, lat2, r2, symb2, fp2)

  call calc_fpd(nat1, natx_sphere, ns, np, fp1, fp2, fp_dist)

  print*, fp_dist

contains
  subroutine help
    implicit none
    print*, ''
    print*, 'This is the fingerprint-distance program of the periodicIO library.&
    & It calculates the fingerprint distance structure of two structures with &
    & both periodic and free boundary conditions. It can read and compare all the&
    & file formats supported by the periodicIO library.'
    print*, ''
    print*, 'Usage: fingerprint-distance file1 fil2'
    print*, ''
    stop
  end subroutine help
end program fingerprint_distance
