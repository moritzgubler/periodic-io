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

  integer, parameter :: natx_sphere=120, ns=1, np=1
  real(8) :: fp_dist
  real(8), allocatable, dimension(:,:) :: fp1, fp2
  real(8), parameter :: width_cutoff = 3.d0

  call get_command_argument(1,file1, status=stat)
  if ( stat/= 0 ) then
    stop "first argument must contain input filename"
  end if
  call get_command_argument(2,file2, status=stat)
  if ( stat /= 0 ) then
    stop "second argument must contain output filename"
  end if

  call get_nat_periodic(file1, nat1)
  call get_nat_periodic(file2, nat2)
  if (nat1 /= nat2) stop 'input files contain different number of atoms'
  allocate(r1(3, nat1), r2(3,nat2), symb1(nat1), symb2(nat2))
  call read_periodic(trim(file1), nat1, r1, lat1, symb1, comment)
  call read_periodic(trim(file2), nat2, r2, lat2, symb2, comment)

  allocate(fp1((ns + 3*np)*natx_sphere, nat1), fp2((ns + 3*np)*natx_sphere, nat2))

  call back2cell(nat1, r1, lat1)
  call fingerprint(nat1, natx_sphere, ns, np, width_cutoff, lat1, r1, symb1, fp1)

  call back2cell(nat2, r2, lat2)
  call fingerprint(nat2, natx_sphere, ns, np, width_cutoff, lat2, r2, symb2, fp2)

  call calc_fpd(nat1, natx_sphere, ns, np, fp1, fp2, fp_dist)

  print*, fp_dist


end program fingerprint_distance
