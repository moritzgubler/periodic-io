program main

  implicit none
  integer :: nat1, nat2
  !! positions in au
  real*8, allocatable, dimension(:,:) :: rxyz
  !! alat in au
  real*8, dimension(3,3) :: alat
  integer, parameter :: natx_sphere=120, ns=1, np=1
  real*8 :: fp_dist
  real*8, allocatable, dimension(:,:) :: fp1, fp2
  character(len=250) :: file1, file2
  character(len=2), allocatable, dimension(:) :: atomnames
  real*8, allocatable, dimension(:) :: rcov
  character(80) :: comment
  integer :: i

  call get_command_argument(1,file1)
  if ( len_trim(file1) == 0 ) then
    stop "first argument must contain first filename"
  end if
  call get_command_argument(2,file2)
  if ( len_trim(file2) == 0 ) then
    stop "second argument must contain second filename"
  end if
  call get_nat_periodic(file1, nat1)
  call get_nat_periodic(file2, nat2)
  if ( nat1 /= nat2 ) then
    stop "both input files must contain the same number of atoms."
  end if
  allocate(rxyz(3,nat1))
  allocate(fp1((ns + 3*np)*natx_sphere, nat1), fp2((ns + 3*np)*natx_sphere, nat1))
  allocate(atomnames(nat1))
  allocate(rcov(nat1))
  call read_periodic(file1, nat1, rxyz, alat, atomnames, comment)
  do i = 1, nat1, 1
    call sym2rcov(atomnames(i), rcov(i))
  end do
  call back2cell(nat1, rxyz, alat)
  call fingerprint_eval(nat1, natx_sphere, ns, np, alat, rxyz, rcov, fp1)
  call read_periodic(file2, nat2, rxyz, alat, atomnames, comment)
  do i = 1, nat1, 1
    call sym2rcov(atomnames(i), rcov(i))
  end do
  call back2cell(nat2, rxyz, alat)
  call fingerprint_eval(nat1, natx_sphere, ns, np, alat, rxyz, rcov, fp2)
  call calc_fpd(nat1, natx_sphere, ns, np, fp1, fp2, fp_dist)
  !print"(a,f10.6)", "the two structures have a fingerprint distance of: ", fp_dist
  print*, fp_dist
end program main
