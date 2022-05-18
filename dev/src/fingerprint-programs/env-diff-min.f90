program main
  !! Calculates the index of the atom of the the two most different atomic environments
  !! and prints it to the standard output.
  implicit none
  integer :: nat1, nat2
  !! positions in au
  real*8, allocatable, dimension(:,:) :: rxyz
  !! alat in au
  real*8, dimension(3,3) :: alat
  integer, parameter :: natx_sphere=400, ns=1, np=1
  real*8 :: env_dist, env_distx
  real*8, allocatable, dimension(:,:) :: fp1, fp2
  character(len=250) :: file1, file2
  character(len=2), allocatable, dimension(:) :: atomnames1, atomnames2
  real*8, allocatable, dimension(:) :: rcov
  character(80) :: comment
  integer :: i, j, ix, jx
  real(8), parameter :: width_cutoff = 3.d0

  env_distx = huge(env_dist)

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

  allocate(rxyz(3,nat1))
  allocate(fp1((ns + 3*np)*natx_sphere, nat1), fp2((ns + 3*np)*natx_sphere, nat2))
  allocate(atomnames1(nat1))
  allocate(rcov(nat1))
  call read_periodic(file1, nat1, rxyz, alat, atomnames1, comment)

  call back2cell(nat1, rxyz, alat)
  call fingerprint(nat1, natx_sphere, ns, np, width_cutoff, alat, rxyz, atomnames1, fp1)

  deallocate(rxyz, rcov)
  allocate(rxyz(3,nat2))
  allocate(atomnames2(nat2))
  allocate(rcov(nat2))

  call read_periodic(file2, nat2, rxyz, alat, atomnames2, comment)
  call back2cell(nat2, rxyz, alat)

  call fingerprint(nat2, natx_sphere, ns, np, width_cutoff, alat, rxyz, atomnames2, fp2)

  do i = 1, nat1, 1
    do j = 1, nat2, 1
      env_dist = norm2(fp1(:,i) - fp2(:,j))
      if ( env_dist < env_distx ) then
        env_distx = env_dist
        ix = i
        jx = j
      end if
    end do
  end do

  write(*,'(a, a, a, i0)') 'Atom-1: ', trim(atomnames1(ix)), ' index: ', ix
  write(*,'(a, a, a, i0)') 'Atom-2: ', trim(atomnames1(jx)), ' index: ', jx
  write(*,'(a, g0.8)') 'Distance: ', env_distx
end program main
