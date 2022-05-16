program transform_periodic
  !! reads a periodic file calculates a supercell and writes it into new file
  implicit none
  integer :: nat, nsuper
  real*8, allocatable, dimension(:,:) :: rxyz, rsuper
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: alat, asuper
  !! lattice vectors (bohr)
  character(len=250) :: file_in, file_out
  !! input/ output file
  character(len=2), allocatable, dimension(:) :: atomnames, atoms_super
  !! chemical name of all the atoms
  character(len=80) :: comment
  !! contents of first comment line

  integer :: stat, i, j, k, iat, dims(3), scounter
  character(len=10) :: superstring


  do i = 1, 3, 1
    call get_command_argument(i, superstring, status=stat)
    if (stat /= 0 ) stop 'error reading superstring'
    read(superstring, *, iostat = stat) dims(i)
    if (stat /= 0) stop 'error parsing superstring'
    if ( dims(i) <= 0 ) then
      stop ' error negative supercell size given'
    end if
  end do

  call get_command_argument(4,file_in, status=stat)
  if ( stat /= 0 ) then
    stop "4th argument must contain input filename"
  end if
  call get_command_argument(5,file_out, status=stat)
  if ( stat /= 0 ) then
    stop "5th argument must contain output filename"
  end if

  ! get positions
  call get_nat_periodic(file_in, nat)
  allocate(rxyz(3,nat), atomnames(nat))
  call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)

  !calculate supercell
  nsuper = product(dims) * nat
  allocate(rsuper(3, nsuper), atoms_super(nsuper))

  asuper(:, 1) = dims(1) * alat(:, 1)
  asuper(:, 2) = dims(2) * alat(:, 2)
  asuper(:, 3) = dims(3) * alat(:, 3)

  scounter = 0

  ! loop over all relevant periodic images
  do i = 0, dims(1) - 1, 1
    do j = 0, dims(2) - 1 , 1
      do k = 0, dims(3) - 1 , 1
        do iat = 1, nat, 1
          rsuper(:,scounter*nat + iat) = rxyz(:, iat) + i * alat(:, 1) &
            + j * alat(:, 2) + k * alat(:, 3)
        end do
        atoms_super(scounter*nat + 1 : (scounter + 1)*nat) = atomnames
        scounter = scounter + 1
      end do
    end do
  end do

  call write_periodic(trim(file_out), nsuper, rsuper, asuper, atoms_super, comment)

end program transform_periodic
