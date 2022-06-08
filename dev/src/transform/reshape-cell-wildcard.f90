program transform_periodic
  !! this program reads a periodic structure, tries to correct its shape
  !! and writes it to the filename given in the second filename
  implicit none
  integer :: nat
  real*8, allocatable, dimension(:,:) :: rxyz
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: alat
  !! lattice vectors (bohr)
  character(len=250) :: file_in
  !! input/ output file
  character(len=2), allocatable, dimension(:) :: atomnames
  !! chemical name of all the atoms
  character(len=80) :: comment
  !! contents of first comment line

  integer :: i, numargs

  numargs = command_argument_count()
  if ( numargs <= 0 ) then
    print*, 'This program reads all the files given as command line arguments, reshapes the cell'
    print*, 'if possible and overwrites each file with the new geometry. Works also with regex.'
    stop ''
  end if
  do i = 1, numargs
    call get_command_argument(i,file_in)
    call get_nat_periodic(file_in, nat)
    allocate(rxyz(3,nat), atomnames(nat))
    call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)

    call reshapecell(nat, alat, rxyz)
    call write_periodic(trim(file_in), nat, rxyz, alat, atomnames , comment)

    deallocate(rxyz, atomnames)
  end do

end program transform_periodic
