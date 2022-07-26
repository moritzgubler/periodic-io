program set_comment
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
  integer :: i

  call get_command_argument(1,file_in, status = i)
  if ( i /= 0 ) then
    stop "first argument must contain input filename"
  end if

  call get_nat_periodic(file_in, nat)
  allocate(rxyz(3,nat), atomnames(nat))
  call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)

  call get_command_argument(2,comment, status=i)
  if ( i/= 0 ) then
    stop "second argument must contain comment"
  end if

  call write_periodic(trim(file_in), nat, rxyz, alat, atomnames , comment)
end program set_comment
