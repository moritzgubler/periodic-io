program periodic2qe
  !! reads a periodic file and converts it to the qe format
  implicit none
  integer :: nat
  real*8, allocatable, dimension(:,:) :: rxyz
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: alat
  !! lattice vectors (bohr)
  character(len=250) :: file_in, file_out
  !! input/ output file
  character(len=2), allocatable, dimension(:) :: atomnames
  !! chemical name of all the atoms
  character(len=80) :: comment
  !! contents of first comment line

  call get_command_argument(1,file_in)
  if ( len_trim(file_in) == 0 ) then
    stop "first argument must contain input filename"
  end if
  call get_command_argument(2,file_out)
  if ( len_trim(file_out) == 0 ) then
    stop "second argument must contain output filename"
  end if
  call get_nat_periodic(file_in, nat)
  allocate(rxyz(3,nat), atomnames(nat))
  call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)
  call write_quantum_espresso(file_out, nat, rxyz, alat, atomnames , comment)
end program periodic2qe