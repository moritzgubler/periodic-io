program cif2ascii
  !! reads an ascii file and converts it to the .in format
  implicit none
  integer :: nat
  real*8, allocatable, dimension(:,:) :: rxyz
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: alat
  character(len=250) :: file_in, file_out
  character(len=2), allocatable, dimension(:) :: atomnames
  character(len=80) :: comment

  call get_command_argument(1,file_in)
  if ( len_trim(file_in) == 0 ) then
    stop "first argument must contain .cif filename"
  end if
  call get_command_argument(2,file_out)
  if ( len_trim(file_out) == 0 ) then
    stop "second argument must contain .ascii filename"
  end if
  call get_nat_cif(file_in, nat)
  allocate(rxyz(3,nat), atomnames(nat))
  call read_cif(trim(file_in), nat, rxyz, alat, atomnames, comment)
  call write_ascii(trim(file_out), nat, rxyz, alat, atomnames, comment)
end program cif2ascii
