program ascii2in
  !! reads an ascii file and converts it to the .in format
  implicit none
  integer :: nat
  real*8, allocatable, dimension(:,:) :: rxyz
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: alat
  character(len=250) :: file_in, file_out
  character(len=2), allocatable, dimension(:) :: atomnames
  character(len=80) :: comment
  character(len=10) :: filetype
  integer :: i, num_args

  num_args = command_argument_count()
  if ( num_args == 0 ) then
    print*, "Converts ascii files to .in files. Command line arguments:"
    print*, "List of .ascii files. Wildcards are allowed or src and dest where dest is a .in file"
    print*, "The new file is generated be replacing .ascii by .in"
    print*, "All files will be written in the current working directory."
    stop
  end if
  if ( num_args == 2 ) then
    call get_command_argument(2,file_out)
    call get_file_type(file_out, filetype)
    if ( trim(filetype) == "in" ) then
      call get_command_argument(1,file_in)
      call get_nat_ascii(file_in, nat)
      allocate(rxyz(3,nat), atomnames(nat))
      call read_ascii(trim(file_in), nat, rxyz, alat, atomnames, comment)
      call write_in(trim(file_out), nat, rxyz, alat, atomnames, comment)
      !print*, trim(file_in)//" "//trim(file_out)
      deallocate(rxyz, atomnames)
      stop
    end if
  end if

  i = 1
  do
    call get_command_argument(i,file_in)
    if ( len_trim(file_in) == 0 ) then
      exit
    end if
    call get_file_type(file_in, filetype)
    if ( trim(filetype) /= "ascii" ) then
      print*, "wrong file ending for file (.ascii expected): "//trim(file_in)
      print*, "skipping"
      i = i + 1
      cycle
    end if
    call get_basename(file_in, file_out)
    file_out = trim(file_out)//".in"
    call get_nat_ascii(file_in, nat)
    allocate(rxyz(3,nat), atomnames(nat))
    call read_ascii(trim(file_in), nat, rxyz, alat, atomnames, comment)
    call write_in(trim(file_out), nat, rxyz, alat, atomnames, comment)
    !print*, trim(file_in)//" "//trim(file_out)
    deallocate(rxyz, atomnames)
    i = i + 1
  end do

end program ascii2in
