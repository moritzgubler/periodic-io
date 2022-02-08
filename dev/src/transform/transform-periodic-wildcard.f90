program transform_periodic_wildcard
  !! reads multiple periodic files and converts it to the output format
  implicit none
  integer :: nat
  !! number of atoms
  real*8, allocatable, dimension(:,:) :: rxyz
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: alat
  !! lattice vectors (bohr)
  character(len=250) :: file_in, file_out, outdir, basename
  character(len=2), allocatable, dimension(:) :: atomnames
  character(len=80) :: comment, output_format
  integer :: numargs, ios, i

  numargs = command_argument_count()

  !! check number of arguments
  if (numargs < 3) then
    call usage
    stop
  end if

  !! check if first argument contains output format
  call get_command_argument(1, output_format, status=ios)
  if (ios /= 0 ) stop "error parsing first argument."
  if ( index(output_format, "-") /= 1 ) then
    call usage
    stop "first argument doesnt start with -"
  end if
  output_format(1:1) = " "
  output_format = adjustl(output_format)

  !! get output directory
  call get_command_argument(numargs, outdir, status=ios)
  if (ios /= 0) stop "error reading outdir command argument."
  call execute_command_line("mkdir -p "//trim(adjustl(outdir)), exitstat=ios)
  if (ios /= 0) stop "outdir not a directory or not writable"
  if (outdir(len_trim(outdir):len_trim(outdir)) /= "/") outdir=trim(outdir)//"/"

  do i = 2, numargs-1, 1
    call get_command_argument(i, file_in, status = ios)
    if (ios /= 0) stop "error getting argument"

    call get_nat_periodic(file_in, nat)
    allocate(rxyz(3,nat), atomnames(nat))
    call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)

    call get_basename(file_in, basename)
    file_out = trim(outdir)//trim(basename)//trim(output_format)
    call write_periodic(trim(file_out), nat, rxyz, alat, atomnames , comment)
    deallocate(rxyz, atomnames)
  end do

contains
  subroutine usage
    print*, "This program  transforms periodic file format"
    print*, "Usage:"
    print*, "transform-periodic-wildcard -output_format list_of_input_files destination"
    print*, "Example:"
    print*, "transform-periodic-wildcard -cif *.ascii *.in data/destination"
    print*, "this reads in all .ascii and .in files from the current working directory,"
    print*, "transforms them to cif files and saves them in the directory data/destination."
  end subroutine usage

end program transform_periodic_wildcard
