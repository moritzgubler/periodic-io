program back2cell_wildcard
  !! reads multiple periodic files and converts it to the output format
  implicit none
  integer :: nat
  !! number of atoms
  real*8, allocatable, dimension(:,:) :: rxyz
  !! atomic positions always in bohr.
  real*8, dimension(3,3) :: alat
  !! lattice vectors (bohr)
  character(len=300) :: file_in, file_out, outdir, basename
  character(len=2), allocatable, dimension(:) :: atomnames
  character(len=80) :: comment
  integer :: numargs, ios, i

  numargs = command_argument_count()

  !! check number of arguments
  if (numargs < 2) then
    call usage
    stop
  end if


  !! get output directory
  call get_command_argument(numargs, outdir, status=ios)
  if (ios /= 0) stop "error reading outdir command argument."
  call execute_command_line("mkdir -p "//trim(adjustl(outdir)), exitstat=ios)
  if (ios /= 0) stop "outdir not a directory or not writable"
  if (outdir(len_trim(outdir):len_trim(outdir)) /= "/") outdir=trim(outdir)//"/"

  do i = 1, numargs-1, 1
    call get_command_argument(i, file_in, status = ios)
    if (ios /= 0) stop "error getting argument"

    call get_nat_periodic(file_in, nat)
    allocate(rxyz(3,nat), atomnames(nat))
    call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)

    call back2cell(nat, rxyz, alat)

    call get_basename(file_in, basename)
    file_out = trim(outdir)//trim(basename)// '-back2cell.ascii'
    call write_periodic(trim(file_out), nat, rxyz, alat, atomnames , comment)
    deallocate(rxyz, atomnames)
  end do

contains
  subroutine usage
    print*, "This program puts atoms back into the unit cell."
    print*, "Usage:"
    print*, "back2cell-wildcard list_of_input_files destination"
    print*, "Example:"
    print*, "back2cell-wildcard *.ascii *.in data/destination"
    print*, "this reads in all .ascii and .in files from the current working directory,"
    print*, "puts the atoms back into the cell and saves them in the directory data/destination."
  end subroutine usage

end program back2cell_wildcard
