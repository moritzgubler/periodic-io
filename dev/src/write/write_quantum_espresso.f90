
subroutine write_quantum_espresso(filename, nat, rxyz, alat, atomnames)
  !! This subroutine appends cell parameters and atomic positions to the file specified in
  !! the filename variable. It is useful to create input files for quantum espresso
  use, intrinsic :: iso_fortran_env
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real(8), dimension(3, nat), intent(in) :: rxyz
  !! atom positions in bohr
  real(8), dimension(3, 3), intent(in) :: alat
  !! lattice vectors
  character(len=2), intent(in) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.

  integer :: ntypes
  character(len=2), dimension(:), allocatable :: atom_types
  integer :: io, ios
  integer :: i

  call get_atom_types

  call set_nat_qe(filename, nat)

  open(file=filename, newunit=io, iostat=ios, position="append")
  if (ios /= 0) stop "error reading file in write_quantum_espresso (libmoleculario)"

  write(io, *) "CELL_PARAMETERS bohr"
  do i = 1, 3, 1
    write(io,*) alat(1, i), alat(2, i), alat(3, i)
  end do

  write(io,*) "ATOMIC_POSITIONS bohr"
  do i = 1, nat, 1
    write(io,*) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
  end do

  close(io)
contains
  !! deprecated but should work. The subroutine write_quantum_espresso does and
  !! should not call this subroutine.
  subroutine get_atom_types
    integer :: ii, jj
    character(len=2), dimension(nat) :: an_copy

    ntypes = 0
    an_copy = atomnames
    do ii = 1, nat, 1
      if (len_trim(an_copy(ii)) /= 0) then ! new atom species found. delete all the following from list
        ntypes = ntypes + 1
        do jj = ii + 1, nat, 1
          if (an_copy(jj) == an_copy(ii)) then
            an_copy(jj) = ""
          end if
        end do
      end if
    end do
    if (allocated(atom_types)) deallocate(atom_types)
    allocate (atom_types(ntypes))
    jj = 0
    do ii = 1, nat, 1
      if (an_copy(ii) /= "") then
        jj = jj + 1
        atom_types(jj) = an_copy(ii)
      end if
    end do
  end subroutine get_atom_types
end subroutine write_quantum_espresso

subroutine set_nat_qe(fname, nat)

  implicit none
  character(len=300), intent(in) :: fname
  integer, intent(in) :: nat
  integer :: i, stat, io, nlines
  character(len=300) :: all_line
  character(len=300), dimension(:), allocatable :: all_file
  logical :: nat_found
  
  nat_found = .FALSE.
  
  nlines = count_lines(fname)
  allocate(all_file(nlines))
  
  open(file=fname, newunit=io, iostat=stat)
  if ( stat /= 0 ) then
    stop 'error opening input file'
  end if
  
  i = 0
  do
    read(io, '(a300)', iostat =stat) all_line
    if ( stat /=0 ) exit
    i = i + 1
    all_file(i) = all_line
    call lower_case(all_line)
    if ( index(all_line, 'nat') > 0 ) then
      nat_found = .TRUE.
      write(all_file(i), '(a, i5)') 'nat = ', nat
    end if
  end do
  
  close(io)
  
  open(file=fname, newunit=io, iostat=stat)
  
  i = 1
  do
    write(io, '(a300)') all_file(i)
    i = i + 1
  end do
  close(io)
  
contains
  
  function count_lines(filename) result(nlinesp)
    implicit none
    character(len=*) :: filename
    integer :: nlinesp, ipers, iop
    open(newunit=ipers,file=filename, status="old", iostat=iop)
    if (iop /= 0 ) then
      print*, trim(filename)
      stop "error opening file in countlines, keysystm"
    end if
    nlinesp = 0
    do
    read(ipers,*,iostat=iop)
    if (iop/=0) exit
    nlinesp = nlinesp + 1
    end do
    close(ipers)
  end function count_lines

end subroutine set_nat_qe
