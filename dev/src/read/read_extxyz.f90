subroutine read_extxyz(filename, nat, rxyz, alat, atomnames, comment)
  !! reads the position, lattice vecors and atomnames from a file. File format is
  !! detected based on the file ending. Supported endings:
  !! .ascii, .cif, .in, .gen, .qe, .vasp and POSCAR (vasp)
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3, nat), intent(out) :: rxyz
  !! atom positions (bohr)
  real*8, dimension(3, 3), intent(out) :: alat
  !! lattice vectors (bohr)
  character(len=2), intent(out) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(out) :: comment
  integer :: io, ios
  logical :: eof_reached

  open(newunit=io, file=filename, iostat=ios)
  if ( ios /= 0 ) then
    stop 'error opening file read_extxyz'
  end if

  call get_partial_nat_extxyz(io, ios, eof_reached)
  if (eof_reached) then
    stop 'eof_reached'
  end if
  if (ios /= nat) then
    print*, ios, nat
    stop 'conflicting number of atoms read_extxyz'
  end if

  call read_partial_extxyz(io, nat, rxyz, alat, atomnames, comment)

  close(io)

end subroutine read_extxyz


subroutine read_partial_extxyz(io, nat, rxyz, alat, atomnames, comment)
    !! reads xyz file. assumes units of file are in angstrom and converts them to a.u.
    implicit none
    integer, parameter :: lineLength = 500
    character(6), parameter :: all_line_format = '(a500)'
    integer, intent(in) :: io
    !! io unit already opened
    integer, intent(in) :: nat
    !! number of atoms
    real(8), intent(out), dimension(3,nat) :: rxyz
    !! positions of the atoms (bohr)
    real(8), dimension(3, 3) :: alat
    !! lattice vectors in bohr
    character(len=2), intent(out), dimension(nat) :: atomnames
    !! chemical name of the atoms
    character(len=80), intent(out) :: comment
    !! contents of the first comment line

    real(8), parameter :: Bohr_Ang = 0.52917721067
    integer :: ios, i
    character(len=lineLength) :: all_line

    ! read comment line
    read(io, all_line_format, iostat=ios) all_line
    if (ios /= 0) stop "error reading commentline in readextxyz."
    do while(iscommentorempty(all_line))
      read(io, all_line_format, iostat=ios) all_line
      if (ios /= 0) stop "error reading commentline in readextxyz."
    end do
    comment = all_line(1:80)

    alat = getLattice(all_line)

    i = 1
    do while (i <= nat)
      read(io, all_line_format, iostat = ios) all_line
      if ( ios/=0 ) then
        print*, trim(all_line)
        print*, "error reading line: ", i
        stop
      end if
      if (iscommentorempty(all_line)) cycle
      read(all_line,*, iostat = ios) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
      if ( ios/=0 ) then
        print*, trim(all_line)
        print*, "error internal reading line: ", i
        stop
      end if
      i = i + 1
    end do

    alat = alat / Bohr_Ang
    rxyz = rxyz / Bohr_Ang

  contains

    function getLattice(line_in) result(lat)
      implicit none
      character(len=lineLength) :: line_in
      character(len=lineLength) :: line_copy
      real(8) :: lat(3,3)
      integer :: ind
      integer :: nwords
      integer, parameter :: len_words = 50
      character(len=len_words), dimension(:), allocatable :: words

      ind = index(line_in, 'Lattice=') + 9
      if (ind <= 0) stop "no lattice tag present in read_extxyz"
      
      line_copy = line_in(ind:)
      if (index(line_copy, '"') <= 0 .and. index(line_copy, "'") <= 0) stop "quotation mark not closed in readextxyz"
      ind = max(index(line_copy, '"'), index(line_copy, "'"))
      line_copy(ind:ind) = ' '

      call get_num_words(line_copy, nwords, lineLength)
      allocate(words(nwords))
      call get_words(line_copy, nwords, lineLength, words, len_words)

      read(words(1), *, iostat=ios) lat(1,1)
      if (ios /= 0) stop "error parsing lattice element 1 in readextxyz"
      read(words(2), *, iostat=ios) lat(2,1)
      if (ios /= 0) stop "error parsing lattice element 2 in readextxyz"
      read(words(3), *, iostat=ios) lat(3,1)
      if (ios /= 0) stop "error parsing lattice element 3 in readextxyz"
      read(words(4), *, iostat=ios) lat(1,2)
      if (ios /= 0) stop "error parsing lattice element 4 in readextxyz"
      read(words(5), *, iostat=ios) lat(2,2)
      if (ios /= 0) stop "error parsing lattice element 5 in readextxyz"
      read(words(6), *, iostat=ios) lat(3,2)
      if (ios /= 0) stop "error parsing lattice element 6 in readextxyz"
      read(words(7), *, iostat=ios) lat(1,3)
      if (ios /= 0) stop "error parsing lattice element 7 in readextxyz"
      read(words(8), *, iostat=ios) lat(2,3)
      if (ios /= 0) stop "error parsing lattice element 8 in readextxyz"
      read(words(9), *, iostat=ios) lat(3,3)
      if (ios /= 0) stop "error parsing lattice element 9 in readextxyz"
    end function getLattice

    function iscommentorempty(line_in) result(bool)
      implicit none
      character(len=lineLength) :: line_in
      logical :: bool
      bool = .FALSE.
      if ( len_trim(line_in) == 0 ) then !line is empty
        bool = .TRUE.
      end if
      if ( line_in(1:1) == "#" ) then
        bool = .TRUE. ! line is comment
      end if
      return
    end function iscommentorempty

  end subroutine read_partial_extxyz
  
  subroutine get_partial_nat_extxyz(io, nat, end_of_file)
    implicit none
    integer, intent(in) :: io
    integer, intent(out) :: nat
    logical :: end_of_file
    integer, parameter :: lineLength = 500
    character(6), parameter :: all_line_format = '(a500)'
    integer :: ios, i
    character(len=lineLength) :: all_line

    !! get number of atoms line
    read(io, all_line_format, iostat=ios) all_line

    if (ios < 0) then
      end_of_file = .true.
      nat = - 1
      return
    end if
    if (ios > 0) stop "error reading line in readextxyz."

    do while(iscommentorempty(all_line))
      read(io, all_line_format, iostat=ios) all_line
      if (ios < 0) then
        end_of_file = .true.
        nat = - 1
        return
      end if
      if (ios > 0) stop "error reading line in readextxyz."
    end do

    read(all_line, *, iostat=ios) nat
    if(ios /= 0) stop 'error parsing nat in readextxyz'
    end_of_file = .false.

    contains
    function iscommentorempty(line_in) result(bool)
      implicit none
      character(len=lineLength) :: line_in
      logical :: bool
      bool = .FALSE.
      if ( len_trim(line_in) == 0 ) then !line is empty
        bool = .TRUE.
      end if
      if ( line_in(1:1) == "#" ) then
        bool = .TRUE. ! line is comment
      end if
      return
    end function iscommentorempty
  end subroutine get_partial_nat_extxyz

  subroutine get_nat_extxyz(filename, nat)
    !! counts the number of atoms in a .xyz file
    implicit none
    character(len=*), intent(in) :: filename
    !! filename that will be parsed
    integer, intent(out) :: nat
    !! number of atoms
    integer :: io, ios
    open(newunit=io,file=filename,iostat=ios, status="old")
    if(ios/=0) stop "error opening input file get_nat_extxyz"
    read(io,*) nat
    close(io)
  end subroutine get_nat_extxyz
  