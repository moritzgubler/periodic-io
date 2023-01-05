subroutine read_extxyz(io, nat, rxyz, alat, atomnames, comment)
    !! reads xyz file. assumes units of file are in angstrom and converts them to a.u.
    implicit none
    integer, parameter :: lineLength = 500
    character(6), parameter :: all_line_format = '(a500)'
    integer, intent(in) :: io
    !! io unit already opened
    integer, intent(inout) :: nat
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


    !! get number of atoms line
    read(io, all_line_format, iostat=ios) all_line
    if (ios /= 0) stop "error reading line in readextxyz."
    do while(.not. iscommentorempty(all_line))
      read(io, all_line_format, iostat=ios) all_line
      if (ios /= 0) stop "error reading line in readextxyz."
    end do
    read(all_line, *, iostat=ios) i
    if(ios /= 0) stop 'error parsing nat in readextxyz'
    if (i /= nat) stop "provided nat and nat from file differ"


    ! read comment line
    read(io, all_line_format, iostat=ios) all_line
    if (ios /= 0) stop "error reading line in readextxyz."
    do while(.not. iscommentorempty(all_line))
      read(io, all_line_format, iostat=ios) all_line
      if (ios /= 0) stop "error reading line in readextxyz."
    end do
    comment = all_line(1:80)

    alat = getLattice(all_line)

    i = 1
    do while (i <= nat)
      read(io, all_line_format, iostat = ios) all_line
      if ( ios/=0 ) then
        print*, "error reading line: ", i
        stop
      end if
      if (iscommentorempty(all_line)) cycle
      read(all_line,*, iostat = ios) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
      if ( ios/=0 ) then
        print*, "error reading line: ", i
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

      ind = index(line_in, 'Lattice=') + 8
      if (ind <= 0) stop "no lattice tag present in read_extxyz"
      
      line_copy = line_in(ind:)
      if (index(line_copy, '"') <= 0 .and. index(line_copy, "'") <= 0) stop "quotation mark not closed in readextxyz"
      ind = max(index(line_copy, '"'), index(line_copy, "'"))
      line_copy(ind:ind) = ' '

      call get_num_words(line_copy, nwords, lineLength)
      allocate(words(nwords))
      call get_words(line_copy, nwords, lineLength, words, len_words)

      read(words(1), *, iostat=ios) alat(1,1)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(2), *, iostat=ios) alat(2,1)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(3), *, iostat=ios) alat(3,1)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(4), *, iostat=ios) alat(1,2)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(5), *, iostat=ios) alat(2,2)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(6), *, iostat=ios) alat(3,2)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(7), *, iostat=ios) alat(1,3)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(8), *, iostat=ios) alat(2,3)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
      read(words(9), *, iostat=ios) alat(3,3)
      if (ios /= 0) stop "error parsing lattice element in readextxyz"
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

  end subroutine read_extxyz
  
  subroutine get_nat_extxyz(filename, nat)
    !! counts the number of atoms in a .xyz file
    implicit none
    character(len=*), intent(in) :: filename
    !! filename that will be parsed
    integer, intent(out) :: nat
    !! number of atoms
    integer :: io, ios
    open(newunit=io,file=filename,iostat=ios, status="old")
    if(ios/=0) stop "error opening input file"
    read(io,*) nat
    close(io)
  end subroutine get_nat_extxyz
  