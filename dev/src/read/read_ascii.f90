subroutine read_ascii(filename, nat, rxyz, alat, atomnames, comment)
  !! reads an ascii file with the specified filename units are assumed to be angstroem.
  !! units are converted to hartree units before they are returned
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is read
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3, nat), intent(out) :: rxyz
  !! atom positions (in bohr)
  real*8, dimension(3, 3), intent(out) :: alat
  !! lattice vectors in bohr
  character(len=2), intent(out) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(out) :: comment
  !! string that was read from comment line (first line in ascii format)
  character(len=250) :: all_line
  !! sting containing entire line
  integer :: i, io, ios
  real*8 :: alat_temp(2, 3)
  real*8 :: Bohr_Ang = 0.52917721067
  !! Bohr to angstrom conversion factor.
  logical :: is_reduced
  logical :: is_in_bohr
  character(len=250):: key_string
  real(8) :: xyzred(3, nat)

  is_reduced = .false.
  is_in_bohr = .false.

  open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status="old")
  if (ios /= 0) then
    print *, "error opening ascii file: "//filename
    stop
  end if
  read (io, "(a80)", iostat=ios) comment
  if (ios /= 0) then
    print *, trim(adjustl(filename)), ios
    stop "error reading file "
  end if
  i = 1
  alat = 0
  do while ( i < 3 )
    read(io,"(a250)", iostat=ios) all_line
    if (ios /= 0) stop "error reading lattice in moleculario"
    all_line = adjustl(all_line)
    if ( len_trim(all_line) == 0 ) cycle
    call parse_keyword
    if ( index(all_line, "#") == 1 ) cycle
    read(all_line, *, iostat=ios) alat_temp(i,:)
    if (ios /= 0) stop "error parsing lattice in moleculario"
    i = i + 1
  end do
  alat = 0.0
  alat(1, 1) = alat_temp(1, 1)
  alat(1, 2) = alat_temp(1, 2)
  alat(2, 2) = alat_temp(1, 3)
  alat(1, 3) = alat_temp(2, 1)
  alat(2, 3) = alat_temp(2, 2)
  alat(3, 3) = alat_temp(2, 3)
  i = 1
  do while (i <= nat)
    read(io,"(a250)", iostat=ios) all_line
    all_line = adjustl(all_line)
    if ( ios/= 0 ) then
      print*, all_line
      stop "error reading line in read_ascii"
    end if
    if ( len_trim(all_line) == 0 ) cycle
    call parse_keyword
    if ( index(all_line, "#") == 1 ) cycle
    read (all_line, *, iostat=ios) rxyz(1, i), rxyz(2, i), rxyz(3, i), atomnames(i)
    if ( ios/= 0 ) then
      print*, all_line
      stop "error parsing line in read_ascii"
    end if
    i = i + 1
  end do
  close (io)

  if ( is_reduced ) then
    xyzred = rxyz
    if ( is_in_bohr ) then
      call frac2cart(nat, alat, xyzred, rxyz)
    else
      alat = alat / Bohr_Ang
      call frac2cart(nat, alat, xyzred, rxyz)
    end if
  elseif (.not. is_in_bohr) then
    alat = alat / Bohr_Ang
    rxyz = rxyz / Bohr_Ang
  end if

  contains
  subroutine parse_keyword()
    if ( index(all_line, '#keyword:') == 1 ) then ! keyword line found
      key_string = all_line(11:)
      if ( trim(key_string) == 'reduced' ) is_reduced = .true.
      if ( trim(key_string) ==  'bohr' ) is_in_bohr = .true.
    end if  
  end subroutine parse_keyword

end subroutine read_ascii

subroutine get_nat_ascii(filename, nat)
  !! counts the number of atoms in an *.ascii file.
  !! nat doesn't need to be printed on first line.
  implicit none
  character(len=*), intent(in) :: filename
  !! name of the ascii file wich will be read.
  integer, intent(out) :: nat
  !! number of atoms
  ! private variables
  integer :: ios, io
  real*8 :: place(3)
  character(len=250) :: all_line
  character(len=2) :: atname

  open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status="old")
  if (ios /= 0) then
    print *, "Error opening ascii file to get number of atoms from. Filename: "//filename
    stop
  end if
  read (io, '(a250)', iostat=ios) all_line
  if (ios /= 0) then
    print *, "ios", ios
    print *, "filename:_", filename
    stop "Error reading file in getnat "!//filename
  end if
  nat = 1
  do while (nat < 3)
    read (io, *, iostat=ios) place(1), place(2), place(3)
    if (ios > 0) then
      cycle
    end if
    if (ios < 0) stop "end of file in get nat ascii"!//filename
    nat = nat + 1
  end do
  nat = 0
  do
    read (io, *, iostat=ios) place(1), place(2), place(3), atname
    if (ios > 0) cycle
    if (ios < 0) exit
    nat = nat + 1
  end do
  close (io)
  !nat = 44
end subroutine get_nat_ascii
