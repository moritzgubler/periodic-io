subroutine read_periodic(filename, nat, rxyz, alat, atomnames, comment)
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
  !! content of comment line
  character(len=10) :: filetype
  integer :: l, vaspind

  vaspind = index(filename, 'POSCAR')
  l = len_trim(filename)
  !if (vaspind > 0) print*, 'vaspind', l -vaspind
  if ( l - vaspind == 5 .and. vaspind > 0) then !! vasp file
    call read_vasp(filename, nat, rxyz, alat, atomnames, comment)
    return
  end if


  call get_file_type(filename, filetype)

  select case (trim(filetype))
  case ("ascii")
    call read_ascii(filename, nat, rxyz, alat, atomnames, comment)
  case("cif")
    call read_cif(filename, nat, rxyz, alat, atomnames, comment)
  case("in")
    call read_in(filename, nat, rxyz, alat, atomnames, comment)
  case("next_step")
    call read_in(filename, nat, rxyz, alat, atomnames, comment)
  case("gen")
    call read_dftb(filename, nat, rxyz, alat, atomnames, comment)
  case("qe")
    call read_quantum_espresso(filename, nat, rxyz, alat, atomnames)
    comment = ''
  case('vasp')
    call read_vasp(filename, nat, rxyz, alat, atomnames, comment)
  case default
    print*, "filetype, ", filetype
    stop "unknown filetype read_periodic"
  end select

end subroutine read_periodic

subroutine get_nat_periodic(filename, nat)
  !! counts the number of atoms of a vasp file
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is read
  integer, intent(out) :: nat
  !! number of atoms
  character(len=10) :: filetype
  integer :: l, vaspind
  vaspind = index(filename, 'POSCAR')
  l = len_trim(filename)
  !if (vaspind > 0) print*, 'vaspind', l -vaspind
  if ( l - vaspind == 5 .and. vaspind > 0) then !! vasp file
    call get_nat_vasp(filename, nat)
    return
  end if

  call get_file_type(filename, filetype)

  select case (trim(filetype))
  case ("ascii")
    call get_nat_ascii(filename, nat)
  case("cif")
    call get_nat_cif(filename, nat)
  case("in")
    call get_nat_in(filename, nat)
  case('qe')
    call get_nat_quantum_espresso(filename, nat)
  case("gen")
    call get_nat_dftb(filename, nat)
  case("vasp")
    call get_nat_vasp(filename, nat)
  case default
    print*, "filetype, ", filetype
    stop "unknown filetype get_nat_periodic"
  end select
end subroutine get_nat_periodic
