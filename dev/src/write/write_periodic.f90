subroutine write_periodic(filename, nat, rxyz, alat, atomnames, comment)
  !! writes a periodic structure to a file specified in filename variable.
  !! format is detected in ending of filename.
  !! allowed values:
  !! .asci, .cif, .gen, .in, .qe, .vasp and POSCAR (vasp)
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3, nat), intent(in) :: rxyz
  !! atom positions
  real*8, dimension(3, 3), intent(in) :: alat
  !! lattice vectors
  character(len=2), intent(in) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(in) :: comment
  !! content that will be written to comment line.
  integer :: l, vaspind
  character(len=10) :: filetype

  vaspind = index(filename, 'POSCAR')
  l = len_trim(filename)

  if ( l - vaspind == 5 .and. vaspind > 0) then !! vasp file
    call write_vasp(filename, nat, rxyz, alat, atomnames, comment)
    return
  end if

  call get_file_type(filename, filetype)

  select case (trim(filetype))
  case ("ascii")
    call write_ascii(filename, nat, rxyz, alat, atomnames, comment)
  case("cif")
    call write_cif(filename, nat, rxyz, alat, atomnames, comment)
  case("in")
    call write_in(filename, nat, rxyz, alat, atomnames, comment)
  case("gen")
    call write_dftb(filename, nat, rxyz, alat, atomnames, comment)
  case("qe")
    call write_quantum_espresso(filename, nat, rxyz, alat, atomnames)
  case('vasp')
    call write_vasp(filename, nat, rxyz, alat, atomnames, comment)
  case('alm')
    call write_alamode(filename, nat, rxyz, alat, atomnames, comment)
  case default
    print*, "filetype,", filetype
    print*, "filename, ", trim(filename)
    stop "unknown filetype write_periodic"
  end select

end subroutine write_periodic
