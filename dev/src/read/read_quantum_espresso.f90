subroutine read_quantum_espresso(filename, nat, rxyz, alat, atomnames)
  !! reads the periodic structure of a quantum espresso file. units are returned in bohr
  !! only ibrav = 0 is currently supported.
  use, intrinsic :: iso_fortran_env
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

  real(real64) :: Bohr_Ang = 0.52917721067
  integer :: io, ios, i, iend
  character(len=200) :: aline
  integer :: ibrav
  logical :: ibravset
  real(real64) :: convert
  logical :: latread, posread

  ibravset = .FALSE.
  latread = .FALSE.
  posread = .FALSE.

  open(newunit=io, file=filename, iostat=ios, status="old")
  if ( ios /= 0 ) then
    print*, "error opening file: ",trim(filename)
    stop "ioerror in read_quantum_espresso (molecularIO)"
  end if

  floop: do
      read(io, "(a200)", iostat = ios) aline
      if (ios /= 0) exit floop
      if ( index(aline, "ibrav") > 0 .or. index(aline, "IBRAV") > 0) then
        i = index(aline, "=") + 1
        iend = index(aline, ",")
        if (iend <=0 ) iend = 200
        ibravset = .TRUE.
        read(aline(i:iend), *, iostat = ios) ibrav
        if ( ibrav /= 0 ) then
          stop "only ibrav = 0 supported at the moment"
        end if
      end if

      if ( index(aline, "CELL_PARAMETERS") > 0 ) then ! cell parameters found
        if(index(aline, "bohr") > 0) then
           convert = 1.0
        else if( index(aline, "angstrom") > 0) then
          convert = 1.0 / Bohr_Ang
        else
          stop "unknown units"
        end if
        do i = 1, 3, 1
          read(io, *, iostat = ios) alat(i, 1), alat(i, 2), alat(i, 3)
          if ( ios /= 0 ) then
            stop "error reading lattice vectors"
          end if
        end do
        alat = alat * convert
        latread = .TRUE.
      end if

      ! atoms found
      if ( index(aline, "ATOMIC_POSITIONS") > 0 .or. index(aline, "atomic_positions") > 0 ) then
        ! get units:
        if ( index(aline, "bohr") > 0 ) then
          convert = 1.0
        else if( index(aline, "angstrom") > 0) then
          convert = 1.0 / Bohr_Ang
        else
          print*, "unkown unit on this line:", trim(aline)
          stop "unknown unkwonsn unit"
        end if
        do i = 1, nat, 1
          read(io, *, iostat=ios) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
          if ( ios /= 0 ) then
            stop "error reading positions in qe read moleculario"
          end if
        end do
        rxyz = rxyz * convert
        posread = .TRUE.
      end if
  end do floop
  if ( .not. ibravset ) then
    print*, 'ibrav not set in quantum espresso file. Proceed with care.'
  end if
  if ( .not. latread ) stop "lattice vectors not set"
  if ( .not. posread ) stop "positions not red"
  close(io)
end subroutine read_quantum_espresso

subroutine get_nat_quantum_espresso(filename, nat)
  !! counts the number of atoms in a quantum espresso file.
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(out) :: nat
  !! number of atoms
  character(len=*) :: filename
  !! filename that will be read.

  integer :: io, ios, i, iend
  character(len=200) :: aline

  open(newunit=io, file=filename, iostat=ios, status="old")
  if ( ios /= 0 ) then
    print*, "error opening file: ",trim(filename)
    stop "ioerror in read_quantum_espresso (molecularIO)"
  end if
  nat = -1

  floop: do
    read(io, "(a200)", iostat = ios) aline
    if (ios /= 0) exit floop
    if ( index(aline, "nat") > 0 .or. index(aline, "NAT") > 0) then
      i = index(aline, "=") + 1
      if (i < 1) stop "no equal sign in nat line qe in file (molecularIO)"
      iend = index(aline, ",")
      if (iend <= 0) iend = 200
        read(aline(i:iend), *, iostat = ios) nat
        if (ios /= 0) stop "error getting nat from quantum espresso file"
      end if
  end do floop
  close(io)
end subroutine get_nat_quantum_espresso
