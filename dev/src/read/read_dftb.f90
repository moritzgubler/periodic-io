
subroutine read_dftb(filename, nat, rxyz, alat, atomnames, comment)
  !! reads the .gen format of the dftb+ package. Cell type S (supercell)
  !! and cell type F (supercel with fractional coordintes) are implemented.
  !! Cell type C (cluster) causes the subroutine to stop.
  !! Documentation of file format in Appendix D:
  !! https://github.com/dftbplus/dftbplus/releases/download/21.1/manual.pdf
  use iso_fortran_env
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real(real64), dimension(3,nat), intent(out) :: rxyz
  !! atom positions (in bohr)
  real(real64), dimension(3,3), intent(out) :: alat
  !! lattice vectors in bohr
  character(len=2), dimension(nat), intent(out) :: atomnames
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(out) :: comment
  !! string that was read from comment line
  integer, parameter :: max_line_length = 250
  character(len=max_line_length) :: all_line
  integer :: io, ios, i, iat
  real(real64) :: Bohr_Ang = 0.52917721067
  character(len=1) :: celltype
  logical :: isfirstcomment
  !! number of different atoms
  integer :: ntypes
  !! Array containing all the different elements. Each element is only listed once.
  !! Only the first ntype entries contain data, the rest is garbage.
  character(len=2), dimension(nat) :: atomtypes
  !! Encodes each element with a number. The element name of iat is accessed the following way:
  !! atomtypes(atomnumbernames(iat))
  integer :: atomnumbernames(nat)
  integer :: natread
  real(real64) :: shift(3)
  real(real64), dimension(:,:), allocatable :: xyzred

  isfirstcomment = .TRUE.
  comment = ""

  open(newunit=io, file=filename, iostat=ios, status="old")
  if ( ios /= 0 ) then
    print*, "error opening file: ",trim(filename)
    stop "ioerror"
  end if
  ! read first line (nat and cell type)
  do
    read(io, "(a250)", iostat=ios) all_line
    if (ios /= 0) stop "error reading line in read_dftb moleculario"
    if (index(adjustl(all_line), "#") == 1) then !! comment line
      if ( isfirstcomment ) then
        isfirstcomment = .FALSE.
        comment = all_line(3:82)
      end if
      cycle
    end if
    read(all_line, *, iostat=ios) natread, celltype
    if (ios /= 0) stop "error parsing line containg dftb nat celltype moleculario"
    exit
  end do

  if ( nat /= natread ) then
    print*, "nat_input, natread :", nat, natread
    stop "nat from argument in read_dftb is not equal to nat read in dftb file. Aborting..."
  end if

  if ( celltype /= "S" .and. celltype /= "F" ) then
    print*, "not accepted cell type: ", celltype
    stop "not recognized cell type in read_dftb moleculario"
  end if

  ! read atomtypes
  do
    read(io, "(a250)", iostat=ios) all_line
    if (ios /= 0) stop "error reading line in read_dftb moleculario"
    if (index(adjustl(all_line), "#") == 1) cycle !! comment line
    all_line = adjustl(all_line)
    call get_num_words(all_line, ntypes, max_line_length)
    call get_words(all_line, ntypes, max_line_length, atomtypes(1:ntypes), 2)
    exit
  end do

  ! read all atoms
  i = 0
  do
    if (nat == i) exit
    read(io, "(a250)", iostat=ios) all_line
    if (ios /= 0) stop "error reading pos in read_dftb moleculario"
    if (index(adjustl(all_line), "#") == 1) cycle !! comment line
    i = i + 1
    read(all_line, *, iostat = ios) iat, atomnumbernames(i), rxyz(1, i), rxyz(2,i), rxyz(3,i)
    if ( ios /= 0 ) then
      stop "error parsing position in read_dftb (moleculario)"
    end if
    if (iat /= i) stop "inconsistent internal state in read_dftb (moleculario)"
    atomnames(i) = atomtypes(atomnumbernames(i))
  end do

  ! read shift of origin from lattice cell
  do
    read(io, "(a250)", iostat=ios) all_line
    if (ios /= 0) stop "error reading shift in read_dftb moleculario"
    if (index(adjustl(all_line), "#") == 1) cycle !! comment line
    read(all_line, *, iostat=ios) shift(1), shift(2), shift(3)
    if (ios /=0 ) stop "error parsing shift in read_dftb (moleculario)"
    exit
  end do

  ! read lattice cell
  i = 0
  do
    if (i == 3) exit
    read(io, "(a250)", iostat=ios) all_line
    if (ios /= 0) stop "error reading lattice in read_dftb moleculario"
    if (index(adjustl(all_line), "#") == 1) cycle !! comment line
    i = i + 1
    read(all_line, *, iostat=ios) alat(1,i), alat(2,i), alat(3,i)
    if (ios /= 0) stop "error parsing lattice in read_dftb moleculario"
  end do
  if (celltype == "F" ) then
    allocate(xyzred(3,nat))
    xyzred = rxyz
    call frac2cart(nat, alat, xyzred, rxyz)
    deallocate(xyzred)
  end if

  !! Positions and lattice vectors are in Angstroem in .gen file format. Convert it
  !! to hartree units.
  rxyz = rxyz / Bohr_Ang
  alat = alat / Bohr_Ang

end subroutine read_dftb

subroutine get_nat_dftb(filename, nat)
  !! counts the number of atoms in an *.gen file.
  implicit none
  integer, intent(out) :: nat
  !! number of atoms
  character(len=*), intent(in) :: filename
    !! name of the gen file wich will be read.
  integer :: io, ios
  character(len=250) :: all_line

  open(newunit=io, iostat=ios, file=filename, status="old")
  if (ios /= 0) stop "error opening dftb file to get number of atoms from (moleculario)"
  ! try and read first line and ignore comments.
  do
    read(io, "(a250)", iostat=ios) all_line
    if (ios /=0) stop "error reading line in get_nat_dftb moleculario"
    if (index(adjustl(all_line), "#") == 1) cycle !! comment line
    read(all_line, *, iostat=ios) nat
    if (ios /= 0) stop "error parsing line containg dftb nat moleculario"
    exit
  end do
  close(io)
end subroutine get_nat_dftb
