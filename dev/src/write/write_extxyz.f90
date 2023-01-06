subroutine write_extxyz(filename, nat, rxyz, alat, atomnames, comment)
  !! reads the position, lattice vecors and atomnames from a file. File format is
  !! detected based on the file ending. Supported endings:
  !! .ascii, .cif, .in, .gen, .qe, .vasp and POSCAR (vasp)
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3, nat), intent(in) :: rxyz
  !! atom positions (bohr)
  real*8, dimension(3, 3), intent(in) :: alat
  !! lattice vectors (bohr)
  character(len=2), intent(in) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(in) :: comment
  integer :: io, ios

  open(newunit=io, file=filename, iostat=ios)
  if ( ios /= 0 ) then
    stop 'error opening file write_extxyz'
  end if

  call write_partial_extxyz(io, nat, rxyz, alat, atomnames)

  close(io)

end subroutine write_extxyz

subroutine write_partial_extxyz(io, nat, rxyz, alat, atomnames)
  !! writes xyz file. assumes units of rxyz are in bohr and writes them to file
  !! in a.u.
  implicit none
  integer, intent(in) :: io
  !! io unit
  integer, intent(in) :: nat
  !! number of atoms
  real(8), intent(in), dimension(3,nat) :: rxyz
  !! atomic positions in a.u
  real(8), intent(in), dimension(3, 3) :: alat
  character(len=2), intent(in), dimension(nat) :: atomnames
  !! chemical symbol of the atoms
  integer :: ios, i
  real(8) :: Bohr_Ang = 0.52917721067
  real(8) :: xyz_convert(3,nat), lat_A(3,3)
  character(len=500) :: ext_comment
  character(len=100) :: propertystring

  xyz_convert = Bohr_Ang * rxyz
  lat_A = Bohr_Ang * alat

  ! create comment
  write(ext_comment, '(a9, 8(g0.7, 1x), g0.7, 1a)') 'Lattice="', lat_A, '"'

  propertystring = 'Properties="species:S:1:pos:R:3"'

  ext_comment = trim(ext_comment) // ' ' // trim(propertystring)

  write(io,*) nat
  write(io,*) trim(ext_comment)
  do i = 1, nat, 1
    write(io,*, iostat = ios) atomnames(i), xyz_convert(1,i), xyz_convert(2,i), xyz_convert(3,i)
    if ( ios/=0 ) then
      print*, "error writing line: ", i
      stop
    end if
  end do

end subroutine write_partial_extxyz
