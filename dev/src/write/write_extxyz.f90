subroutine write_xyz(io, nat, rxyz, alat, atomnames, comment, energy, forces, stress)
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
  character(len=80), intent(in) :: comment
  !! content that will be written to the comment line
  real(8), optional :: energy
  real(8), optional, dimension(3, nat) :: forces
  real(8), optional, dimension(3, 3) :: stress
  integer :: ios, i
  real(8) :: Bohr_Ang = 0.52917721067
  real(8) :: xyz_convert(3,nat)
  character(len=500) :: ext_comment
  character(len=100) :: propertystring

  xyz_convert = Bohr_Ang * rxyz

  ! create comment
  write(ext_comment, '(a9, 8(g0.7, 1x), g0.7, 1a)') 'Lattice="', alat(1,1), alat(2,1), alat(3,1), alat(1,2), alat(2,2), alat(3,2)&
    , alat(1,3), alat(2,3), alat(3,3), '"'

  if ( present(energy) ) then
    write(ext_comment, '(a, a, g0.10)') trim(ext_comment), ' energy=', energy
  end if
  if ( present(stress) ) then
    write(ext_comment, '(a7, 8(g0.7, 1x), g0.7, 1a)') 'stress="', stress, '"'
  end if
  if (present(forces)) then
    propertystring = 'Properties="species:S:1:pos:R:3:forces:R:3"'
  else
    propertystring = 'Properties="species:S:1:pos:R:3:forces:R:3"'
  end if

  if(ios/=0) stop "error opening output file"
  write(io,*) nat
  write(io,*) trim(comment)
  do i = 1, nat, 1
    write(io,*, iostat = ios) atomnames(i), xyz_convert(1,i), xyz_convert(2,i), xyz_convert(3,i)
    if ( ios/=0 ) then
      print*, "error writing line: ", i
      stop
    end if
  end do

end subroutine write_xyz
