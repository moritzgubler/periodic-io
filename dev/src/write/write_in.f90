subroutine write_in(filename, nat, rxyz, alat, atomnames, comment)
  !! writes a periodic structure to a . in file
  implicit none
  character(len=*) :: filename
  !! output file
  integer, intent(in) :: nat
  !! number of atoms
  real*8, intent(in), dimension(3,nat) :: rxyz
  !! atomic positions in bohr
  character(len=2), intent(in), dimension(nat) :: atomnames
  !! atomnames
  real*8, dimension(3,3), intent(in) :: alat
  !! lattice vectors in bohr
  character(len=80), intent(in) :: comment
  !! content that will be written to comment line
  integer :: io, ios, i
  real*8 :: Bohr_Ang = 0.52917721067
  real*8 :: alat_convert(3,3), xyz_convert(3,nat)

  alat_convert = alat * Bohr_Ang
  xyz_convert = rxyz * Bohr_Ang

  open(newunit=io,file=filename,iostat=ios)
  if(ios /= 0) stop "error opening output file"
  if ( len_trim(comment) > 0 ) then
    write(io, *) "# "//trim(comment)
  end if
  write(io,*) "lattice_vector", alat_convert(1,1), alat_convert(2,1), alat_convert(3,1)
  write(io,*) "lattice_vector", alat_convert(1,2), alat_convert(2,2), alat_convert(3,2)
  write(io,*) "lattice_vector", alat_convert(1,3), alat_convert(2,3), alat_convert(3,3)
  do i = 1, nat, 1
    write(io,*) "atom", xyz_convert(1,i), xyz_convert(2,i), xyz_convert(3,i), atomnames(i)
  end do
  close(io)
end subroutine write_in
