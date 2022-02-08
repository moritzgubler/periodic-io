subroutine write_xyz(filename, nat, rxyz, atomnames, comment)
  !! writes xyz file. assumes units of rxyz are in bohr and writes them to file
  !! in a.u.
  implicit none
  character(len=*), intent(in) :: filename
  !! name of the output file
  integer, intent(in) :: nat
  !! number of atoms
  real*8, intent(in), dimension(3,nat) :: rxyz
  !! atomic positions in a.u
  character(len=2), intent(in), dimension(nat) :: atomnames
  !! chemical symbol of the atoms
  character(len=80), intent(in) :: comment
  !! content that will be written to the comment line
  integer :: io, ios, i
  real*8 :: Bohr_Ang = 0.52917721067
  real*8 :: xyz_convert(3,nat)

  xyz_convert = Bohr_Ang * rxyz
  open(newunit=io,file=filename,iostat=ios)
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
  close(io)
end subroutine write_xyz
