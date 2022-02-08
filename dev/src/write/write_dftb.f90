subroutine write_dftb(filename, nat, rxyz, alat, atomnames, comment)
  !! writes periodic structure to a file in the dftb+ (.gen) format
  use iso_fortran_env
  implicit none
  !! filename of the file which is created
  character(len=*), intent(in) :: filename
  !! number of atoms
  integer, intent(in) :: nat
  !! atom positions in bohr
  real(real64), dimension(3, nat), intent(in) :: rxyz
  !! lattice vectors
  real(real64), dimension(3, 3), intent(in) :: alat
  !! String containing the chemical symbol of each atom.
  character(len=2), intent(in) :: atomnames(nat)
  !! comment string that will be written to file
  character(len=80), intent(in) :: comment
  !! number of different atoms
  integer :: ntypes
  !! Array containing all the different elements. Each element is only listed once.
  !! Only the first ntype entries contain data, the rest is garbage.
  character(len=2), dimension(nat) :: atomtypes
  !! Encodes each element with a number. The element name of iat is accessed the following way:
  !! atomtypes(atomnumbernames(iat))
  integer :: atomnumbernames(nat)

  integer :: u, ityp, iat, ios, i, j

  real(real64), parameter :: Bohr_Ang = 0.52917721067
  real(real64) :: pos_ang(3,nat), lat_ang(3,3)

  pos_ang = rxyz * Bohr_Ang
  lat_ang = alat * Bohr_Ang

  call back2cell(nat, pos_ang, lat_ang)

  ! calculate ntypes and atomtypes
  ntypes = 1
  atomtypes(1) = atomnames(1)
  do i = 2, nat, 1
    do j = 1, ntypes, 1
      if ( atomnames(i) == atomtypes(j) ) then
        goto 123
      end if
    end do
    ntypes = ntypes + 1
    atomtypes(ntypes) = atomnames(i)
    123 continue
  end do

  ! calculate atomnumbernames
  do i = 1, nat, 1
    do j = 1, ntypes, 1
      if ( atomnames(i) == atomtypes(j) ) then
        atomnumbernames(i) = j
      end if
    end do
  end do

  open(newunit=u, file=filename, iostat=ios)
  if ( ios /= 0 ) stop "error opening dftb output file in moleculario"

  if ( len_trim(comment) > 0 ) then
      write(u, "(a1, 1x, a)", iostat=ios) "#", trim(comment)
      if (ios /= 0) stop "error writing comment in write dftb in moleculario"
  end if

  ! first non comment line contains number of atoms and format type. Here 'S' is
  ! used for a periodic supercell
  write(u,'(i5.5,a)', iostat=ios) nat, " S"
  if (ios /= 0) stop "error writing to dftb  output file in moleculario"

  write(u,*, iostat=ios) (trim(adjustl(atomtypes(ityp)))//" ", ityp=1,ntypes)
  if (ios /= 0) stop "error writing atomtypes to dftb  output file in moleculario"

  do iat = 1, nat
      write(u,'(i5,1x,i5,3(1x,es25.15))', iostat=ios) iat, atomnumbernames(iat), pos_ang(:,iat)
      if (ios /= 0) stop "error writing positions to dftb  output file in moleculario"
  end do

  ! write origin of coordinate system
  write(u,'(3(1x,es25.15))', iostat=ios) 0.d0,0.d0,0.d0
  if (ios /= 0) stop "error writing zeros to dftb output file in moleculario"

  write(u,'(3(1x,es25.15))', iostat=ios) lat_ang(:, 1)
  if (ios /= 0) stop "error writing lattice vector 1 to dftb output file in moleculario"

  write(u,'(3(1x,es25.15))', iostat=ios) lat_ang(:, 2)
  if (ios /= 0) stop "error writing lattice vector 2 to dftb output file in moleculario"

  write(u,'(3(1x,es25.15))', iostat=ios) lat_ang(:, 3)
  if (ios /= 0) stop "error writing lattice vector 3 to dftb output file in moleculario"

  close(u)
end subroutine write_dftb
