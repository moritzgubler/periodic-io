subroutine write_vasp(filename, nat, rxyz, alat, atomnames, comment)
  !! writes periodic geometry  in  vasp format
  !! Multiple instances of an atom type must be contigouos. This subroutine
  !! takes care of this. ( C C C C H H Si Si ...)
  !! https://www.vasp.at/wiki/index.php/POSCAR
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real(8), dimension(3, nat), intent(in) :: rxyz
  !! atom positions in bohr
  real(8), dimension(3, 3), intent(in) :: alat
  !! lattice vectors (bohr)
  character(len=2), intent(in) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(in) :: comment
  !! content that will be written to the comment line

  ! private variables:
  integer :: io, ios, i
  real(8) :: Bohr_Ang = 0.52917721067
  real(8) :: p_out(3,nat), a_out(3,3)
  !! number of different atoms
  integer :: ntypes
  !! Array containing all the different elements. Each element is only listed once.
  !! Only the first ntype entries contain data, the rest is garbage.
  character(len=2), dimension(nat) :: atomtypes
  !! Encodes each element with a number. The element name of iat is accessed the following way:
  !! atomtypes(atomnumbernames(iat))
  integer :: atomnumbernames(nat), itype
  integer, allocatable, dimension(:) :: typepos, typecount

  if ( nat <= 0 ) then
    stop 'number of atoms is non positive in write_vasp'
  end if

  ! calculate ntypes and atomtypes
  ntypes = 1
  atomtypes(1) = atomnames(1)
  do i = 2, nat, 1
    do itype = 1, ntypes, 1
      if ( atomnames(i) == atomtypes(itype) ) then
        goto 123
      end if
    end do
    ntypes = ntypes + 1
    atomtypes(ntypes) = atomnames(i)
    123 continue
  end do

  ! calculate atomnumbernames
  allocate(typepos(ntypes), typecount(ntypes))
  typecount = 0
  do i = 1, nat, 1
    do itype = 1, ntypes, 1
      if ( atomnames(i) == atomtypes(itype) ) then
        atomnumbernames(i) = itype
        typecount(itype) = typecount(itype) + 1
      end if
    end do
  end do
  typepos(1) = 1
  do i = 2 ,ntypes
    typepos(i) = typecount(i-1) + typepos(i - 1)
  end do

  !! sort array
  do i = 1, nat, 1
    p_out(:, typepos(atomnumbernames(i)) ) = rxyz(:,i)
    typepos(atomnumbernames(i)) = typepos(atomnumbernames(i)) + 1
  end do

  deallocate(typepos)


  ! open file
  open (newunit=io, file=trim(adjustl(filename)), iostat=ios)
  if (ios /= 0) then
    print *, "Error opening file"//trim(filename)
    stop "Error opening file"
  end if

  !! write comment to file
  write(io, '(a)', iostat=ios) trim(comment)
  if(ios /= 0) then
    print*, 'error writing to file in write_vasp moleculario'
    print*, trim(filename)
    stop 'error writing to file'
  end if

  ! convert to angstroem
  p_out = p_out * Bohr_Ang
  a_out = alat * Bohr_Ang

  write(io, '(a6)', iostat=ios) '  1.00'
  if(ios/=0) stop 'error in write_vasp'

  !! write lattice
  do i = 1, 3, 1
    write(io, *, iostat=ios) a_out(1, i), a_out(2, i), a_out(3, i)
    if(ios /= 0) stop 'error writing lattice in write_vasp'
  end do

  ! write atomtypes
  write(io, "( *(a, 1x) )", iostat=ios) atomtypes(1:ntypes)
  if (ios /= 0 ) stop ' error writing atomtypes in moleculario'

  ! write atomcounts
  write(io, *, iostat=ios) typecount
  if (ios /= 0) stop 'error writing typecount in moleculario'

  ! write cartesian flag
  write(io, '(a)', iostat= ios) 'Cartesian'
  if (ios /= 0) stop 'error writing cartesian flag in moleculario'

  ! write positions
  do i = 1, nat, 1
    write(io, *, iostat  = ios) p_out(1,i), p_out(2,i), p_out(3,i)
    if ( ios /= 0) stop 'error writing position vasp moleculario'
  end do

  close(io, iostat=ios)
  if (ios /= 0 ) stop 'error closing vasp file moleculario'

end subroutine write_vasp
