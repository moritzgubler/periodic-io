subroutine write_alamode(filename, nat, rxyz, alat, atomnames, comment)
  !! writes a periodic structure to a .alm file.
  !!IMPORTANT:  puts atoms back into cell.
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
  integer :: io, ios, i, j
  real*8 :: xyzred(3,nat)
  character(len=2), dimension(nat) :: an_copy
  integer :: ntypes
  character(len=2), allocatable, dimension(:) :: atom_types
  character(len=30) :: fstring

  ntypes = 0
  an_copy = atomnames
  do i = 1, nat, 1
    if (len_trim(an_copy(i)) /= 0) then ! new atom species found. delete all the following from list
      ntypes = ntypes + 1
      do j = i + 1, nat, 1
        if (an_copy(j) == an_copy(i)) then
          an_copy(j) = ""
        end if
      end do
    end if
  end do
  if (allocated(atom_types)) then
    deallocate (atom_types)
  end if
  allocate (atom_types(ntypes))
  j = 0
  do i = 1, nat, 1
    if (an_copy(i) /= "") then
      j = j + 1
      atom_types(j) = an_copy(i)
    end if
  end do

  call cart2frac(nat, alat, rxyz, xyzred)
  xyzred = modulo(xyzred, 1.d0)
  open(newunit=io,file=filename,iostat=ios)
  if(ios /= 0) stop "error opening output file"

  ! write general tag
  write(io, '(a)') '&general'
  write(io, '(a, i5)', iostat=ios) '  NAT = ', nat
  if (ios /= 0) stop 'error writing nat '
  write(io, '(a, i5)', iostat=ios) '  NKD = ', ntypes
  if(ios /= 0) stop 'error writing ntypes alamode'
  write(fstring, '(a, i5.5, a)', iostat=ios) '(a,', ntypes,'a )'
  if( ios /= 0) stop 'error in internal write ntypes alamode'
  write(io, fstring, iostat=ios) '  KD = ', ( atom_types(i)//' ', i=1,ntypes )
  if (ios /=0) stop 'error writing atomtypes in alamode'
  write(io, '(a)') '/'

  ! write atomic positions
  write(io, '(a)') '&cell'
  write(io, '(a)') '  1.0'
  do i = 1, 3, 1
    write(io, *, iostat=ios) alat(i, 1), alat(i, 2), alat(i, 3)
    if (ios /=0) stop 'error writing lattice vectors in alamode'
  end do
  write(io, '(a)') '/'

  ! write positions
  write(io, '(a)') '&position'
  do i = 1, nat, 1
    do j = 1, ntypes, 1
      if ( atomnames(i) == atom_types(j) ) exit
    end do
    write(io, '(2x, i5, 3(1x, f20.15))', iostat=ios) j, xyzred(1, i), xyzred(2, i), xyzred(3, i)
    if (ios /= 0) stop 'error writing alamode positions'
  end do

  write(io, '(a)') '/'

  close(io)
end subroutine write_alamode
