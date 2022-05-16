subroutine read_vasp(filename, nat, rxyz, alat, atomnames, comment)
  implicit none
  !! reads a vasp file with the specified filename
  !! units are converted to hartree units before they are returned
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real(8), dimension(3, nat) :: rxyz
  !! atom positions
  real(8), dimension(3, 3) :: alat
  !! lattice vectors
  character(len=2) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(out) :: comment
  !! content of the comment line
  character(len=300) :: all_line
  integer :: io, ios, i
  real(8) :: Bohr_Ang = 0.52917721067
  real(8) :: conver
  integer :: ntypes
  character(len=2), allocatable, dimension(:) :: typenames
  integer, allocatable, dimension(:) :: typecount, typesum
  logical :: is_cartesian
  real(8), allocatable, dimension(:, :) :: xyzred

  ! open file
  open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status='old')
  if (ios /= 0) then
    print*, 'io', io
    print*, 'ios', ios
    print *, "Error opening file "//trim(filename)
    stop "Error opening file"
  end if

  read(io, '(a80)', iostat=ios) comment
  if (ios /= 0 ) stop 'error reading comment in readvasp moleculario'

  read(io, *, iostat = ios) conver
  if (ios /= 0 ) stop ' error reading scaling factor readvasp moleculario'

  do i = 1, 3, 1
    read(io, *, iostat=ios) alat(1,i), alat(2,i), alat(3,i)
    if(ios /=0 ) stop 'error reading lattice readvasp moleculario'
  end do

  read(io, '(a300)', iostat =ios) all_line
  if (ios /=  0 ) stop ' error reading all_line (for atomtypes) in readvasp moleculario'

  call get_num_words(all_line, ntypes, 300)
  allocate(typenames(ntypes), typecount(ntypes), typesum(ntypes))
  call get_words(all_line, ntypes, 300, typenames, 2)

  read(io, *, iostat=ios) typecount
  if(ios /= 0 ) stop 'error reading typcount readvasp moleculario'

  !! check coordinate type (line could contain selective dynamics check that)
  read(io, '(a300)', iostat=ios) all_line
  if (ios /=0 ) stop 'error reading coordinate type readvasp moleculario'
  if ( index(all_line, 'Selective dynamics') > 0 &
      .or. index(all_line, 'selective dynamics') > 0) then
    read(io, '(a300)', iostat=ios) all_line
    if (ios /=0 ) stop 'error reading coordinate type readvasp moleculario'
  end if

  is_cartesian = calc_is_cartesian()

  do i = 1, nat, 1
    read(io, *, iostat = ios) rxyz(1,i), rxyz(2,i), rxyz(3,i)
    if (ios /= 0) stop 'error reading coordinate line in readvasp moleculario'
  end do

  alat = conver*alat / Bohr_Ang

  if ( .not. is_cartesian) then
    allocate(xyzred(3,nat))
    xyzred = rxyz
    call frac2cart(nat, alat, xyzred, rxyz)
    deallocate(xyzred)
  else
    rxyz = conver * rxyz / Bohr_Ang
  end if

  typesum(1) = 1
  do i = 2 ,ntypes
    typesum(i) = typecount(i-1) + typesum(i - 1)
  end do

  do i = 1, ntypes -1, 1
    atomnames(typesum(i): typesum(i+1) -1) = typenames(i)
  end do
  atomnames(typesum(ntypes):nat) = typenames(ntypes)


contains
  logical function calc_is_cartesian()
    all_line = adjustl(all_line)
    if (all_line(1:1) == 'c' .or. all_line(1:1) =='C' &
      .or. all_line(1:1) == 'k' .or. all_line(1:1) =='K') then
      calc_is_cartesian = .TRUE.
      return
    end if
    if ( all_line(1:1) =='d' .or. all_line(1:1) == 'D' ) then
      calc_is_cartesian = .FALSE.
      return
    end if
    stop 'error reading coordinate type line in readvasp moleculario'
  end function calc_is_cartesian

end subroutine read_vasp


subroutine get_nat_vasp(filename, nat)
  !! counts the number of atoms in a vasp1 file.
  implicit none
  character(len=*), intent(in) :: filename
  !! name of the ascii file wich will be read.
  integer, intent(out) :: nat
  ! private variables
  integer :: ios, io, i, ntypes
  real*8 :: place(3)
  character(len=250) :: all_line
  integer, allocatable, dimension(:) :: typecount

  ! open file
  open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status='old')
  if (ios /= 0) then
    print *, "Error opening file"//trim(filename)
    stop "Error opening file"
  end if

  read(io, '(a250)', iostat=ios) all_line
  if (ios /= 0 ) stop 'error reading comment in getnatvasp moleculario'

  read(io, *, iostat = ios) place(1)
  if (ios /= 0 ) stop ' error reading scaling factor getnatvasp moleculario'

  do i = 1, 3, 1
    read(io, *, iostat=ios) place(1), place(2), place(3)
    if(ios /=0 ) stop 'error reading lattice getnatvasp moleculario'
  end do

  read(io, '(a250)', iostat=ios) all_line
  if(ios /=0) stop 'error reading typpenames in getnatvasp moleculario'

  call get_num_words(all_line, ntypes, 250)
  allocate(typecount(ntypes))

  read(io, *, iostat=ios) typecount
  if(ios /=0) stop 'error reading typecount getnatvasp moleculario'

  nat = sum(typecount)
end subroutine get_nat_vasp
