subroutine read_extxyz(filename, nat, rxyz, alat, atomnames, comment)
    !! reads xyz file. assumes units of file are in angstrom and converts them to a.u.
    implicit none
    integer, parameter :: lineLength = 500
    character(len=*), intent(in) :: filename
    !! name of the file that will be read
    integer, intent(inout) :: nat
    !! number of atoms
    real(8), intent(out), dimension(3,nat) :: rxyz
    !! positions of the atoms (bohr)
    real(8), dimension(3, 3) :: alat
    !! lattice vectors in bohr
    character(len=2), intent(out), dimension(nat) :: atomnames
    !! chemical name of the atoms
    character(len=80), intent(out) :: comment
    !! contents of the first comment line

    real(8), parameter :: Bohr_Ang = 0.52917721067
    integer :: io, ios, i
    character(len=lineLength) :: all_line


    open(newunit=io, file=filename, iostat=ios, status="old")
    if(ios/=0) stop "error opening input file"
    read(io,*) nat
    read(io,"(a80)") comment
    do i = 1, nat, 1
      read(io,*, iostat = ios) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
      if ( ios/=0 ) then
        print*, "error reading line: ", i
        stop
      end if
    end do
    close(io)
    rxyz = rxyz / Bohr_Ang

  contains
    function iscommentorempty(line_in) result(bool)
      implicit none
      character(len=lineLength) :: line_in
      logical :: bool
      bool = .FALSE.
      if ( len_trim(line_in) == 0 ) then !line is empty
        bool = .TRUE.
      end if
      if ( line_in(1:1) == "#" ) then
        bool = .TRUE. ! line is comment
      end if
      return
    end function iscommentorempty

  end subroutine read_extxyz
  
  subroutine get_nat_extxyz(filename, nat)
    !! counts the number of atoms in a .xyz file
    implicit none
    character(len=*), intent(in) :: filename
    !! filename that will be parsed
    integer, intent(out) :: nat
    !! number of atoms
    integer :: io, ios
    open(newunit=io,file=filename,iostat=ios, status="old")
    if(ios/=0) stop "error opening input file"
    read(io,*) nat
    close(io)
  end subroutine get_nat_extxyz
  