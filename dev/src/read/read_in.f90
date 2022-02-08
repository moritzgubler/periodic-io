subroutine read_in(filename, nat, rxyz, alat, atomnames, comment)
  !! reads a *.in file with the specified filename.
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is read
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3,nat), intent(out) :: rxyz
  !! atom positions (in bohr)
  real*8, dimension(3,3), intent(out) :: alat
  !! lattice vectors in bohr
  character(len=2), dimension(nat), intent(out) :: atomnames
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(out) :: comment
  !! string that was read from comment line (first line in ascii format)
  integer, parameter :: max_line_length = 250
  character(len=max_line_length) :: all_line
  character(len=max_line_length) :: dat_line
  integer :: io, ios, i
  real*8 :: Bohr_Ang = 0.52917721067
  integer :: comment_count

  comment_count = 0

  open(newunit=io, file=filename, iostat=ios, status="old")
  if ( ios /= 0 ) then
    print*, "error opening file: ",trim(filename)
    stop "ioerror"
  end if
  ! read lattice vectors
  i = 1
  lattice_loop: do while ( i <= 3 )
    read(io, "(a250)") all_line
    if ( ios /= 0 ) then
      print*, "error reading lattice vectors in file: ", trim(filename)
      stop "ioerror"
    end if
    all_line = adjustl(all_line)
    if(iscommentorempty(all_line)) then
      if ( index(all_line, "#") == 1 ) then !! comment, read it
        comment_count = comment_count + 1
        if ( comment_count == 1 ) then
          comment = all_line(2:81)
          comment = adjustl(comment)
        end if
      end if
      cycle
    end if

    if ( all_line(1:14) == "lattice_vector" ) then ! lattice vector found
      dat_line = adjustl(all_line(15:max_line_length))
      read(dat_line,*, iostat = ios) alat(1,i), alat(2,i), alat(3,i)
      if ( ios /= 0 ) then
        print*, "error in line containing:"
        print*, trim(all_line)
        stop "error reading lattice vectors"
      end if
      i = i + 1
      cycle
    end if

    !check if an atom tag is found
    if ( all_line(1:4) == "atom" ) then
      print*, "atom found before all lattice vectors were given"
      stop "format error in .in file"
    end if
    print*, "unable to format line: "//trim(all_line)
    stop "format error"
  end do lattice_loop
  ! read atomic positions
  i = 1
  do while ( i <= nat )
    read(io, "(a250)", iostat=ios) all_line
    if ( ios /= 0 ) then
      print*, "error reading atom in file: ", trim(filename)
      stop "ioerror"
    end if
    all_line = adjustl(all_line)
    if(iscommentorempty(all_line)) cycle

    if ( all_line(1:4) == "atom" ) then ! atomic position found
      dat_line = adjustl(all_line(5:max_line_length))
      read(dat_line,*, iostat = ios) rxyz(1,i), rxyz(2,i), rxyz(3,i), atomnames(i)
      if ( ios /= 0 ) then
        print*, "error in line containing:"
        print*, trim(all_line)
        stop "error reading atom"
      end if
      i = i + 1
      cycle
    end if
    if ( all_line(1:9) == "atom_frac" ) then ! reduced atomic position found
      dat_line = adjustl(all_line(5:max_line_length))
      read(dat_line,*, iostat = ios) rxyz(1,i), rxyz(2,i), rxyz(3,i), atomnames(i)
      if ( ios /= 0 ) then
        print*, "error in line containing:"
        print*, trim(all_line)
        stop "error reading atom"
      end if
      rxyz(:,i) = matmul(alat, rxyz(:,i))
      i = i + 1
      cycle
    end if

  end do
  close(io)
  rxyz = rxyz / Bohr_Ang
  alat = alat / Bohr_Ang
  if ( comment_count == 0 ) then
    comment = ""
  end if
contains
  function iscommentorempty(line_in) result(bool)
    implicit none
    character(len=250) :: line_in
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
end subroutine read_in

subroutine get_nat_in(filename, nat)
  !! counts the number of atoms in an *.in file.
  implicit none
  integer, intent(out) :: nat
  !! number of atoms
  character(len=*), intent(in) :: filename
  !! name of the in file wich will be read.
  integer :: io, ios
  character(len=250) :: all_line
  nat = 0
  open(newunit=io, iostat=ios, file=filename, status="old")
  if ( ios /= 0 ) then
    nat = 0
    return
  end if
  ios=0
  do while ( ios == 0 )
    read(io,"(a250)",iostat=ios) all_line
    if ( ios /= 0 ) exit
    all_line = adjustl(all_line)
    if ( all_line(1:4) == "atom" ) then
      nat = nat + 1
    end if
  end do
  if ( ios > 0 ) then !! io error occured
    print*, filename
    stop "io error in getting number of atoms"
  end if
  close(io)
end subroutine get_nat_in
