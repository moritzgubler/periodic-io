subroutine read_cif(filename,nat , rxyz, alat, atomnames, comment)
  !! This reader will not be able to read everything. It does not understand the full
  !! cif syntax. But it should understand most of them.
  implicit none
  integer :: nat
  !! number of atoms
  character(len=*) :: filename
  !! filename of the file which is created
  real*8, intent(out), dimension(3, nat) :: rxyz
  !! atom positions (in bohr)
  real*8, dimension(3,3), intent(out) :: alat
  !! lattice vectors in bohr
  character(len=2), dimension(nat), intent(out) :: atomnames
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(out) :: comment
  !! string that was read from first comment line
  integer :: io, ios, lc, k
  character(len=250) :: all_line, num_string
  real*8 :: length_a, length_b, length_c
  real*8 :: alpha, beta, gamma
  logical :: aset = .FALSE., bset = .FALSE., cset = .FALSE.
  logical :: aaset = .FALSE., bbset = .FALSE., ccset = .FALSE.
  integer :: atom_site_count, atom_label_pos, x_pos, y_pos, z_pos
  integer :: nat_read
  integer, parameter :: len_words = 30
  character(len=len_words), dimension(:), allocatable :: words
  integer :: nwords, num_words
  real*8 :: xyzred(3, nat)
  real*8 :: Bohr_Ang = 0.52917721067
  integer :: comment_count
  nat_read = 0
  atom_site_count = 0
  x_pos = -1
  y_pos = -1
  z_pos = -1
  comment_count = 0

  open(newunit=io, file=filename, iostat=ios, status="old")
  if ( ios /= 0 ) stop "Error opening file in read_cif"
  lc = 0
  ios = 0
  do
    read(io, fmt="(a250)", iostat=ios) all_line
    if ( ios < 0 ) then
      exit
    end if
    lc = lc + 1
    if ( ios > 0 ) then
      print*, "line", lc, "was formatted wrong and is ignored."
      cycle
    end if
    all_line = adjustl(all_line)
    !check if comment line or ompty
    if ( (index(all_line, "#") == 1) .or. len_trim(all_line) == 0 ) then
      if ( (index(all_line, "#") == 1) ) then ! comment line. read it.
        comment_count = comment_count + 1
        if ( comment_count == 1 ) then
          comment = all_line(2:81)
          comment = adjustl(comment)
        end if
      end if
      cycle
    end if
    k = index(all_line, "_cell_length_a")
    if ( k == 1) then !! line contains cell length a
      num_string = all_line((k+15):len_trim(all_line))
      call remove_uncert_from_numstring(num_string)
      read(num_string, *, iostat=ios) length_a
      if ( ios /= 0 ) then
        print*, 'error reading cell length a'
        print*, trim(all_line)
        stop 'error reading _cell_length_a in read_cif periodicIO'
      end if
      aset = .TRUE.
      cycle
    end if
    k = index(all_line, "_cell_length_b")
    if ( k == 1) then !! line contains cell length b
      num_string = all_line((k+15):len_trim(all_line))
      call remove_uncert_from_numstring(num_string)
      read(num_string, *, iostat=ios) length_b
      if(ios /=0) then
        print*, 'error reading cell length b'
        print*, trim(all_line)
        stop 'error reading _cell_length_b in read_cif periodicIO'
      end if
      bset = .TRUE.
      cycle
    end if
    k = index(all_line, "_cell_length_c")
    if ( k == 1) then !! line contains cell length c
      num_string = all_line((k+15):len_trim(all_line))
      call remove_uncert_from_numstring(num_string)
      read(num_string, *, iostat=ios) length_c
      if(ios /=0) then
        print*, 'error reading cell length c'
        print*, trim(all_line)
        stop 'error reading _cell_length_c in read_cif periodicIO'
      end if
      cset = .TRUE.
      cycle
    end if
    k = index(all_line, "_cell_angle_alpha")
    if ( k == 1) then !! line contains angle alpha
      num_string = all_line((k+18):len_trim(all_line))
      call remove_uncert_from_numstring(num_string)
      read(num_string, *, iostat=ios) alpha
      if ( ios /= 0 ) then
        stop 'error reading cell_angle_alpha in readcif'
      end if
      aaset = .TRUE.
      cycle
    end if
    k = index(all_line, "_cell_angle_beta")
    if ( k == 1) then !! line contains angle beta
      num_string = all_line((k+17):len_trim(all_line))
      call remove_uncert_from_numstring(num_string)
      read(num_string, *, iostat=ios) beta
      if ( ios /= 0 ) then
        stop 'error reading cell_angle_beta in readcif'
      end if
      bbset = .TRUE.
      cycle
    end if
    k = index(all_line, "_cell_angle_gamma")
    if ( k == 1) then !! line contains angle gamma
      num_string = all_line((k+18):len_trim(all_line))
      call remove_uncert_from_numstring(num_string)
      read(num_string, *, iostat=ios) gamma
      if ( ios /= 0 ) then
        stop 'error reading cell_angle_gamma in readcif'
      end if
      ccset = .TRUE.
      cycle
    end if

    if ( index(all_line, "_atom_site_") == 1 ) then
      atom_site_count = atom_site_count + 1
      if(index(all_line, "_atom_site_type_symbol") == 1) atom_label_pos = atom_site_count
      if(index(all_line, "_atom_site_fract_x") == 1) x_pos = atom_site_count
      if(index(all_line, "_atom_site_fract_y") == 1) y_pos = atom_site_count
      if(index(all_line, "_atom_site_fract_z") == 1) z_pos = atom_site_count
      cycle
    end if

    if ( (x_pos > 0) .and. (y_pos > 0) .and. (z_pos > 0) ) then ! atom line expected after _atom_* lines have been given.
      nat_read = nat_read + 1
      if ( nat_read == 1 ) then
        nwords = atom_site_count
        allocate(words(nwords))
      end if
      if ( nwords /= atom_site_count ) then !! unexpected _atom_site_ line
        stop "nwords /= atom_site_count"
      end if
      call get_num_words(all_line, num_words, 250)
      if ( num_words /= nwords ) then
        print*, trim(all_line)
        stop "numwords /= nwords"
      end if
      call get_words(all_line, nwords, 250, words, len_words)
      if ( len_trim(words(atom_label_pos)) > 2 ) then
        print*, trim(words(atom_label_pos))
        stop "expected atom symbol but got line above"
      end if
      atomnames(nat_read) = trim(words(atom_label_pos))
      read(words(x_pos), *, iostat=ios) xyzred(1, nat_read)
      if ( ios/= 0 ) stop "error reading atom position"
      read(words(y_pos), *, iostat=ios) xyzred(2, nat_read)
      if ( ios/= 0 ) stop "error reading atom position"
      read(words(z_pos), *, iostat=ios) xyzred(3, nat_read)
      if ( ios/= 0 ) stop "error reading atom position"
      if ( nat_read == nat ) then
        exit
      end if
    end if
  end do
  if ( nat /= nat_read ) then
    print*, nat, nat_read
    stop "expected and actual number of atoms do not match in read_cif"
  end if
  length_a = length_a / Bohr_Ang
  length_b = length_b / Bohr_Ang
  length_c = length_c / Bohr_Ang
  call calcAlat
  call frac2cart(nat, alat, xyzred, rxyz)
  if ( comment_count == 0 ) then
    comment = ""
  end if
contains
  subroutine calcAlat
    implicit none
    real*8 :: conv = 0.01745329251994329576 ! converts degree to radians
    alat = 0
    alat(1,1) = length_a
    alat(1,2) = length_b*cos(gamma*conv)
    alat(2,2) = sqrt(length_b*length_b - alat(1,2)*alat(1,2))
    alat(1,3) = length_c*cos(beta*conv)
    alat(2,3) = (length_b * length_c * cos(alpha*conv) &
                  - alat(1,2) * alat(1,3)) / alat(2,2)
    if(length_c**2 - alat(1,3)**2 - alat(2,3)**2 < 0) then
      stop "square root of negativa number has to be calcuted in read cif."
    end if
    alat(3,3) = sqrt(length_c**2 - alat(1,3)**2 - alat(2,3)**2)
  end subroutine calcAlat

  subroutine remove_uncert_from_numstring(numstr)
    !! in some .cif file numbers are written with brackets in the end: 3.141(5)
    !! This subroutine removes these brackets and the number in between them.
    implicit none
    character(len=*), intent(inout) :: numstr
    !! string that may contain brackets that will be removed if present.
    integer :: len_word, numwords, brack_index

    len_word = len(numstr)
    call get_num_words(numstr, numwords, len_word)
    if ( numwords /= 1 ) then
      print*, trim(numstr)
      print*, numwords
      stop 'multiple words supplied in remove_encert_from_numstr'
    end if

    brack_index = index(numstr, '(')
    if (brack_index <= 0) return
    numstr(brack_index:len_word) = ''
  end subroutine remove_uncert_from_numstring

end subroutine read_cif

subroutine get_nat_cif(filename, nat)
  !! counts the number of atoms in an *.cif file.
  implicit none
  integer, intent(out) :: nat
  !! number of atoms
  character(len=*) :: filename
  !! name of the cif file which will be read
  integer :: io, ios
  character(len=250) :: all_line
  logical :: look_for_atoms, isEle
  integer :: atom_site_pos, nwords
  integer, parameter :: len_words = 30
  character(len=len_words), dimension(:), allocatable :: words
  atom_site_pos = 1
  look_for_atoms = .FALSE.
  open(newunit=io, file=filename, iostat=ios, status="old")
  if ( ios /= 0 ) then
    print*, filename
    stop "error opening cif file to read number of atoms."
  end if
  nat = 0
  do
    read(io, "(a)", iostat=ios) all_line
    if(ios < 0) exit
    if ( ios > 0 ) then
      print*, "reading line in cif file", filename
      stop "error reading line"
    end if
    all_line = adjustl(all_line)
    !call downcase(all_line)
    if (iscommentorempty(all_line)) cycle
    if(index(all_line, "_atom_site_type_symbol") == 1) then
      look_for_atoms = .TRUE.
      allocate(words(atom_site_pos))
      cycle
    end if

    if ( .not. look_for_atoms ) then !! look for _atom_site labels
      if ( index(all_line, "_atom_site_") == 1) then
        atom_site_pos = atom_site_pos + 1
        cycle
      end if
    else ! look for atoms.
      call get_num_words(all_line, nwords, 250)
      if ( nwords >= atom_site_pos ) then !! posible atom
        call get_words(all_line, atom_site_pos, 250, words, len_words)
        call isElement(words(atom_site_pos)(1:2), isEle)
        if(isEle) then !! atom found
          nat = nat + 1
        end if
      else !! no atom
        cycle
      end if
    end if
  end do

  close(io)
contains
  function iscommentorempty(line_in) result(bool)
    implicit none
    character(len=*), intent(in) :: line_in
    logical :: bool
    bool = .FALSE.
    if ( len_trim(line_in) == 0 ) then !line is empty
      bool = .TRUE.
    end if
    if ( adjustl(line_in(1:1)) == "#" ) then
      bool = .TRUE. ! line is comment
    end if
    return
  end function iscommentorempty

end subroutine get_nat_cif
