subroutine get_num_words(string, nwords, len_string)
  !! gets the number of words (characters separated by one or multiple spaces) in a string
  implicit none
  integer, intent(out) :: nwords
  !! number of words
  integer, intent(in) :: len_string
  !! length of input string
  character(len=len_string), intent(in) :: string
  !! input string were number of words should be counted.
  integer :: ichar
  !! iteration variable of characters of string
  logical :: in_word

  if ( len_trim(string) == 0 ) then
    nwords = 0
    return
  end if
  ichar = 1
  ! advance until first word starts
  do while ( string(ichar:ichar) == " " .and. ichar <= len_string)
    ichar = ichar + 1
  end do
  nwords = 0
  ichar = 1
  in_word = .FALSE.
  do while ( ichar <= len_string ) !! find the atom site column
    if ( (string(ichar:ichar) == " ") .and. in_word ) then
      in_word = .FALSE.
      nwords = nwords + 1
    else
      if(string(ichar:ichar) /= " ") in_word = .TRUE.
    end if
    ichar = ichar + 1
  end do
end subroutine get_num_words


subroutine get_words(string, nwords, len_string, words, len_words)
  !! splits string seperated by spaces into words.
  implicit none
  !! number of words that are in string (use getnumword to calculate this)
  integer, intent(in) :: nwords
  !! length of the string
  integer, intent(in) :: len_string
  !! string that shoul be split
  character(len=len_string), intent(in) :: string
  !! max length of each word
  integer, intent(in) :: len_words
  !! string array containing all the words
  character(len=len_words), dimension(nwords) :: words
  integer :: ichar, iword, wordend
  ichar = 1
  words = ""
  do iword = 1, nwords, 1
    do while ( string(ichar:ichar) == " " .and. ichar <= len_string)
      ichar = ichar + 1
    end do
    wordend = index(string(ichar:len_string), " ")
    if ( wordend <= 0 ) then
      words(iword) = string(ichar:len_string)
      return
    end if
    words(iword) = string(ichar:(ichar + wordend-1))
    ichar = ichar + wordend - 1
  end do
end subroutine get_words

subroutine get_file_type(filename, filetype)
  implicit none
  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: filetype
  integer :: sep_pos, len_file
  len_file = len_trim(filename)
  do sep_pos = len_file, 1, -1
    if ( filename(sep_pos:sep_pos) == "." ) then
      exit
    end if
  end do
  if ( len_file == 0 ) then
    stop "empy filename in get_file_type"
  end if
  if ( sep_pos <= 0 ) then
    filetype = ""
    return
  end if
  if ( sep_pos == len_file ) then
    filetype = ""
    return
  end if
  filetype = filename((sep_pos + 1):len_file)
end subroutine get_file_type

subroutine get_basename(filename, bname)
  implicit none
  !! filename where basename should be extracted from
  character(len=*), intent(in) :: filename
  !! basename of filename (stripped of ending and prepending path)
  character(len=*), intent(out) :: bname
  integer :: lenfile, i
  !! position of point in filename
  integer :: point_pos
  !! position of last "/".
  integer :: slash_pos
  lenfile = len_trim(filename)
  if ( lenfile == 0) then
    stop "empty string in get basename"
  end if
  do point_pos = lenfile, 1, -1
    if ( filename(point_pos:point_pos) == "." ) then
      exit
    end if
  end do
  if ( point_pos <= 0 ) then
    stop "filename does not contain a dot in moleculario library."
  end if
  slash_pos = 0
  do i = 1, point_pos, 1
    if ( filename(i:i) == "/" ) then
      slash_pos = i
    end if
  end do

  bname = filename(slash_pos+1:point_pos - 1)
end subroutine get_basename

!program test
!  implicit none
!  character(len=100) :: string
!  integer :: len_string=100, nwords
!  character(len=40), allocatable, dimension(:) :: words
!  integer :: i
!  print*, "write test to console!"
!  read(*,"(a100)") string
!  call get_num_words(string, nwords, len_string)
!  print*, nwords
!  allocate(words(nwords))
!  call get_words(string, nwords, len_string, words, 40)
!  do i = 1, nwords, 1
!    print*, words(i)
!  end do
!end program test
