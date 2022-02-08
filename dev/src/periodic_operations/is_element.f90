subroutine isElement(atom_in, is_element)
  !! returns wether chemical symbol is an element or not
  implicit none
  character(len=2), intent(in) :: atom_in
  !! String containing chemical symbol
  logical, intent(out) :: is_element
  !! logical, true if string is element, false if not
  select case (trim(atom_in))
  case('H')
     is_element = .true.
  case('He')
     is_element = .true.
  case('Li')
     is_element = .true.
  case('Be')
     is_element = .true.
  case('B')
     is_element = .true.
  case('C')
     is_element = .true.
  case('N')
     is_element = .true.
  case('O')
     is_element = .true.
  case('F')
     is_element = .true.
  case('Ne')
     is_element = .true.
  case('Na')
     is_element = .true.
  case('Mg')
     is_element = .true.
  case('Al')
     is_element = .true.
  case('Si')
     is_element = .true.
  case('P')
     is_element = .true.
  case('S')
     is_element = .true.
  case('Cl')
     is_element = .true.
  case('Ar')
     is_element = .true.
  case('K')
     is_element = .true.
  case('Ca')
     is_element = .true.
  case('Sc')
     is_element = .true.
  case('Ti')
     is_element = .true.
  case('V')
     is_element = .true.
  case('Cr')
     is_element = .true.
  case('Mn')
     is_element = .true.
  case('Fe')
     is_element = .true.
  case('Co')
     is_element = .true.
  case('Ni')
     is_element = .true.
  case('Cu')
     is_element = .true.
  case('Zn')
     is_element = .true.
  case('Ga')
     is_element = .true.
  case('Ge')
     is_element = .true.
  case('As')
     is_element = .true.
  case('Se')
     is_element = .true.
  case('Br')
     is_element = .true.
  case('Kr')
     is_element = .true.
  case('Rb')
     is_element = .true.
  case('Sr')
     is_element = .true.
  case('Y')
     is_element = .true.
  case('Zr')
     is_element = .true.
  case('Nb')
     is_element = .true.
  case('Mo')
     is_element = .true.
  case('Tc')
     is_element = .true.
  case('Ru')
     is_element = .true.
  case('Rh')
     is_element = .true.
  case('Pd')
     is_element = .true.
  case('Ag')
     is_element = .true.
  case('Cd')
     is_element = .true.
  case('In')
     is_element = .true.
  case('Sn')
     is_element = .true.
  case('Sb')
     is_element = .true.
  case('Te')
     is_element = .true.
  case('I')
     is_element = .true.
  case('Xe')
     is_element = .true.
  case('Cs')
     is_element = .true.
  case('Ba')
     is_element = .true.
  case('La')
     is_element = .true.
     !     case('Ce')
     !     case('Pr')
     !     case('Nd')
     !     case('Pm')
     !     case('Sm')
     !     case('Eu')
     !     case('Gd')
     !     case('Tb')
     !     case('Dy')
     !     case('Ho')
     !     case('Er')
     !     case('Tm')
     !     case('Yb')
  case('Lu')
     is_element = .true.
  case('Hf')
     is_element = .true.
  case('Ta')
     is_element = .true.
  case('W')
     is_element = .true.
  case('Re')
     is_element = .true.
  case('Os')
     is_element = .true.
  case('Ir')
     is_element = .true.
  case('Pt')
     is_element = .true.
  case('Au')
     is_element = .true.
  case('Hg')
     is_element = .true.
  case('Tl')
     is_element = .true.
  case('Pb')
     is_element = .true.
  case('Bi')
     is_element = .true.
     !     case('Po')
     !     case('At')
  case('Rn')
     is_element = .true.
    case default
      is_element = .FALSE.
  end select
end subroutine isElement
