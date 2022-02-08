subroutine sym2rcov(sym, rcov)
! returns the covalent radius of atom with chemical symbol sym
  implicit none
  real*8  :: rcov
  character(len=2) :: sym  ! chemical symbol
  select case (trim(sym))
    ! covalet radius in Angstrom taken from WebElements: http://www.webelements.com/periodicity/covalent_radius/
  case ('H')
    rcov = 0.37d0
  case ('He')
    rcov = 0.32d0
  case ('Li')
    rcov = 1.34d0
  case ('Be')
    rcov = 0.90d0
  case ('B')
    rcov = 0.82d0
  case ('C')
    rcov = 0.77d0
  case ('N')
    rcov = 0.75d0
  case ('O')
    rcov = 0.73d0
  case ('F')
    rcov = 0.71d0
  case ('Ne')
    rcov = 0.69d0
  case ('Na')
    rcov = 1.54d0
  case ('Mg')
    rcov = 1.30d0
  case ('Al')
    rcov = 1.18d0
  case ('Si')
    rcov = 1.11d0
  case ('P')
    rcov = 1.06d0
  case ('S')
    rcov = 1.02d0
  case ('Cl')
    rcov = 0.99d0
  case ('Ar')
    rcov = 0.97d0
  case ('K')
    rcov = 1.96d0
  case ('Ca')
    rcov = 1.74d0
  case ('Sc')
    rcov = 1.44d0
  case ('Ti')
    rcov = 1.36d0
  case ('V')
    rcov = 1.25d0
  case ('Cr')
    rcov = 1.27d0
  case ('Mn')
    rcov = 1.39d0
  case ('Fe')
    rcov = 1.25d0
  case ('Co')
    rcov = 1.26d0
  case ('Ni')
    rcov = 1.21d0
  case ('Cu')
    rcov = 1.38d0
  case ('Zn')
    rcov = 1.31d0
  case ('Ga')
    rcov = 1.26d0
  case ('Ge')
    rcov = 1.22d0
  case ('As')
    rcov = 1.19d0
  case ('Se')
    rcov = 1.16d0
  case ('Br')
    rcov = 1.14d0
  case ('Kr')
    rcov = 1.10d0
  case ('Rb')
    rcov = 2.11d0
  case ('Sr')
    rcov = 1.92d0
  case ('Y')
    rcov = 1.62d0
  case ('Zr')
    rcov = 1.48d0
  case ('Nb')
    rcov = 1.37d0
  case ('Mo')
    rcov = 1.45d0
  case ('Tc')
    rcov = 1.56d0
  case ('Ru')
    rcov = 1.26d0
  case ('Rh')
    rcov = 1.35d0
  case ('Pd')
    rcov = 1.31d0
  case ('Ag')
    rcov = 1.53d0
  case ('Cd')
    rcov = 1.48d0
  case ('In')
    rcov = 1.44d0
  case ('Sn')
    rcov = 1.41d0
  case ('Sb')
    rcov = 1.38d0
  case ('Te')
    rcov = 1.35d0
  case ('I')
    rcov = 1.33d0
  case ('Xe')
    rcov = 1.30d0
  case ('Cs')
    rcov = 2.25d0
  case ('Ba')
    rcov = 1.98d0
  case ('La')
    rcov = 1.69d0
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
  case ('Lu')
    rcov = 1.60d0
  case ('Hf')
    rcov = 1.50d0
  case ('Ta')
    rcov = 1.38d0
  case ('W')
    rcov = 1.46d0
  case ('Re')
    rcov = 1.59d0
  case ('Os')
    rcov = 1.28d0
  case ('Ir')
    rcov = 1.37d0
  case ('Pt')
    rcov = 1.28d0
  case ('Au')
    rcov = 1.44d0
  case ('Hg')
    rcov = 1.49d0
  case ('Tl')
    rcov = 1.48d0
  case ('Pb')
    rcov = 1.47d0
  case ('Bi')
    rcov = 1.46d0
    !     case('Po')
    !     case('At')
  case ('Rn')
    rcov = 1.45d0
  case ('LJ')   ! Lennard Jones atom
    rcov = 0.25d0   ! Assuming distances are about 1
  case ('LA')   ! Lennard Jones atom
    rcov = 1.122462048309373d0
  case ('LB')  ! Lennard Jones atom
    rcov = 0.9877666025122482d0
    !     case('Fr')
    !     case('Ra')
    !     case('Ac')
    !     case('Th')
    !     case('Pa')
    !     case('U')
    !     case('Np')
    !     case('Pu')
    !     case('Am')
    !     case('Cm')
  case default
    print *, " Not recognized atomic type "//sym; stop
  end select

  rcov = rcov/0.52917720859d0   ! convert to atomic units

!  write(*,*) rcov

end subroutine sym2rcov
