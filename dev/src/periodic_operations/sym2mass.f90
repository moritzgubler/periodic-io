subroutine sym2amu(sym, amu)
  !! returns the atomic mass of atom with chemical symbol sym
  !! mass in amu taken from https://www.angelo.edu/faculty/kboudrea/periodic/structure_numbers.htm
  implicit none
  real*8  :: amu
  !! atomic mass unit
  character(len=2) :: sym
  !! chemical symbol
  select case (trim(sym))
  case ('H')
    amu = 1.00797
  case ('He')
    amu = 4.00260
  case ('Li')
    amu = 6.941
  case ('Be')
    amu = 9.01218
  case ('B')
    amu = 10.81
  case ('C')
    amu = 12.011
  case ('N')
    amu = 14.0067
  case ('O')
    amu = 15.9994
  case ('F')
    amu = 18.998403
  case ('Ne')
    amu = 20.179
  case ('Na')
    amu = 22.98977
  case ('Mg')
    amu = 24.305
  case ('Al')
    amu = 26.98154
  case ('Si')
    amu = 28.0855
  case ('P')
    amu = 30.97376
  case ('S')
    amu = 32.06
  case ('Cl')
    amu = 35.453
  case ('Ar')
    amu = 39.948
  case ('K')
    amu = 39.0983
  case ('Ca')
    amu = 40.08
  case ('Sc')
    amu = 44.9559
  case ('Ti')
    amu = 47.90
  case ('V')
    amu = 50.9415
  case ('Cr')
    amu = 51.996
  case ('Mn')
    amu = 54.9380
  case ('Fe')
    amu = 55.847
  case ('Co')
    amu = 58.9332
  case ('Ni')
    amu = 58.70
  case ('Cu')
    amu = 63.546
  case ('Zn')
    amu = 65.38
  case ('Ga')
    amu = 69.72
  case ('Ge')
    amu = 72.59
  case ('As')
    amu = 74.9216
  case ('Se')
    amu = 78.96
  case ('Br')
    amu = 79.904
  case ('Kr')
    amu = 83.80
  case ('Rb')
    amu = 85.4678
  case ('Sr')
    amu = 87.62
  case ('Y')
    amu = 88.9059
  case ('Zr')
    amu = 91.22
  case ('Nb')
    amu = 92.9064
  case ('Mo')
    amu = 95.94
  case ('Tc')
    amu = 98
  case ('Ru')
    amu = 101.07
  case ('Rh')
    amu = 102.9055
  case ('Pd')
    amu = 106.4
  case ('Ag')
    amu = 107.868
  case ('Cd')
    amu = 112.41
  case ('In')
    amu = 114.82
  case ('Sn')
    amu = 118.69
  case ('Sb')
    amu = 121.75
  case ('Te')
    amu = 127.60
  case ('I')
    amu = 126.9045
  case ('Xe')
    amu = 131.30
  case ('Cs')
    amu = 132.9054
  case ('Ba')
    amu = 137.33
  case ('La')
    amu = 138.9055
  case('Ce')
    amu = 140.12
  case('Pr')
    amu = 140.9077
  case('Nd')
    amu = 144.24
  case('Pm')
    amu = 145
  case('Sm')
    amu = 150.4
  case('Eu')
    amu = 151.96
  case('Gd')
    amu = 157.25
  case('Tb')
    amu = 158.9254
  case('Dy')
    amu = 162.50
  case('Ho')
    amu = 164.9304
  case('Er')
    amu = 167.26
  case('Tm')
    amu = 168.9342
  case('Yb')
    amu = 173.04
  case ('Lu')
    amu = 174.967
  case ('Hf')
    amu = 178.49
  case ('Ta')
    amu = 180.9479
  case ('W')
    amu = 183.85
  case ('Re')
    amu = 186.207
  case ('Os')
    amu = 190.2
  case ('Ir')
    amu = 192.22
  case ('Pt')
    amu = 195.09
  case ('Au')
    amu = 196.9665
  case ('Hg')
    amu = 200.59
  case ('Tl')
    amu = 204.37
  case ('Pb')
    amu = 207.2
  case ('Bi')
    amu = 208.9804
  case('Po')
    amu = 209
  case('At')
    amu = 210
  case ('Rn')
    amu = 222
  case('Fr')
    amu = 223
  case('Ra')
    amu = 226.0254
  case('Ac')
    amu = 227.0278
  case('Th')
    amu = 232.0381
  case('Pa')
    amu = 231.0359
  case('U')
    amu = 238.029
  case('Np')
    amu = 237.0482
  case('Pu')
    amu = 242
  case('Am')
    amu = 243
  case('Cm')
    amu = 247
  case default
    print *, " Not recognized atomic type "//sym
    stop
  end select

end subroutine sym2amu
