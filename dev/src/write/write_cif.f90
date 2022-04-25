subroutine write_cif(filename, nat, rxyz, alat, atomnames, comment)
  !! writes a peridic structure to a file with the cif standard.
  implicit none
  integer :: nat
  !! number of atoms
  character(len=*) :: filename
  !! name of the file
  real*8, intent(in), dimension(3, nat) :: rxyz
  !! positions of the atoms (units in bohr)
  real*8, dimension(3,3), intent(in) :: alat
  !! lattice vectors (units in bohr)
  character(len=2), dimension(nat), intent(in) :: atomnames
  !! chemical names of the atoms
  character(len=80), intent(in) :: comment
  !! content that will be written to the comment line
  integer :: io, ios, i
  real*8 :: length_a, length_b, length_c
  real*8 :: alpha, beta, gamma
  real*8 :: xyzred(3,nat)
  real*8 :: Bohr_Ang = 0.52917721067
  real*8 :: alat_convert(3,3)

  alat_convert = Bohr_Ang * alat
  call cif_angles()
  call cart2frac(nat, alat, rxyz, xyzred)

  open(newunit=io, file=filename, iostat=ios)
  if ( ios /= 0 ) then
    print*, "cif filename: ", trim(filename)
    stop "error opening file in write_cif"
  end if
  if(len_trim(comment) > 0) then
    write(io, *) "# "//trim(comment)
  end if
  write(io,*) "data_cif"
  write(io, "(a)") "_space_group_name_H-M_alt 'P 1'"
  write(io,"(a, 1x, f20.15)") "_cell_length_a", length_a
  write(io,"(a, 1x, f20.15)") "_cell_length_b", length_b
  write(io,"(a, 1x, f20.15)") "_cell_length_c", length_c
  write(io,"(a, 1x, f20.15)") "_cell_angle_alpha", alpha
  write(io,"(a, 1x, f20.15)") "_cell_angle_beta", beta
  write(io,"(a, 1x, f20.15)") "_cell_angle_gamma", gamma

  write(io, "(a)") "loop_"
  write(io, "(a)") " _atom_site_type_symbol"
  write(io, "(a)") " _atom_site_label"
  write(io, "(a)") " _atom_site_fract_x"
  write(io, "(a)") " _atom_site_fract_y"
  write(io, "(a)") " _atom_site_fract_z"

  do i = 1, nat, 1
    write(io, "(2x, a2, 1x, a, i3.3, 1x, 3f18.14)") atomnames(i)&
          , trim(atomnames(i)), i, xyzred(1,i), xyzred(2,i), xyzred(3,i)
  end do
  close(io)
contains
  subroutine cif_angles
    implicit none
    real*8, parameter :: converts = 57.295779513
    !! converts radian to degree
    length_a = norm2(alat_convert(:,1))
    length_b = norm2(alat_convert(:,2))
    length_c = norm2(alat_convert(:,3))
    alpha = acos(dot_product(alat_convert(:,2), alat_convert(:,3)) / (length_b * length_c))
    beta = acos(dot_product(alat_convert(:,1), alat_convert(:,3)) / (length_a * length_c))
    gamma = acos(dot_product(alat_convert(:,1), alat_convert(:,2)) / (length_a * length_b))
    alpha = converts * alpha
    beta = converts * beta
    gamma = converts * gamma
  end subroutine cif_angles
end subroutine write_cif
