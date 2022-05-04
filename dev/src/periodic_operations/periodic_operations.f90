
subroutine cart2frac(nat, alat, rxyz, xyzred)
  !! converts cartesian coordinates rxyz to reduced coordinates xyzred
  implicit none
  integer, intent(in) :: nat
  !! Number of Atoms
  real*8, intent(in), dimension(3, 3) :: alat
  !! Lattice Vectors.
  real*8, intent(in), dimension(3, nat) :: rxyz
  !! Position of the Atoms in cartesian coorinates.
  real*8, intent(out), dimension(3, nat) :: xyzred
  !! Position of the Atoms in reduced coordinates.

  !private variables
  real*8, dimension(3, 3) :: alatinv
  !! inverse of the lattice matrix
  real*8 :: div
  !! inverse volume of the lattice matrix

  div = alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
        alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
        alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1)
  if ( abs(div) < 1.d-7 ) then
    stop "cell is singular in cart2frac"
  end if
  div = 1.d0/div
  alatinv(1, 1) = (alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2))
  alatinv(1, 2) = -(alat(1, 2)*alat(3, 3) - alat(1, 3)*alat(3, 2))
  alatinv(1, 3) = (alat(1, 2)*alat(2, 3) - alat(1, 3)*alat(2, 2))
  alatinv(2, 1) = -(alat(2, 1)*alat(3, 3) - alat(2, 3)*alat(3, 1))
  alatinv(2, 2) = (alat(1, 1)*alat(3, 3) - alat(1, 3)*alat(3, 1))
  alatinv(2, 3) = -(alat(1, 1)*alat(2, 3) - alat(1, 3)*alat(2, 1))
  alatinv(3, 1) = (alat(2, 1)*alat(3, 2) - alat(2, 2)*alat(3, 1))
  alatinv(3, 2) = -(alat(1, 1)*alat(3, 2) - alat(1, 2)*alat(3, 1))
  alatinv(3, 3) = (alat(1, 1)*alat(2, 2) - alat(1, 2)*alat(2, 1))
  alatinv = alatinv * div
  xyzred = matmul(alatinv, rxyz)
end subroutine cart2frac

subroutine frac2cart(nat, alat, xyzred, rxyz)
  !! Converts reduced coordinates xyzred to cartesian coordinates rxyz
  implicit none
  integer, intent(in) :: nat
  !! Number of atoms.
  real*8, intent(in), dimension(3, 3) :: alat
  !! Lattice Vecors
  real*8, dimension(3, nat), intent(in) :: xyzred
  !! Position of the atoms in reduced coordinates.
  real*8, dimension(3, nat), intent(out) :: rxyz
  !! Position of the atoms in cartesian coordinates.
  rxyz = matmul(alat, xyzred)
end subroutine frac2cart

subroutine back2cell(nat, rxyz, alat)
  !! Translates atoms outside the cell back into the cell.
  implicit none
  integer, intent(in) :: nat
  !! Number of atoms
  real*8, dimension(3, nat), intent(inout) :: rxyz
  !! Positions of the atoms.
  real*8, dimension(3, 3), intent(in) :: alat
  !! Lattice vectors of the atoms.

  ! private variables
  real*8, dimension(3, nat) :: xyzred
  !! Position of the atoms in cartesian coordinates.

  call cart2frac(nat, alat, rxyz, xyzred)
  xyzred = modulo(xyzred, 1.d0)
  call frac2cart(nat, alat, xyzred, rxyz)
end subroutine back2cell
