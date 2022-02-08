
subroutine write_ascii(filename, nat, rxyz, alat, atomnames, comment)
  !! writes the position and lattice vectors to a file in ascii format.
  implicit none
  character(len=*), intent(in) :: filename
  !! filename of the file which is created
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3, nat), intent(in) :: rxyz
  !! atom positions in bohr
  real*8, dimension(3, 3), intent(in) :: alat
  !! lattice vectors
  character(len=2), intent(in) :: atomnames(nat)
  !! String containing the chemical symbol of each atom.
  character(len=80), intent(in) :: comment
  !! content that will be written to the comment line

  ! private variables:
  integer :: io, ios, i
  real*8, dimension(3, nat) :: rxyz_copy
  real*8, dimension(3, 3) :: alat_copy
  real*8 :: Bohr_Ang = 0.52917721067

  alat_copy = alat
  rxyz_copy = rxyz
  call alat2ascii(nat, rxyz_copy, alat_copy)

  open (newunit=io, file=trim(adjustl(filename)), iostat=ios)
  if (ios /= 0) then
    print *, "Error opening file"//trim(filename)
    stop "Error opening file"
  end if

  !Transofrm output to angstroem
  alat_copy = alat_copy*Bohr_Ang
  rxyz_copy = rxyz_copy*Bohr_Ang

  write (io, *) trim(comment)
  write (unit=io, fmt=*) alat_copy(1, 1), alat_copy(1, 2), alat_copy(2, 2)
  write (unit=io, fmt=*) alat_copy(1, 3), alat_copy(2, 3), alat_copy(3, 3)
  do i = 1, nat, 1
    write (unit=io, fmt=*) rxyz_copy(:, i), atomnames(i)
  end do
  close (io)

contains
  subroutine alat2ascii(nat1, rxyz1, alat_io)
    !! Converts the lattice vectors in ascii format.
   implicit none
   integer, intent(in) :: nat1
   !! Number of atoms
   real*8, dimension(3, nat1), intent(inout) :: rxyz1
   !! Position of atoms.
   real*8, intent(inout), dimension(3, 3) :: alat_io
   !! Lattice vectors that should to be transformed.

   ! private variables:
   real*8, dimension(3, 3) :: alat1, t, alat_in_inv
   real*8 :: r1, r2, r3
   integer :: it
   alat1 = 0
   r1 = norm2(alat_io(:, 1))
   r2 = norm2(alat_io(:, 2))
   r3 = norm2(alat_io(:, 3))
   alat1(1, 1) = r1
   alat1(1, 2) = dot_product(alat_io(:, 1), alat_io(:, 2))/r1
   alat1(1, 3) = dot_product(alat_io(:, 1), alat_io(:, 3))/r1
   alat1(2, 2) = sqrt(r2**2 - alat1(1, 2)**2)
   alat1(2, 3) = (dot_product(alat_io(:, 2), alat_io(:, 3)) - alat1(1, 2)*alat1(1, 3)) &
                /alat1(2, 2)
   alat1(3, 3) = sqrt(r3**2 - alat1(1, 3)**2 - alat1(2, 3)**2)
   call matinv3(alat_io, alat_in_inv)
   t = matmul(alat1, alat_in_inv)
   do it = 1, nat1, 1
     rxyz1(:, it) = matmul(t, rxyz1(:, it))
   end do
   alat_io = alat1
   !call back2cell(nat1, rxyz1, alat1)
  end subroutine alat2ascii

  subroutine matinv3(A, B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    implicit none
    real*8, intent(in), dimension(3, 3) :: A
    !! Input Matrix
    real*8, intent(out), dimension(3, 3) :: B
    !! Inverse matrix
    real*8 :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) &
                - A(1, 2)*A(2, 1)*A(3, 3) + A(1, 2)*A(2, 3)*A(3, 1) &
                + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1))

    ! Calculate the inverse of the matrix
    B(1, 1) = +detinv*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))
    B(2, 1) = -detinv*(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))
    B(3, 1) = +detinv*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
    B(1, 2) = -detinv*(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))
    B(2, 2) = +detinv*(A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))
    B(3, 2) = -detinv*(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))
    B(1, 3) = +detinv*(A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))
    B(2, 3) = -detinv*(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))
    B(3, 3) = +detinv*(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
  end subroutine matinv3
end subroutine write_ascii
