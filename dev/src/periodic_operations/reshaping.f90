subroutine reshapecell(nat, alat0, rxyz)
  !! reshapes to cell to a more convenient shape if possible
  implicit none
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3, nat), intent(inout) :: rxyz
  !! cartesian coordinates of atoms
  real*8, dimension(3, 3), intent(inout) :: alat0
  !! lattice matrix
  real(8) :: p(3, 3), pos_temp(3, nat), pinv(3, 3)
  integer :: imax
  logical :: success

  call alat2ascii(nat, rxyz, alat0)
  imax = 6

  ! try reshaping in original upper triangular lattice matrix
  call reshape_ascii(alat0, imax, success)

  ! search in first permutation of lattice vectors:
  ! 1 -> 3
  ! 2 -> 1
  ! 3 -> 2
  p = 0.d0
  p(1, 2) = 1.d0
  p(2, 3) = 1.d0
  p(3, 1) = 1.d0
  call matinv3(p, pinv)
  alat0 = matmul(alat0, p)
  call alat2ascii(nat, rxyz, alat0)
  call reshape_ascii(alat0, imax, success)
  if ( success ) then
    alat0 = matmul(alat0, pinv)
  end if
  call back2cell(nat, rxyz, alat0)

  ! search in second permuation of lattice vectors:
  ! 1 -> 2
  ! 2 -> 3
  ! 3 -> 1
  pos_temp = rxyz 
  p = 0.d0
  p(1, 3) = 1.d0
  p(2, 1) = 1.d0
  p(3, 2) = 1.d0
  call matinv3(p, pinv)
  alat0 = matmul(alat0, p)
  call alat2ascii(nat, pos_temp, alat0)
  call reshape_ascii(alat0, imax, success)
  if ( success ) then
    rxyz = pos_temp
    alat0 = matmul(alat0, pinv)
    call back2cell(nat, rxyz, alat0)
  end if

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

end subroutine reshapecell

subroutine reshape_ascii(alat0, imax, new_lattice_found)
  implicit none
  real(8), intent(inout) :: alat0(3, 3)
  integer, intent(in) :: imax
  logical, intent(out) :: new_lattice_found
  real(8) :: vol0, alat(3, 3), vol, area, areamin
  real(8) :: area1, area2, area3, alatmin(3, 3)
  integer :: iba, ica, icb
  new_lattice_found = .false.
  vol0 = alat0(1,1) * alat0(2,2) * alat0(3, 3)
  call cellsurface(alat0(1, 1), alat0(1, 2), area1)
  call cellsurface(alat0(1, 1), alat0(1, 3), area2)
  call cellsurface(alat0(1, 2), alat0(1, 3), area3)
  areamin = area1 + area2 + area3
  alat(:, 1) = alat0(:, 1)
  do iba = - imax, imax
    alat(:, 2) = alat0(:, 2) + iba * alat0(:,1)
    do ica = -imax, imax
      do icb = -imax, imax
        alat(:, 3) = alat0(:, 3) + ica * alat0(:, 1) + icb * alat0(:, 2)
        vol = alat(1, 1) * alat(2, 2) * alat(3,3)
        call cellsurface(alat(1, 1), alat(1, 2), area1)
        call cellsurface(alat(1, 1), alat(1, 3), area2)
        call cellsurface(alat(1, 2), alat(1, 3), area3)
        area = area1 + area2 + area3
        if (area .lt. areamin) then
          new_lattice_found = .true.
          areamin = area
          alatmin = alat 
        end if
      end do
    end do  
  end do
  alat0 = alatmin 

contains

  subroutine cellsurface(a, b, area_in)
    !! calculates surface spanned by two 3d vectors
    implicit none
    real*8, dimension(3) :: a
    real*8, dimension(3) :: b
    real*8, dimension(3) :: c
    real*8 :: area_in

    c(1) = a(2)*b(3) - b(2)*a(3)
    c(2) = a(3)*b(1) - b(3)*a(1)
    c(3) = a(1)*b(2) - b(1)*a(2)
    area_in = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
  end subroutine cellsurface

end subroutine reshape_ascii
