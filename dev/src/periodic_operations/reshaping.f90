
subroutine reshapecell(nat, alat0, rxyz)
  !! reshapes to cell to a more convenient shape if possible
  implicit none
  integer, intent(in) :: nat
  !! number of atoms
  real*8, dimension(3, nat), intent(inout) :: rxyz
  !! cartesian coordinates of atoms
  real*8, dimension(3, 3), intent(out) :: alat0
  !! lattice matrix
  real*8, dimension(3, 3) :: alat
  !! lattice vectors
  real*8, dimension(3, 3) :: alatmin
  real*8 :: area, area1, area2, area3, areamin
  integer :: i, i12, i13, i21, i23, i31, i32, j, imax
  real*8 :: vol, vol0
  logical :: alatminset = .FALSE.

  imax = 6

  vol0 = abs(alat0(1, 1)*alat0(2, 2)*alat0(3, 3) - alat0(1, 1)*alat0(2, 3)*alat0(3, 2) - &
             alat0(1, 2)*alat0(2, 1)*alat0(3, 3) + alat0(1, 2)*alat0(2, 3)*alat0(3, 1) + &
             alat0(1, 3)*alat0(2, 1)*alat0(3, 2) - alat0(1, 3)*alat0(2, 2)*alat0(3, 1))

  areamin = 1.d100

  do i13 = -imax, imax
    do i12 = -imax, imax
      alat(1, 1) = alat0(1, 1) + i12*alat0(1, 2) + i13*alat0(1, 3)
      alat(2, 1) = alat0(2, 1) + i12*alat0(2, 2) + i13*alat0(2, 3)
      alat(3, 1) = alat0(3, 1) + i12*alat0(3, 2) + i13*alat0(3, 3)
      do i21 = -imax, imax
        do i23 = -imax, imax
          alat(1, 2) = alat0(1, 2) + i21*alat0(1, 1) + i23*alat0(1, 3)
          alat(2, 2) = alat0(2, 2) + i21*alat0(2, 1) + i23*alat0(2, 3)
          alat(3, 2) = alat0(3, 2) + i21*alat0(3, 1) + i23*alat0(3, 3)
          do i31 = -imax, imax
            do i32 = -imax, imax
              alat(1, 3) = alat0(1, 3) + i31*alat0(1, 1) + i32*alat0(1, 2)
              alat(2, 3) = alat0(2, 3) + i31*alat0(2, 1) + i32*alat0(2, 2)
              alat(3, 3) = alat0(3, 3) + i31*alat0(3, 1) + i32*alat0(3, 2)

              vol = abs(alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
                        alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
                        alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1))

              if (abs(vol - vol0) .lt. 1.d-8) then
                call cellsurface(alat(1, 1), alat(1, 2), area1)
                call cellsurface(alat(1, 1), alat(1, 3), area2)
                call cellsurface(alat(1, 2), alat(1, 3), area3)
                area = area1 + area2 + area3
                if (area .lt. areamin) then
                  alatminset = .TRUE.
                  areamin = area
                  do j = 1, 3; do i = 1, 3
                      alatmin(i, j) = alat(i, j)
                  end do; end do
                end if
              end if

            end do
          end do
        end do
      end do
    end do
  end do
  if ( alatminset ) then
    alat0 = alatmin
  else
    write(*, *) 'reshape cell would not have worked.'
    write(666, *) 'reshape cell would not have worked.'
    flush(666)
  end if

  call back2cell(nat, rxyz, alat0)

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

end subroutine reshapecell
