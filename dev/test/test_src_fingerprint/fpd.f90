!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Subroutines for FingerPrint !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine atoms_sphere(width_cutoff, nex_cutoff, lat, llat, ixyzmax, nat, natx_sphere, &
                        nat_sphere, alat, rxyz, rxyz_sphere, rcov, rcov_sphere, indat, amplitude, deramplitude)
  implicit none
  real*8 :: width_cutoff
  integer :: nex_cutoff
  integer :: lat, llat, ixyzmax, nat, natx_sphere
  real*8, dimension(3, 3) :: alat
  real*8, dimension(3, nat) :: rxyz
  real*8, dimension(nat) :: rcov
  real*8, dimension(3, natx_sphere) :: rxyz_sphere
  real*8, dimension(natx_sphere) :: rcov_sphere
  real*8, dimension(natx_sphere) :: amplitude
  real*8, dimension(natx_sphere) :: deramplitude
  integer, dimension(natx_sphere) :: indat
  integer :: nat_sphere

  real*8 :: dist2, factor_cutoff
  integer :: ix, iy, iz, jat
  real*8 :: radius_cutoff, radius_cutoff2, temp, xj, yj, zj

  radius_cutoff = sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2 = radius_cutoff**2
  factor_cutoff = 1.d0/(2.d0*nex_cutoff*width_cutoff**2)
  ! write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff

  nat_sphere = 0
  do jat = 1, nat
    do ix = -ixyzmax, ixyzmax
      do iy = -ixyzmax, ixyzmax
        do iz = -ixyzmax, ixyzmax
          xj = rxyz(1, jat) + ix*alat(1, 1) + iy*alat(1, 2) + iz*alat(1, 3)
          yj = rxyz(2, jat) + ix*alat(2, 1) + iy*alat(2, 2) + iz*alat(2, 3)
          zj = rxyz(3, jat) + ix*alat(3, 1) + iy*alat(3, 2) + iz*alat(3, 3)
          dist2 = (xj - rxyz(1, lat))**2 + (yj - rxyz(2, lat))**2 + (zj - rxyz(3, lat))**2
          ! write(*,*) xj,rxyz(1, lat),yj,rxyz(2, lat),zj,rxyz(3,lat)

          if (dist2 .le. radius_cutoff2) then
            nat_sphere = nat_sphere + 1
            if (nat_sphere .gt. natx_sphere) then
              print*, nat_sphere, natx_sphere
              stop 'enlarge natx_sphere h'
            end if
            !amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
            temp = (1.d0 - dist2*factor_cutoff)**(nex_cutoff - 1)
            amplitude(nat_sphere) = temp*(1.d0 - dist2*factor_cutoff)
            deramplitude(nat_sphere) = -2.d0*factor_cutoff*nex_cutoff*temp

            rxyz_sphere(1, nat_sphere) = xj
            rxyz_sphere(2, nat_sphere) = yj
            rxyz_sphere(3, nat_sphere) = zj
            rcov_sphere(nat_sphere) = rcov(jat)
            indat(nat_sphere) = jat
            if (dist2 .eq. 0.d0) llat = nat_sphere
          end if
        end do
      end do
    end do
  end do
end subroutine atoms_sphere

subroutine crtovrlp(nat, rxyz, alpha, cs, cp, ns, np, ovrlp)
  implicit none
  integer :: nat
  integer :: ns
  integer :: np
  real*8 :: rxyz(3, nat)
  real*8 :: ovrlp(nat*(ns + 3*np), nat*(ns + 3*np))
  real*8 :: alpha(nat), cs(10), cp(10)

  real*8 :: ai, aj
  integer :: iat, iorb, ip, is, jat, jp, js, jorb
  real*8 :: r2, sij, t1, t2, t3, t4, t5, xi, xij, xj, yi, yij, yj, zi, zij, zj

  if (ns > 10 .or. np > 10) stop 'ns > 10   .or.  np > 10  !'

  ! 1- setup the overlap matrix

  !  <s|s>
  do jat = 1, nat
    do js = 1, ns
      jorb = (jat - 1)*(ns + 3*np) + js
      aj = alpha(jat)/cs(js)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do is = 1, ns
          !!iorb=iat+(is-1)*nat
          iorb = (iat - 1)*(ns + 3*np) + is
          ai = alpha(iat)/cs(is)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2
          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          ovrlp(iorb, jorb) = sij

        end do
      end do
    end do
  end do

  !  <pi|sj>
  do jat = 1, nat
    do js = 1, ns

      jorb = (jat - 1)*(ns + 3*np) + js
      aj = alpha(jat)/cs(js)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do ip = 1, np
          !!iorb=1+(iat-1)*3+ns*nat + (ip-1)*3*nat
          iorb = (iat - 1)*(ns + 3*np) + ns + ip
          ai = alpha(iat)/cp(ip)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2

          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t3 = -2.d0*sqrt(ai)*aj/t2
          ovrlp(iorb, jorb) = t3*xij*sij
          ovrlp(iorb + 1, jorb) = t3*yij*sij
          ovrlp(iorb + 2, jorb) = t3*zij*sij

        end do
      end do
    end do
  end do

  !  <si|pj>
  do jat = 1, nat
    do jp = 1, np

      jorb = (jat - 1)*(ns + 3*np) + ns + jp
      aj = alpha(jat)/cp(jp)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do is = 1, ns
          !!iorb=iat+(is-1)*nat
          iorb = (iat - 1)*(ns + 3*np) + is
          ai = alpha(iat)/cs(is)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2

          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t3 = +2.d0*sqrt(aj)*ai/t2
          ovrlp(iorb, jorb) = t3*xij*sij
          ovrlp(iorb, jorb + 1) = t3*yij*sij
          ovrlp(iorb, jorb + 2) = t3*zij*sij

        end do
      end do
    end do
  end do

  !  <p|p>
  do jat = 1, nat
    do jp = 1, np

      jorb = (jat - 1)*(ns + 3*np) + ns + jp
      aj = alpha(jat)/cp(jp)
      xj = rxyz(1, jat); yj = rxyz(2, jat); zj = rxyz(3, jat)

      do iat = 1, nat
        do ip = 1, np
          iorb = (iat - 1)*(ns + 3*np) + ns + ip
          ai = alpha(iat)/cp(ip)
          xi = rxyz(1, iat); yi = rxyz(2, iat); zi = rxyz(3, iat)

          xij = xi - xj; yij = yi - yj; zij = zi - zj
          r2 = xij**2 + yij**2 + zij**2
          t1 = ai*aj
          t2 = ai + aj

          ! normalized GTOs:
          sij = sqrt(2.d0*sqrt(t1)/t2)**3*exp(-t1/t2*r2)
          t4 = 2.d0*sqrt(t1)/t2
          t5 = -2.d0*t1/t2

          ovrlp(iorb, jorb) = t4*(1.d0 + t5*xij*xij)*sij
          ovrlp(iorb + 1, jorb) = t4*(t5*yij*xij)*sij
          ovrlp(iorb + 2, jorb) = t4*(t5*zij*xij)*sij
          ovrlp(iorb, jorb + 1) = t4*(t5*xij*yij)*sij
          ovrlp(iorb + 1, jorb + 1) = t4*(1.d0 + t5*yij*yij)*sij
          ovrlp(iorb + 2, jorb + 1) = t4*(t5*zij*yij)*sij
          ovrlp(iorb, jorb + 2) = t4*(t5*xij*zij)*sij
          ovrlp(iorb + 1, jorb + 2) = t4*(t5*yij*zij)*sij
          ovrlp(iorb + 2, jorb + 2) = t4*(1.d0 + t5*zij*zij)*sij

        end do
      end do
    end do
  end do
end subroutine crtovrlp

subroutine multampspd(nat, ovrlp, amplitude, norb, ns, np, nd, ovrla)
  implicit none
  integer :: nat
  integer :: norb
  integer :: ns
  integer :: np
  integer :: nd
  real*8 :: amplitude(nat), ovrlp(norb, norb), ovrla(norb, norb)

  integer :: i, iat, iorb, j, jat, jorb

  do jat = 1, nat
    do j = 1, (ns + 3*np + 5*nd)
      jorb = (jat - 1)*(ns + 3*np + 5*nd) + j
      do iat = 1, nat
        do i = 1, (ns + 3*np + 5*nd)
          iorb = (iat - 1)*(ns + 3*np + 5*nd) + i
          ovrla(iorb, jorb) = ovrlp(iorb, jorb)*amplitude(iat)*amplitude(jat)
        end do
      end do
    end do
  end do

end subroutine multampspd

subroutine multampoff(nat, ovrlp, amplitude, norb, ns, np, ovrla)
  implicit none
  integer :: nat
  integer :: norb
  integer :: ns
  integer :: np
  real*8 :: amplitude(nat), ovrlp(norb, norb), ovrla(norb, norb)
  integer :: i, iat, iorb, j, jat, jorb

  do jat = 1, nat
    do j = 1, (ns + 3*np)
      jorb = (jat - 1)*(ns + 3*np) + j
      do iat = 1, nat
        do i = 1, (ns + 3*np)
          iorb = (iat - 1)*(ns + 3*np) + i
          if (iat .eq. jat) then
            ovrla(iorb, jorb) = ovrlp(iorb, jorb)
          else
            ovrla(iorb, jorb) = ovrlp(iorb, jorb)*amplitude(iat)*amplitude(jat)
          end if
        end do
      end do
    end do
  end do

end subroutine multampoff

subroutine multamp(nat, ovrlp, amplitude, norb, ns, np, ovrla)
  implicit none
  integer :: nat
  integer :: norb
  integer :: ns
  integer :: np
  real*8 :: amplitude(nat), ovrlp(norb, norb), ovrla(norb, norb)
  integer :: i, iat, iorb, j, jat, jorb

  do jat = 1, nat
    do j = 1, (ns + 3*np)
      jorb = (jat - 1)*(ns + 3*np) + j
      do iat = 1, nat
        do i = 1, (ns + 3*np)
          iorb = (iat - 1)*(ns + 3*np) + i
          ovrla(iorb, jorb) = ovrlp(iorb, jorb)*amplitude(iat)*amplitude(jat)
        end do
      end do
    end do
  end do

end subroutine multamp

subroutine calc_fpd(nat, natx_sphere, ns, np, fp1, fp2, fpd)
  implicit none
  integer :: nat, natx_sphere
  integer :: ns, np
  real*8 :: fp1((ns + 3*np)*natx_sphere, nat), fp2((ns + 3*np)*natx_sphere, nat)
  real*8 :: tt, cost1(nat, nat)
  real *8 :: fpd
  integer :: iat, jat, iassign(nat), l


  do iat = 1, nat
    do jat = 1, nat
      !tt = 0.d0
      !do l = 1, (ns + 3*np)*natx_sphere
      !  tt = tt + (fp1(l, iat) - fp2(l, jat))**2
      !end do
      !tt = sqrt(tt)
      tt = norm2(fp1(:,iat) -fp2(:,jat))
      cost1(iat, jat) = tt
    end do
  end do
  call apc(nat, cost1, iassign, fpd)

end subroutine calc_fpd

subroutine fingerprint_eval(nat, natx_sphere, ns, np, alat, rxyz, rcov, fp)
  implicit real*8(a - h, o - z)
  parameter(nwork=100)
  dimension workalat(nwork)
  dimension rxyz_sphere(3, natx_sphere), rcov_sphere(natx_sphere), indat(natx_sphere)
  dimension fp((ns + 3*np)*natx_sphere, nat), amplitude(natx_sphere), deramplitude(natx_sphere)
  dimension rxyz(3, nat), rcov(nat)!,eval((ns+3*np)*natx_sphere)
  dimension alat(3, 3), alatalat(3, 3), eigalat(3)
  dimension alpha(natx_sphere), cs(10), cp(10)
  real*8, allocatable ::  ovrlp(:, :), ovrla(:, :), eval(:)

! parameters for cutoff function: width_cutoff is the width of the Gauusian
! approximated by a polynomial with exponent nex_cutoff
  width_cutoff = 3.d0
  nex_cutoff = 2
! The following line has to be idential to the corresponding line in xyz2devaldr
  radius_cutoff = sqrt(2.d0*nex_cutoff)*width_cutoff  ! radius where the polynomial is zero

  if (alat(1, 1)*alat(2, 2)*alat(3, 3) .eq. 0.d0) then
    ixyzmax = 0
  else

    do i = 1, 3
    do j = 1, 3
      alatalat(i, j) = alat(1, i)*alat(1, j) + alat(2, i)*alat(2, j) + alat(3, i)*alat(3, j)
    end do
    end do
    call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
    !  write(*,*) 'alat EVals',eigalat
    !  write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
    ixyzmax = int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1
  end if
!    write(*,*) 'ixyzmax ',ixyzmax

  do i = 1, 10
    cs(i) = sqrt(2.d0)**(i - 1)
    cp(i) = sqrt(2.d0)**(i - 1)
  end do

! loop over all center atoms
  !natsmax = 0
  !natsmin = 1000000
  !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(FP)
  !$omp do
  do lat = 1, nat
    call atoms_sphere(width_cutoff, nex_cutoff, lat, llat, ixyzmax, nat, natx_sphere, nat_sphere, alat, rxyz, rxyz_sphere, &
                      rcov, rcov_sphere, indat, amplitude, deramplitude)
    ! write(*,*) lat,rcov(lat),nat_sphere
    !natsmax = max(natsmax, nat_sphere)
    !natsmin = min(natsmin, nat_sphere)
!       call xyz2eval(nat_sphere,rxyz_sphere,rcov_sphere,ns,np,amplitude,fp(1,lat))

    norb = nat_sphere*(ns + np*3)
    allocate (ovrlp(norb, norb), ovrla(norb, norb), eval(norb))

    !do iat = 1, nat_sphere
    !  alpha(iat) = .5d0/(1.0d0*rcov_sphere(iat))**2
    !end do
    alpha(1:nat_sphere) = 0.5/(rcov_sphere(1:nat_sphere))**2

    call crtovrlp(nat_sphere, rxyz_sphere, alpha, cs, cp, ns, np, ovrlp)
    call multamp(nat_sphere, ovrlp, amplitude, norb, ns, np, ovrla)
    !trace = 0.d0
    !do i = 1, norb
    !  trace = trace + ovrla(i, i)
    !end do
    !tracein = 1.d0/trace

    call diagonalizematrix(norb, ovrla, eval)
! eigenvalues in decreasing order
    !evals = 0.d0
    !do i = 1, norb
    !  evals = evals + eval(norb - i + 1)
      !fp(i,lat)=sqrt(eval(norb-i+1)*tracein)
    !  fp(i, lat) = eval(norb - i + 1)
    !end do
    fp(1:norb, lat) = eval(norb:1:-1)
    norbx = natx_sphere*(ns + np*3)
    fp( (norb+1):norbx, lat ) = 0.0
    !do i = norb + 1, norbx
    !  fp(i, lat) = 0.d0
    !end do
    deallocate (ovrlp, ovrla, eval)
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  !write (*, *) 'min, max number of atoms in sphere ', natsmin, natsmax
end subroutine fingerprint_eval

subroutine diagonalizeMatrix(n, aa, eval)
  implicit real*8(a - h, o - z)

  ! Calling arguments
  real*8, dimension(n, n) :: aa
  real*8, dimension(n) :: eval

  ! Local variables
  real*8, dimension(:), allocatable:: work

  lwork = 100*n
  allocate (work(lwork))
  call dsyev('v', 'l', n, aa, n, eval, work, lwork, info)
  if (info /= 0) stop ' ERROR in dsyev'
  deallocate (work)

end subroutine diagonalizeMatrix
