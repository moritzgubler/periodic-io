! Copyright (C) 2021 Marco Krummenacher
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.



!------------------------------------------------------------------------------
! OM Fingerprint Package
!------------------------------------------------------------------------------
!
! MODULE: Fingerprint
!
!> Marco Krummenacher
!> University of Basel
!
! DESCRIPTION:
!> In this module the Overlap Matrix Fingerprint (OMFP) and its derivative are
!> be calculated in an efficient way.
!
!------------------------------------------------------------------------------

! DESCRIPTION:
!> This routine calculates the OM finperprint vector.
!
!> @param[in] nat number of atoms in molecule or unit cell
!> @param[in] nat_sphere_max maximal number of atoms in sphere
!> @param[in] ns if 0 NO s-orbitals; if 1 s-orbitals
!> @param[in] np if 0 NO p-orbitals; if 1 px, py, pz -orbitals
!> @param[in] width_cutoff width of the cutoff sphere
!> @param[in] alat unit cell vectors; if cluster alat = 0.0
!> @param[in] rxyz xyz-positions of the atoms
!> @param[in] symb atomic symbols
!> @param[out] fp fingerprint vectors
!---------------------------------------------------------------------------
  SUBROUTINE fingerprint(nat, nat_sphere_max, ns, np, width_cutoff, alat,rxyz, symb, fp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nat      !number of atoms in system
    INTEGER, INTENT(IN) :: ns       !s-orbitals
    INTEGER, INTENT(IN) :: np       !p-orbitals
    INTEGER, INTENT(IN) :: nat_sphere_max  !max. number of atoms in sphere
    REAL(8), INTENT(IN) :: width_cutoff  !cutoff radius
    REAL(8), DIMENSION(3,nat), INTENT(IN) :: rxyz  !coorinates
    REAL(8), DIMENSION(3,3), INTENT(IN) :: alat  !unit cell vectors, if a cluster then alat=0.d0
    REAL(8), DIMENSION(nat_sphere_max*(ns+3*np),nat), INTENT(OUT) :: fp  !fingerprint vectors

    INTEGER, PARAMETER :: nex_cutoff = 2
    INTEGER, PARAMETER :: nwork = 100
    INTEGER :: nenv
    INTEGER :: len_fp  !length of fingerprint
    INTEGER :: ixyzmax
    INTEGER :: ienv
    INTEGER :: llat
    INTEGER :: iorb
    INTEGER :: i
    INTEGER :: j
    INTEGER :: norb
    INTEGER :: nat_sphere
    INTEGER :: iat
    INTEGER :: info
    INTEGER, DIMENSION(nat_sphere_max) :: indat

    REAL(8) :: radius_cutoff
    REAL(8) :: t1
    REAL(8) :: t2
    REAL(8) :: volume

    REAL(8), DIMENSION(3,nat_sphere_max) :: rxyz_sphere
    REAL(8), DIMENSION(nat) :: rcov
    REAL(8), DIMENSION(nat_sphere_max) :: rcov_sphere
    REAL(8), DIMENSION(3,3) :: alatalat
    REAL(8), DIMENSION(3) :: eigalat
    REAL(8), DIMENSION(nwork) :: workalat
    REAL(8), DIMENSION(nat_sphere_max) :: amplitude_tmp
    REAL(8), DIMENSION(nat_sphere_max) :: deramplitude_tmp
    REAL(8), DIMENSION(nat_sphere_max) :: alpha
    REAL(8), DIMENSION(10) :: cs
    REAL(8), DIMENSION(10) :: cp

    REAL(8), ALLOCATABLE :: ovrlp(:,:)
    REAL(8), ALLOCATABLE :: ampovrlp(:,:)
    REAL(8), ALLOCATABLE :: eval(:)
    REAL(8), ALLOCATABLE :: amplitude(:)
    REAL(8), ALLOCATABLE :: deramplitude(:)

    CHARACTER(len=2), DIMENSION(nat) :: symb

    len_fp = nat_sphere_max*(ns+3*np)
    nenv = nat

    !calculate covalent radii corresponding to atom type
    DO iat = 1, nat
      CALL sym2rcov_om_fp(symb(iat), rcov(iat))
    ENDDO
    fp = 0.d0

    IF (ns .gt. 3) THEN
      WRITE(*,*) "WARNING: NUMBER OF S-ORBITALS GREATER THAN 4"
      WRITE(*,*) "WARNING: NOT TESTED!"
    ENDIF

    IF (np .gt. 1) THEN
      WRITE(*,*) "WARNING: NUMBER OF P-ORBITALS GREATER THAN 1"
      WRITE(*,*) "WARNING: NOT TESTED!"
    ENDIF
    !loop over all atoms to get a fingerprint vector for each environment

    radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff

    volume = abs(alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
                 alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
                 alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1))


    !$omp parallel private(ienv, ixyzmax, alatalat, amplitude,&
    !$omp                  llat, nat_sphere, rxyz_sphere, rcov_sphere, indat,&
    !$omp                  deramplitude, alpha, cs, cp, norb, ovrlp, ampovrlp,&
    !$omp                  eval,eigalat,workalat, amplitude_tmp, deramplitude_tmp)
    !$omp do schedule(dynamic)

    DO ienv=1, nenv

      IF (volume.lt.0.d0 ) THEN ! no periodic boundary condition
          ixyzmax=0
      ELSE  ! periodic boundary conditions
        DO i=1,3
          DO j=1,3
            alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
          ENDDO
        ENDDO
        CALL dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
        ! ixyzmax determines over how many periodiv images one has to search to fill the sphere with atoms
        ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1
      ENDIF

      amplitude_tmp = 0.d0
      deramplitude_tmp = 0.d0
      !determine the atoms which are within the sphere
      CALL atoms_sphere(width_cutoff,nex_cutoff,ienv,llat,ixyzmax,nat,nat_sphere_max,nat_sphere,alat,rxyz,rxyz_sphere, &
                              rcov,rcov_sphere,indat,amplitude_tmp,deramplitude_tmp)

      ALLOCATE(amplitude(nat_sphere))
      ALLOCATE(deramplitude(nat_sphere))
      DO iat=1,nat_sphere
        alpha(iat)=.5d0/rcov_sphere(iat)**2
        amplitude(iat) = amplitude_tmp(iat)
        deramplitude(iat) = deramplitude_tmp(iat)
      ENDDO
      ! Specify the width of the Gaussians if several Gaussians per l-channel are used


      DO i=1,10
        cs(i)=sqrt(2.d0)**(float(i)-((float(ns)+1.d0)/2.d0))
        cp(i)=sqrt(2.d0)**(i-1)
      ENDDO
      norb = (ns+3*np)*nat_sphere


      ALLOCATE(ovrlp(norb,norb))
      ALLOCATE(ampovrlp(norb,norb))
      ALLOCATE(eval(norb))
      ovrlp = 0.d0
      ampovrlp = 0.d0
      eval = 0.d0

      ! calculate overlapmatrix
      CALL crtovrlp(nat_sphere,rxyz_sphere,alpha,cs,cp,ns,np,ovrlp)
      ! multiply cutoff function to overlapmatrix
      CALL multamp(nat_sphere,ovrlp,amplitude,norb,ns,np,ampovrlp)
      ! generate eigenvalues
      CALL diagonalizeMatrix(norb, ampovrlp, eval)

      ! eigenvalues in descending order
      DO i=1,norb/2
         t1=eval(i)
         t2=eval(norb-i+1)
         eval(i)=t2
         eval(norb-i+1)=t1
      ENDDO

      DO iorb=1,norb
        fp(iorb,ienv) = eval(iorb)
      ENDDO

      DEALLOCATE(amplitude)
      DEALLOCATE(deramplitude)
      DEALLOCATE(ovrlp)
      DEALLOCATE(ampovrlp)
      DEALLOCATE(eval)
      !write(*,*) "Fingerprint of environment", ienv, " : DONE"
    ENDDO
    !$omp end do
    !$omp end parallel

  END SUBROUTINE fingerprint

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

! subroutines used for calculations
!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> Selecting the atoms which are within the sphere of the central atom.
!> @param[in] width_cutoff width of the cutoff sphere
!> @param[in] nex_cutoff ??
!> @param[in] lat array position of central atom in rxyz
!> @param[out] llat array position of central atom in rxyz_sphere
!> @param[in] ixyzmax number of neighboring cells to build the sphere
!> @param[in] nat number of atoms in molecule or unit cell
!> @param[in] natx_sphere maximal number of atoms in sphere
!> @param[out] nat_sphere number of atoms in the sphere
!> @param[in] alat unit cell vectors; if cluster alat = 0.0
!> @param[in] rxyz xyz-positions of the atoms
!> @param[out] rxyz_sphere xyz-positions of the atoms in the sphere
!> @param[in] rcov covalent radii of the atoms corresponding to rxyz
!> @param[out] rcov_sphere covalent radii of the atoms corresponding to rxyz_sphere
!> @param[out] indat connecting rxyz and rxyz_sphere
!> @param[out] amplitude cutoff function values
!> @param[out] deramplitude derivative of the cutoff function
!---------------------------------------------------------------------------



  subroutine atoms_sphere(width_cutoff,nex_cutoff,lat,llat,ixyzmax,nat,natx_sphere,nat_sphere,alat,rxyz,rxyz_sphere, &
                          rcov,rcov_sphere,indat,amplitude,deramplitude)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nex_cutoff
    INTEGER, INTENT(IN) :: lat
    INTEGER, INTENT(IN) :: ixyzmax
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: natx_sphere
    INTEGER, INTENT(OUT) :: llat
    INTEGER, INTENT(OUT) :: nat_sphere
    INTEGER :: ix
    INTEGER :: iy
    INTEGER :: iz
    INTEGER :: jat
    REAL(8) :: xj
    REAL(8) :: yj
    REAL(8) :: zj
    REAL(8) :: radius_cutoff
    REAL(8) :: radius_cutoff2
    REAL(8) :: temp
    REAL(8), INTENT(IN) :: width_cutoff
    REAL(8) :: dist2
    REAL(8) :: factor_cutoff
    REAL(8), DIMENSION(3,natx_sphere), INTENT(OUT) :: rxyz_sphere
    REAL(8), DIMENSION(natx_sphere), INTENT(OUT) :: rcov_sphere
    REAL(8), DIMENSION(natx_sphere), INTENT(OUT) :: amplitude
    REAL(8), DIMENSION(natx_sphere), INTENT(OUT) :: deramplitude
    REAL(8), DIMENSION(3,nat), INTENT(IN) :: rxyz
    REAL(8), DIMENSION(nat), INTENT(IN) :: rcov
    REAL(8), DIMENSION(3,3), INTENT(IN) :: alat
    INTEGER, DIMENSION(natx_sphere), INTENT(OUT) :: indat


    radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff
    radius_cutoff2=radius_cutoff**2
    factor_cutoff=1.d0/(2.d0*nex_cutoff*width_cutoff**2)
  !!  write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff
       nat_sphere=0
       DO jat = 1, nat
           DO ix = -ixyzmax,ixyzmax
             DO iy = -ixyzmax,ixyzmax
               DO iz = -ixyzmax,ixyzmax
                  xj = rxyz(1, jat) + ix*alat(1,1)+iy*alat(1,2)+iz*alat(1,3)
                  yj = rxyz(2, jat) + ix*alat(2,1)+iy*alat(2,2)+iz*alat(2,3)
                  zj = rxyz(3, jat) + ix*alat(3,1)+iy*alat(3,2)+iz*alat(3,3)
                  dist2 = (xj-rxyz(1, lat))**2+(yj-rxyz(2, lat))**2+(zj-rxyz(3, lat))**2
                  !write(*,*) xj,rxyz(1, lat),yj,rxyz(2, lat),zj,rxyz(3,lat)


                  IF (dist2.le.radius_cutoff2) THEN
                      nat_sphere=nat_sphere+1
                      IF (nat_sphere.gt.natx_sphere) STOP 'enlarge natx_sphere'
                      !amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
                      temp = (1.d0 - dist2*factor_cutoff)**(nex_cutoff-1)
                      amplitude(nat_sphere) = temp*(1.d0 - dist2*factor_cutoff)
                      deramplitude(nat_sphere) = -2.d0*factor_cutoff*nex_cutoff*temp

                      rxyz_sphere(1,nat_sphere)=xj
                      rxyz_sphere(2,nat_sphere)=yj
                      rxyz_sphere(3,nat_sphere)=zj
                      rcov_sphere(nat_sphere)=rcov(jat)
                      indat(nat_sphere)=jat
                      IF (dist2.eq.0.d0) llat=nat_sphere
                  ENDIF
               ENDDO
             ENDDO
           ENDDO
       ENDDO
  end subroutine atoms_sphere

  !---------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !> Construction of the overlap matrix (OM)
  !> @param[in] nat number of atoms in molecule or unit cell
  !> @param[in] rxyz xyz-positions of the atoms to cunstruct the OM
  !> @param[in] alpha ??
  !> @param[in] cs ??
  !> @param[in] cp ??
  !> @param[in] ns if 0 NO s-orbitals; if 1 s-orbitals
  !> @param[in] np if 0 NO p-orbitals; if 1 px, py, pz -orbitals
  !> @param[out] ovrlp overlap matrix
  !---------------------------------------------------------------------------
  subroutine crtovrlp(nat,rxyz,alpha,cs,cp,ns,np,ovrlp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: ns
    INTEGER, INTENT(IN) :: np
    INTEGER :: iat
    INTEGER :: jat
    INTEGER :: iorb
    INTEGER :: jorb
    INTEGER :: is
    INTEGER :: ip
    INTEGER :: js
    INTEGER :: jp
    REAL(8) :: r2
    REAL(8) :: sij
    REAL(8) :: ai
    REAL(8) :: aj
    REAL(8) :: t1
    REAL(8) :: t2
    REAL(8) :: t3
    REAL(8) :: t4
    REAL(8) :: t5
    REAL(8) :: xi
    REAL(8) :: yi
    REAL(8) :: zi
    REAL(8) :: xj
    REAL(8) :: yj
    REAL(8) :: zj
    REAL(8) :: xij
    REAL(8) :: yij
    REAL(8) :: zij
    REAL(8), DIMENSION(3,nat), INTENT(IN) :: rxyz
    REAL(8), DIMENSION(nat*(ns+3*np),nat*(ns+3*np)), INTENT(OUT) :: ovrlp
    REAL(8), DIMENSION(nat), INTENT(IN) :: alpha
    REAL(8), DIMENSION(10), INTENT(IN) :: cs
    REAL(8), DIMENSION(10), INTENT(IN) :: cp


    IF(ns>10 .or. np > 10) STOP 'ns > 10   .or.  np > 10  !'


   ! 1- setup the overlap matrix

    !  <s|s>
    DO jat=1,nat
      DO js=1,ns
        jorb=(jat-1)*(ns+3*np)+js
        aj=alpha(jat)/cs(js)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO is=1,ns
            !!iorb=iat+(is-1)*nat
            iorb=(iat-1)*(ns+3*np)+is
            ai= alpha(iat)/cs(is)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2
            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            ovrlp(iorb,jorb)=sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <pi|sj>
    DO jat=1,nat
      DO js=1,ns

        jorb=(jat-1)*(ns+3*np)+js
        aj=alpha(jat)/cs(js)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO ip=1,np
            !!iorb=1+(iat-1)*3+ns*nat + (ip-1)*3*nat
            iorb=(iat-1)*(ns+3*np)+ns+ip
            ai= alpha(iat)/cp(ip)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2

            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t3=-2.d0*sqrt(ai)*aj/t2
            ovrlp(iorb  ,jorb  )= t3 * xij *sij
            ovrlp(iorb+1,jorb  )= t3 * yij *sij
            ovrlp(iorb+2,jorb  )= t3 * zij *sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <si|pj>
    DO jat=1,nat
      DO jp=1,np

        jorb=(jat-1)*(ns+3*np)+ns+jp
        aj=alpha(jat)/cp(jp)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO is=1,ns
            !!iorb=iat+(is-1)*nat
            iorb=(iat-1)*(ns+3*np)+is
            ai= alpha(iat)/cs(is)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2

            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t3=+2.d0*sqrt(aj)*ai/t2
            ovrlp(iorb,jorb  )= t3 * xij *sij
            ovrlp(iorb,jorb+1)= t3 * yij *sij
            ovrlp(iorb,jorb+2)= t3 * zij *sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !  <p|p>
    DO jat=1,nat
      DO jp=1,np

        jorb=(jat-1)*(ns+3*np)+ns+jp
        aj=alpha(jat)/cp(jp)
        xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

        DO iat=1,nat
          DO ip=1,np
            iorb=(iat-1)*(ns+3*np)+ns+ip
            ai= alpha(iat)/cp(ip)
            xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

            xij=xi-xj; yij=yi-yj; zij=zi-zj
            r2=xij**2 + yij**2 + zij**2
            t1=ai*aj
            t2=ai+aj

            ! normalized GTOs:
            sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
            t4= 2.d0*sqrt(t1)/t2
            t5=-2.d0*t1/t2

            ovrlp(iorb  ,jorb  )= t4 *(1.d0 + t5* xij* xij)  * sij
            ovrlp(iorb+1,jorb  )= t4 *(       t5* yij* xij)  * sij
            ovrlp(iorb+2,jorb  )= t4 *(       t5* zij* xij)  * sij
            ovrlp(iorb  ,jorb+1)= t4 *(       t5* xij* yij)  * sij
            ovrlp(iorb+1,jorb+1)= t4 *(1.d0+  t5* yij* yij)  * sij
            ovrlp(iorb+2,jorb+1)= t4 *(       t5* zij* yij)  * sij
            ovrlp(iorb  ,jorb+2)= t4 *(       t5* xij* zij)  * sij
            ovrlp(iorb+1,jorb+2)= t4 *(       t5* yij* zij)  * sij
            ovrlp(iorb+2,jorb+2)= t4 *(1.d0+  t5* zij* zij)  * sij

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  end subroutine crtovrlp

  !---------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !> Multiplying the cutoff function to the overlap matrix (OM)
  !> @param[in] ovrlp overlap matrix
  !> @param[in] amplitude cutoff function values
  !> @param[in] norb number of orbitals ((ns+3*np)*nat_sphere)
  !> @param[in] ns if 0 NO s-orbitals; if 1 s-orbitals
  !> @param[in] np if 0 NO p-orbitals; if 1 px, py, pz -orbitals
  !> @param[out] ovrla overlap matrix with including cutoff function
  !---------------------------------------------------------------------------
  subroutine multamp(nat,ovrlp,amplitude,norb,ns,np,ovrla)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: norb
    INTEGER, INTENT(IN) :: ns
    INTEGER, INTENT(IN) :: np
    INTEGER :: i
    INTEGER :: j
    INTEGER :: iat
    INTEGER :: jat
    INTEGER :: iorb
    INTEGER :: jorb
    REAL(8), DIMENSION(nat), INTENT(IN) :: amplitude
    REAL(8), DIMENSION(norb, norb), INTENT(IN) :: ovrlp
    REAL(8), DIMENSION(norb, norb), INTENT(OUT) :: ovrla

    DO jat = 1,nat
      DO j = 1,(ns+3*np)
        jorb = (jat -1)*(ns+3*np) + j
        DO iat = 1,nat
          DO i = 1,(ns+3*np)
            iorb = (iat-1)*(ns+3*np) +i
            ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  end subroutine multamp



  !---------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !> Diagonalization of a matrix
  !> @param[in] n size of the n by n matrix
  !> @param[inout] aa n by n matrix / output are the eigenvectors
  !> @param[out] eval eigenvalues of the n by n matrix
  !---------------------------------------------------------------------------

  subroutine diagonalizeMatrix(n, aa, eval)
    IMPLICIT NONE
    ! Calling arguments
    INTEGER, INTENT(IN) :: n
    REAL(8),DIMENSION(n,n), INTENT(INOUT) :: aa
    REAL(8),DIMENSION(n), INTENT(IN) :: eval
    ! Local variables
    INTEGER :: info
    INTEGER :: lwork
    REAL(8),DIMENSION(:),ALLOCATABLE:: work

    lwork=100*n
    ALLOCATE(work(lwork))
    CALL dsyev('v','l', n, aa, n, eval, work, lwork, info)
    IF(info/=0) STOP ' ERROR in dsyev'
    DEALLOCATE(work)

  end subroutine diagonalizeMatrix

  !---------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !> Converting symbols to covalent radii
  !> @param[in] sym atomic symbol
  !> @param[out] rcov covalent radius
  !---------------------------------------------------------------------------

  subroutine sym2rcov_om_fp(sym,rcov)
  ! returns the covalent radius of atom with chemical symbol sym
    REAL(8), INTENT(OUT)  :: rcov
    CHARACTER(len=2), INTENT(IN) :: sym  ! chemical symbol
    SELECT CASE (adjustl(trim(sym)))
       ! covalet radius in Angstrom taken from WebElements: http://www.webelements.com/periodicity/covalent_radius/
    CASE('H')
       rcov= 0.37d0
    CASE('He')
       rcov= 0.32d0
    CASE('Li')
       rcov= 1.34d0
    CASE('Be')
       rcov= 0.90d0
    CASE('B')
       rcov= 0.82d0
    CASE('C')
       rcov= 0.77d0
    CASE('N')
       rcov= 0.75d0
    CASE('O')
       rcov= 0.73d0
    CASE('F')
       rcov= 0.71d0
    CASE('Ne')
       rcov= 0.69d0
    CASE('Na')
       rcov= 1.54d0
    CASE('Mg')
       rcov= 1.30d0
    CASE('Al')
       rcov= 1.18d0
    CASE('Si')
       rcov= 1.11d0
    CASE('P')
       rcov= 1.06d0
    CASE('S')
       rcov= 1.02d0
    CASE('Cl')
       rcov= 0.99d0
    CASE('Ar')
       rcov= 0.97d0
    CASE('K')
       rcov= 1.96d0
    CASE('Ca')
       rcov= 1.74d0
    CASE('Sc')
       rcov= 1.44d0
    CASE('Ti')
       rcov= 1.36d0
    CASE('V')
       rcov= 1.25d0
    CASE('Cr')
       rcov= 1.27d0
    CASE('Mn')
       rcov= 1.39d0
    CASE('Fe')
       rcov= 1.25d0
    CASE('Co')
       rcov= 1.26d0
    CASE('Ni')
       rcov= 1.21d0
    CASE('Cu')
       rcov= 1.38d0
    CASE('Zn')
       rcov= 1.31d0
    CASE('Ga')
       rcov= 1.26d0
    CASE('Ge')
       rcov= 1.22d0
    CASE('As')
       rcov= 1.19d0
    CASE('Se')
       rcov= 1.16d0
    CASE('Br')
       rcov= 1.14d0
    CASE('Kr')
       rcov= 1.10d0
    CASE('Rb')
       rcov= 2.11d0
    CASE('Sr')
       rcov= 1.92d0
    CASE('Y')
       rcov= 1.62d0
    CASE('Zr')
       rcov= 1.48d0
    CASE('Nb')
       rcov= 1.37d0
    CASE('Mo')
       rcov= 1.45d0
    CASE('Tc')
       rcov= 1.56d0
    CASE('Ru')
       rcov= 1.26d0
    CASE('Rh')
       rcov= 1.35d0
    CASE('Pd')
       rcov= 1.31d0
    CASE('Ag')
       rcov= 1.53d0
    CASE('Cd')
       rcov= 1.48d0
    CASE('In')
       rcov= 1.44d0
    CASE('Sn')
       rcov= 1.41d0
    CASE('Sb')
       rcov= 1.38d0
    CASE('Te')
       rcov= 1.35d0
    CASE('I')
       rcov= 1.33d0
    CASE('Xe')
       rcov= 1.30d0
    CASE('Cs')
       rcov= 2.25d0
    CASE('Ba')
       rcov= 1.98d0
    CASE('La')
       rcov= 1.69d0
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
    CASE('Lu')
       rcov= 1.60d0
    CASE('Hf')
       rcov= 1.50d0
    CASE('Ta')
       rcov= 1.38d0
    CASE('W')
       rcov= 1.46d0
    CASE('Re')
       rcov= 1.59d0
    CASE('Os')
       rcov= 1.28d0
    CASE('Ir')
       rcov= 1.37d0
    CASE('Pt')
       rcov= 1.28d0
    CASE('Au')
       rcov= 1.44d0
    CASE('Hg')
       rcov= 1.49d0
    CASE('Tl')
       rcov= 1.48d0
    CASE('Pb')
       rcov= 1.47d0
    CASE('Bi')
       rcov= 1.46d0
       !     case('Po')
       !     case('At')
    CASE('Rn')
       rcov= 1.45d0
    CASE('LA')   ! Lennard Jones atom
       rcov= 1.122462048309373d0
    CASE('LB')  ! Lennard Jones atom
       rcov= 0.9877666025122482d0
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
    CASE DEFAULT
       PRINT*, " Not recognized atomic type "//sym ; STOP
    ENDSELECT

    rcov = rcov /  0.52917720859d0   ! convert to atomic units
    !rcov = 1.25d0 * rcov

  endsubroutine sym2rcov_om_fp
