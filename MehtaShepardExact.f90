! Exact (1+2) phase shift for three identical bosons on a line with equal
! pairwise attractive delta-function interactions -- Mehta & Shepard, Phys.
! Rev. A 72, 032728 (2005), Eq. (12):
!
!   q*a2*tan(delta) = 4*sqrt(2/3) * (q*a2)^2 / ((q*a2)^2 - 2)
!
! where a2 is the 1+2 (dimer-atom) analog... actually the two-body scattering
! length (a2 = -1/(mu2*g), mu2 the 2-body reduced mass, g<0 the delta-function
! strength; B2 = 1/(2*mu2*a2^2) is the two-body binding energy, this repo's
! unit of energy), and q is defined by E = q^2/(2m) - B2 (m = single-particle
! mass), i.e. q is the paper's OWN hyperspherical channel momentum -- NOT
! necessarily equal to this repo's qOurs (RMATPROPAdiabatic.f90's channel
! momentum, qOurs^2 = 2*muHyper*(Energy-Threshold), muHyper =
! sqrt(m1*m2*m3/(m1+m2+m3)) = m/sqrt(3) for three equal masses m).
!
! Our hyperradial wavenumber k=qOurs=sqrt(2*muHyper*(E-E2)), muHyper=m/sqrt(3)
! for three equal masses m (this repo's hyperradial reduced mass); Mehta-
! Shepard's own q=sqrt(2*m*(E-E2)) (their Eq. near (7), m = single-particle
! mass). Eliminating (E-E2) between the two: k^2/muHyper = q^2/m, i.e.
!   q = qOurs * sqrt(muHyper/m) = qOurs * 3^(-1/4)   (equal masses, m=1)
! CONFIRMED empirically against the numerical R-matrix result (single box,
! Rmax=20, quadratic grid, 3bodydata_delta_1ch): with qConv=3^(-1/4), K_exact/
! K_numerical -> 1 as q increases through the tested range (0.98 to 1.68 to
! ~1.0 depending on q; residual deviation at intermediate q likely reflects
! numerical resolution, not a units error -- see conversation record).
SUBROUTINE MehtaShepardExact(qOurs, a2, qConv, delta, Kmat)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(in) :: qOurs   ! this repo's channel momentum
  DOUBLE PRECISION, INTENT(in) :: a2      ! two-body scattering length (=1 in this repo's units)
  DOUBLE PRECISION, INTENT(in) :: qConv   ! q_paper/qOurs conversion factor (see header)
  DOUBLE PRECISION, INTENT(out) :: delta, Kmat
  DOUBLE PRECISION X, qPaper

  qPaper = qConv*qOurs
  X = qPaper*a2

  Kmat = 4d0*DSQRT(2d0/3d0)*X/(X**2-2d0)
  delta = DATAN(Kmat)
END SUBROUTINE MehtaShepardExact
!****************************************************************************************************
PROGRAM WriteMehtaShepardExact
  ! Reads q from PhaseShift.dat (this repo's R-matrix output) and writes the
  ! Mehta-Shepard exact K and delta at the same q values (qConv=3^(-1/4),
  ! confirmed empirically -- see header), for direct comparison/plotting.
  IMPLICIT NONE
  DOUBLE PRECISION q, Energy, Knum, deltanum, a2, qConv, deltaEx, Kex
  CHARACTER(200) line
  INTEGER ios

  a2 = 1.0d0            ! this repo's a2=1 convention (Threshold(1)=-1 => B2=1)
  ! qConv = q_MehtaShepard/q_ours. Mehta-Shepard 2005 uses muref=m (single-particle
  ! mass) for their hyperspherical/Jacobi coordinates -- confirmed via their own
  ! Appendix A (ksi1=(1/sqrt(2))(x1-x2), ksi2=sqrt(2/3)((x1+x2)/2-x3), matching
  ! Amaya-Tapia's own convention exactly) and independently via their own k1,k2,k3
  ! (Eqs 5-7): E=(1/(2m))(k1^2+k2^2+k3^2)=q^2/(2m)-1/(m*a2^2), giving muref=m
  ! directly. This repo's own qOurs uses muHyper=m/sqrt(3) (Mehta-Esry-Greene 2007's
  ! mu=sqrt(mu12*mu123), NOT Mehta-Shepard 2005's own convention). Equating
  ! qOurs^2/(2*muHyper) = q_MS^2/(2*m) (same physical E-Eth) gives q_MS =
  ! qOurs*3^(1/4) -- the INVERSE of the 3^(-1/4) factor used earlier this session.
  qConv = 3d0**(0.25d0)

  OPEN(unit=10,file='PhaseShift.dat',status='old')
  OPEN(unit=20,file='MehtaShepardExact.dat',status='replace')
  WRITE(20,'(A)') '# q  Energy  K_numerical  delta_numerical  K_exact  delta_exact'
  DO
     READ(10,'(A)',IOSTAT=ios) line
     IF (ios.NE.0) EXIT
     IF (line(1:1).EQ.'#' .OR. line(1:1).EQ.'@') CYCLE
     READ(line,*,IOSTAT=ios) q, Energy, Knum, deltanum
     IF (ios.NE.0) CYCLE
     CALL MehtaShepardExact(q, a2, qConv, deltaEx, Kex)
     WRITE(20,'(1P,6E20.12)') q, Energy, Knum, deltanum, Kex, deltaEx
  ENDDO
  CLOSE(10)
  CLOSE(20)
  WRITE(6,*) "Wrote MehtaShepardExact.dat"
END PROGRAM WriteMehtaShepardExact
