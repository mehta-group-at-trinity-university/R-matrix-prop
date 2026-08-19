!****************************************************************************************************
! Concrete AdiabaticInterfaces.f90 plugins for the equal-mass three-boson sech^2-well problem used
! to generate 3bodydata_sech2_5ch/3bodydata_sech2_5ch_fine (Adiabatic-Scattering-BoundStates'
! DriveAdiabaticSolver1D.x, physics from DriveAdiabaticSolver1D_5ch.par: m1=m2=m3=1, V2Depth=2.0,
! PotRange=1.0, alpha=1.0, Order/Left/Right=5/1/0). Used to independently re-diagonalize the SAME
! physics near the origin (via GenericSVDChannelBasis.f90's ARPACK path) instead of trusting the
! tabulated/interpolated Uad.dat/Pmat.dat/Qmat.dat baseline (RMATPROPAdiabatic.x) there -- a
! genuinely independent numerical cross-check, mirroring what DeltaFunctionPlugins.f90 already did
! for the delta-function benchmark.
!
! Unlike the delta-function problem, this potential is NOT a boundary condition (no Robin jump) --
! it is an ordinary smooth (finite) potential evaluated at every angular quadrature point, and its
! R-dependence enters only through the pairwise distances r12/r13/r23 = f(R,phi), exactly as
! Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90's own CalcHamiltonian does. SechGridMaker
! below is an exact replica of that repo's GridMakerHHL (r0=PotRange=1.0) so the angular grid used
! here matches the one used to generate the tabulated baseline -- a genuinely apples-to-apples
! comparison, not just "some reasonable grid, some reasonable answer".
!****************************************************************************************************
MODULE SechFunctionPlugins
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: SechV2Depth = 2.0d0
  DOUBLE PRECISION, PARAMETER :: SechPotRange = 1.0d0
  DOUBLE PRECISION, PARAMETER :: SechAlpha = 1.0d0
CONTAINS

  ! Two-body sech^2 well and its r-derivative (Adiabatic-Scattering-BoundStates/lib/Potential.f90's
  ! V2body/DV2body, inlined here rather than linked -- this repo doesn't otherwise depend on that
  ! file, and the form is a two-line closed formula).
  DOUBLE PRECISION FUNCTION SechV2body(r)
    DOUBLE PRECISION, INTENT(in) :: r
    SechV2body = -SechV2Depth/DCOSH(r/SechPotRange)**2
  END FUNCTION SechV2body

  DOUBLE PRECISION FUNCTION SechDV2body(r)
    DOUBLE PRECISION, INTENT(in) :: r
    SechDV2body = (2d0*SechV2Depth/SechPotRange)*DTANH(r/SechPotRange)/DCOSH(r/SechPotRange)**2
  END FUNCTION SechDV2body

  !--- (1) PotentialProc: pairwise-sum sech^2 potential at every (LegPoints x xNumPoints-1)
  !    quadrature point, and its R-derivative. r12/r13/r23 = sqrt(2/sqrt(3))*R*|cos(phi-offset)|
  !    for offsets {0,-pi/3,+pi/3} -- equal-mass 1D 3-body coordinates, matching
  !    adiabaticSolver1D.f90 lines 405-427 exactly (phi=x0 here is already the mapped quadrature
  !    point, per AdiabaticInterfaces.f90's own header). dPot uses sumpairwisepotDR3's
  !    (alpha/R)*sum(r*DV2body(r)) form -- note R=0 is never actually reached here regardless of
  !    Box1GridType ('Radau' or 'LobattoDirichlet' both structurally avoid evaluating a1, see
  !    SVDChannelBasis.f90's own BuildSVDBox header), so the 1/R factor is never a problem.
  SUBROUTINE SechPotential(R,LegPoints,xNumPoints,x0,Pot,dPot)
    INTEGER, INTENT(in) :: LegPoints,xNumPoints
    DOUBLE PRECISION, INTENT(in) :: R
    DOUBLE PRECISION, INTENT(in) :: x0(LegPoints,xNumPoints)
    DOUBLE PRECISION, INTENT(out) :: Pot(LegPoints,xNumPoints), dPot(LegPoints,xNumPoints)
    DOUBLE PRECISION, PARAMETER :: RPair = 1.0745699318235419d0  ! sqrt(2/sqrt(3)) -- NOT 2/sqrt(3)
    DOUBLE PRECISION, PARAMETER :: PiThird = 1.0471975511965976d0 ! pi/3
    DOUBLE PRECISION :: r12,r13,r23,phi
    INTEGER :: lx,kx
    Pot = 0d0
    dPot = 0d0
    DO kx = 1,xNumPoints-1
       DO lx = 1,LegPoints
          phi = x0(lx,kx)
          r12 = RPair*R*DABS(DCOS(phi))
          r13 = RPair*R*DABS(DCOS(phi-PiThird))
          r23 = RPair*R*DABS(DCOS(phi+PiThird))
          Pot(lx,kx) = SechAlpha*(SechV2body(r12)+SechV2body(r13)+SechV2body(r23))
          dPot(lx,kx) = (SechAlpha/R)*(r12*SechDV2body(r12)+r13*SechDV2body(r13)+r23*SechDV2body(r23))
       ENDDO
    ENDDO
  END SUBROUTINE SechPotential

  !--- (2) GridMakerProcI: exact replica of Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90's
  !    GridMakerHHL (r0=PotRange=1.0), so the angular grid used to re-diagonalize here matches the
  !    one used to generate the tabulated 3bodydata_sech2_5ch/_fine baseline. Genuinely R-dependent
  !    (deltax=r0/R narrows with R) -- caller must set GridRebuildEveryR=.TRUE., unlike
  !    DeltaWedgeGridMaker's fixed grid.
  SUBROUTINE SechGridMaker(R,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
    INTEGER, INTENT(in) :: xNumPoints
    INTEGER, INTENT(out) :: CalcNewBasisFunc
    DOUBLE PRECISION, INTENT(in) :: R,xMin,xMax
    DOUBLE PRECISION, INTENT(out) :: xPoints(xNumPoints)
    DOUBLE PRECISION, PARAMETER :: Pi = 3.1415926535897932385d0
    DOUBLE PRECISION :: phi23,deltax,xRswitch,x0,x1,x2,x3,x4,xDelt
    INTEGER :: i,k
    phi23 = Pi/6.0d0
    deltax = SechPotRange/R
    xRswitch = 10.0d0*SechPotRange/Pi
    IF (R.GT.xRswitch .AND. xMax.LE.phi23) THEN
       ! Symmetry-reduced domain (e.g. [0,pi/6] for 3 identical bosons) whose right edge
       ! IS the r23=0 coalescence itself, not strictly interior to it -- the ordinary
       ! 4-segment scheme below assumes phi23<xMax strictly, and produces a REVERSED
       ! (negative-width) segment otherwise (confirmed: this exact case crashed
       ! ARPACK/LAPACK -- "DLASCL parameter 4 illegal value" -- when first tried against
       ! this domain). Same fix as adiabaticSolver1D.f90's GridMakerHHL (both repo
       ! copies): fall back to a plain coarse + refined-near-xMax 2-segment split.
       x0 = xMin
       x2 = xMax
       x1 = xMax - deltax
       IF (x1.LE.x0) x1 = x0 + (x2-x0)*0.5d0
       k = 1
       xDelt = (x1-x0)/DBLE(xNumPoints/2)
       DO i = 1,xNumPoints/2
          xPoints(k) = x0 + (i-1)*xDelt
          k = k+1
       ENDDO
       xDelt = (x2-x1)/DBLE(xNumPoints-xNumPoints/2-1)
       DO i = 1,xNumPoints-xNumPoints/2
          xPoints(k) = x1 + (i-1)*xDelt
          k = k+1
       ENDDO
    ELSE IF (R.GT.xRswitch) THEN
       x0 = xMin
       x1 = phi23 - deltax
       x2 = phi23 + deltax
       x3 = xMax  - deltax
       x4 = xMax
       IF (x1.LE.x0) x1 = x0 + (phi23-x0)*0.5d0
       IF (x2.GE.x3) x2 = x3 - (x3-phi23)*0.5d0
       k = 1
       xDelt = (x1-x0)/DBLE(xNumPoints/4)
       DO i = 1,xNumPoints/4
          xPoints(k) = x0 + (i-1)*xDelt
          k = k+1
       ENDDO
       xDelt = (x2-x1)/DBLE(xNumPoints/4)
       DO i = 1,xNumPoints/4
          xPoints(k) = x1 + (i-1)*xDelt
          k = k+1
       ENDDO
       xDelt = (x3-x2)/DBLE(xNumPoints/4)
       DO i = 1,xNumPoints/4
          xPoints(k) = x2 + (i-1)*xDelt
          k = k+1
       ENDDO
       xDelt = (x4-x3)/DBLE(xNumPoints/4-1)
       DO i = 1,xNumPoints/4
          xPoints(k) = x3 + (i-1)*xDelt
          k = k+1
       ENDDO
    ELSE
       x0 = xMin
       xDelt = (xMax-x0)/DBLE(xNumPoints-1)
       DO i = 1,xNumPoints
          xPoints(i) = x0 + (i-1)*xDelt
       ENDDO
    ENDIF
    CalcNewBasisFunc = 1
  END SUBROUTINE SechGridMaker

  ! (3) RobinCoeffProc: Left=1 (Neumann at phi=0)/Right=0 (Dirichlet at phi=pi/2) are both ordinary,
  ! non-Robin BC codes -- placeholder sentinel values only, never read by CalcBasisFuncsBP's
  ! non-Robin branches (matches this codebase's existing "unused" convention, see
  ! AdiabaticInterfaces.f90's own header).
  SUBROUTINE SechRobinCoeff(R,kLeft,kRight)
    DOUBLE PRECISION, INTENT(in) :: R
    DOUBLE PRECISION, INTENT(out) :: kLeft,kRight
    kLeft = 1d20
    kRight = 1d20
  END SUBROUTINE SechRobinCoeff
END MODULE SechFunctionPlugins
