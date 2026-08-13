!****************************************************************************************************
! Concrete AdiabaticInterfaces.f90 plugins for the equal-mass three-boson delta-function problem
! (Mehta, Esry & Greene, PRA 76, 022711 (2007)), the same system KMSFormulas.f/RobinSVDChannelBasis.f90
! already solve exactly/numerically. Used to validate OneDimChannelsGeneric/CalcHamiltonianGeneric
! (AdiabaticSolverGeneric.f90) against that already-validated benchmark -- a genuinely independent
! code path (ARPACK shift-invert Lanczos via MyDsband, not RobinSVDChannelBasis.f90's own dense
! DSYGV) for the SAME physics.
!
! Between phi=0 and phi=pi/6 the delta-function problem has NO potential at all (its entire effect
! is the Robin log-derivative jump condition at phi=pi/6) -- so DeltaZeroPotential is trivially zero
! everywhere, and DeltaWedgeGridMaker's grid never depends on R (only the Robin coefficient
! kRight=sqrt(2*mu)*R, from DeltaRobinCoeff, does) -- see AdiabaticInterfaces.f90's own header for
! why grid cadence and Robin-coefficient update are independent axes; this problem is the case that
! motivated separating them (fixed grid, R-dependent kRight).
!****************************************************************************************************
MODULE DeltaFunctionPlugins
  IMPLICIT NONE
CONTAINS
  SUBROUTINE DeltaZeroPotential(R,LegPoints,xNumPoints,x0,Pot,dPot)
    INTEGER, INTENT(in) :: LegPoints,xNumPoints
    DOUBLE PRECISION, INTENT(in) :: R
    DOUBLE PRECISION, INTENT(in) :: x0(LegPoints,xNumPoints)
    DOUBLE PRECISION, INTENT(out) :: Pot(LegPoints,xNumPoints), dPot(LegPoints,xNumPoints)
    Pot = 0d0
    dPot = 0d0
  END SUBROUTINE DeltaZeroPotential

  ! Fixed uniform grid over phi in [0,pi/6], matching RobinSVDChannelBasis.f90's own BPDPhi%xPoints
  ! exactly -- ignores R entirely (this problem's angular domain never moves), so this is only ever
  ! called once (GridRebuildEveryR=.FALSE. in the test driver).
  SUBROUTINE DeltaWedgeGridMaker(R,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
    USE GlobalVars, ONLY: Pi
    INTEGER, INTENT(in) :: xNumPoints
    INTEGER, INTENT(out) :: CalcNewBasisFunc
    DOUBLE PRECISION, INTENT(in) :: R,xMin,xMax
    DOUBLE PRECISION, INTENT(out) :: xPoints(xNumPoints)
    INTEGER :: iphi
    DO iphi = 1,xNumPoints
       xPoints(iphi) = (Pi/6.0d0)*DBLE(iphi-1)/DBLE(xNumPoints-1)
    ENDDO
    CalcNewBasisFunc = 1
  END SUBROUTINE DeltaWedgeGridMaker

  ! kLeft unused (Left=1, Neumann, matching Phi'(0)=0 -- see KMSFormulas.f's own header on the
  ! even-parity boundary condition); kRight=sqrt(2*mu)*R exactly, the delta-function's log-derivative
  ! jump condition (verified channel-independent to ~1e-13-1e-16 against SolveQ earlier this session).
  SUBROUTINE DeltaRobinCoeff(R,kLeft,kRight)
    USE GlobalVars, ONLY: reducedmass
    DOUBLE PRECISION, INTENT(in) :: R
    DOUBLE PRECISION, INTENT(out) :: kLeft,kRight
    kLeft = 1d20
    kRight = DSQRT(2d0*reducedmass)*R
  END SUBROUTINE DeltaRobinCoeff
END MODULE DeltaFunctionPlugins
