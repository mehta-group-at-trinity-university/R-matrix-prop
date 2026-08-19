!****************************************************************************************************
! Direct test of the user's factorization hypothesis: at E=-0.999 (deep below the 3-body breakup
! threshold), channels 2-5 are strongly closed at every R in this box (their thresholds are all
! near 0, so kappa_i*R >> 1 throughout), and CouplingFlag=0 means BuildSVDBoxGeneric never computes
! an explicit P/Q coupling term at all -- any channel-1-to-2..5 interaction can only come from the
! DVR's own channel-overlap matrix O. If the physics genuinely factors as F_1(R)*Phi_1(R,phi) (the
! excited channels contributing negligibly), then embedding channel 1 in a NumStates=5 SVD box
! construction should reproduce the SAME channel-1 R-matrix as building it ALONE (NumStates=1,
! already validated). If it does NOT, that isolates the bug to how BuildSVDBoxGeneric handles
! NumStates>1 specifically (its O/Gam0/Overlap assembly), not to genuine multichannel coupling.
!
! Single box spanning [0,xEnd] (no chaining), SechFunctionPlugins physics, Left1D=Right1D=1,
! GridType='Radau', AlphaFactor=0 -- otherwise identical setup to the already-validated
! CompareRmatrixSVDvsAdiabatic.f90, just at E=-0.999 and with both NumStates=1 and NumStates=5
! SVD box constructions directly compared against each other.
!****************************************************************************************************
PROGRAM CompareRmatrixMultichannel
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE AdiabaticPotential
  USE AdiabaticInterfaces
  USE SechFunctionPlugins
  USE GenericSVDChannelBasis, ONLY: BuildSVDBoxGeneric
  IMPLICIT NONE

  CHARACTER(LEN=64) DataDir
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim, EnergyLocal, xEndLocal
  INTEGER NumDataPoints
  DOUBLE PRECISION RmatOne, RmatFive
  DOUBLE PRECISION xEndBase, kChan1
  ! Sensitivity-hypothesis test (user request): if the 175%/17.7% Rmat disagreements are an
  ! artifact of xEnd landing near an antinode of u (bf=-u'/u accidentally near zero, since
  ! Rmat=Zf^2/bf diverges there), shifting xEnd by a quarter wavelength pi/(2k) of channel 1's
  ! asymptotic oscillation should move bf away from zero and shrink the disagreement.
  LOGICAL, PARAMETER :: UseShiftedXEnd = .TRUE.

  DataDir = '3bodydata_sech2_3boson'   ! full 5-channel dataset, for real Threshold/Leff(2:5)
  xEndBase = 40.0d0
  EnergyLocal = -0.999d0

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = 10
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  CALL ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alpha,ddim, &
       Threshold,Leff,UInterp,PInterp,QInterp)
  reducedmass = muLocal
  AlphaFactor = 0d0
  EffDim = ddim
  Threshold = 2d0*muLocal*Threshold   ! CalcK/production convention

  WRITE(6,*) "Thresholds (rescaled) = ",Threshold
  WRITE(6,*) "2*mu*E = ",2d0*muLocal*EnergyLocal

  ! Channel 1 is open: k1 = sqrt(2*mu*E - Threshold(1)) (both already in the "2*mu*" rescaled
  ! convention). Shift xEnd by a quarter wavelength pi/(2*k1) to move bf away from wherever it
  ! landed near zero at xEndBase.
  kChan1 = DSQRT(2d0*muLocal*EnergyLocal - Threshold(1))
  WRITE(6,*) "k_channel1 = ",kChan1," wavelength = ",2d0*Pi/kChan1
  IF (UseShiftedXEnd) THEN
     xEndLocal = xEndBase + Pi/(2d0*kChan1)
  ELSE
     xEndLocal = xEndBase
  ENDIF
  WRITE(6,*) "xEndBase = ",xEndBase," xEndLocal (used) = ",xEndLocal

  WRITE(6,*) "kappa_i*xEnd for channels 2..5 (should be >>1, deep closed):"
  BLOCK
    INTEGER ii
    DOUBLE PRECISION kappa
    DO ii = 2,NumChannels
       kappa = DSQRT(Threshold(ii)-2d0*muLocal*EnergyLocal)
       WRITE(6,*) "  channel ",ii," kappa=",kappa," kappa*xEnd=",kappa*xEndLocal
    ENDDO
  END BLOCK

  CALL RunSVDSingleBox(1,RmatOne)
  WRITE(6,*) "=== NumStates=1 (channel 1 alone) ==="
  WRITE(6,*) "Rmat(1,1) = ",RmatOne

  CALL RunSVDSingleBox(NumChannels,RmatFive)
  WRITE(6,*) "=== NumStates=",NumChannels," (channel 1 embedded, 2..",NumChannels," deep closed) ==="
  WRITE(6,*) "Rmat(1,1) = ",RmatFive

  WRITE(6,*)
  WRITE(6,*) "=== Comparison ==="
  WRITE(6,*) "Rmat(1,1), NumStates=1        = ",RmatOne
  WRITE(6,*) "Rmat(1,1), NumStates=",NumChannels," = ",RmatFive
  WRITE(6,*) "relative diff                  = ",DABS(RmatFive-RmatOne)/DABS(RmatOne)

CONTAINS

  SUBROUTINE RunSVDSingleBox(NStates,Rmat11)
    INTEGER, INTENT(in) :: NStates
    DOUBLE PRECISION, INTENT(out) :: Rmat11
    TYPE(BPData) :: BPDSVD
    TYPE(GenEigVal) :: EIGSVD
    TYPE(BoxData) :: Bnull, Box1
    DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
    DOUBLE PRECISION, ALLOCATABLE :: ODiagL(:), ODiagR(:)
    DOUBLE PRECISION :: wsq1Out, wsqLOut, ShiftLocal
    CHARACTER*64 :: LegendreFile64, datadirLoc
    INTEGER, PARAMETER :: LSVD = 24
    DOUBLE PRECISION :: Zf11, Zfp11, bf11

    LegendreFile64 = 'Legendre.dat'
    datadirLoc = '.'
    ShiftLocal = -5.0d0
    ALLOCATE(ODiagL(NStates),ODiagR(NStates))

    ! NumOpenR passed to BuildSVDBoxGeneric/Box1 is the count of PHYSICALLY open channels
    ! (=1 at E=-0.999, regardless of NStates) -- matching production's "last box" convention
    ! (Boxes(NumBoxes)%NumOpenR=NumOpenChannels, not NumChannels). Channels 2..NStates, even
    ! when constructed, must be eliminated at the boundary too (closed channels never carry
    ! retained boundary DOF) -- getting this wrong here (as an earlier version of this file
    ! did, using NumOpenR=NStates) would silently keep closed-channel DOF "open" instead of
    ! correctly eliminating them, which is NOT what the production driver does.
    CALL BuildSVDBoxGeneric(0.0d0,xEndLocal,LSVD,NStates,5,1,1, &
         200,0d0,dacos(-1d0)/6.0d0,LegendreFile64,10,ShiftLocal,datadirLoc, &
         .FALSE.,SechPotential,SechGridMaker,SechRobinCoeff, &
         0,1,'Radau',BPDSVD,EIGSVD,ODiagL,ODiagR,wsq1Out,wsqLOut)

    CALL CombineGam(EIGSVD,EnergyLocal)

    Box1%NumOpenL = 0
    Box1%NumOpenR = 1   ! only channel 1 is physically open at E=-0.999
    CALL AllocateBox(Box1)
    CALL InitZeroBox(Box1)

    ALLOCATE(evalRed(Box1%betaMax),evecRed(Box1%betaMax,Box1%betaMax),LamooDiag(Box1%betaMax))
    CALL PartitionAndEliminate(BPDSVD,EIGSVD,Box1%NumOpenL,Box1%NumOpenR,evalRed,evecRed,LamooDiag)

    Bnull%NumOpenL = 0
    Bnull%NumOpenR = 0
    CALL AllocateBox(Bnull)
    CALL InitZeroBox(Bnull)
    CALL BoxMatch(Bnull, Box1, BPDSVD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)

    Zf11 = Box1%Zf(1,1)
    Zfp11 = Box1%Zfp(1,1)
    bf11 = -Zfp11/Zf11
    Rmat11 = Zf11**2/bf11

    DEALLOCATE(evalRed,evecRed,LamooDiag,ODiagL,ODiagR)
    CALL DeAllocateBox(Box1)
    CALL DeAllocateBox(Bnull)
    CALL DeAllocateEIG(EIGSVD)
  END SUBROUTINE RunSVDSingleBox

END PROGRAM CompareRmatrixMultichannel
