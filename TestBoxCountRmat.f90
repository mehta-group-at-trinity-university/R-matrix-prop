!****************************************************************************************************
! Sanity check requested directly: does the final propagated R-matrix (Rmat=Zf^2/bf at R=xEnd)
! depend on how many boxes [0,xEnd] is split into? It should NOT -- BoxMatch chaining is supposed
! to be exact (up to basis truncation within each box), so 1/2/3-box configurations spanning the
! SAME total range at the SAME energy should agree closely. Tests the adiabatic (B-spline) path
! first (believed working), then the SVD/direct path (LSVD held at 120 for every box in every
! configuration -- deliberately over-resolved even for the widest single-box case, so this isolates
! box-COUNT/chaining consistency from the resolution question already settled separately).
!****************************************************************************************************
PROGRAM TestBoxCountRmat
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
  DOUBLE PRECISION muLocal, alpha, ddim, EnergyLocal, xMinLocal, xEndLocal, xStartLocal
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  INTEGER NumDataPoints, OrderLocal, LegPointsLocal, xNumPoints
  DOUBLE PRECISION RmatOut

  DataDir = '3bodydata_sech2_3boson_ch1'
  OrderLocal = 5
  xNumPoints = 40
  LegPointsLocal = 10
  xMinLocal = 1.0d-4
  xStartLocal = 0.0d0
  xEndLocal = 40.0d0
  EnergyLocal = 0.5d0

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = LegPointsLocal
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  CALL ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alpha,ddim, &
       Threshold,Leff,UInterp,PInterp,QInterp)
  reducedmass = muLocal
  AlphaFactor = 0d0
  EffDim = ddim

  WRITE(6,*) "=== Path A: adiabatic (B-spline) box-count sanity check ==="
  CALL RunAdiabaticChain(1,RmatOut); WRITE(6,*) "NBoxes=1  Rmat = ",RmatOut
  CALL RunAdiabaticChain(2,RmatOut); WRITE(6,*) "NBoxes=2  Rmat = ",RmatOut
  CALL RunAdiabaticChain(3,RmatOut); WRITE(6,*) "NBoxes=3  Rmat = ",RmatOut

  WRITE(6,*)
  WRITE(6,*) "=== Path B: SVD/direct box-count sanity check (LSVD=120 per box, every config) ==="
  CALL RunSVDChain(1,RmatOut); WRITE(6,*) "NBoxes=1  Rmat = ",RmatOut
  CALL RunSVDChain(2,RmatOut); WRITE(6,*) "NBoxes=2  Rmat = ",RmatOut
  CALL RunSVDChain(3,RmatOut); WRITE(6,*) "NBoxes=3  Rmat = ",RmatOut

CONTAINS

  SUBROUTINE RunAdiabaticChain(NBoxesTest,RmatFinal)
    INTEGER, INTENT(in) :: NBoxesTest
    DOUBLE PRECISION, INTENT(out) :: RmatFinal
    TYPE(BPData) BPD
    TYPE(GenEigVal) EIG
    TYPE(BoxData), ALLOCATABLE :: Boxes(:)
    TYPE(BoxData) Bnull
    DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
    DOUBLE PRECISION, ALLOCATABLE :: CoulombC(:)
    DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
    INTEGER iBox
    DOUBLE PRECISION ZfFinal, ZfpFinal, bfFinal

    ALLOCATE(CoulombC(NumChannels))
    CoulombC = 0d0
    ALLOCATE(Boxes(NBoxesTest))
    Bnull%NumOpenL = 0
    Bnull%NumOpenR = 0
    CALL AllocateBox(Bnull)
    CALL InitZeroBox(Bnull)

    DO iBox = 1,NBoxesTest
       BPD%NumChannels = NumChannels
       BPD%Order = OrderLocal
       BPD%Left = 2
       BPD%Right = 2
       BPD%xNumPoints = xNumPoints
       BPD%kl = 1d20
       BPD%kr = 1d20
       CALL AllocateBPD(BPD)
       BPD%xl = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox-1)/DBLE(NBoxesTest))**BoxSpacingPower
       BPD%xr = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox)/DBLE(NBoxesTest))**BoxSpacingPower
       IF (iBox.EQ.1 .AND. BPD%xl.EQ.0d0) BPD%xl = xMinLocal
       CALL GridMaker(BPD%xPoints,BPD%xNumPoints,BPD%xl,BPD%xr,"quadratic")
       CALL MakeBasis(BPD)
       BPD%lam(1:NumChannels) = Leff(1:NumChannels)
       CALL SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp,CoulombC)

       EIG%MatrixDim = BPD%MatrixDim
       CALL AllocateEIG(EIG)
       Boxes(iBox)%NumOpenL = MERGE(0,1,iBox.EQ.1)
       Boxes(iBox)%NumOpenR = 1
       CALL CalcGamLam(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
       CALL CombineGam(EIG,EnergyLocal)
       CALL AllocateBox(Boxes(iBox))
       CALL InitZeroBox(Boxes(iBox))

       ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
       ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
       CALL PartitionAndEliminate(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR,evalRed,evecRed,LamooDiag)

       IF (iBox.EQ.1) THEN
          CALL BoxMatch(Bnull, Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
       ELSE
          CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
       ENDIF

       DEALLOCATE(evalRed,evecRed,LamooDiag)
       CALL DeAllocateEIG(EIG)
       CALL DeAllocateBPD(BPD)
    ENDDO

    ZfFinal = Boxes(NBoxesTest)%Zf(1,1)
    ZfpFinal = Boxes(NBoxesTest)%Zfp(1,1)
    bfFinal = -ZfpFinal/ZfFinal
    RmatFinal = ZfFinal**2/bfFinal

    DO iBox = 1,NBoxesTest
       CALL DeAllocateBox(Boxes(iBox))
    ENDDO
    CALL DeAllocateBox(Bnull)
    DEALLOCATE(Boxes,CoulombC)
  END SUBROUTINE RunAdiabaticChain

  SUBROUTINE RunSVDChain(NBoxesTest,RmatFinal)
    INTEGER, INTENT(in) :: NBoxesTest
    DOUBLE PRECISION, INTENT(out) :: RmatFinal
    TYPE(BPData) BPDSVD
    TYPE(GenEigVal) EIGSVD
    TYPE(BoxData), ALLOCATABLE :: BoxesB(:)
    TYPE(BoxData) BnullB
    DOUBLE PRECISION, ALLOCATABLE :: evalRedB(:), evecRedB(:,:), LamooDiagB(:)
    DOUBLE PRECISION, ALLOCATABLE :: ODiagL(:), ODiagR(:)
    DOUBLE PRECISION :: wsq1Out, wsqLOut, ShiftLocal, a1, a2
    CHARACTER*64 :: LegendreFile64, datadirLoc
    CHARACTER(LEN=16) :: GType
    DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
    INTEGER, PARAMETER :: LSVD = 120
    INTEGER iBox
    DOUBLE PRECISION ZfFinal, ZfpFinal, bfFinal

    LegendreFile64 = 'Legendre.dat'
    datadirLoc = '.'
    ShiftLocal = -5.0d0
    ALLOCATE(ODiagL(1),ODiagR(1))
    ALLOCATE(BoxesB(NBoxesTest))
    BnullB%NumOpenL = 0
    BnullB%NumOpenR = 0
    CALL AllocateBox(BnullB)
    CALL InitZeroBox(BnullB)

    DO iBox = 1,NBoxesTest
       a1 = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox-1)/DBLE(NBoxesTest))**BoxSpacingPower
       a2 = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox)/DBLE(NBoxesTest))**BoxSpacingPower
       GType = MERGE('Radau   ','Lobatto ',iBox.EQ.1)

       IF (ALLOCATED(BPDSVD%lam)) DEALLOCATE(BPDSVD%lam)
       CALL BuildSVDBoxGeneric(a1,a2,LSVD,1,5,1,1, &
            200,0d0,dacos(-1d0)/6.0d0,LegendreFile64,10,ShiftLocal,datadirLoc, &
            .FALSE.,SechPotential,SechGridMaker,SechRobinCoeff, &
            MERGE(0,1,iBox.EQ.1),1,TRIM(GType),BPDSVD,EIGSVD,ODiagL,ODiagR,wsq1Out,wsqLOut)

       CALL CombineGam(EIGSVD,EnergyLocal)

       BoxesB(iBox)%NumOpenL = MERGE(0,1,iBox.EQ.1)
       BoxesB(iBox)%NumOpenR = 1
       CALL AllocateBox(BoxesB(iBox))
       CALL InitZeroBox(BoxesB(iBox))

       ALLOCATE(evalRedB(BoxesB(iBox)%betaMax),evecRedB(BoxesB(iBox)%betaMax,BoxesB(iBox)%betaMax))
       ALLOCATE(LamooDiagB(BoxesB(iBox)%betaMax))
       CALL PartitionAndEliminate(BPDSVD,EIGSVD,BoxesB(iBox)%NumOpenL,BoxesB(iBox)%NumOpenR, &
            evalRedB,evecRedB,LamooDiagB)

       IF (iBox.EQ.1) THEN
          CALL BoxMatch(BnullB, BoxesB(iBox), BPDSVD, evalRedB, evecRedB, LamooDiagB, EffDim, AlphaFactor)
       ELSE
          CALL BoxMatch(BoxesB(iBox-1), BoxesB(iBox), BPDSVD, evalRedB, evecRedB, LamooDiagB, EffDim, AlphaFactor)
       ENDIF

       DEALLOCATE(evalRedB,evecRedB,LamooDiagB)
       CALL DeAllocateEIG(EIGSVD)
    ENDDO

    ZfFinal = BoxesB(NBoxesTest)%Zf(1,1)
    ZfpFinal = BoxesB(NBoxesTest)%Zfp(1,1)
    bfFinal = -ZfpFinal/ZfFinal
    RmatFinal = ZfFinal**2/bfFinal

    DO iBox = 1,NBoxesTest
       CALL DeAllocateBox(BoxesB(iBox))
    ENDDO
    CALL DeAllocateBox(BnullB)
    DEALLOCATE(BoxesB,ODiagL,ODiagR)
  END SUBROUTINE RunSVDChain

END PROGRAM TestBoxCountRmat
