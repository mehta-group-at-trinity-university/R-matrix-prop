!****************************************************************************************************
! Standalone A/B comparator (Validation Playbook technique, doc/UserGuide.tex): single channel,
! single box, two independent code paths for the SAME physics (equal-mass 3-identical-boson sech^2
! problem, ground/lowest adiabatic curve), diffing the propagated R-matrix AND the reconstructed
! wavefunction u(R) directly -- isolates box-construction/propagation from CalcK's asymptotic
! Bessel-matching step, and from the 50-box/xNumPointsPhi=200 cost of the full comparison.
!
! Path A (adiabatic/tabulated): exactly TestAlphaZero.f90's already-validated recipe --
! 3bodydata_sech2_3boson_ch1 (tabulated Uad/Pmat/Qmat), Left=2/NumOpenL=0, AlphaFactor=0.
! Wavefunction reconstruction via PartitionAndEliminateFull + BPD1%u basis-function sum, same
! pattern as TestAlphaZero.f90/PlotWavefunctionFree.f90.
!
! Path B (SVD/direct): BuildSVDBoxGeneric with SechFunctionPlugins' SechPotential/SechGridMaker
! (re-diagonalizes the true sech^2 potential via ARPACK at every DVR node), domain [0,pi/6],
! Left1D=Right1D=1 (Neumann-Neumann, 3-identical-boson even parity), GridType='Radau'
! (AlphaFactor=0's SVD-side mechanism), NumStates=1 (just the lowest/ground curve, matching
! Path A's single-channel dataset). Wavefunction reconstruction via PartitionAndEliminateFull +
! the DVR Lagrange-basis Kronecker-delta property u(R_j)=cfull_j/sqrt(w_j) -- same pattern as
! PlotWavefunctionSVD.f90 (that file's GridType='LobattoDirichlet'; here GridType='Radau', so
! QOffset=0 and LQuad=L exactly, no "+1" node-count offset needed).
!
! Both paths use a SINGLE box spanning [0,xEnd] (no box-chaining at all), so any disagreement is
! attributable to box CONSTRUCTION (the physics provider), not to BoxMatch/box-to-box propagation.
! Compares the physical R-matrix Rmat=Zf^2/bf (bf=-Zfp/Zf, this repo's BoxData%Zf/%Zfp/%bf
! convention -- see memory wavefn-reduction-term-bug), not raw Zfp/Zf, since that's the
! convention-independent quantity (RmatPropCore.f90's own CalcK does the same conversion).
!****************************************************************************************************
PROGRAM CompareRmatrixSVDvsAdiabatic
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE AdiabaticPotential
  USE AdiabaticInterfaces
  USE SechFunctionPlugins
  USE GenericSVDChannelBasis, ONLY: BuildSVDBoxGeneric
  USE SVDChannelBasis, ONLY: GetGaussRadauFactors
  IMPLICIT NONE

  CHARACTER(LEN=64) DataDir
  TYPE(BPData) BPD1
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull, Box1
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim, EnergyLocal, xMinLocal, xEndLocal
  INTEGER NumDataPoints, xNumPoints, LegPointsLocal
  INTEGER OrderLocal
  DOUBLE PRECISION bfA, RmatA, ZfA, ZfpA
  DOUBLE PRECISION bfB, RmatB, ZfB, ZfpB
  DOUBLE PRECISION, ALLOCATABLE :: Xout(:,:), copen(:), cclosed(:), cfull(:)
  INTEGER, ALLOCATABLE :: OpenIdxOut(:), ClosedIdxOut(:)
  INTEGER NumOpen, Ncc, ix, lx, kx
  DOUBLE PRECISION xEndBase, kChan1, ThresholdRescaled1
  LOGICAL, PARAMETER :: UseShiftedXEnd = .TRUE.

  ! ------------------------------------------------------------------
  ! Common setup
  ! ------------------------------------------------------------------
  DataDir = '3bodydata_sech2_3boson_ch1'
  OrderLocal = 5
  xNumPoints = 400
  LegPointsLocal = 10
  xMinLocal = 1.0d-4
  xEndBase = 40.0d0
  EnergyLocal = -0.999d0

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

  ! Sensitivity-hypothesis test (user request): shift xEnd by a quarter wavelength of channel 1's
  ! asymptotic oscillation to move bf=-Zfp/Zf away from wherever it landed near zero at xEnd=40.
  ThresholdRescaled1 = 2d0*muLocal*Threshold(1)
  kChan1 = DSQRT(2d0*muLocal*EnergyLocal - ThresholdRescaled1)
  WRITE(6,*) "k_channel1 = ",kChan1," wavelength = ",2d0*Pi/kChan1
  IF (UseShiftedXEnd) THEN
     xEndLocal = xEndBase + Pi/(2d0*kChan1)
     ! Path A's tabulated dataset (3bodydata_sech2_3boson_ch1) only extends to R=50 -- cap there
     ! (full quarter-wavelength shift would need R=86.4, out of range for the tabulated data).
     xEndLocal = DMIN1(xEndLocal, 49.9d0)
  ELSE
     xEndLocal = xEndBase
  ENDIF

  WRITE(6,*) "reducedmass=",reducedmass," EffDim=",EffDim," Energy=",EnergyLocal, &
       " xEndBase=",xEndBase," xEndLocal (used)=",xEndLocal

  ! ------------------------------------------------------------------
  ! Path A: tabulated/adiabatic, single B-spline box, Left=2/NumOpenL=0 (TestAlphaZero.f90 recipe)
  ! ------------------------------------------------------------------
  BLOCK
    ALLOCATE(CoulombC(NumChannels))
    CoulombC = 0d0

    BPD1%NumChannels = NumChannels
    BPD1%Order = OrderLocal
    BPD1%Left = 2
    BPD1%Right = 2
    BPD1%xNumPoints = xNumPoints
    BPD1%kl = 1d20
    BPD1%kr = 1d20
    CALL AllocateBPD(BPD1)
    BPD1%xl = xMinLocal
    BPD1%xr = xEndLocal
    CALL GridMaker(BPD1%xPoints,BPD1%xNumPoints,BPD1%xl,BPD1%xr,"quadratic")
    CALL MakeBasis(BPD1)
    BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
    CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)

    EIG%MatrixDim = BPD1%MatrixDim
    CALL AllocateEIG(EIG)
    CALL CalcGamLam(BPD1,EIG,0,1)
    CALL CombineGam(EIG,EnergyLocal)

    Box1%NumOpenL = 0
    Box1%NumOpenR = 1
    CALL AllocateBox(Box1)
    CALL InitZeroBox(Box1)

    NumOpen = Box1%betaMax
    Ncc = BPD1%MatrixDim - NumOpen
    ALLOCATE(evalRed(NumOpen),evecRed(NumOpen,NumOpen),LamooDiag(NumOpen))
    ALLOCATE(Xout(Ncc,NumOpen),OpenIdxOut(NumOpen),ClosedIdxOut(Ncc))
    CALL PartitionAndEliminateFull(BPD1,EIG,Box1%NumOpenL,Box1%NumOpenR,evalRed,evecRed,LamooDiag, &
         Xout,OpenIdxOut,ClosedIdxOut)

    Bnull%NumOpenL = 0
    Bnull%NumOpenR = 0
    CALL AllocateBox(Bnull)
    CALL InitZeroBox(Bnull)
    CALL BoxMatch(Bnull, Box1, BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)

    ZfA = Box1%Zf(1,1)
    ZfpA = Box1%Zfp(1,1)
    bfA = -ZfpA/ZfA
    RmatA = ZfA**2/bfA

    WRITE(6,*) "=== Path A (tabulated/adiabatic, single B-spline box) ==="
    WRITE(6,*) "Zf  = ",ZfA
    WRITE(6,*) "Zfp = ",ZfpA
    WRITE(6,*) "bf  = ",bfA
    WRITE(6,*) "Rmat= ",RmatA

    ALLOCATE(copen(NumOpen),cclosed(Ncc),cfull(BPD1%MatrixDim))
    copen = evecRed(:,1)
    cclosed = -MATMUL(Xout,copen)
    DO ix = 1,NumOpen
       cfull(OpenIdxOut(ix)) = copen(ix)
    ENDDO
    DO ix = 1,Ncc
       cfull(ClosedIdxOut(ix)) = cclosed(ix)
    ENDDO

    OPEN(unit=91,file='wavefunction_adiabatic_singlebox.dat',status='replace')
    DO kx = 1,BPD1%xNumPoints-1
       DO lx = 1,LegPoints
          BLOCK
            DOUBLE PRECISION Rval, uval
            Rval = BPD1%x(lx,kx)
            uval = 0d0
            DO ix = 1,BPD1%xDim
               uval = uval + cfull((ix-1)*NumChannels+1)*BPD1%u(lx,kx,ix)
            ENDDO
            WRITE(91,*) Rval, uval
          END BLOCK
       ENDDO
    ENDDO
    CLOSE(91)
    WRITE(6,*) "Wrote wavefunction_adiabatic_singlebox.dat"

    DEALLOCATE(copen,cclosed,cfull,Xout,OpenIdxOut,ClosedIdxOut)
    CALL DeAllocateBox(Box1)
    CALL DeAllocateBox(Bnull)
    DEALLOCATE(evalRed,evecRed,LamooDiag)
    CALL DeAllocateEIG(EIG)
    CALL DeAllocateBPD(BPD1)
  END BLOCK

  ! ------------------------------------------------------------------
  ! Path B: SVD/direct, single box, NumStates=1, domain [0,pi/6], Left1D=Right1D=1, GridType='Radau'
  ! ------------------------------------------------------------------
  BLOCK
    TYPE(BPData) :: BPDSVD
    TYPE(GenEigVal) :: EIGSVD
    TYPE(BoxData) :: Box1SVD, BnullSVD
    DOUBLE PRECISION, ALLOCATABLE :: evalRedB(:), evecRedB(:,:), LamooDiagB(:)
    DOUBLE PRECISION, ALLOCATABLE :: XoutB(:,:), copenB(:), cclosedB(:), cfullB(:)
    INTEGER, ALLOCATABLE :: OpenIdxOutB(:), ClosedIdxOutB(:)
    DOUBLE PRECISION :: ODiagL(1), ODiagR(1), wsq1Out, wsqLOut
    CHARACTER*64 :: LegendreFile64, datadirLoc
    INTEGER, PARAMETER :: LSVD = 24
    ! NOT a literal here -- BuildSVDBoxGeneric's Shift dummy is declared INTENT(in) but is
    ! passed straight through into OneDimChannelsGeneric's own Shift argument, which WRITES
    ! to it (the adaptive per-R ARPACK shift update). RunHybridSVDGeneric always passes a
    ! genuine local variable, so this INTENT mismatch never bit it -- but a literal constant
    ! argument here gets placed in protected memory by the compiler, and the in-place write
    ! then SIGBUSes (KERN_PROTECTION_FAILURE). See memory/commit notes on this bug.
    DOUBLE PRECISION :: ShiftLocal
    DOUBLE PRECISION, ALLOCATABLE :: nodesR(:), weightsR(:)
    INTEGER :: NumOpenB, NccB, j

    LegendreFile64 = 'Legendre.dat'
    datadirLoc = '.'
    ShiftLocal = -5.0d0

    CALL BuildSVDBoxGeneric(0.0d0,xEndLocal,LSVD,1,5,1,1, &
         200,0d0,dacos(-1d0)/6.0d0,LegendreFile64,10,ShiftLocal,datadirLoc, &
         .FALSE.,SechPotential,SechGridMaker,SechRobinCoeff, &
         0,1,'Radau',BPDSVD,EIGSVD,ODiagL,ODiagR,wsq1Out,wsqLOut)

    CALL CombineGam(EIGSVD,EnergyLocal)

    Box1SVD%NumOpenL = 0
    Box1SVD%NumOpenR = 1
    CALL AllocateBox(Box1SVD)
    CALL InitZeroBox(Box1SVD)

    NumOpenB = Box1SVD%betaMax
    NccB = BPDSVD%MatrixDim - NumOpenB
    ALLOCATE(evalRedB(NumOpenB),evecRedB(NumOpenB,NumOpenB),LamooDiagB(NumOpenB))
    ALLOCATE(XoutB(NccB,NumOpenB),OpenIdxOutB(NumOpenB),ClosedIdxOutB(NccB))
    CALL PartitionAndEliminateFull(BPDSVD,EIGSVD,Box1SVD%NumOpenL,Box1SVD%NumOpenR,evalRedB,evecRedB,LamooDiagB, &
         XoutB,OpenIdxOutB,ClosedIdxOutB)

    BnullSVD%NumOpenL = 0
    BnullSVD%NumOpenR = 0
    CALL AllocateBox(BnullSVD)
    CALL InitZeroBox(BnullSVD)
    CALL BoxMatch(BnullSVD, Box1SVD, BPDSVD, evalRedB, evecRedB, LamooDiagB, EffDim, AlphaFactor)

    ZfB = Box1SVD%Zf(1,1)
    ZfpB = Box1SVD%Zfp(1,1)
    bfB = -ZfpB/ZfB
    RmatB = ZfB**2/bfB

    WRITE(6,*) "=== Path B (SVD/direct, single box, NumStates=1) ==="
    WRITE(6,*) "Zf  = ",ZfB
    WRITE(6,*) "Zfp = ",ZfpB
    WRITE(6,*) "bf  = ",bfB
    WRITE(6,*) "Rmat= ",RmatB

    ALLOCATE(copenB(NumOpenB),cclosedB(NccB),cfullB(BPDSVD%MatrixDim))
    copenB = evecRedB(:,1)
    cclosedB = -MATMUL(XoutB,copenB)
    DO j = 1,NumOpenB
       cfullB(OpenIdxOutB(j)) = copenB(j)
    ENDDO
    DO j = 1,NccB
       cfullB(ClosedIdxOutB(j)) = cclosedB(j)
    ENDDO

    ! Independently rebuild the L-point Radau grid (QOffset=0, LQuad=L exactly for Radau --
    ! see BuildSVDBoxGeneric's own header/GridType branch), matching its construction exactly,
    ! to recover the DVR nodes' own R-values/weights (BuildSVDBoxGeneric doesn't export them).
    ! Lagrange DVR basis pi_j(R)=L_j(R)/sqrt(w_j) has pi_j(R_j)=1/sqrt(w_j), pi_m(R_j)=0 for
    ! m/=j (Kronecker-delta property) -- so u(R_j)=cfullB_j/sqrt(w_j) directly, same as
    ! PlotWavefunctionSVD.f90's own reconstruction.
    ALLOCATE(nodesR(LSVD),weightsR(LSVD))
    CALL GetGaussRadauFactors(LSVD,nodesR,weightsR)
    nodesR = 0.5d0*(0.0d0+xEndLocal) + 0.5d0*(xEndLocal-0.0d0)*nodesR
    weightsR = 0.5d0*(xEndLocal-0.0d0)*weightsR

    OPEN(unit=92,file='wavefunction_svd_singlebox.dat',status='replace')
    DO j = 1,LSVD
       WRITE(92,*) nodesR(j), cfullB(j)/DSQRT(weightsR(j))
    ENDDO
    CLOSE(92)
    WRITE(6,*) "Wrote wavefunction_svd_singlebox.dat"
  END BLOCK

  WRITE(6,*)
  WRITE(6,*) "=== Comparison ==="
  WRITE(6,*) "Rmat_A (adiabatic) = ",RmatA
  WRITE(6,*) "Rmat_B (SVD)        = ",RmatB
  WRITE(6,*) "relative diff       = ",DABS(RmatB-RmatA)/DABS(RmatA)

END PROGRAM CompareRmatrixSVDvsAdiabatic
