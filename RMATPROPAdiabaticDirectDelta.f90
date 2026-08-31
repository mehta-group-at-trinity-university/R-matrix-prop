!****************************************************************************************************
! No-interpolation, ordinary-adiabatic-box (NOT SVD/DVR) analog of RMATPROPAdiabatic.f90 for the
! equal-mass 3-boson delta-function problem, mirroring RMATPROPAdiabaticDirect3Boson.f90 exactly but
! using DeltaAdiabaticPotentialDirect.f90 (DeltaFunctionPlugins physics: Robin BC at phi=pi/6,
! kRight=sqrt(2*mu)*R) instead of SechAdiabaticPotentialDirect.f90, and hardcoded analytic
! Threshold/Leff (same formula as RMATPROPHybridSVDGenericDelta.f90 -- channel 1 Threshold=-1,
! Leff=0.5; channels i>=2 Threshold=0, Leff=6i-9) instead of reading Fit.data/FitLeff.data, since
! this driver needs no DataDir at all.
!
! Purpose: this session's final diagnostic. RMATPROPHybridSVDGenericDelta.x (SVD/DVR boxes, generic
! ARPACK-based per-R diagonalization via GenericSVDChannelBasis.f90) gives badly wrong results (both
! with and without a B-spline tail) while RMATPROPHybridSVDAnalytic.x/RMATPROPHybridSVDRobin.x/
! RMATPROPAdiabatic.x (tabulated-interpolated) all agree with the exact Mehta-Shepard result. This
! driver isolates whether OneDimChannelsGeneric's ARPACK+FixPhase diagonalization itself is at fault,
! or whether the fault is specific to GenericSVDChannelBasis.f90's shared-single-R-overlap-matrix
! cross-R construction (this session's working diagnosis) -- by running the SAME ARPACK/FixPhase
! engine (OneDimChannelsGeneric) through the ORDINARY adiabatic box propagator (MakeBasis/CalcGamLam,
! needing only per-quadrature-point Uad/P/Q via Hellmann-Feynman -- no cross-R overlap matrix at
! all), on the SAME box grid RMATPROPHybridSVDGenericDelta.x's own no-tail R=100 test used.
!****************************************************************************************************
PROGRAM main
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE DeltaAdiabaticPotentialDirect

  IMPLICIT NONE
  CHARACTER(LEN=16) EnergyGridType
  INTEGER ios
  TYPE(BPData) BPD,BPD0,BPD1,BPDLast
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: xprim(:)
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Boxes(:,:,:), OverlapBoxes(:,:,:), LamBoxes(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Box1(:,:), OverlapBox1(:,:), LamBox1(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  DOUBLE PRECISION muLocal
  DOUBLE PRECISION Emin, Emax, qOurs, Kelem, rawDelta, delta, prevDelta
  INTEGER NumChannelsLoc, NumOpenChannels, xNumPoints, LegPointsLocal
  INTEGER i, j, n, iBox, lx, kx, NumEnergies, ie, nBranch, PhaseFile, KMatFile, RecombFile
  LOGICAL havePrevDelta
  COMPLEX*16, ALLOCATABLE :: Amat(:,:), Smat(:,:)
  DOUBLE PRECISION recombProb
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION, PARAMETER :: DeltaMu = 0.57735026918962584d0
  DOUBLE PRECISION, PARAMETER :: DeltaAlpha = 0d0
  DOUBLE PRECISION, PARAMETER :: DeltaDdim = 2d0

  !----------------------------------------------------------------------
  ! Same .inp format as RMATPROPHybridSVDGenericDelta.inp's own first 5 blocks -- no DataDir,
  ! Threshold/Leff generated analytically below instead.
  !----------------------------------------------------------------------
  OPEN(unit=7,file='RMATPROPAdiabaticDirectDelta.inp',status='old')
  READ(7,*)
  READ(7,*) NumChannelsLoc
  READ(7,*)
  READ(7,*)
  READ(7,*) Order, xNumPoints, LegPointsLocal
  READ(7,*)
  READ(7,*)
  READ(7,*) NumBoxes, xStart, xEnd
  READ(7,*)
  READ(7,*)
  READ(7,*) Emin, Emax, NumEnergies
  READ(7,*,IOSTAT=ios) EnergyGridType
  IF (ios.NE.0) EnergyGridType = "linear"
  CLOSE(7)

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = LegPointsLocal
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  NumChannels = NumChannelsLoc
  muLocal = DeltaMu
  AlphaFactor = DeltaAlpha
  EffDim = DeltaDdim
  reducedmass = muLocal

  ! Equal-mass 3-boson delta-function problem's exact hyperspherical channel tower (Mehta, Esry &
  ! Greene, PRA 76, 022711 (2007)): channel 1 (atom-dimer, Threshold=-1) has Leff=1/2; excited
  ! (3-body-threshold, Threshold=0) channels have Leff=6i-9 -- same formula
  ! RMATPROPHybridSVDGenericDelta.f90 uses.
  ALLOCATE(Threshold(NumChannels),Leff(NumChannels))
  Threshold(1) = -1d0
  Leff(1) = 0.5d0
  DO i = 2,NumChannels
     Threshold(i) = 0d0
     Leff(i) = 6d0*DBLE(i) - 9d0
  ENDDO
  Threshold = 2d0*muLocal*Threshold

  WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," ddim = ",EffDim," mu = ",reducedmass
  WRITE(6,*) "Energy scan: Emin = ",Emin," Emax = ",Emax," NumEnergies = ",NumEnergies

  PhaseFile = 50
  OPEN(unit=PhaseFile,file='PhaseShift.dat',status='replace')
  WRITE(PhaseFile,'(A)') '# q  Energy  K  delta(rad, unwrapped, no anchor/offset applied)'
  WRITE(PhaseFile,'(A)') '@    xaxis label "q"'
  WRITE(PhaseFile,'(A)') '@    yaxis label "delta (rad)"'

  KMatFile = 51
  OPEN(unit=KMatFile,file='KMatrixFull.dat',status='replace')

  RecombFile = 52
  OPEN(unit=RecombFile,file='Recombination.dat',status='replace')
  WRITE(RecombFile,'(A)') '# Energy  NumOpenChannels  1-|S(1,1)|^2  (McGuire exact-unitarity test: should -> 0)'
  WRITE(RecombFile,'(A)') '@    xaxis label "Energy"'
  WRITE(RecombFile,'(A)') '@    yaxis label "1 - |S11|^2"'

  !---------------------------------------------------------------------
  ! ONE-TIME setup (energy-independent): box 1, interior boxes 2..NumBoxes-1, and the last box,
  ! all built and their real potential evaluated via SetAdiabaticPotentialDirect exactly once.
  !---------------------------------------------------------------------
  BPD1%NumChannels = NumChannels
  BPD1%Order = Order
  BPD1%Left = 2
  BPD1%Right = 2
  BPD1%xNumPoints = xNumPoints
  BPD1%kl = 1d20
  BPD1%kr = 1d20

  BPD0%NumChannels = NumChannels
  BPD0%Order = Order
  BPD0%Left = 2
  BPD0%Right = 2
  BPD0%xNumPoints = xNumPoints
  BPD0%kl = 1d20
  BPD0%kr = 1d20
  CALL AllocateBPD(BPD0)
  ALLOCATE(xprim(xNumPoints))
  CALL GridMakerLinear(BPD0%xNumPoints,0d0,1d0,xprim)
  BPD0%xPoints = xprim
  CALL MakeBasis(BPD0)
  BPD0%lam(1:NumChannels) = Leff(1:NumChannels)

  BPD = BPD0
  BPDLast = BPD0

  EIG%MatrixDim = BPD%MatrixDim
  CALL AllocateEIG(EIG)

  WRITE(6,*) "Building interior boxes (SetAdiabaticPotentialDirect, one ARPACK diagonalization per quadrature point)..."
  IF (NumBoxes.GE.3) THEN
     ALLOCATE(Gam0Boxes(BPD%MatrixDim,BPD%MatrixDim,2:NumBoxes-1))
     ALLOCATE(OverlapBoxes(BPD%MatrixDim,BPD%MatrixDim,2:NumBoxes-1))
     ALLOCATE(LamBoxes(BPD%MatrixDim,BPD%MatrixDim,2:NumBoxes-1))
     DO iBox = 2,NumBoxes-1
        BPD%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xPoints = BPD%xl + (BPD%xr-BPD%xl)*xprim
        DO kx = 1,BPD%xNumPoints
           DO lx = 1,LegPoints
              BPD%ux(lx,kx,1:BPD%xDim) = BPD0%ux(lx,kx,1:BPD%xDim)/(BPD%xr-BPD%xl)
           ENDDO
        ENDDO
        CALL SetAdiabaticPotentialDirect(BPD,muLocal)
        CALL CalcGamLam(BPD,EIG,NumChannels,NumChannels)
        Gam0Boxes(:,:,iBox) = EIG%Gam0
        OverlapBoxes(:,:,iBox) = EIG%Overlap
        LamBoxes(:,:,iBox) = EIG%Lam
        WRITE(6,'(A,I5,A,I5)') "  interior box ",iBox," / ",NumBoxes-1
     ENDDO
  ENDIF

  WRITE(6,*) "Building box 1..."
  IF (NumBoxes.GE.2) THEN
     ALLOCATE(Gam0Box1(BPD%MatrixDim,BPD%MatrixDim))
     ALLOCATE(OverlapBox1(BPD%MatrixDim,BPD%MatrixDim))
     ALLOCATE(LamBox1(BPD%MatrixDim,BPD%MatrixDim))
     CALL AllocateBPD(BPD1)
     BPD1%xl = xStart
     BPD1%xr = xStart + (xEnd-xStart)*(1d0/DBLE(NumBoxes))**BoxSpacingPower
     CALL GridMaker(BPD1%xPoints,BPD1%xNumPoints,BPD1%xl,BPD1%xr,"quadratic")
     CALL MakeBasis(BPD1)
     BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
     CALL SetAdiabaticPotentialDirect(BPD1,muLocal)
     CALL CalcGamLam(BPD1,EIG,0,NumChannels)
     Gam0Box1 = EIG%Gam0
     OverlapBox1 = EIG%Overlap
     LamBox1 = EIG%Lam
  ENDIF

  WRITE(6,*) "Building last box..."
  BPDLast%xl = xStart + (xEnd-xStart)*(DBLE(NumBoxes-1)/DBLE(NumBoxes))**BoxSpacingPower
  BPDLast%xr = xStart + (xEnd-xStart)*(DBLE(NumBoxes)/DBLE(NumBoxes))**BoxSpacingPower
  BPDLast%xPoints = BPDLast%xl + (BPDLast%xr-BPDLast%xl)*xprim
  DO kx = 1,BPDLast%xNumPoints
     DO lx = 1,LegPoints
        BPDLast%ux(lx,kx,1:BPDLast%xDim) = BPD0%ux(lx,kx,1:BPDLast%xDim)/(BPDLast%xr-BPDLast%xl)
     ENDDO
  ENDDO
  CALL SetAdiabaticPotentialDirect(BPDLast,muLocal)
  WRITE(6,*) "All boxes built. Starting energy scan..."

  havePrevDelta = .FALSE.

  DO ie = 1,NumEnergies
     IF (NumEnergies.GT.1) THEN
        IF (EnergyGridType.EQ."quadratic") THEN
           Energy = Emin + (Emax-Emin)*(DBLE(ie-1)/DBLE(NumEnergies-1))**2
        ELSE
           Energy = Emin + (Emax-Emin)*DBLE(ie-1)/DBLE(NumEnergies-1)
        ENDIF
     ELSE
        Energy = Emin
     ENDIF

     NumOpenChannels = 0
     DO i = 1,NumChannels
        IF (2d0*reducedmass*Energy.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
     ENDDO

     ALLOCATE(Boxes(NumBoxes))
     DO i = 1,NumBoxes
        Boxes(i)%NumOpenL = MERGE(0,NumChannels,i.EQ.1)
        IF (i.EQ.NumBoxes) THEN
           Boxes(i)%NumOpenR = NumOpenChannels
        ELSE
           Boxes(i)%NumOpenR = NumChannels
        ENDIF
        IF (i.EQ.1) THEN
           Boxes(i)%xl = xStart
        ELSE
           Boxes(i)%xl = Boxes(i-1)%xr
        ENDIF
        Boxes(i)%xr = xStart + (xEnd-xStart)*(DBLE(i)/DBLE(NumBoxes))**BoxSpacingPower
        CALL AllocateBox(Boxes(i))
        CALL InitZeroBox(Boxes(i))
     ENDDO

     Bnull%NumOpenL = 0
     Bnull%NumOpenR = 0
     CALL AllocateBox(Bnull)
     CALL InitZeroBox(Bnull)

     CALL AllocateScat(SD,NumChannels)

     IF (NumBoxes.GE.2) THEN
        EIG%Gam0 = Gam0Box1
        EIG%Overlap = OverlapBox1
        EIG%Lam = LamBox1
        CALL CombineGam(EIG,Energy)
        ALLOCATE(evalRed(Boxes(1)%betaMax),evecRed(Boxes(1)%betaMax,Boxes(1)%betaMax))
        ALLOCATE(LamooDiag(Boxes(1)%betaMax))
        CALL PartitionAndEliminate(BPD1,EIG,Boxes(1)%NumOpenL,Boxes(1)%NumOpenR,evalRed,evecRed,LamooDiag)
        CALL BoxMatch(Bnull, Boxes(1), BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        DEALLOCATE(evalRed,evecRed,LamooDiag)
     ENDIF

     DO iBox = 2,NumBoxes-1
        BPD%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
        EIG%Gam0 = Gam0Boxes(:,:,iBox)
        EIG%Overlap = OverlapBoxes(:,:,iBox)
        EIG%Lam = LamBoxes(:,:,iBox)
        CALL CombineGam(EIG,Energy)
        ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
        ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
        CALL PartitionAndEliminate(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR,evalRed,evecRed,LamooDiag)
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        DEALLOCATE(evalRed,evecRed,LamooDiag)
     ENDDO

     iBox = NumBoxes
     CALL CalcGamLam(BPDLast,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
     CALL CombineGam(EIG,Energy)
     ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
     ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
     CALL PartitionAndEliminate(BPDLast,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR,evalRed,evecRed,LamooDiag)
     IF (NumBoxes.GE.2) THEN
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPDLast, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
     ELSE
        CALL BoxMatch(Bnull, Boxes(iBox), BPDLast, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
     ENDIF
     DEALLOCATE(evalRed,evecRed,LamooDiag)

     CALL CalcK(Boxes(NumBoxes),BPDLast,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)

     IF (NumOpenChannels.EQ.NumChannels .AND. NumChannels.GT.1) THEN
        WRITE(KMatFile,'(A,1P,E20.12)') '# Energy = ',Energy
        DO i = 1,NumOpenChannels
           WRITE(KMatFile,'(1P,20E20.12)') SD%K(i,1:NumOpenChannels)
        ENDDO
     ENDIF

     ALLOCATE(Amat(NumOpenChannels,NumOpenChannels),Smat(NumOpenChannels,NumOpenChannels))
     Amat = (0d0,0d0)
     Smat = (0d0,0d0)
     DO i = 1,NumOpenChannels
        Amat(i,1:NumOpenChannels) = -(0d0,1d0)*SD%K(i,1:NumOpenChannels)
        Smat(i,1:NumOpenChannels) = (0d0,1d0)*SD%K(i,1:NumOpenChannels)
        Amat(i,i) = Amat(i,i) + (1d0,0d0)
        Smat(i,i) = Smat(i,i) + (1d0,0d0)
     ENDDO
     CALL CompSqrMatInv(Amat,NumOpenChannels)
     Smat = MATMUL(Amat,Smat)
     recombProb = 1d0 - ABS(Smat(1,1))**2
     WRITE(RecombFile,'(1P,E20.12,I5,E20.12)') Energy, NumOpenChannels, recombProb
     DEALLOCATE(Amat,Smat)

     IF (NumOpenChannels.EQ.1) THEN
        Kelem = SD%K(1,1)
        rawDelta = DATAN(Kelem)
        IF (havePrevDelta) THEN
           nBranch = NINT((prevDelta-rawDelta)/Pi)
           delta = rawDelta + DBLE(nBranch)*Pi
        ELSE
           delta = rawDelta
        ENDIF
        prevDelta = delta
        havePrevDelta = .TRUE.
        qOurs = DSQRT(2d0*reducedmass*Energy-Threshold(1))
        WRITE(PhaseFile,'(1P,4E20.12)') qOurs, Energy, Kelem, delta
        WRITE(6,'(A,I5,A,1P,E14.6,A,E14.6,A,E14.6)') " ie=",ie," Energy=",Energy," q=",qOurs," delta=",delta
     ELSE
        WRITE(6,'(A,I5,A,1P,E14.6,A,I3,A)') " ie=",ie," Energy=",Energy," NumOpenChannels=",NumOpenChannels, &
             " (skipped -- phase shift only written for exactly 1 open channel)"
     ENDIF

     DO i = 1,NumBoxes
        CALL DeAllocateBox(Boxes(i))
     ENDDO
     DEALLOCATE(Boxes)
     CALL DeAllocateBox(Bnull)
     CALL DeAllocateScat(SD)
  ENDDO

  CALL DeAllocateBPD(BPD0)
  IF (NumBoxes.GE.2) THEN
     CALL DeAllocateBPD(BPD1)
     DEALLOCATE(Gam0Box1,OverlapBox1,LamBox1)
  ENDIF
  CALL DeAllocateBPD(BPDLast)
  DEALLOCATE(xprim)
  CALL DeAllocateEIG(EIG)
  IF (ALLOCATED(Gam0Boxes)) DEALLOCATE(Gam0Boxes,OverlapBoxes,LamBoxes)

  CLOSE(PhaseFile)
  CLOSE(KMatFile)
  CLOSE(RecombFile)
  WRITE(6,*)
  WRITE(6,*) "Phase shift written to PhaseShift.dat"

END PROGRAM main
