!****************************************************************************************************
! Hybrid SVD / B-spline R-matrix propagator. Structure:
!   box 1            : SVD box, Gauss-RADAU DVR grid (no node at R=0 -- see
!                       SVDChannelBasis.f90's GetGaussRadauFactors header note),
!                       NumOpenL=0 (regularity: no left-boundary DOF exists to retain).
!   boxes 2..NumSVDBoxes  : SVD boxes, ordinary Gauss-LOBATTO DVR grid (both edges are
!                       nodes, all channels retained at both edges -- interior boxes).
!   boxes NumSVDBoxes+1..NumBoxes : ordinary B-spline / smooth-adiabatic boxes, reading
!                       tabulated Uad/Pmat/Qmat via AdiabaticPotential's interpolants
!                       (unchanged from RMATPROPAdiabatic.f90).
! Rswitch = xStart + (xEnd-xStart)*NumSVDBoxes/NumBoxes falls out of the existing
! uniform box-edge spacing formula (BoxSpacingPower=1).
!
! Every box, regardless of type, is fed through the SAME CombineGam/PartitionAndEliminate
! /BoxMatch pipeline (RMatPropCore.f90) -- confirmed basis-agnostic. All boxes here are
! energy-independent (SVD boxes: OneDimChannels' per-R diagonalization doesn't depend on
! scattering energy; interior B-spline boxes: same as RMATPROPAdiabatic.f90) except the
! last box, whose NumOpenR tracks NumOpenChannels(Energy) and is rebuilt every energy.
!
! No CoulombC/AnalyticComboNearOrigin near-origin patch is needed here (unlike
! RMATPROPAdiabatic.f90's box 1): that machinery existed specifically to handle
! tabulated-data interpolation's Runge-ringing/domain-edge behavior near R=0. Box 1 here
! computes fresh ARPACK adiabatic curves directly at its own (Radau, R>0) DVR points --
! no interpolation, no near-origin patch needed, PROVIDED the SVD region stays below
! GridMakerHHL's xRswitch=10*PotRange/Pi (~3.18 for PotRange=1), where it uses a plain
! uniform angular grid with no r0/R term at all. This is expected to hold for the
! short-range/avoided-crossing region SVD targets; if the SVD region needs to extend
! beyond xRswitch, GridMakerHHL's r0/R term would need revisiting for R very close to 0.
!****************************************************************************************************
PROGRAM main
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE AdiabaticPotential
  USE SVDChannelBasis
  USE potential, ONLY: V2Depth, PotRange, alpha

  IMPLICIT NONE
  CHARACTER(LEN=64) DataDir
  CHARACTER(LEN=16) EnergyGridType
  INTEGER ios
  TYPE(BPData) BPD,BPD0
  TYPE(BPData), ALLOCATABLE :: BPDSVDBoxes(:)
  TYPE(GenEigVal) EIG,EIGSVD
  TYPE(BoxData) Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: xprim(:)
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Boxes(:,:,:), OverlapBoxes(:,:,:), LamBoxes(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Gam0SVDBoxes(:,:,:), OverlapSVDBoxes(:,:,:), LamSVDBoxes(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: ODiagLSVDBoxes(:,:), ODiagRSVDBoxes(:,:), wsq1SVDBoxes(:), wsqLSVDBoxes(:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alphaGlobal, EffDimLocal
  DOUBLE PRECISION Emin, Emax, qOurs, Kelem, rawDelta, delta, prevDelta
  INTEGER NumDataPoints, NumOpenChannels, xNumPoints, LegPointsLocal
  INTEGER i, iBox, lx, kx, NumEnergies, ie, nBranch, PhaseFile, KMatFile, RecombFile, UnitFile
  LOGICAL havePrevDelta
  COMPLEX*16, ALLOCATABLE :: Amat(:,:), Smat(:,:)
  DOUBLE PRECISION recombProb
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION Rswitch

  INTEGER NumSVDBoxes, LSVD, Order1D, Left1D, Right1D, xNumPoints1D, LegPoints1D
  DOUBLE PRECISION xMin1D, xMax1D, Shift1D
  INTEGER MatrixDimSVD

  OPEN(unit=7,file='RMATPROPHybridSVD.inp',status='old')
  READ(7,*)
  READ(7,*) DataDir
  READ(7,*)
  READ(7,*)
  READ(7,*) Order, xNumPoints, LegPoints
  READ(7,*)
  READ(7,*)
  READ(7,*) NumBoxes, xStart, xEnd
  READ(7,*)
  READ(7,*)
  READ(7,*) Emin, Emax, NumEnergies
  READ(7,*,IOSTAT=ios) EnergyGridType
  IF (ios.NE.0) EnergyGridType = "linear"
  READ(7,*)
  READ(7,*)
  READ(7,*) NumSVDBoxes, LSVD, Order1D, Left1D, Right1D
  READ(7,*)
  READ(7,*)
  READ(7,*) xNumPoints1D, xMin1D, xMax1D, LegPoints1D
  READ(7,*)
  READ(7,*)
  READ(7,*) V2Depth, PotRange, alpha
  CLOSE(unit=7)
  ! Was -20.d0, back when adiabaticSolver1D.f90 unconditionally overwrote whatever
  ! Shift the caller passed in (a bug fixed during the free-particle diagnostic --
  ! see that file's header note). Now that Shift genuinely seeds ARPACK's shift-
  ! invert search, -20 is badly mismatched for a purely-repulsive (all-positive-
  ! eigenvalue) spectrum: -1 stays safely on the correct side of zero without being
  ! absurdly far from the target eigenvalues, matching what converged cleanly for
  ! both the sech^2 and free-particle tabulated-dataset generation runs.
  Shift1D = -1.d0

  IF (NumSVDBoxes.LT.1 .OR. NumSVDBoxes.GT.NumBoxes) THEN
     WRITE(6,*) "RMATPROPHybridSVD requires 1 <= NumSVDBoxes <= NumBoxes"
     STOP
  ENDIF

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  CALL ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alphaGlobal,EffDimLocal, &
       Threshold,Leff,UInterp,PInterp,QInterp)
  reducedmass = muLocal
  AlphaFactor = alphaGlobal
  EffDim = EffDimLocal

  ! No near-origin Coulomb patch needed (see file header) -- CoulombC stays zero,
  ! only used by the B-spline-tail region's SetAdiabaticPotential calls (which never
  ! reach R=0 here since Rswitch>0), where it's an ordinary no-op default.
  ALLOCATE(CoulombC(NumChannels))
  CoulombC = 0d0

  Threshold = 2d0*reducedmass*Threshold
  Rswitch = xStart + (xEnd-xStart)*DBLE(NumSVDBoxes)/DBLE(NumBoxes)
  WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," EffDim = ",EffDim," mu = ",reducedmass
  WRITE(6,*) "Energy scan: Emin = ",Emin," Emax = ",Emax," NumEnergies = ",NumEnergies
  WRITE(6,*) "NumSVDBoxes = ",NumSVDBoxes," Rswitch = ",Rswitch," NumBoxes = ",NumBoxes," xEnd = ",xEnd
  WRITE(6,*) "SVD box: L = ",LSVD," V2Depth = ",V2Depth," PotRange = ",PotRange," potAlpha = ",alpha

  PhaseFile = 50
  OPEN(unit=PhaseFile,file='PhaseShift.dat',status='replace')
  WRITE(PhaseFile,'(A)') '# q  Energy  K  delta(rad, unwrapped, no anchor/offset applied)'
  WRITE(PhaseFile,'(A)') '@    xaxis label "q"'
  WRITE(PhaseFile,'(A)') '@    yaxis label "delta (rad)"'

  KMatFile = 51
  OPEN(unit=KMatFile,file='KMatrixFull.dat',status='replace')

  RecombFile = 52
  OPEN(unit=RecombFile,file='Recombination.dat',status='replace')
  WRITE(RecombFile,'(A)') '# Energy  NumOpenChannels  ||S^dagger S - I||_max  (general unitarity check)'
  WRITE(RecombFile,'(A)') '@    xaxis label "Energy"'
  WRITE(RecombFile,'(A)') '@    yaxis label "||S^dagger S - I||_max"'

  !---------------------------------------------------------------------
  ! ONE-TIME setup: SVD boxes 1..NumSVDBoxes (box 1 Radau, rest Lobatto), cached
  ! Gam0/Overlap/Lam for B-spline interior boxes NumSVDBoxes+1..NumBoxes-1.
  !---------------------------------------------------------------------
  ALLOCATE(BPDSVDBoxes(NumSVDBoxes))
  MatrixDimSVD = LSVD*NumChannels
  ALLOCATE(Gam0SVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
  ALLOCATE(OverlapSVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
  ALLOCATE(LamSVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
  ! Stashed per-box boundary data (see SVDChannelBasis.f90's RebuildSVDBoxLam) --
  ! only actually needed for whichever box ends up being "the last box"
  ! (index NumSVDBoxes, when NumSVDBoxes==NumBoxes: no B-spline tail at all, so the
  ! usual "rebuild the last box fresh every energy" trick has to apply to an SVD box
  ! instead), but cheap enough to just capture for every SVD box uniformly.
  ALLOCATE(ODiagLSVDBoxes(NumChannels,NumSVDBoxes),ODiagRSVDBoxes(NumChannels,NumSVDBoxes))
  ALLOCATE(wsq1SVDBoxes(NumSVDBoxes),wsqLSVDBoxes(NumSVDBoxes))

  DO iBox = 1,NumSVDBoxes
     BLOCK
       DOUBLE PRECISION :: a1,a2
       INTEGER :: NOpL,NOpR
       CHARACTER(LEN=16) :: GType
       a1 = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
       a2 = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
       IF (iBox.EQ.1) THEN
          ! True Dirichlet u(0)=0 at AlphaFactor=(EffDim-1)/2 -- 'Radau' avoided handing
          ! R=0 to OneDimChannels but left the box unconstrained/open there, which is the
          ! wrong regularity condition for this alpha (see SVDChannelBasis.f90's grid-
          ! setup comment and delta_function_validation.md).
          GType = 'LobattoDirichlet'
          NOpL = 0
       ELSE
          GType = 'Lobatto'
          NOpL = NumChannels
       ENDIF
       NOpR = NumChannels   ! cache with all channels retained -- Gam0/Overlap don't
                             ! depend on NumOpenR anyway; the true (possibly smaller,
                             ! energy-dependent) NumOpenR for a last-box SVD box is
                             ! applied per-energy below via RebuildSVDBoxLam.
       UnitFile = 200+iBox
       ! BuildSVDBox allocates EIGSVD internally (EIG%MatrixDim=BPD%MatrixDim;
       ! CALL AllocateEIG(EIG)) -- deallocate after each use so the next box's call
       ! doesn't hit "already allocated".
       CALL BuildSVDBox(a1,a2,LSVD,NumChannels,Order1D,Left1D,Right1D, &
            xNumPoints1D,xMin1D,xMax1D,LegendreFile,LegPoints1D,Shift1D,DataDir, &
            NOpL,NOpR,GType,BPDSVDBoxes(iBox),EIGSVD, &
            ODiagLSVDBoxes(:,iBox),ODiagRSVDBoxes(:,iBox),wsq1SVDBoxes(iBox),wsqLSVDBoxes(iBox))
       Gam0SVDBoxes(:,:,iBox) = EIGSVD%Gam0
       OverlapSVDBoxes(:,:,iBox) = EIGSVD%Overlap
       LamSVDBoxes(:,:,iBox) = EIGSVD%Lam
       CALL DeAllocateEIG(EIGSVD)
       WRITE(6,'(A,I4,A,1P,E14.6,A,E14.6,A,A)') " SVD box ",iBox," xl=",a1," xr=",a2," grid=",TRIM(GType)
     END BLOCK
  ENDDO

  ! If the SVD region reaches all the way to the outer edge (no B-spline tail), the
  ! last SVD box also plays the role CalcK's asymptotic Bessel matching needs Leff
  ! for -- BuildSVDBox otherwise leaves BPD%lam at its placeholder zero.
  IF (NumSVDBoxes.EQ.NumBoxes) THEN
     BPDSVDBoxes(NumSVDBoxes)%lam(1:NumChannels) = Leff(1:NumChannels)
  ENDIF

  ! Re-allocate EIGSVD once more as the scratch structure reused every energy in the
  ! main loop below (CombineGam/PartitionAndEliminate) -- same MatrixDim for every SVD
  ! box, so one allocation serves all of them.
  EIGSVD%MatrixDim = MatrixDimSVD
  CALL AllocateEIG(EIGSVD)

  ! No B-spline tail at all when NumSVDBoxes==NumBoxes -- skip BPD0/BPD/EIG entirely
  ! in that case (nothing downstream references them unless a B-spline box exists).
  IF (NumSVDBoxes.LT.NumBoxes) THEN
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

     EIG%MatrixDim = BPD%MatrixDim
     CALL AllocateEIG(EIG)

     IF (NumBoxes.GE.NumSVDBoxes+2) THEN
        ALLOCATE(Gam0Boxes(BPD%MatrixDim,BPD%MatrixDim,NumSVDBoxes+1:NumBoxes-1))
        ALLOCATE(OverlapBoxes(BPD%MatrixDim,BPD%MatrixDim,NumSVDBoxes+1:NumBoxes-1))
        ALLOCATE(LamBoxes(BPD%MatrixDim,BPD%MatrixDim,NumSVDBoxes+1:NumBoxes-1))
        DO iBox = NumSVDBoxes+1,NumBoxes-1
           BPD%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
           BPD%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
           BPD%xPoints = BPD%xl + (BPD%xr-BPD%xl)*xprim
           DO kx = 1,BPD%xNumPoints
              DO lx = 1,LegPoints
                 BPD%ux(lx,kx,1:BPD%xDim) = BPD0%ux(lx,kx,1:BPD%xDim)/(BPD%xr-BPD%xl)
              ENDDO
           ENDDO
           CALL SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp,CoulombC)
           CALL CalcGamLam(BPD,EIG,NumChannels,NumChannels)
           Gam0Boxes(:,:,iBox) = EIG%Gam0
           OverlapBoxes(:,:,iBox) = EIG%Overlap
           LamBoxes(:,:,iBox) = EIG%Lam
        ENDDO
     ENDIF
  ENDIF

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

     !--- SVD boxes 1..NumSVDBoxes ---
     DO iBox = 1,NumSVDBoxes
        EIGSVD%Gam0 = Gam0SVDBoxes(:,:,iBox)
        EIGSVD%Overlap = OverlapSVDBoxes(:,:,iBox)
        IF (iBox.EQ.NumBoxes) THEN
           ! This SVD box IS the terminal box (only possible when NumSVDBoxes==
           ! NumBoxes, i.e. no B-spline tail at all) -- its NumOpenR tracks
           ! NumOpenChannels(Energy), so the cached Lam (built with NumOpenR=
           ! NumChannels at setup time) isn't valid here; rebuild cheaply instead
           ! of redoing the box's ARPACK diagonalization every energy.
           CALL RebuildSVDBoxLam(EIGSVD,ODiagLSVDBoxes(:,iBox),ODiagRSVDBoxes(:,iBox), &
                wsq1SVDBoxes(iBox),wsqLSVDBoxes(iBox),BPDSVDBoxes(iBox)%xl,BPDSVDBoxes(iBox)%xr, &
                NumChannels,LSVD,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
        ELSE
           EIGSVD%Lam = LamSVDBoxes(:,:,iBox)
        ENDIF
        CALL CombineGam(EIGSVD,Energy)
        ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
        ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
        CALL PartitionAndEliminate(BPDSVDBoxes(iBox),EIGSVD,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR, &
             evalRed,evecRed,LamooDiag)
        IF (iBox.EQ.1) THEN
           CALL BoxMatch(Bnull, Boxes(1), BPDSVDBoxes(1), evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        ELSE
           CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPDSVDBoxes(iBox), evalRed, evecRed, LamooDiag, &
                EffDim, AlphaFactor)
        ENDIF
        DEALLOCATE(evalRed,evecRed,LamooDiag)
     ENDDO

     IF (NumSVDBoxes.LT.NumBoxes) THEN
        !--- interior B-spline boxes NumSVDBoxes+1..NumBoxes-1 ---
        DO iBox = NumSVDBoxes+1,NumBoxes-1
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

        !--- last box (rebuilt fresh every energy: NumOpenR can change with energy) ---
        iBox = NumBoxes
        BPD%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xPoints = BPD%xl + (BPD%xr-BPD%xl)*xprim
        DO kx = 1,BPD%xNumPoints
           DO lx = 1,LegPoints
              BPD%ux(lx,kx,1:BPD%xDim) = BPD0%ux(lx,kx,1:BPD%xDim)/(BPD%xr-BPD%xl)
           ENDDO
        ENDDO
        CALL SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp,CoulombC)
        CALL CalcGamLam(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
        CALL CombineGam(EIG,Energy)
        ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
        ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
        CALL PartitionAndEliminate(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR,evalRed,evecRed,LamooDiag)
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        DEALLOCATE(evalRed,evecRed,LamooDiag)

        CALL CalcK(Boxes(NumBoxes),BPD,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
     ELSE
        ! Every box is SVD -- the terminal box's own BPDSVDBoxes descriptor already
        ! carries the correct xl/xr/lam (Leff, set once above).
        CALL CalcK(Boxes(NumBoxes),BPDSVDBoxes(NumBoxes),SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
     ENDIF

     IF (NumOpenChannels.EQ.NumChannels .AND. NumChannels.GT.1) THEN
        WRITE(KMatFile,'(A,1P,E20.12)') '# Energy = ',Energy
        DO i = 1,NumOpenChannels
           WRITE(KMatFile,'(1P,20E20.12)') SD%K(i,1:NumOpenChannels)
        ENDDO
     ENDIF

     !---------------------------------------------------------------------
     ! GENERAL multichannel unitarity check: ||S^dagger S - I||_max over ALL open
     ! channels (not the trivial single-channel Cayley identity |S11|=1, which holds
     ! for ANY real 1x1 K regardless of physics -- see conversation). This is a
     ! nontrivial check whenever NumOpenChannels>1.
     !---------------------------------------------------------------------
     BLOCK
       COMPLEX*16, ALLOCATABLE :: SdagS(:,:)
       DOUBLE PRECISION :: maxdev
       INTEGER :: ii,jj
       ALLOCATE(Amat(NumOpenChannels,NumOpenChannels),Smat(NumOpenChannels,NumOpenChannels))
       ALLOCATE(SdagS(NumOpenChannels,NumOpenChannels))
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
       SdagS = MATMUL(CONJG(TRANSPOSE(Smat)),Smat)
       maxdev = 0d0
       DO ii = 1,NumOpenChannels
          DO jj = 1,NumOpenChannels
             maxdev = MAX(maxdev,ABS(SdagS(ii,jj)-MERGE((1d0,0d0),(0d0,0d0),ii.EQ.jj)))
          ENDDO
       ENDDO
       recombProb = 1d0 - ABS(Smat(1,1))**2
       WRITE(RecombFile,'(1P,E20.12,I5,E20.12,A,E20.12)') Energy, NumOpenChannels, maxdev, &
            '  (1-|S11|^2=',recombProb
       DEALLOCATE(Amat,Smat,SdagS)
     END BLOCK

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
             " (phase shift skipped -- see Recombination.dat for the multichannel unitarity check)"
     ENDIF

     DO i = 1,NumBoxes
        CALL DeAllocateBox(Boxes(i))
     ENDDO
     DEALLOCATE(Boxes)
     CALL DeAllocateBox(Bnull)
     CALL DeAllocateScat(SD)
  ENDDO

  DO iBox = 1,NumSVDBoxes
     CALL DeAllocateSVDBPD(BPDSVDBoxes(iBox))
  ENDDO
  DEALLOCATE(BPDSVDBoxes,Gam0SVDBoxes,OverlapSVDBoxes,LamSVDBoxes)
  DEALLOCATE(ODiagLSVDBoxes,ODiagRSVDBoxes,wsq1SVDBoxes,wsqLSVDBoxes)
  CALL DeAllocateEIG(EIGSVD)
  IF (NumSVDBoxes.LT.NumBoxes) THEN
     CALL DeAllocateBPD(BPD0)
     DEALLOCATE(xprim)
     CALL DeAllocateEIG(EIG)
     IF (ALLOCATED(Gam0Boxes)) DEALLOCATE(Gam0Boxes,OverlapBoxes,LamBoxes)
  ENDIF

  CLOSE(PhaseFile)
  CLOSE(KMatFile)
  CLOSE(RecombFile)
  WRITE(6,*)
  WRITE(6,*) "Phase shift written to PhaseShift.dat (xmgrace: xmgrace PhaseShift.dat)"

END PROGRAM main
