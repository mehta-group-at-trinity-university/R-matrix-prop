!****************************************************************************************************
! DIAGNOSTIC VARIANT of RMATPROPAdiabatic.f90: box 1 (R in [0,xr1], the region where the
! interpolated-tabulated-data pipeline was found to break down badly for the non-Coulomb
! channels -- see AdiabaticPotential.f90's SetAdiabaticPotential and the surrounding
! near-origin discussion) is NOT built from interpolated data at all here. Instead, its own
! R-matrix (all NumChannels retained, at its own right edge) is read in from
! 'Box1Rmatrix.dat' -- written by Adiabatic-Scattering-BoundStates/DeltaScat/DeltaScat.f90,
! which evaluates the exact KMS closed-form potential at every quadrature point, no
! interpolation, no near-origin domain limit -- one block per energy, diagonalized (DSYEV)
! to reconstruct an equivalent Zf/bf eigenchannel representation, and used to seed
! Boxes(1) directly. Boxes 2..NumBoxes (interior + last box) are built EXACTLY as in
! RMATPROPAdiabatic.f90, unchanged -- interpolated data, same CalcGamLam/
! PartitionAndEliminate/BoxMatch/CalcK pipeline. Requires NumBoxes>=2 (box 1 must be
! distinct from the last box) and the exact same Emin/Emax/NumEnergies/EnergyGridType as
! the DeltaScat run that produced Box1Rmatrix.dat (so energies line up block-for-block).
!
! Purpose: isolate whether the ~50x K-matrix discrepancy against DeltaScat found with the
! ordinary interpolated pipeline is entirely a box-1/near-origin effect. If substituting
! DeltaScat's exact box-1 R-matrix here reproduces DeltaScat's own K/delta closely, that
! confirms (a) the near-origin interpolation breakdown IS the culprit, and (b) everything
! downstream of box 1 (interior boxes, last box, CalcK) in RMATPROPAdiabatic's own pipeline
! is correct.
!****************************************************************************************************
PROGRAM main
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE AdiabaticPotential

  IMPLICIT NONE
  CHARACTER(LEN=64) DataDir
  CHARACTER(LEN=16) EnergyGridType
  INTEGER ios
  TYPE(BPData) BPD,BPD0
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: xprim(:)
  ! Cached energy-independent Gam0/Overlap/Lam for interior boxes 2..NumBoxes-1 (built
  ! once, reused every energy via CombineGam). The last box is rebuilt fresh each energy,
  ! since its NumOpenR=NumOpenChannels can change with energy.
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Boxes(:,:,:), OverlapBoxes(:,:,:), LamBoxes(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  ! Box 1 seed: read from Box1Rmatrix.dat (DeltaScat's exact box-1 R-matrix), one
  ! NumChannels x NumChannels block per energy, diagonalized via DSYEV to reconstruct an
  ! equivalent Zf/bf eigenchannel representation for Boxes(1).
  INTEGER Box1RFile, jSeed, betaSeed, dsyevInfo, dsyevLwork
  DOUBLE PRECISION seedEnergy, seedXr
  DOUBLE PRECISION, ALLOCATABLE :: Rmat1(:,:), dsyevWork(:), Rmat1Eval(:)
  CHARACTER(LEN=200) seedHeaderLine
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim
  DOUBLE PRECISION Emin, Emax, qOurs, Kelem, rawDelta, delta, prevDelta
  INTEGER NumDataPoints, NumOpenChannels, xNumPoints, LegPointsLocal
  INTEGER i, iBox, lx, kx, NumEnergies, ie, nBranch, PhaseFile, KMatFile, RecombFile
  LOGICAL havePrevDelta
  ! Full open-channel S-matrix (Cayley transform of SD%K, S=(1-iK)^-1(1+iK)), used to test
  ! McGuire's exact-elastic-unitarity result for the delta-function 3-body problem:
  ! 1-|S(1,1)|^2 (total probability the incident dimer+atom channel scatters into any other
  ! open channel -- i.e. 3-body recombination/breakup) must vanish EXACTLY for the true
  ! solution at every energy, for any number of channels. A nonzero numerical value is a
  ! finite-channel-count/finite-RMax truncation artifact that should shrink as NumChannels
  ! (and RMax, for small Energy) increase -- see AnalyticComboNearOrigin's neighbor note and
  ! McGuire, J. Math. Phys. 5, 622 (1964).
  COMPLEX*16, ALLOCATABLE :: Amat(:,:), Smat(:,:)
  DOUBLE PRECISION recombProb
  ! Box-boundary spacing power: 1.0 = linear (uniform box widths), 2.0 = quadratic
  ! (clusters more/narrower boxes near xStart). BoxEdge(i) = xStart + (xEnd-xStart)*(i/NumBoxes)^BoxSpacingPower.
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0

  !----------------------------------------------------------------------
  ! Read the run's numerical parameters, and the physics (NumChannels,
  ! alpha, ddim, muLocal, per-channel Threshold/Leff, U/P/Qtilde
  ! interpolants) from DataDir via AdiabaticPotential's ReadAdiabaticData.
  !----------------------------------------------------------------------
  OPEN(unit=7,file='RMATPROPAdiabatic.inp',status='old')
  READ(7,*)
  READ(7,*) DataDir
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

  CALL ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alpha,ddim, &
       Threshold,Leff,UInterp,PInterp,QInterp)
  reducedmass = muLocal
  AlphaFactor = alpha
  EffDim = ddim

  ! Bound-pair (two-body-threshold) channels have a genuine Coulomb-like -C/R
  ! divergence in U+Q/2mu near the origin (Mehta, Esry & Greene, PRA 76, 022711
  ! (2007), Eqs. 51-53: U(R)->gamma0/(2*mu*a*R), gamma0=-12/(3^(1/4) 2^(1/2) pi),
  ! for equal-mass three-body systems -- see kLeft_phase_shift_derivation notes).
  ! Threshold here is still the RAW (not yet 2*mu-rescaled) two-body binding
  ! energy B2=-Threshold(i); a = 1/sqrt(2*mu2*B2), mu2=0.5 for two unit masses.
  !
  ! *** NOT GENERAL -- delta-function/equal-mass-specific, needs revisiting for any
  ! other problem: *** this block silently ASSUMES any channel with Threshold<0 is
  ! this specific kind of genuine Coulomb-diverging bound-pair channel, and applies
  ! the equal-mass gamma0 formula above unconditionally. A different potential (or
  ! unequal masses) with a genuine negative-threshold channel that does NOT have
  ! this particular near-origin -C/R form would get the WRONG CoulombC here, which
  ! then silently substitutes the wrong near-origin physics via
  ! AnalyticComboNearOrigin (RMatPropCore.f90) inside SetAdiabaticPotential's
  ! R<RcutoffAnalytic override (AdiabaticPotential.f90). If porting to a new
  ! problem: either zero out CoulombC for channels that don't have this exact
  ! divergence, or replace this whole formula with whatever near-origin physics
  ! actually applies there. Left as-is here to preserve this file's validated
  ! delta-function benchmark unchanged (see RMATPROPAdiabatic.inp, DataDir=
  ! 3bodydata_delta_5ch) -- confirmed to reproduce Amaya-Tapia, Larsen & Popiel
  ! (Few-Body Systems 23, 87-109, 1997)'s own coupled-channel phase-shift figure.
  ALLOCATE(CoulombC(NumChannels))
  CoulombC = 0d0
  DO i = 1,NumChannels
     IF (Threshold(i).LT.0d0) THEN
        BLOCK
          DOUBLE PRECISION gamma0, B2, a2
          gamma0 = -12d0/(3d0**0.25d0*DSQRT(2d0)*Pi)
          B2 = -Threshold(i)
          a2 = 1d0/DSQRT(2d0*0.5d0*B2)
          CoulombC(i) = -gamma0/(2d0*reducedmass*a2)
        END BLOCK
     ENDIF
  ENDDO

  ! ReadAdiabaticData/FitLeff.data give raw threshold energies (matching Uad.dat's raw
  ! U(R) convention). CalcK's k(i)^2=2*mu*EE-Eth(i) formula (RMatPropCore.f90:493-494)
  ! expects Eth already pre-multiplied by 2*mu -- rescale once here so both the
  ! NumOpenChannels count below and the CalcK call use a consistent convention. This was
  ! silently correct for the Baluja benchmark (reducedmass=1, so 2*mu*threshold=threshold)
  ! but wrong for any reducedmass != 1.
  Threshold = 2d0*reducedmass*Threshold
  WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," ddim = ",EffDim," mu = ",reducedmass
  WRITE(6,*) "Energy scan: Emin = ",Emin," Emax = ",Emax," NumEnergies = ",NumEnergies

  !----------------------------------------------------------------------
  ! Plain-text, xmgrace-readable phase-shift output. Columns: q (this
  ! channel's own wavenumber, k(1)=sqrt(2*mu*E-Threshold(1)); NOT rescaled
  ! to any other paper's units/mass convention -- e.g. comparing to
  ! Amaya-Tapia et al.'s q requires an extra factor, see DeltaBoundValidation
  ! notes), Energy, K (the single-channel K-matrix element), delta (radians,
  ! unwrapped for continuity across the scan; no Levinson/comparison-specific
  ! offset applied -- add by hand if needed).
  !----------------------------------------------------------------------
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

  Box1RFile = 53
  OPEN(unit=Box1RFile,file='Box1Rmatrix.dat',status='old')
  READ(Box1RFile,'(A)') seedHeaderLine
  READ(Box1RFile,'(A)') seedHeaderLine
  ALLOCATE(Rmat1(NumChannels,NumChannels),Rmat1Eval(NumChannels))
  dsyevLwork = MAX(1,3*NumChannels-1)
  ALLOCATE(dsyevWork(dsyevLwork))

  !---------------------------------------------------------------------
  ! ONE-TIME setup (energy-independent): primitive basis BPD0 on [0,1], and cached
  ! Gam0/Overlap/Lam for interior boxes 2..NumBoxes-1, whose NumOpenL=NumOpenR=NumChannels
  ! never changes with energy so their matrix elements never need to be recomputed after
  ! this. Box 1 is read from Box1Rmatrix.dat fresh each energy (see the energy loop); the
  ! last box is handled fresh inside the energy loop as usual.
  !---------------------------------------------------------------------
  ! Quadratic box-boundary spacing (not just quadratic grid WITHIN box 1): boundary(i) =
  ! xStart + (xEnd-xStart)*(i/NumBoxes)^2, clustering more/narrower boxes near xStart so
  ! resolution is concentrated where the delta-function channel's U/Q curves vary fastest.
  ! (BoxEdge(i), defined below at each use site, replaces the old uniform xStart+i*xDelt.)

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

  BPD = BPD0  ! workhorse set actually used for boxes 2..NumBoxes

  EIG%MatrixDim = BPD%MatrixDim
  CALL AllocateEIG(EIG)

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
        CALL SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp,CoulombC)
        ! interior box-to-box boundary: all channels retained, always -- energy-independent
        CALL CalcGamLam(BPD,EIG,NumChannels,NumChannels)
        Gam0Boxes(:,:,iBox) = EIG%Gam0
        OverlapBoxes(:,:,iBox) = EIG%Overlap
        LamBoxes(:,:,iBox) = EIG%Lam
     ENDDO
  ENDIF

  IF (NumBoxes.LT.2) THEN
     WRITE(6,*) "RMATPROPAdiabatic_seeded: requires NumBoxes>=2 (box 1 must be distinct ", &
          "from the last box) -- got NumBoxes=",NumBoxes
     STOP
  ENDIF

  havePrevDelta = .FALSE.

  DO ie = 1,NumEnergies
     IF (NumEnergies.GT.1) THEN
        IF (EnergyGridType.EQ."quadratic") THEN
           ! Clusters points near Emin (i.e. near threshold, q~0) -- same
           ! formula as GridMaker's "quadratic" R-grid scale, applied here to
           ! Energy directly.
           Energy = Emin + (Emax-Emin)*(DBLE(ie-1)/DBLE(NumEnergies-1))**2
        ELSE
           Energy = Emin + (Emax-Emin)*DBLE(ie-1)/DBLE(NumEnergies-1)
        ENDIF
     ELSE
        Energy = Emin
     ENDIF

     !---------------------------------------------------------------------
     ! Fixed open-channel count for THIS energy (mirrors CABA.EXT.f's
     ! NumOpenChannels/RMATPROP2016.f90's Baluja main: chosen once per
     ! energy, not per-box). Channels are assumed ordered open-before-closed
     ! (same convention as CABA/FitLeff.data). Recomputed fresh every energy
     ! since it can change which channels the last box's NumOpenR retains.
     !---------------------------------------------------------------------
     NumOpenChannels = 0
     DO i = 1,NumChannels
        IF (2d0*reducedmass*Energy.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
     ENDDO

     !---------------------------------------------------------------------
     ! Regularity start at xMin: box 1 has NumOpenL=0 (no left-boundary DOF
     ! retained), so BoxMatch's loops over "the previous box" are all
     ! vacuous for it -- Bnull is a trivial, never-dereferenced placeholder.
     ! Boxes(:) itself is cheap (no integration) so it's still rebuilt fresh
     ! every energy -- only the last box's NumOpenR can actually change.
     !---------------------------------------------------------------------
     ALLOCATE(Boxes(NumBoxes))
     DO i = 1,NumBoxes
        Boxes(i)%NumOpenL = MERGE(0,NumChannels,i.EQ.1)
        IF (i.EQ.NumBoxes) THEN
           Boxes(i)%NumOpenR = NumOpenChannels  ! final box: only physically-open channels retained
        ELSE
           Boxes(i)%NumOpenR = NumChannels      ! interior box-to-box boundary: all channels retained
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

     !---------------------------------------------------------------------
     ! Box 1: read DeltaScat's exact R-matrix for THIS energy from Box1Rmatrix.dat
     ! (written block-by-block in the same energy order this scan uses), diagonalize it
     ! (R = V*diag(Lambda)*V^T via DSYEV), and reconstruct an equivalent Zf/bf eigenchannel
     ! representation: bf(beta)=1/Lambda(beta), Zf(i,beta)=V(i,beta), Zfp=-Zf*bf (BoxMatch's
     ! own normalization convention -- see RMatPropCore.f90). This exactly reproduces
     ! Rmat1(i,j)=sum_beta Zf(i,beta)*Zf(j,beta)/bf(beta), so Boxes(1) is a valid substitute
     ! for whatever BoxMatch(Bnull,Boxes(1),...) would have produced, without ever touching
     ! the interpolated near-origin potential.
     !---------------------------------------------------------------------
     READ(Box1RFile,'(A)') seedHeaderLine
     READ(seedHeaderLine(12:31),*) seedEnergy
     READ(seedHeaderLine(39:58),*) seedXr
     IF (DABS(seedEnergy-Energy).GT.1d-9*MAX(1d0,DABS(Energy))) THEN
        WRITE(6,*) "RMATPROPAdiabatic_seeded: energy grid mismatch -- seed file has ", &
             seedEnergy," this run wants ",Energy,". Was Box1Rmatrix.dat generated with ", &
             "the same Emin/Emax/NumEnergies/EnergyGridType as RMATPROPAdiabatic.inp?"
        STOP
     ENDIF
     IF (DABS(seedXr-Boxes(1)%xr).GT.1d-9*MAX(1d0,DABS(Boxes(1)%xr))) THEN
        WRITE(6,*) "RMATPROPAdiabatic_seeded: box-1 edge mismatch -- seed file xr=", &
             seedXr," this run's Boxes(1)%xr=",Boxes(1)%xr,". Was Box1Rmatrix.dat generated ", &
             "with the same NumBoxes/xStart/xEnd/BoxSpacingPower?"
        STOP
     ENDIF
     DO i = 1,NumChannels
        READ(Box1RFile,*) Rmat1(i,1:NumChannels)
     ENDDO

     CALL DSYEV('V','U',NumChannels,Rmat1,NumChannels,Rmat1Eval,dsyevWork,dsyevLwork,dsyevInfo)
     IF (dsyevInfo.NE.0) THEN
        WRITE(6,*) "RMATPROPAdiabatic_seeded: DSYEV failed, info=",dsyevInfo," at Energy=",Energy
        STOP
     ENDIF
     ! Rmat1 now holds the eigenvectors V (DSYEV overwrites it in place); Rmat1Eval holds Lambda.
     DO betaSeed = 1,NumChannels
        Boxes(1)%bf(betaSeed) = 1d0/Rmat1Eval(betaSeed)
        DO jSeed = 1,NumChannels
           Boxes(1)%Zf(jSeed,betaSeed) = Rmat1(jSeed,betaSeed)
           Boxes(1)%Zfp(jSeed,betaSeed) = -Boxes(1)%Zf(jSeed,betaSeed)*Boxes(1)%bf(betaSeed)
        ENDDO
     ENDDO

     !---------------------------------------------------------------------
     ! Interior boxes 2..NumBoxes-1: reuse the cached Gam0/Overlap/Lam built
     ! once above -- just the cheap CombineGam recombination plus the linear
     ! algebra (PartitionAndEliminate/BoxMatch), no re-integration.
     !---------------------------------------------------------------------
     DO iBox = 2,NumBoxes-1
        ! BoxMatch uses BPD%xl/xr directly (R-power normalization) -- refresh to this
        ! box's own edges every iteration (see the matching fix/comment in
        ! RMATPROPAdiabatic.f90; numerically inert for this dataset but fixes a real bug).
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

     !---------------------------------------------------------------------
     ! Last box (iBox=NumBoxes, only when NumBoxes>=2): rebuilt fresh every
     ! energy, since its NumOpenR=NumOpenChannels can change with energy.
     !---------------------------------------------------------------------
     IF (NumBoxes.GE.2) THEN
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
     ENDIF

     CALL CalcK(Boxes(NumBoxes),BPD,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)

     !---------------------------------------------------------------------
     ! Diagnostic: full K-matrix dump for the all-channels-open diagonal-
     ! asymptotic test (no cross-channel coupling by construction, so K
     ! should come out exactly diagonal with each diagonal eigenphase
     ! atan(K_nn) equal to zero).
     !---------------------------------------------------------------------
     IF (NumOpenChannels.EQ.NumChannels .AND. NumChannels.GT.1) THEN
        WRITE(KMatFile,'(A,1P,E20.12)') '# Energy = ',Energy
        DO i = 1,NumOpenChannels
           WRITE(KMatFile,'(1P,20E20.12)') SD%K(i,1:NumOpenChannels)
        ENDDO
     ENDIF

     !---------------------------------------------------------------------
     ! McGuire elastic-unitarity test: S = (1-iK)^-1(1+iK) (Cayley transform of the
     ! real, energy-normalized, symmetric open-channel K-matrix SD%K, dimension
     ! NumOpenChannels x NumOpenChannels -- see CalcK), then 1-|S(1,1)|^2. Channel 1 is
     ! always the dimer+atom elastic channel here, so this is exactly the incident-
     ! channel recombination/breakup probability. Computed at every energy (trivially
     ! ~0 to machine precision when NumOpenChannels=1, which is itself a useful sanity
     ! check on this Cayley-transform code before trusting it in the genuinely
     ! multi-open-channel regime above the 3-body threshold).
     !---------------------------------------------------------------------
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

     !---------------------------------------------------------------------
     ! Phase shift only has a plain scalar meaning with exactly one open
     ! channel (our physical regime here: energy between the 2-body and
     ! 3-body thresholds). Other energies in the scan are still fully
     ! propagated above (so the R-matrix/K-matrix printout is always valid),
     ! just not written to PhaseShift.dat.
     !---------------------------------------------------------------------
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

     !---------------------------------------------------------------------
     ! Per-energy teardown (Fortran runtime bounds-checking will abort on a
     ! re-ALLOCATE of an already-allocated variable, so anything reallocated
     ! fresh next iteration -- Boxes, Bnull, SD -- must be freed here).
     ! BPD0/BPD/EIG/xprim/Gam0Boxes/OverlapBoxes/LamBoxes now live for the
     ! whole run (built once before the energy loop) and are torn down once
     ! after it below.
     !---------------------------------------------------------------------
     DO i = 1,NumBoxes
        CALL DeAllocateBox(Boxes(i))
     ENDDO
     DEALLOCATE(Boxes)
     CALL DeAllocateBox(Bnull)
     CALL DeAllocateScat(SD)
  ENDDO

  CALL DeAllocateBPD(BPD0)
  DEALLOCATE(xprim)
  CALL DeAllocateEIG(EIG)
  IF (ALLOCATED(Gam0Boxes)) DEALLOCATE(Gam0Boxes,OverlapBoxes,LamBoxes)

  CLOSE(PhaseFile)
  CLOSE(KMatFile)
  CLOSE(RecombFile)
  CLOSE(Box1RFile)
  WRITE(6,*)
  WRITE(6,*) "Phase shift written to PhaseShift.dat (xmgrace: xmgrace PhaseShift.dat)"

END PROGRAM main
