!****************************************************************************************************
! R-matrix propagation driven by real adiabatic hyperspherical data (P/Q-tilde/U from
! Adiabatic-Scattering-BoundStates' adiabatic solvers) instead of the synthetic Baluja
! test potential in RMATPROP2016.f90. Reuses RMatPropCore.f90's CalcGamLam/
! PartitionAndEliminate/BoxMatch/CalcK unchanged, and AdiabaticPotential.f90 for data
! import + interpolation. Structurally mirrors RMATPROP2016.f90's Baluja-benchmark main
! (box loop, PartitionAndEliminate, BoxMatch, final CalcK), but starts propagation at
! xMin with a regularity condition (Boxes(1)%NumOpenL=0, no retained left DOF, no
! literature-R-matrix seed needed) instead of a diagonalized literature-R-matrix seed --
! this pattern (and the BPD1/Left=0-vs-BPD0/Left=2 basis split it requires) is carried
! over from an earlier prototype, ~/Documents/GitHub/Adiabatic-R-Mat-Prop/RMATPROP2016.f90,
! which validates the same regularity-start approach; unlike that prototype, Leff/
! Threshold here come from FitLeff.data rather than lam=0, so CalcK's asymptotic Bessel
! matching uses the correct generalized angular momentum.
!
! Runs an ENERGY SCAN from Emin to Emax (NumEnergies points, linearly spaced), redoing the
! full box/basis setup and propagation at each energy (the basis itself doesn't depend on
! Energy, but redoing it is simple and correct, and the per-point cost is small). For any
! energy where exactly one channel is open, computes the single-channel phase shift
! delta=atan(K), unwrapped for continuity across the scan (branch chosen closest to the
! previous written point; no absolute anchor/offset is imposed -- the physical bound-state-
! count-based Levinson anchor, or any comparison-specific convention offset, is left to the
! user/analysis step), and writes it to a plain-text, xmgrace-readable file.
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
  TYPE(BPData) BPD,BPD0,BPD1
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: xprim(:)
  ! Cached energy-independent Gam0/Overlap/Lam for interior boxes 2..NumBoxes-1 (built
  ! once, reused every energy via CombineGam). Box 1 and the last box are rebuilt fresh
  ! each energy: box 1 because Part B's log-derivative BC will be energy-dependent, the
  ! last box because NumOpenR=NumOpenChannels can change with energy.
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Boxes(:,:,:), OverlapBoxes(:,:,:), LamBoxes(:,:,:)
  ! Box 1's per-channel, energy-dependent log-derivative-combined boundary basis (see
  ! BuildBox1CombinedBasis/CalcGamLamBox1 in RMatPropCore.f90).
  DOUBLE PRECISION, ALLOCATABLE :: u1Box1(:,:,:), ux1Box1(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: kLeftBox1(:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim
  DOUBLE PRECISION Emin, Emax, qOurs, Kelem, rawDelta, delta, prevDelta
  INTEGER NumDataPoints, NumOpenChannels, xNumPoints, LegPointsLocal
  INTEGER i, iBox, lx, kx, NumEnergies, ie, nBranch, PhaseFile, KMatFile
  LOGICAL havePrevDelta
  ! TEMP diagnostic toggle: use the simple, energy-independent Neumann-sum box-1 basis
  ! (mirroring 2006/2007 adrmatprop.f's Left=1) instead of the Robin/Coulomb-matched one.
  ! Intended for use with a very small xMin (near the smallest tabulated R).
  LOGICAL, PARAMETER :: UseNeumannBox1 = .TRUE.
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

  !---------------------------------------------------------------------
  ! ONE-TIME setup (energy-independent): primitive basis BPD0 on [0,1],
  ! BPD1's allocation (its geometry/basis is rebuilt fresh each energy inside
  ! the loop -- see note there), and cached Gam0/Overlap/Lam for interior
  ! boxes 2..NumBoxes-1, whose NumOpenL=NumOpenR=NumChannels never changes
  ! with energy so their matrix elements never need to be recomputed after
  ! this. Box 1 and the last box are handled fresh inside the energy loop.
  !---------------------------------------------------------------------
  ! Quadratic box-boundary spacing (not just quadratic grid WITHIN box 1): boundary(i) =
  ! xStart + (xEnd-xStart)*(i/NumBoxes)^2, clustering more/narrower boxes near xStart so
  ! resolution is concentrated where the delta-function channel's U/Q curves vary fastest.
  ! (BoxEdge(i), defined below at each use site, replaces the old uniform xStart+i*xDelt.)

  BPD1%NumChannels = NumChannels
  BPD1%Order = Order
  ! TEST (per user request): Left=2 (both B1,B2 independently retained at xMin),
  ! NumOpenL still 0 (Boxes(1)%NumOpenL=0 below) -- letting the ordinary
  ! CalcGamLam/PartitionAndEliminate elimination mechanism (already exercised for
  ! every other box's closed channels) handle box 1's boundary instead of the
  ! specialized BuildBox1CombinedBasis/CalcGamLamBox1 Robin-matched basis. Earlier
  ! note claimed this is a no-op (basis change within an already-eliminated
  ! subspace) -- re-verifying now that the near-origin potential/BC pipeline has
  ! been substantially fixed since that claim was made.
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
     ! Box 1: rebuilt fully fresh every energy (regularity basis at xMin).
     ! Its left-boundary (ix=1) basis function is, per channel, an exact
     ! log-derivative (Robin) combination of Left=2's two independent boundary
     ! functions -- BuildBox1CombinedBasis/CalcGamLamBox1 -- matched analytically
     ! to the true regular solution at R=xMin, which is genuinely energy-dependent
     ! through k(channel) and so must be rebuilt every energy.
     !---------------------------------------------------------------------
     CALL AllocateBPD(BPD1)
     BPD1%xl = Boxes(1)%xl
     BPD1%xr = Boxes(1)%xr
     CALL GridMaker(BPD1%xPoints,BPD1%xNumPoints,BPD1%xl,BPD1%xr,"quadratic")
     CALL MakeBasis(BPD1)
     BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
     CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)
     ! TEST: ordinary Left=2 basis, NumOpenL=0 elimination via the generic
     ! CalcGamLam/PartitionAndEliminate mechanism -- no special box-1 Robin/Neumann
     ! basis construction at all (see note at BPD1%Left assignment above).
     CALL CalcGamLam(BPD1,EIG,Boxes(1)%NumOpenL,Boxes(1)%NumOpenR)
     CALL CombineGam(EIG,Energy)
     ALLOCATE(evalRed(Boxes(1)%betaMax),evecRed(Boxes(1)%betaMax,Boxes(1)%betaMax))
     ALLOCATE(LamooDiag(Boxes(1)%betaMax))
     CALL PartitionAndEliminate(BPD1,EIG,Boxes(1)%NumOpenL,Boxes(1)%NumOpenR,evalRed,evecRed,LamooDiag)
     CALL BoxMatch(Bnull, Boxes(1), BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
     DEALLOCATE(evalRed,evecRed,LamooDiag)
     ! BPD1 stays allocated past here -- CalcK still needs BPD1%xr/%lam below when
     ! NumBoxes==1 (see note there); deallocated after that call instead.

     !---------------------------------------------------------------------
     ! Interior boxes 2..NumBoxes-1: reuse the cached Gam0/Overlap/Lam built
     ! once above -- just the cheap CombineGam recombination plus the linear
     ! algebra (PartitionAndEliminate/BoxMatch), no re-integration.
     !---------------------------------------------------------------------
     DO iBox = 2,NumBoxes-1
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

     ! NumBoxes=1: there is no separate "last box" (box 1 itself is the only box),
     ! so BPD%xr is never set to the physical box edge -- Box 1's R-matrix was actually
     ! produced using BPD1 (regularity basis), so CalcK must be matched using BPD1's own
     ! xr, not BPD's (which would silently be 0, from BPD0's own [0,1] construction,
     ! feeding hyperrjry an x=k*0=0 argument and crashing bessjy).
     IF (NumBoxes.EQ.1) THEN
        CALL CalcK(Boxes(NumBoxes),BPD1,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
     ELSE
        CALL CalcK(Boxes(NumBoxes),BPD,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
     ENDIF
     CALL DeAllocateBPD(BPD1)

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
     ! after it below; BPD1 is already deallocated above (rebuilt every
     ! energy, per the box-1 comment there).
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
  WRITE(6,*)
  WRITE(6,*) "Phase shift written to PhaseShift.dat (xmgrace: xmgrace PhaseShift.dat)"

END PROGRAM main
