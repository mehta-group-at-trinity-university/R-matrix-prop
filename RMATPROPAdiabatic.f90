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
  ! Gates the equal-mass delta-function-specific CoulombC/gamma0 near-origin correction
  ! below (see that block's own header comment for the full "NOT GENERAL" rationale) --
  ! defaults to .FALSE. (optional trailing .inp line, absent = off) so any dataset whose
  ! Threshold(i)<0 does NOT mean "genuine Coulomb-diverging bound-pair channel" (e.g. a
  ! short-range sech^2-well bound-dimer channel, where the near-origin U+Q/2mu behavior
  ! is regular, not -C/R) doesn't silently get the wrong near-origin physics. Only
  ! RMATPROPAdiabatic.inp (3bodydata_delta_5ch) sets this .TRUE., to keep that validated
  ! benchmark unchanged.
  LOGICAL ApplyCoulombC
  TYPE(BPData) BPD,BPD0,BPD1
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: xprim(:)
  ! Cached energy-independent Gam0/Overlap/Lam for interior boxes 2..NumBoxes-1 (built
  ! once, reused every energy via CombineGam). The last box is rebuilt fresh each energy,
  ! since its NumOpenR=NumOpenChannels can change with energy.
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Boxes(:,:,:), OverlapBoxes(:,:,:), LamBoxes(:,:,:)
  ! Cached energy-independent Gam0/Overlap/Lam for box 1 (only when NumBoxes>=2, i.e. box 1
  ! is NOT simultaneously the last box -- see the caching block before the energy loop and
  ! its use inside the loop for why NumBoxes==1 can't use this and must rebuild box 1 fresh
  ! every energy instead).
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Box1(:,:), OverlapBox1(:,:), LamBox1(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, EffDimLocal
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
  ! alpha, EffDimLocal, muLocal, per-channel Threshold/Leff, U/P/Qtilde
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
  READ(7,*,IOSTAT=ios) ApplyCoulombC
  IF (ios.NE.0) ApplyCoulombC = .FALSE.
  CLOSE(7)

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = LegPointsLocal
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  ! Builds B-spline interpolants for U_n(R) (diagonal), P_mn(R) (first-derivative coupling),
  ! Qtilde_mn(R) (second-derivative coupling) from the tabulated Uad/Pmat/Qmat.dat.
  CALL ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alpha,EffDimLocal, &
       Threshold,Leff,UInterp,PInterp,QInterp)
  reducedmass = muLocal
  ! AlphaFactor is forced to 0 rather than taking the dataset's own alpha (from Fit.data,
  ! typically (EffDim-1)/2=0.5 for EffDim=2): the wavefunction-reduction term this produces,
  ! AlphaFactor*(AlphaFactor-EffDim+2)/(2mu R^2), sits exactly at the critical/marginal
  ! inverse-square coupling (T=1/4, Frobenius index nu=0) and is mishandled somewhere
  ! downstream of provably-exact matrix elements (see memory wavefn-reduction-term-bug --
  ! root cause still open, localized to PartitionAndEliminate/BoxMatch). AlphaFactor=0 makes
  ! the term vanish identically regardless of EffDim and was validated to ~1e-7 relative
  ! against independent NDSolve, both single-box (TestAlphaZero.f90) and through this
  ! driver's full box-chained pipeline. Every dataset in this repo either already has
  ! alpha=0 in Fit.data (3bodydata_delta_5ch -- never exposed to the bug, which is why it
  ! validates against Amaya-Tapia without needing this override) or alpha=(EffDim-1)/2 at
  ! the exact critical coefficient (every other dataset) -- there is currently no dataset
  ! for which alpha=(EffDim-1)/2 is the correct choice, so this is unconditional, not a flag.
  IF (alpha.NE.0d0) THEN
     WRITE(6,*) "NOTE: overriding dataset's own alpha=",alpha," with AlphaFactor=0 -- see", &
          " this line's own comment / memory wavefn-reduction-term-bug."
  ENDIF
  AlphaFactor = 0d0
  EffDim = EffDimLocal

  ! Bound-pair (two-body-threshold) channels have a genuine Coulomb-like -C/R
  ! divergence in U+Q/2mu near the origin (Mehta, Esry & Greene, PRA 76, 022711
  ! (2007), Eqs. 51-53: U(R)->gamma0/(2*mu*a*R), gamma0=-12/(3^(1/4) 2^(1/2) pi),
  ! for equal-mass three-body systems -- see kLeft_phase_shift_derivation notes).
  ! Threshold here is still the RAW (not yet 2*mu-rescaled) two-body binding
  ! energy B2=-Threshold(i); a = 1/sqrt(2*mu2*B2), mu2=0.5 for two unit masses.
  !
  ! *** NOT GENERAL -- delta-function/equal-mass-specific. *** Naively keying this off
  ! Threshold<0 would silently ASSUME any negative-threshold channel is this specific
  ! kind of genuine Coulomb-diverging bound-pair channel, and apply the equal-mass
  ! gamma0 formula above unconditionally -- WRONG for a different potential (or unequal
  ! masses) whose negative-threshold channel does NOT have this particular near-origin
  ! -C/R form (e.g. a short-range sech^2-well bound-dimer channel: Threshold<0 there too,
  ! but U+Q/2mu is regular at the origin, not -C/R divergent). This block would then
  ! silently substitute the wrong near-origin physics via AnalyticComboNearOrigin
  ! (RMatPropCore.f90) inside SetAdiabaticPotential's R<RcutoffAnalytic override
  ! (AdiabaticPotential.f90). Gated behind ApplyCoulombC (read from the .inp file,
  ! defaults to .FALSE.) instead of inferring it from Threshold's sign -- only
  ! RMATPROPAdiabatic.inp (DataDir=3bodydata_delta_5ch) sets it .TRUE., to preserve this
  ! file's validated delta-function benchmark unchanged -- confirmed to reproduce
  ! Amaya-Tapia, Larsen & Popiel (Few-Body Systems 23, 87-109, 1997)'s own coupled-channel
  ! phase-shift figure. New datasets with a genuine negative-threshold channel that ISN'T
  ! this exact divergence (e.g. 3bodydata_sech2_5ch) simply leave ApplyCoulombC unset.
  ALLOCATE(CoulombC(NumChannels))
  CoulombC = 0d0
  IF (ApplyCoulombC) THEN
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
  ENDIF

  ! ReadAdiabaticData/FitLeff.data give raw threshold energies (matching Uad.dat's raw
  ! U(R) convention). CalcK's k(i)^2=2*mu*EE-Eth(i) formula (RMatPropCore.f90:493-494)
  ! expects Eth already pre-multiplied by 2*mu -- rescale once here so both the
  ! NumOpenChannels count below and the CalcK call use a consistent convention. This was
  ! silently correct for the Baluja benchmark (reducedmass=1, so 2*mu*threshold=threshold)
  ! but wrong for any reducedmass != 1.
  Threshold = 2d0*reducedmass*Threshold
  WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," EffDim = ",EffDim," mu = ",reducedmass
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


  !---------------------------------------------------------------------
  ! ONE-TIME setup (energy-independent): primitive basis BPD0 on [0,1],
  ! BPD1's allocation (its geometry/basis is rebuilt fresh each energy inside
  ! the loop -- see note there), and cached Gam0/Overlap/Lam for interior
  ! boxes 2..NumBoxes-1, whose NumOpenL=NumOpenR=NumChannels never changes
  ! with energy so their matrix elements never need to be recomputed after
  ! this. Box 1 and the last box are handled fresh inside the energy loop.
  !---------------------------------------------------------------------
  ! Box boundaries are boundary(i) = xStart + (xEnd-xStart)*(i/NumBoxes)**BoxSpacingPower --
  ! uniform at the BoxSpacingPower=1 set above; change BoxSpacingPower for a different box
  ! grid (>1 clusters narrower boxes near xStart).  Distinct from the quadratic grid used
  ! WITHIN box 1 by GridMaker.

  BPD1%NumChannels = NumChannels
  BPD1%Order = Order
  ! Left=2 (both B1,B2 independently retained at xMin), NumOpenL still 0
  ! (Boxes(1)%NumOpenL=0 below): the ordinary CalcGamLam/PartitionAndEliminate
  ! elimination mechanism (already exercised for every other box's closed
  ! channels) handles box 1's boundary directly, in place of the earlier
  ! specialized BuildBox1CombinedBasis/CalcGamLamBox1 Robin-matched basis.
  ! Confirmed equivalent to that older approach (and to a plain Neumann box-1
  ! basis) once the near-origin potential/BC pipeline was fixed -- this is now
  ! the permanent implementation, not a variant under test.
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
        ! Fills Pot(m,n)=Qtilde_mn(R)/(2*mu)+delta_mn*U_m(R) and Pcoup=P(R) at this box's
        ! quadrature points, from the interpolants ReadAdiabaticData built above.
        CALL SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp,CoulombC)
        ! interior box-to-box boundary: all channels retained, always -- energy-independent.
        ! Builds Gam0_ij=Int[B_i'B_j' - 2*mu*Pot*B_iB_j - alpha(alpha-D+2)/(2*mu*R^2)*B_iB_j]
        ! R^(D-1-2alpha)dR (+ Pcoup coupling term) and Overlap_ij=Int[B_iB_j]R^(D-1-2alpha)dR.
        CALL CalcGamLam(BPD,EIG,NumChannels,NumChannels)
        Gam0Boxes(:,:,iBox) = EIG%Gam0
        OverlapBoxes(:,:,iBox) = EIG%Overlap
        LamBoxes(:,:,iBox) = EIG%Lam
     ENDDO
  ENDIF

  !---------------------------------------------------------------------
  ! Box 1's basis (geometry, MakeBasis, SetAdiabaticPotential, CalcGamLam) is
  ! energy-independent -- same as the interior boxes above -- now that it uses the
  ! ordinary Left=2/NumOpenL=0 elimination mechanism instead of the old energy-dependent
  ! Robin-matched basis. The only exception is NumBoxes==1, where box 1 doubles as the
  ! last box and its NumOpenR=NumOpenChannels genuinely changes with energy (see the
  ! NumBoxes==1 branch inside the energy loop, which rebuilds box 1 fresh every energy
  ! exactly as before).
  !---------------------------------------------------------------------
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
     CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)
     CALL CalcGamLam(BPD1,EIG,0,NumChannels)
     Gam0Box1 = EIG%Gam0
     OverlapBox1 = EIG%Overlap
     LamBox1 = EIG%Lam
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
     ! Box 1: NumBoxes>=2 reuses the Gam0Box1/OverlapBox1/LamBox1 cached once before the
     ! energy loop (BPD1 itself is also already allocated/built there and stays live for
     ! the whole run) -- just the cheap CombineGam recombination plus the linear algebra
     ! (PartitionAndEliminate/BoxMatch), exactly mirroring the interior-box pattern below.
     ! NumBoxes==1 is the one case that can't be cached: box 1 doubles as the last box, so
     ! its NumOpenR=NumOpenChannels genuinely changes with energy and it must be rebuilt
     ! fully fresh every time (matching the last-box branch's own NumOpenR-changes-with-
     ! energy reasoning).
     !---------------------------------------------------------------------
     IF (NumBoxes.GE.2) THEN
        EIG%Gam0 = Gam0Box1
        EIG%Overlap = OverlapBox1
        EIG%Lam = LamBox1
        ! Gam(E) = Gam0 + E*Overlap -- the actual energy-dependent generalized-eigenvalue
        ! matrix PartitionAndEliminate needs, from the energy-independent Gam0/Overlap above.
        CALL CombineGam(EIG,Energy)
        ALLOCATE(evalRed(Boxes(1)%betaMax),evecRed(Boxes(1)%betaMax,Boxes(1)%betaMax))
        ALLOCATE(LamooDiag(Boxes(1)%betaMax))
        ! Schur-eliminates the interior/closed DOF (Omega_oo = Gam_oo - Gam_oc*Gam_cc^-1*Gam_co),
        ! then solves the small generalized eigenproblem Omega_oo*c = beta*Lam_oo*c on the
        ! boundary-retained ("open") DOF only -- evalRed=beta, evecRed=c.
        CALL PartitionAndEliminate(BPD1,EIG,Boxes(1)%NumOpenL,Boxes(1)%NumOpenR,evalRed,evecRed,LamooDiag)
        ! Builds this box's normalized log-derivative amplitudes Z,Z'=-beta*Z and matches them
        ! to the previous box (here: the trivial regularity start) via a generalized eigenvalue
        ! match enforcing R-matrix/log-derivative continuity (Burke thesis Eqs. 3.13-3.14).
        CALL BoxMatch(Bnull, Boxes(1), BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        DEALLOCATE(evalRed,evecRed,LamooDiag)
     ELSE
        CALL AllocateBPD(BPD1)
        BPD1%xl = Boxes(1)%xl
        BPD1%xr = Boxes(1)%xr
        CALL GridMaker(BPD1%xPoints,BPD1%xNumPoints,BPD1%xl,BPD1%xr,"quadratic")
        CALL MakeBasis(BPD1)
        BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
        CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)
        CALL CalcGamLam(BPD1,EIG,Boxes(1)%NumOpenL,Boxes(1)%NumOpenR)
        CALL CombineGam(EIG,Energy)
        ALLOCATE(evalRed(Boxes(1)%betaMax),evecRed(Boxes(1)%betaMax,Boxes(1)%betaMax))
        ALLOCATE(LamooDiag(Boxes(1)%betaMax))
        CALL PartitionAndEliminate(BPD1,EIG,Boxes(1)%NumOpenL,Boxes(1)%NumOpenR,evalRed,evecRed,LamooDiag)
        CALL BoxMatch(Bnull, Boxes(1), BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        DEALLOCATE(evalRed,evecRed,LamooDiag)
     ENDIF
     ! For NumBoxes==1, BPD1 stays allocated past here -- CalcK still needs BPD1%xr/%lam
     ! below (see note there); deallocated after that call instead. For NumBoxes>=2, BPD1
     ! is persistent (built once before the loop) and is NOT deallocated per-energy at all.

     !---------------------------------------------------------------------
     ! Interior boxes 2..NumBoxes-1: reuse the cached Gam0/Overlap/Lam built
     ! once above -- just the cheap CombineGam recombination plus the linear
     ! algebra (PartitionAndEliminate/BoxMatch), no re-integration.
     !---------------------------------------------------------------------
     DO iBox = 2,NumBoxes-1
        ! BoxMatch uses BPD%xl/xr directly (R-power normalization) -- BPD itself is only
        ! reused here for its energy-independent MatrixDim/Order/NumChannels structure
        ! (PartitionAndEliminate doesn't touch %xl/%xr at all), but BoxMatch does, so these
        ! must be refreshed to THIS box's own edges every iteration -- cheap (two scalar
        ! assignments), fixes a bug where they were left stuck at whatever box the
        ! one-time pre-loop cache-build above last visited (iBox=NumBoxes-1), so every
        ! interior box's BoxMatch call was normalizing against the wrong box edges.
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

     ! NumBoxes=1: there is no separate "last box" (box 1 itself is the only box),
     ! so BPD%xr is never set to the physical box edge -- Box 1's R-matrix was actually
     ! produced using BPD1 (regularity basis), so CalcK must be matched using BPD1's own
     ! xr, not BPD's (which would silently be 0, from BPD0's own [0,1] construction,
     ! feeding hyperrjry an x=k*0=0 argument and crashing bessjy).
     ! Matches Rmat(i,j)=sum_beta Zf(i,beta)Zf(j,beta)/bf(beta) at xEnd to asymptotic
     ! Riccati-Bessel-like reference functions s,c (order=Leff): K=(c+cp*Rmat)*(s+sp*Rmat)^-1.
     IF (NumBoxes.EQ.1) THEN
        CALL CalcK(Boxes(NumBoxes),BPD1,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
        CALL DeAllocateBPD(BPD1)  ! rebuilt fresh every energy in this branch -- see above
     ELSE
        CALL CalcK(Boxes(NumBoxes),BPD,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
        ! BPD1 is persistent here (built once before the loop) -- torn down once after the
        ! loop instead, alongside BPD0/EIG/Gam0Boxes.
     ENDIF

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
  IF (NumBoxes.GE.2) THEN
     CALL DeAllocateBPD(BPD1)
     DEALLOCATE(Gam0Box1,OverlapBox1,LamBox1)
  ENDIF
  DEALLOCATE(xprim)
  CALL DeAllocateEIG(EIG)
  IF (ALLOCATED(Gam0Boxes)) DEALLOCATE(Gam0Boxes,OverlapBoxes,LamBoxes)

  CLOSE(PhaseFile)
  CLOSE(KMatFile)
  CLOSE(RecombFile)
  WRITE(6,*)
  WRITE(6,*) "Phase shift written to PhaseShift.dat (xmgrace: xmgrace PhaseShift.dat)"

END PROGRAM main
