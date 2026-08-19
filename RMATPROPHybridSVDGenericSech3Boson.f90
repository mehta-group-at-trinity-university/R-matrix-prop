!****************************************************************************************************
! All-SVD hybrid driver for the equal-mass THREE-IDENTICAL-BOSON sech^2-well problem
! (3bodydata_sech2_3boson / _10ch / _20ch), the genuinely independent numerical cross-check for
! RMATPROPAdiabatic.x's tabulated-data baseline on that dataset -- mirrors what
! RMATPROPHybridSVDGenericDelta.f90 already did for the delta-function benchmark, and what
! RMATPROPHybridSVDGenericSech.f90 does for the DIFFERENT (2-fermion+1-distinguishable, domain
! [0,pi/2]) sech^2 case -- do not confuse the two. Every box is SVD, re-diagonalizing the true
! sech^2 potential at every quadrature R via ARPACK (AdiabaticSolverGeneric.f90) instead of
! interpolating a tabulated Uad(R)/P(R)/Q(R).
!
! Domain [0,pi/6] (3 identical bosons' fundamental domain), Left1D=Right1D=1 (Neumann-Neumann --
! matches 3bodydata_sech2_3boson*/DriveAdiabaticSolver1D_3boson.par.used's own Order/Left/Right=
! 5/1/1 exactly, the same angular BC used to generate the tabulated baseline this run is compared
! against). SechGridMaker (SechFunctionPlugins.f90) already handles a domain whose right edge
! coincides with the r23=0 coalescence point phi23=pi/6 (its xMax.LE.phi23 branch) -- no separate
! "wedge" grid maker needed.
!
! AlphaFactor is hardcoded to 0, NOT read from the dataset's own Fit.data (alpha=0.5=(EffDim-1)/2)
! -- see memory wavefn-reduction-term-bug and RMATPROPAdiabatic.f90's own comment at the same
! override: AlphaFactor=(EffDim-1)/2 sits exactly at the critical/marginal inverse-square coupling
! for the wavefunction-reduction term and is mishandled downstream of provably-exact matrix
! elements. Box1GridType='Radau' is the SVD-side mechanism matching AlphaFactor=0 (open/regularity,
! not Dirichlet) -- see RMATPROPHybridSVDGeneric.f90's own header table. This mirrors
! RMATPROPHybridSVDGenericDelta.f90, which already used AlphaFactor=0/Radau for its own (unrelated)
! reduced-coordinate reasons -- same mechanism, different reason here.
!
! NumChannels/mu/Threshold/Leff/UInterp/PInterp/QInterp are read directly from DataDir via
! AdiabaticPotential.f90's ReadAdiabaticData (the SAME routine RMATPROPAdiabatic.x itself uses),
! so this driver stays in sync with whichever 3-boson dataset (5/10/20-channel) it's pointed at,
! and gets a genuine tabulated-potential B-spline tail for free.
!
! GENUINE HYBRID: when NumBoxes (read from the .inp) exceeds NumSVDBoxes, boxes
! NumSVDBoxes+1..NumBoxes are ordinary B-spline boxes reading the SAME tabulated Uad/Pmat/Qmat
! RMATPROPAdiabatic.x uses for this dataset -- exactly RMATPROPHybridSVD.f90's original box
! structure (SVD short range, B-spline tail beyond Rswitch), now built on the pluggable SVD
! engine. Set NumBoxes=NumSVDBoxes (and leave the OrderTail/xNumPointsTail/LegPointsTail trailer
! at any values, e.g. 0) for the previous all-SVD-to-xEnd behavior. Since the tail needs a real
! Uad.dat/Pmat.dat/Qmat.dat extending out to xEnd, point DataDir at an "_rmax200"-style dataset
! (fit/tabulated all the way to R~200) when using a tail that reaches that far -- the short-range
! (R=40-50) datasets are fine for an all-SVD (no-tail) run but not for a genuine hybrid one.
!****************************************************************************************************
PROGRAM main
  USE GlobalVars
  USE RMATPROPHybridSVDGenericMod, ONLY: RunHybridSVDGeneric
  USE SechFunctionPlugins
  USE AdiabaticPotential
  IMPLICIT NONE
  CHARACTER(LEN=64) :: DataDir
  CHARACTER(LEN=16) EnergyGridType
  CHARACTER*64 :: LegendreFile64
  CHARACTER*64 :: datadirLoc
  INTEGER ios, NumChannelsLoc, NumDataPointsLoc, NumSVDBoxes, NumBoxesLoc, NumEnergies
  INTEGER LSVD, OrderPhi, xNumPointsPhi, OrderTail, xNumPointsTail, LegPointsTail
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION, PARAMETER :: AlphaFactorFixed = 0d0
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION :: xStartLoc, xEndLoc, Emin, Emax, muLoc, AlphaFactorDataset, EffDimLoc

  OPEN(unit=7,file='RMATPROPHybridSVDGenericSech3Boson.inp',status='old')
  READ(7,*)
  READ(7,*) DataDir
  READ(7,*)
  READ(7,*)
  READ(7,*) NumSVDBoxes, NumBoxesLoc, xStartLoc, xEndLoc
  READ(7,*)
  READ(7,*)
  READ(7,*) Emin, Emax, NumEnergies
  READ(7,*,IOSTAT=ios) EnergyGridType
  IF (ios.NE.0) EnergyGridType = "linear"
  READ(7,*)
  READ(7,*)
  READ(7,*) LSVD
  READ(7,*)
  READ(7,*)
  READ(7,*) OrderPhi, xNumPointsPhi
  ! Optional trailer: B-spline tail's own Order/xNumPoints/LegPoints (radial quadrature knobs,
  ! same meaning as RMATPROPAdiabatic.x's own Order/xNumPoints/LegPoints -- NOT OrderPhi/
  ! xNumPointsPhi above, which are the SVD region's angular resolution). Absent (older .inp,
  ! or a no-tail run) defaults to 0, which is fine since RunHybridSVDGeneric never touches
  ! these when NumBoxesLoc==NumSVDBoxes.
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios) OrderTail, xNumPointsTail, LegPointsTail
  IF (ios.NE.0) THEN
     OrderTail = 0; xNumPointsTail = 0; LegPointsTail = 0
  ENDIF
  CLOSE(7)

  !--- Fit.data/masses.dat/FitLeff.data/Uad.dat/Pmat.dat/Qmat.dat, all in one call -- replaces
  !    this driver's former hand-rolled Fit.data/FitLeff.data-only parsing now that a genuine
  !    B-spline tail needs the tabulated potential interpolants too. AlphaFactorDataset (the
  !    dataset's own alpha=0.5) is read but deliberately NOT used -- see this file's own header. ---
  CALL ReadAdiabaticData(DataDir,NumChannelsLoc,NumDataPointsLoc,muLoc,AlphaFactorDataset,EffDimLoc, &
       Threshold,Leff,UInterp,PInterp,QInterp)
  IF (AlphaFactorDataset.NE.0d0) THEN
     WRITE(6,*) "NOTE: overriding dataset's own alpha=",AlphaFactorDataset," with AlphaFactor=0 -- ", &
          "see this file's own header / memory wavefn-reduction-term-bug."
  ENDIF

  ! No near-origin Coulomb patch for this sech^2 dataset (see RMATPROPAdiabatic.f90's own
  ! ApplyCoulombC gotcha -- a short-range bound-dimer channel's negative threshold does NOT mean
  ! a genuine Coulomb-like -C/R divergence). SVD boxes never reach R=0 with that assumption
  ! anyway, and the tail region only starts at Rswitch>0.
  ALLOCATE(CoulombC(NumChannelsLoc))
  CoulombC = 0d0

  ! FitLeff.data gives the raw threshold energy; CalcK's k(i)^2=2*mu*EE-Eth(i) formula expects
  ! Eth already pre-multiplied by 2*mu -- same rescaling RMATPROPAdiabatic.f90 applies.
  Threshold = 2d0*muLoc*Threshold

  LegendreFile64 = 'Legendre.dat'
  datadirLoc = '.'

  ! Left1D/Right1D = 1/1 (Neumann at phi=0 AND phi=pi/6 -- 3 identical bosons' even-parity BC,
  ! matching DriveAdiabaticSolver1D_3boson.par's own Order/Left/Right=5/1/1), xMax1D=pi/6.
  ! ShiftInit=-5.0d0 matches the .par's own adShift, so ARPACK targets the same low-lying states
  ! used to generate the tabulated baseline. GridRebuildEveryR=.FALSE. is required for SVD/DVR box
  ! construction regardless of SechGridMaker's own R-dependence (see RMATPROPHybridSVDGeneric.f90's
  ! GenericSVDChannelBasis.f90 gotcha) -- box construction anchors the angular basis once per box.
  CALL RunHybridSVDGeneric(NumChannelsLoc,muLoc,AlphaFactorFixed,EffDimLoc,Threshold,Leff, &
       NumSVDBoxes,NumBoxesLoc,xStartLoc,xEndLoc,BoxSpacingPower, &
       Emin,Emax,NumEnergies,EnergyGridType, &
       LSVD,OrderPhi,1,1,xNumPointsPhi,0d0,dacos(-1d0)/6.0d0, &
       LegendreFile64,10,-5.0d0,.FALSE., &
       SechPotential,SechGridMaker,SechRobinCoeff,datadirLoc, &
       'PhaseShiftSechSVD.dat','KMatrixFullSechSVD.dat','RecombinationSechSVD.dat','Radau', &
       OrderTail,xNumPointsTail,LegPointsTail,UInterp,PInterp,QInterp,CoulombC)

END PROGRAM main
