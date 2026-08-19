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
! NumChannels/mu/Threshold/Leff are read directly from DataDir's own Fit.data/masses.dat/
! FitLeff.data (same files/format RMATPROPAdiabatic.x's ReadAdiabaticData uses), so this driver
! stays in sync with whichever 3-boson dataset (5/10/20-channel) it's pointed at.
!****************************************************************************************************
PROGRAM main
  USE GlobalVars
  USE RMATPROPHybridSVDGenericMod, ONLY: RunHybridSVDGeneric
  USE SechFunctionPlugins
  IMPLICIT NONE
  CHARACTER(LEN=64) :: DataDir
  CHARACTER(LEN=16) EnergyGridType
  CHARACTER*64 :: LegendreFile64
  CHARACTER*64 :: datadirLoc
  INTEGER ios, NumChannelsLoc, NumDataPointsLoc, NumSVDBoxes, NumEnergies, LSVD, OrderPhi, xNumPointsPhi
  INTEGER i, j
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION, PARAMETER :: AlphaFactorFixed = 0d0
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  DOUBLE PRECISION :: xStartLoc, xEndLoc, Emin, Emax, muLoc, AlphaFactorDataset, EffDimLoc, dum

  OPEN(unit=7,file='RMATPROPHybridSVDGenericSech3Boson.inp',status='old')
  READ(7,*)
  READ(7,*) DataDir
  READ(7,*)
  READ(7,*)
  READ(7,*) NumSVDBoxes, xStartLoc, xEndLoc
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
  CLOSE(7)

  !--- Fit.data header: NumChannels, NumDataPoints, alpha, ddim (AdiabaticPotential.f90's own
  !    ReadAdiabaticData format). AlphaFactorDataset (the dataset's own alpha=0.5) is read but
  !    deliberately NOT used -- see this file's own header. ---
  OPEN(unit=21,file=TRIM(DataDir)//'/Fit.data',status='old')
  READ(21,*) NumChannelsLoc, NumDataPointsLoc, AlphaFactorDataset, EffDimLoc
  CLOSE(21)
  IF (AlphaFactorDataset.NE.0d0) THEN
     WRITE(6,*) "NOTE: overriding dataset's own alpha=",AlphaFactorDataset," with AlphaFactor=0 -- ", &
          "see this file's own header / memory wavefn-reduction-term-bug."
  ENDIF

  !--- masses.dat: comment line, then mu (hyperradial reduced mass) ---
  OPEN(unit=22,file=TRIM(DataDir)//'/masses.dat',status='old')
  READ(22,*)
  READ(22,*) muLoc
  CLOSE(22)

  !--- FitLeff.data: 2 header lines + blank, then either the 5-column (j,Threshold,residual,Leff,
  !    FitRangeMin) or 4-column (j,Threshold,Leff,FitRangeMin) format -- same dual-format fallback
  !    as AdiabaticPotential.f90's ReadAdiabaticData (see that routine's own comment). ---
  ALLOCATE(Threshold(NumChannelsLoc),Leff(NumChannelsLoc))
  OPEN(unit=26,file=TRIM(DataDir)//'/FitLeff.data',status='old')
  READ(26,*)
  READ(26,*)
  READ(26,*)
  BLOCK
    CHARACTER(LEN=200) :: fitLeffLine
    DOUBLE PRECISION :: fitRangeMinDum
    DO i = 1,NumChannelsLoc
       READ(26,'(A)') fitLeffLine
       READ(fitLeffLine,*,IOSTAT=ios) j, Threshold(i), dum, Leff(i), fitRangeMinDum
       IF (ios.NE.0) THEN
          READ(fitLeffLine,*,IOSTAT=ios) j, Threshold(i), Leff(i), fitRangeMinDum
       ENDIF
       IF (ios.NE.0) THEN
          WRITE(6,*) "RMATPROPHybridSVDGenericSech3Boson: ",TRIM(DataDir)//"/FitLeff.data channel ",i, &
               " matches neither the 5-column nor 4-column FitLeff.data format."
          STOP
       ENDIF
    ENDDO
  END BLOCK
  CLOSE(26)
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
       NumSVDBoxes,xStartLoc,xEndLoc,BoxSpacingPower, &
       Emin,Emax,NumEnergies,EnergyGridType, &
       LSVD,OrderPhi,1,1,xNumPointsPhi,0d0,dacos(-1d0)/6.0d0, &
       LegendreFile64,10,-5.0d0,.FALSE., &
       SechPotential,SechGridMaker,SechRobinCoeff,datadirLoc, &
       'PhaseShiftSechSVD.dat','KMatrixFullSechSVD.dat','RecombinationSechSVD.dat','Radau')

END PROGRAM main
