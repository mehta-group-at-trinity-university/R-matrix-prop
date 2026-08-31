!****************************************************************************************************
! All-SVD hybrid driver for the equal-mass three-boson sech^2-well problem (3bodydata_sech2_5ch /
! 3bodydata_sech2_5ch_fine), reusing RMATPROPHybridSVDGeneric.f90's RunHybridSVDGeneric with
! SechFunctionPlugins.f90's SechPotential/SechGridMaker/SechRobinCoeff plugins -- the genuinely
! independent, non-interpolation numerical cross-check for RMATPROPAdiabatic.x's tabulated-data
! baseline on the same dataset (mirrors what RMATPROPHybridSVDGenericDelta.f90 already did for the
! delta-function benchmark, itself validated against RMATPROPSVDRobin.x). Every box is SVD,
! re-diagonalizing the true sech^2 potential at every quadrature R via ARPACK (AdiabaticSolverGeneric.f90)
! -- unlike the baseline, this never interpolates a tabulated Uad(R)/P(R)/Q(R), and never needs the
! delta-function-specific CoulombC near-origin substitution (RMATPROPAdiabatic.f90) at all, since
! there is no asymptotic near-origin form assumed -- box 1 just diagonalizes the real potential.
!
! Box1GridType='LobattoDirichlet' (true Dirichlet u(a1)=0), matching this repo's ordinary
! AlphaFactor=(EffDim-1)/2 convention (Fit.data: alpha=0.5, ddim=2.0 for this dataset) -- NOT
! 'Radau', which is only correct for AlphaFactor=0 (see RMATPROPHybridSVDGeneric.f90's own header).
!
! NumChannels/AlphaFactor/EffDim/mu/Threshold/Leff are read directly from DataDir's own
! Fit.data/masses.dat/FitLeff.data (same files/format RMATPROPAdiabatic.x's ReadAdiabaticData uses,
! AdiabaticPotential.f90) rather than hardcoded, so this driver stays in sync with whichever
! sech^2 dataset it's pointed at.
!****************************************************************************************************
PROGRAM main
  USE GlobalVars
  USE RMATPROPHybridSVDGenericMod, ONLY: RunHybridSVDGeneric
  USE SechFunctionPlugins
  USE AdiabaticPotential, ONLY: InterpolatingFunction, InterpolatingMatrix, ReadAdiabaticData
  IMPLICIT NONE
  CHARACTER(LEN=64) :: DataDir
  CHARACTER(LEN=16) EnergyGridType
  ! CHARACTER*64, not a bare literal -- see RMATPROPHybridSVDGenericDelta.f90's own comment on
  ! GetGaussFactors/OneDimChannelsGeneric's fixed-length dummies needing genuine blank-padded
  ! buffers, not bare literals threaded through assumed-length hops.
  CHARACTER*64 :: LegendreFile64
  CHARACTER*64 :: datadirLoc
  INTEGER ios, NumChannelsLoc, NumDataPointsLoc, NumSVDBoxes, NumBoxesLoc, NumEnergies, LSVD
  INTEGER OrderPhi, xNumPointsPhi, OrderTail, xNumPointsTail, LegPointsTail
  INTEGER NumChannelsDatasetDum, NumDataPointsDum
  INTEGER i, j
  LOGICAL :: haveTail
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  DOUBLE PRECISION, ALLOCATABLE :: ThresholdDum(:), LeffDum(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION :: muDum, alphaDum, ddimDum
  DOUBLE PRECISION :: xStartLoc, xEndLoc, Emin, Emax, muLoc, AlphaFactorLoc, EffDimLoc, dum

  OPEN(unit=7,file='RMATPROPHybridSVDGenericSech.inp',status='old')
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
  ! Optional hybrid trailer -- NumBoxes (total; >NumSVDBoxes requests a genuine B-spline tail) and
  ! the tail's B-spline radial quadrature. Absent entirely (IOSTAT/=0) means no tail, i.e. exactly
  ! the previous all-SVD behaviour, so existing .inp files keep working unchanged. DataDir is
  ! already read at the top of this file and supplies the tail's tabulated Uad/Pmat/Qmat.
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios) NumBoxesLoc
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios) OrderTail, xNumPointsTail, LegPointsTail
  haveTail = (ios.EQ.0)
  IF (.NOT.haveTail) THEN
     NumBoxesLoc = NumSVDBoxes
     OrderTail = 0; xNumPointsTail = 0; LegPointsTail = 0
  ENDIF
  CLOSE(7)

  !--- Fit.data header: NumChannels, NumDataPoints, alpha, ddim (AdiabaticPotential.f90's own
  !    ReadAdiabaticData format -- read directly here rather than duplicated/hardcoded, so this
  !    driver can't silently drift from whatever DataDir actually contains). ---
  OPEN(unit=21,file=TRIM(DataDir)//'/Fit.data',status='old')
  READ(21,*) NumChannelsLoc, NumDataPointsLoc, AlphaFactorLoc, EffDimLoc
  CLOSE(21)

  !--- masses.dat: comment line, then mu (hyperradial reduced mass) ---
  OPEN(unit=22,file=TRIM(DataDir)//'/masses.dat',status='old')
  READ(22,*)
  READ(22,*) muLoc
  CLOSE(22)

  !--- FitLeff.data: 2 header lines + blank, then j, Threshold(j), residual(unused), Leff(j),
  !    FitRangeMin(j) -- same format/column order as AdiabaticPotential.f90's own read. ---
  ALLOCATE(Threshold(NumChannelsLoc),Leff(NumChannelsLoc))
  OPEN(unit=26,file=TRIM(DataDir)//'/FitLeff.data',status='old')
  READ(26,*)
  READ(26,*)
  READ(26,*)
  DO i = 1,NumChannelsLoc
     READ(26,*) j, Threshold(i), dum, Leff(i)
  ENDDO
  CLOSE(26)
  ! FitLeff.data gives the raw threshold energy; CalcK's k(i)^2=2*mu*EE-Eth(i) formula expects
  ! Eth already pre-multiplied by 2*mu -- same rescaling RMATPROPAdiabatic.f90 applies.
  Threshold = 2d0*muLoc*Threshold

  LegendreFile64 = 'Legendre.dat'
  datadirLoc = '.'

  ! Left1D/Right1D = 1/0 (Neumann at phi=0, Dirichlet at phi=pi/2), matching
  ! DriveAdiabaticSolver1D_5ch.par exactly -- same angular BC used to generate the tabulated
  ! baseline this run is being compared against. OrderPhi should likewise be set to 5 in the
  ! .inp file to match that same .par's B-spline order.
  ! With the trailer present, boxes NumSVDBoxes+1..NumBoxesLoc are ordinary B-spline boxes reading
  ! DataDir's own tabulated Uad/Pmat/Qmat -- the same interpolants RMATPROPAdiabatic.x uses -- so
  ! this driver genuinely switches box type mid-chain. The tail's own Threshold/Leff are discarded
  ! (ThresholdDum/LeffDum): the values read from FitLeff.data above are kept unchanged so a no-tail
  ! run stays bit-for-bit identical to its previous behaviour. CoulombC=0: the sech^2 wells have no
  ! near-origin -C/R divergence (see RMATPROPAdiabatic.f90's ApplyCoulombC discussion), and the
  ! tail never reaches the origin anyway.
  ! SechPotential sums three pairwise sech^2 wells, V(R,phi)=SechAlpha*sum_{pairs}V2body(r_pair),
  ! r_pair=RPair*R*|cos(phi-offset)|.
  ALLOCATE(CoulombC(NumChannelsLoc))
  CoulombC = 0d0
  IF (haveTail) THEN
     CALL ReadAdiabaticData(DataDir,NumChannelsDatasetDum,NumDataPointsDum,muDum,alphaDum,ddimDum, &
          ThresholdDum,LeffDum,UInterp,PInterp,QInterp)
     CALL RunHybridSVDGeneric(NumChannelsLoc,muLoc,AlphaFactorLoc,EffDimLoc,Threshold,Leff, &
          NumSVDBoxes,NumBoxesLoc,xStartLoc,xEndLoc,BoxSpacingPower, &
          Emin,Emax,NumEnergies,EnergyGridType, &
          LSVD,OrderPhi,1,0,xNumPointsPhi,0d0,dacos(-1d0)/2.0d0, &
          LegendreFile64,10,-1.0d0,.FALSE., &
          SechPotential,SechGridMaker,SechRobinCoeff,datadirLoc, &
          'PhaseShiftSech.dat','KMatrixFullSech.dat','RecombinationSech.dat','LobattoDirichlet', &
          OrderTail,xNumPointsTail,LegPointsTail,UInterp,PInterp,QInterp,CoulombC)
  ELSE
     CALL RunHybridSVDGeneric(NumChannelsLoc,muLoc,AlphaFactorLoc,EffDimLoc,Threshold,Leff, &
          NumSVDBoxes,NumBoxesLoc,xStartLoc,xEndLoc,BoxSpacingPower, &
          Emin,Emax,NumEnergies,EnergyGridType, &
          LSVD,OrderPhi,1,0,xNumPointsPhi,0d0,dacos(-1d0)/2.0d0, &
          LegendreFile64,10,-1.0d0,.FALSE., &
          SechPotential,SechGridMaker,SechRobinCoeff,datadirLoc, &
          'PhaseShiftSech.dat','KMatrixFullSech.dat','RecombinationSech.dat','LobattoDirichlet')
  ENDIF

END PROGRAM main
