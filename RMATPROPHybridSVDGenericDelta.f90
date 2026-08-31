!****************************************************************************************************
! Regression test for RMATPROPHybridSVDGeneric.f90's RunHybridSVDGeneric: a thin driver, reusing
! the SAME physics/plugins as RMATPROPHybridSVDRobin.f90 (equal-mass three-boson delta-function
! problem, DeltaFunctionPlugins.f90's DeltaZeroPotential/DeltaWedgeGridMaker/DeltaRobinCoeff), but
! routed through the genuinely reusable RunHybridSVDGeneric subroutine instead of RobinSVDChannel-
! Basis.f90's bespoke box-builder. Must reproduce RMATPROPHybridSVDRobin.x's own validated results
! (matches KMSFormulas.f's exact SolveQ to machine/basis-truncation precision) -- proves the
! genericized propagator itself is correct before it's trusted with new (atom-ion) physics.
! Same RMATPROPHybridSVDRobin.inp input format; reads RMATPROPHybridSVDGenericDelta.inp.
!
! GENUINE HYBRID (optional): an OrderTail/xNumPointsTail/LegPointsTail/DataDir trailer in the
! .inp (all four present) requests a B-spline tail beyond NumSVDBoxes, reading
! 3bodydata_delta_5ch's own tabulated Uad/Pmat/Qmat for boxes NumSVDBoxes+1..NumBoxes -- exactly
! what RMATPROPHybridSVDGenericSech3Boson.f90 does. Absent trailer (or NumBoxes==NumSVDBoxes)
! reproduces this driver's original all-SVD regression-test behavior exactly.
!****************************************************************************************************
PROGRAM main
  USE GlobalVars
  USE RMATPROPHybridSVDGenericMod, ONLY: RunHybridSVDGeneric
  USE DeltaFunctionPlugins
  USE AdiabaticPotential, ONLY: InterpolatingFunction, InterpolatingMatrix, ReadAdiabaticData
  IMPLICIT NONE
  CHARACTER(LEN=16) EnergyGridType
  CHARACTER(LEN=64) :: DataDir
  ! CHARACTER*64, not a bare literal in the CALL below: GetGaussFactors (Quadrature.f90)
  ! trims its filename dummy at the first blank (File(1:INDEX(File,' ')-1)), which only
  ! works if the string arriving there is a genuinely blank-padded 64-byte buffer -- a
  ! bare literal threaded through RunHybridSVDGeneric's assumed-length LegendreFile1D
  ! dummy carries no such padding, so INDEX(...,' ') finds no blank at all and the
  ! resulting File(1:-1) substring opens an empty filename (confirmed: this exact crash,
  ! "Cannot open file ''", before fixing it this way).
  CHARACTER*64 :: LegendreFile64
  ! Same reasoning as LegendreFile64: OneDimChannelsGeneric's own datadir dummy is fixed
  ! CHARACTER(LEN=64), not assumed-length, so the actual argument reaching it must
  ! already be a genuine 64-byte buffer, not a bare literal threaded through several
  ! assumed-length hops.
  CHARACTER*64 :: datadirLoc
  INTEGER ios, NumChannelsLoc, NumSVDBoxes, NumBoxesLoc, NumEnergies, LSVD, OrderPhi, xNumPointsPhi, i
  INTEGER OrderTail, xNumPointsTail, LegPointsTail, NumChannelsDatasetDum, NumDataPointsDum
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION, PARAMETER :: DeltaMu = 0.57735026918962584d0
  DOUBLE PRECISION, PARAMETER :: DeltaAlpha = 0d0
  DOUBLE PRECISION, PARAMETER :: DeltaDdim = 2d0
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  DOUBLE PRECISION, ALLOCATABLE :: ThresholdDum(:), LeffDum(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION :: xStartLoc, xEndLoc, Emin, Emax, muDum, alphaDum, ddimDum
  LOGICAL :: haveTail

  OPEN(unit=7,file='RMATPROPHybridSVDGenericDelta.inp',status='old')
  READ(7,*)
  READ(7,*) NumChannelsLoc
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
  ! Optional hybrid trailer: NumBoxes (total, >NumSVDBoxes for a genuine tail),
  ! OrderTail/xNumPointsTail/LegPointsTail (B-spline radial quadrature), DataDir (dataset
  ! providing the tail's tabulated Uad/Pmat/Qmat). Absent entirely -- IOSTAT/=0 -- means no tail.
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios) NumBoxesLoc
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios) OrderTail, xNumPointsTail, LegPointsTail
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios)
  READ(7,*,IOSTAT=ios) DataDir
  haveTail = (ios.EQ.0)
  IF (.NOT.haveTail) NumBoxesLoc = NumSVDBoxes
  CLOSE(7)

  ! Equal-mass 3-boson delta-function problem's exact hyperspherical channel functions:
  ! channel 1 (atom-dimer, Threshold=-1) has Leff=1/2; excited (3-body-threshold, Threshold=0)
  ! channels have Leff=6i-9 (Mehta, Esry & Greene, PRA 76, 022711 (2007)).
  ALLOCATE(Threshold(NumChannelsLoc),Leff(NumChannelsLoc))
  Threshold(1) = -1d0
  Leff(1) = 0.5d0
  DO i = 2,NumChannelsLoc
     Threshold(i) = 0d0
     Leff(i) = 6d0*DBLE(i) - 9d0
  ENDDO
  Threshold = 2d0*DeltaMu*Threshold
  LegendreFile64 = 'Legendre.dat'
  datadirLoc = '.'

  ! Tail region needs the tabulated Uad/Pmat/Qmat interpolants (its own Threshold/Leff are
  ! discarded -- ThresholdDum/LeffDum -- since 3bodydata_delta_5ch's FitLeff.data already
  ! matches the exact analytic Threshold/Leff above to machine precision; keeping the
  ! hardcoded values above unchanged preserves this driver's own regression-test identity).
  ! No near-origin Coulomb patch needed for the tail: that near-origin divergence is handled
  ! exactly by DeltaRobinCoeff's Robin BC in the SVD region (Rswitch>0), never seen by the tail.
  IF (haveTail) THEN
     CALL ReadAdiabaticData(DataDir,NumChannelsDatasetDum,NumDataPointsDum,muDum,alphaDum,ddimDum, &
          ThresholdDum,LeffDum,UInterp,PInterp,QInterp)
  ENDIF
  ALLOCATE(CoulombC(NumChannelsLoc))
  CoulombC = 0d0

  ! NumBoxesIn=NumSVDBoxes (haveTail=.FALSE.) -- no B-spline tail (all-SVD, exactly the
  ! previously validated behavior); RunHybridSVDGeneric's tail arguments are OPTIONAL and
  ! simply omitted then. DeltaZeroPotential: V=0 identically (the contact interaction enters
  ! only through DeltaRobinCoeff's Robin BC, kRight=sqrt(2*mu)*R, not an ordinary potential).
  IF (haveTail) THEN
     CALL RunHybridSVDGeneric(NumChannelsLoc,DeltaMu,DeltaAlpha,DeltaDdim,Threshold,Leff, &
          NumSVDBoxes,NumBoxesLoc,xStartLoc,xEndLoc,BoxSpacingPower, &
          Emin,Emax,NumEnergies,EnergyGridType, &
          LSVD,OrderPhi,1,3,xNumPointsPhi,0d0,dacos(-1d0)/6.0d0, &
          LegendreFile64,10,-2.0d0,.FALSE., &
          DeltaZeroPotential,DeltaWedgeGridMaker,DeltaRobinCoeff,datadirLoc, &
          'PhaseShift.dat','KMatrixFull.dat','Recombination.dat','Radau', &
          OrderTail,xNumPointsTail,LegPointsTail,UInterp,PInterp,QInterp,CoulombC)
  ELSE
     CALL RunHybridSVDGeneric(NumChannelsLoc,DeltaMu,DeltaAlpha,DeltaDdim,Threshold,Leff, &
          NumSVDBoxes,NumBoxesLoc,xStartLoc,xEndLoc,BoxSpacingPower, &
          Emin,Emax,NumEnergies,EnergyGridType, &
          LSVD,OrderPhi,1,3,xNumPointsPhi,0d0,dacos(-1d0)/6.0d0, &
          LegendreFile64,10,-2.0d0,.FALSE., &
          DeltaZeroPotential,DeltaWedgeGridMaker,DeltaRobinCoeff,datadirLoc, &
          'PhaseShift.dat','KMatrixFull.dat','Recombination.dat','Radau')
  ENDIF

END PROGRAM main
