!****************************************************************************************************
! Regression test for RMATPROPHybridSVDGeneric.f90's RunHybridSVDGeneric: a thin driver, reusing
! the SAME physics/plugins as RMATPROPHybridSVDRobin.f90 (equal-mass three-boson delta-function
! problem, DeltaFunctionPlugins.f90's DeltaZeroPotential/DeltaWedgeGridMaker/DeltaRobinCoeff), but
! routed through the genuinely reusable RunHybridSVDGeneric subroutine instead of RobinSVDChannel-
! Basis.f90's bespoke box-builder. Must reproduce RMATPROPHybridSVDRobin.x's own validated results
! (matches KMSFormulas.f's exact SolveQ to machine/basis-truncation precision) -- proves the
! genericized propagator itself is correct before it's trusted with new (atom-ion) physics.
! Same RMATPROPHybridSVDRobin.inp input format; reads RMATPROPHybridSVDGenericDelta.inp.
!****************************************************************************************************
PROGRAM main
  USE GlobalVars
  USE RMATPROPHybridSVDGenericMod, ONLY: RunHybridSVDGeneric
  USE DeltaFunctionPlugins
  IMPLICIT NONE
  CHARACTER(LEN=16) EnergyGridType
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
  INTEGER ios, NumChannelsLoc, NumSVDBoxes, NumEnergies, LSVD, OrderPhi, xNumPointsPhi, i
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION, PARAMETER :: DeltaMu = 0.57735026918962584d0
  DOUBLE PRECISION, PARAMETER :: DeltaAlpha = 0d0
  DOUBLE PRECISION, PARAMETER :: DeltaDdim = 2d0
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  DOUBLE PRECISION :: xStartLoc, xEndLoc, Emin, Emax

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
  CLOSE(7)

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

  CALL RunHybridSVDGeneric(NumChannelsLoc,DeltaMu,DeltaAlpha,DeltaDdim,Threshold,Leff, &
       NumSVDBoxes,xStartLoc,xEndLoc,BoxSpacingPower, &
       Emin,Emax,NumEnergies,EnergyGridType, &
       LSVD,OrderPhi,1,3,xNumPointsPhi,0d0,dacos(-1d0)/6.0d0, &
       LegendreFile64,10,-2.0d0,.FALSE., &
       DeltaZeroPotential,DeltaWedgeGridMaker,DeltaRobinCoeff,datadirLoc, &
       'PhaseShift.dat','KMatrixFull.dat','Recombination.dat','Radau')

END PROGRAM main
