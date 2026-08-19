!****************************************************************************************************
! Diagnostic-only driver: dumps the sech^2 3-body adiabatic potential curves U_n(R) from the
! "fully SVD" (re-diagonalization via SechFunctionPlugins.f90/AdiabaticSolverGeneric.f90) pathway
! across a plain R grid, in the same column format as Adiabatic-Scattering-BoundStates's own
! Uad.dat -- so it can be directly diffed/plotted against the "fully adiabatic" (tabulated)
! reference produced by DriveAdiabaticSolver1D.x for the same physics. Not part of the R-matrix
! box-chaining machinery (no GridRebuildEveryR=.FALSE. box-sharing constraint here) -- this just
! calls OneDimChannelsGeneric once, directly, over the whole R range, exactly mirroring how
! OneDimChannels itself is driven by DriveAdiabaticSolver1D.f90.
!****************************************************************************************************
PROGRAM SechCurvesSVD
  USE potential
  USE SechFunctionPlugins
  USE AdiabaticSolverGeneric, ONLY: OneDimChannelsGeneric
  IMPLICIT NONE
  INTEGER, PARAMETER :: NumStates=5, Order=5, Left=1, Right=1
  INTEGER, PARAMETER :: xNumPoints=200, LegPoints=10, RSteps=400
  DOUBLE PRECISION, PARAMETER :: mu=0.57735026918962573d0
  DOUBLE PRECISION, PARAMETER :: RMin=1.d-6, RMax=50.d0
  DOUBLE PRECISION :: R(RSteps), Shift
  INTEGER :: PsiDim, i
  DOUBLE PRECISION, ALLOCATABLE :: Uad(:,:,:), Psi(:,:,:), S_out(:,:)
  CHARACTER*64 :: LegendreFile64, datadirLoc
  INTEGER, PARAMETER :: Uunit=10

  V2Depth = 2.0d0
  PotRange = 1.0d0
  alpha = 1.0d0
  LegendreFile64 = 'Legendre.dat'
  datadirLoc = '.'
  PsiDim = xNumPoints + Order - 3

  DO i = 1,RSteps
     R(i) = RMin + (RMax-RMin)*DBLE(i-1)/DBLE(RSteps-1)
  ENDDO

  ALLOCATE(Uad(RSteps,NumStates,2),Psi(RSteps,PsiDim,NumStates),S_out(Order+1,PsiDim))

  Shift = -5.0d0
  OPEN(unit=Uunit,file='3bodydata_sech2_5ch_bosonic/SVDCurvesRaw.log',status='replace')
  CALL OneDimChannelsGeneric(NumStates,PsiDim,0,0,LegendreFile64,LegPoints,Shift,Order,Left,Right,mu, &
       xNumPoints,0d0,dacos(-1d0)/2.0d0,R,RSteps,.TRUE.,0d0, &
       SechPotential,SechGridMaker,SechRobinCoeff, &
       Uunit,72,73,74,75,Uad,Psi,S_out,datadirLoc)
  CLOSE(Uunit)

  OPEN(unit=20,file='3bodydata_sech2_5ch_bosonic/SVDCurves.dat',status='replace')
  WRITE(20,'(A,2I8)') '#', NumStates, RSteps
  DO i = 1,RSteps
     WRITE(20,'(1P,6E22.14)') R(i), Uad(i,1:NumStates,1)
  ENDDO
  CLOSE(20)

END PROGRAM SechCurvesSVD
