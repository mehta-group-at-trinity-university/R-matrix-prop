!****************************************************************************************************
! SVD analog of PlotWavefunctionFree.f90: reconstructs the actual DVR wavefunction u(R) for the
! LOWEST channel (Leff=1) of the '3bodydata_free_5ch' dataset, single SVD box (NumBoxes=1,
! GridType='LobattoDirichlet' so a1=0 is a true node with its DOF eliminated -- Dirichlet u(0)=0,
! matching AlphaFactor=(EffDim-1)/2's regularity requirement), at the same positive scattering
! energy where all 5 channels are open. Written to overlay directly against
! PlotWavefunctionFree.f90's B-spline reconstruction and the exact sqrt(kR)*J_1(kR) solution.
!
! Plotted directly at the box's own DVR nodes -- no interpolation needed: for the
! normalized DVR basis pi_j(R)=L_j(R)/sqrt(w_j), pi_j(R_j)=1/sqrt(w_j) and pi_m(R_j)=0
! for m/=j exactly (Lagrange's own Kronecker-delta property), so u(R_j)=cfull_j/sqrt(w_j)
! directly. The node/weight arrays themselves are independently rebuilt here via the
! same GetGaussLobattoFactors call BuildSVDBox uses internally (deterministic given
! L,a1,a2) since BuildSVDBox doesn't export them.
!****************************************************************************************************
PROGRAM PlotWavefunctionSVD
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE AdiabaticPotential
  USE InterpType
  USE SVDChannelBasis
  IMPLICIT NONE

  CHARACTER(LEN=64) DataDir
  TYPE(BPData) BPD1
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull, Box1
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim, EnergyLocal, a1, a2, Shift1D
  DOUBLE PRECISION xMin1D, xMax1D
  INTEGER NumDataPoints, NumOpenChannels, i, L, LQuad, Order1D, Left1D, Right1D
  INTEGER xNumPoints1D, LegPoints1D

  DOUBLE PRECISION, ALLOCATABLE :: Xout(:,:), copen(:), cclosed(:), cfull(:)
  INTEGER, ALLOCATABLE :: OpenIdxOut(:), ClosedIdxOut(:)
  INTEGER NumOpen, Ncc, betaCh1
  DOUBLE PRECISION, ALLOCATABLE :: ODiagL(:), ODiagR(:)
  DOUBLE PRECISION :: wsq1Out, wsqLOut

  ! -- reconstruction grid state (independently rebuilt, see header note) --
  DOUBLE PRECISION, ALLOCATABLE :: nodesQ(:), weightsQ(:)
  CHARACTER*64 :: LobattoFile64
  INTEGER :: j
  DOUBLE PRECISION :: Rval, uval

  DataDir = '3bodydata_free_5ch'
  a1 = 0.0d0
  a2 = 10.0d0
  L = 63           ! LQuad=L+1=64, tabulated in GaussLobatto.dat
  Order1D = 5
  Left1D = 1
  Right1D = 0
  xNumPoints1D = 200
  xMin1D = 0.d0
  xMax1D = 1.570796326794897d0
  LegPoints1D = 10
  Shift1D = -1.d0
  EnergyLocal = 1.0d0

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = LegPoints1D
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  CALL ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alpha,ddim, &
       Threshold,Leff,UInterp,PInterp,QInterp)
  reducedmass = muLocal
  AlphaFactor = alpha
  EffDim = ddim

  Threshold = 2d0*reducedmass*Threshold
  NumOpenChannels = 0
  DO i = 1,NumChannels
     IF (2d0*reducedmass*EnergyLocal.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
  ENDDO
  WRITE(6,*) "NumChannels=",NumChannels," NumOpenChannels=",NumOpenChannels

  ALLOCATE(ODiagL(NumChannels),ODiagR(NumChannels))
  CALL BuildSVDBox(a1,a2,L,NumChannels,Order1D,Left1D,Right1D, &
       xNumPoints1D,xMin1D,xMax1D,LegendreFile,LegPoints1D,Shift1D,DataDir, &
       0,NumOpenChannels,'LobattoDirichlet',BPD1,EIG,ODiagL,ODiagR,wsq1Out,wsqLOut)
  WRITE(6,*) "CHECKPOINT: BuildSVDBox done, MatrixDim=",BPD1%MatrixDim

  ! BuildSVDBox leaves BPD%lam at its placeholder zero (see its own header note) --
  ! the caller must fill it, exactly as RMATPROPHybridSVD.f90 does for the last SVD
  ! box, or CalcK's outer hyperrjry matching silently uses Leff=0 for every channel.
  BPD1%lam(1:NumChannels) = Leff(1:NumChannels)

  CALL CombineGam(EIG,EnergyLocal)

  Box1%NumOpenL = 0
  Box1%NumOpenR = NumOpenChannels
  CALL AllocateBox(Box1)
  CALL InitZeroBox(Box1)

  ALLOCATE(evalRed(Box1%betaMax),evecRed(Box1%betaMax,Box1%betaMax),LamooDiag(Box1%betaMax))
  NumOpen = Box1%betaMax
  Ncc = BPD1%MatrixDim - NumOpen
  ALLOCATE(Xout(Ncc,NumOpen),OpenIdxOut(NumOpen),ClosedIdxOut(Ncc))
  CALL PartitionAndEliminateFull(BPD1,EIG,Box1%NumOpenL,Box1%NumOpenR,evalRed,evecRed,LamooDiag, &
       Xout,OpenIdxOut,ClosedIdxOut)
  WRITE(6,*) "evalRed=",evalRed

  betaCh1 = 1
  DO i = 1,NumOpen
     IF (ABS(evecRed(1,i)).GT.ABS(evecRed(1,betaCh1))) betaCh1 = i
  ENDDO
  WRITE(6,*) "Selected beta (channel-1-pure mode) = ",betaCh1
  WRITE(6,*) "evecRed(:,betaCh1) = ",evecRed(:,betaCh1)

  ALLOCATE(copen(NumOpen),cclosed(Ncc),cfull(BPD1%MatrixDim))
  copen = evecRed(:,betaCh1)
  cclosed = -MATMUL(Xout,copen)
  DO i = 1,NumOpen
     cfull(OpenIdxOut(i)) = copen(i)
  ENDDO
  DO i = 1,Ncc
     cfull(ClosedIdxOut(i)) = cclosed(i)
  ENDDO

  ! Cross-check against the normal BoxMatch+CalcK pathway.
  Bnull%NumOpenL = 0
  Bnull%NumOpenR = 0
  CALL AllocateBox(Bnull)
  CALL InitZeroBox(Bnull)
  CALL BoxMatch(Bnull, Box1, BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
  CALL AllocateScat(SD,NumChannels)
  CALL CalcK(Box1,BPD1,SD,reducedmass,EffDim,AlphaFactor,EnergyLocal,Threshold)
  WRITE(6,*) "Full K diagonal = ", (SD%K(i,i), i=1,NumChannels)

  ! -- independently rebuild the FULL (L+1)-point Lobatto grid, matching BuildSVDBox's
  !    own GridType='LobattoDirichlet' construction exactly (see header note), just to
  !    get the RETAINED nodes' own R-values and weights back -- no interpolation needed:
  !    for the normalized DVR basis pi_j(R)=L_j(R)/sqrt(w_j), pi_j(R_j)=1/sqrt(w_j) and
  !    pi_m(R_j)=0 for m/=j exactly (Lagrange's own Kronecker-delta property), so
  !    u(R_j) = cfull_j/sqrt(w_j) directly -- no need to evaluate the Lagrange product
  !    at all, just plot at the box's own DVR points. ---
  LQuad = L+1
  ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
  LobattoFile64 = 'GaussLobatto.dat'
  CALL GetGaussLobattoFactors(LobattoFile64,LQuad,nodesQ,weightsQ)
  nodesQ = 0.5d0*(a1+a2) + 0.5d0*(a2-a1)*nodesQ
  weightsQ = 0.5d0*(a2-a1)*weightsQ

  OPEN(unit=99,file='SVD_wavefunction_free.dat',status='replace')
  DO j = 1,L
     Rval = nodesQ(j+1)
     uval = cfull((j-1)*NumChannels+1)/DSQRT(weightsQ(j+1))
     WRITE(99,*) Rval, uval
  ENDDO
  CLOSE(99)
  WRITE(6,*) "Wrote SVD_wavefunction_free.dat"

END PROGRAM PlotWavefunctionSVD
