!****************************************************************************************************
! Single-channel (3bodydata_sech2_3boson_ch1) diagnostic: box 1 (the ONLY box, [xMin,xMax]) built
! with Left=3 (Robin BC) instead of the driver's usual Left=2/PartitionAndEliminate regularity
! mechanism, with kLeft fixed to the TRUE near-R=0 log-derivative from the leading-order Frobenius
! solution -- NOT BuildBox1CombinedBasis's own default (hyperrjry order=Leff, valid only for the
! FAR-FIELD asymptotic fit, not the near-origin behavior of this channel, whose U(R) is finite at
! R=0 so the only near-origin singularity is the universal Langer term -alpha(alpha-D+2)/R^2).
! Indicial equation for that term (D=2,alpha=0.5) has a DOUBLE root s=1/2, so u(R)~sqrt(R) as
! R->0 and L(R)=u'/u -> 1/(2R) (leading order -- same approximation level already validated as
! "sufficient" for the independent NDSolve shooting-method check, see check_ch1_phaseshift.wls).
!
! Reconstructs the FEM wavefunction via PartitionAndEliminateFull's Xout/OpenIdxOut/ClosedIdxOut
! (the "phirecon" equivalent already used by PlotWavefunction.f90), and reports K/delta at a single
! hard-coded energy for direct comparison against the independent NDSolve result there.
!****************************************************************************************************
PROGRAM PlotWavefunctionCh1Robin
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE AdiabaticPotential
  USE InterpType
  IMPLICIT NONE

  CHARACTER(LEN=64) DataDir
  TYPE(BPData) BPD1
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull, Box1
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: u1Box1(:,:,:), ux1Box1(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: uAux(:,:,:), uxAux(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim, EnergyLocal, xMinLocal, xMaxLocal
  DOUBLE PRECISION, EXTERNAL :: MyBSpline
  DOUBLE PRECISION kLeft, constLeft, qWave
  INTEGER NumDataPoints, NumOpenChannels, xNumPoints, LegPointsLocal, i, kx, lx
  INTEGER OrderLocal, xDim2
  INTEGER, ALLOCATABLE :: xBounds2(:)

  DOUBLE PRECISION, ALLOCATABLE :: Xout(:,:), copen(:), cclosed(:), cfull(:)
  INTEGER, ALLOCATABLE :: OpenIdxOut(:), ClosedIdxOut(:)
  INTEGER NumOpen, Ncc, ix

  DataDir = '3bodydata_sech2_3boson_ch1'
  OrderLocal = 5
  xNumPoints = 400
  LegPointsLocal = 10
  xMinLocal = 1.0d-4
  xMaxLocal = 40.0d0
  EnergyLocal = 0.5d0

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
  Order = OrderLocal

  ALLOCATE(CoulombC(NumChannels))
  CoulombC = 0d0   ! this dataset's near-origin U(R) is regular, not Coulomb-divergent

  Threshold = 2d0*reducedmass*Threshold

  NumOpenChannels = 0
  DO i = 1,NumChannels
     IF (2d0*reducedmass*EnergyLocal.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
  ENDDO
  WRITE(6,*) "NumOpenChannels=",NumOpenChannels

  BPD1%NumChannels = NumChannels
  BPD1%Order = OrderLocal
  BPD1%Left = 3
  BPD1%Right = 2
  BPD1%xNumPoints = xNumPoints
  BPD1%kl = 1d0
  BPD1%kr = 1d20
  CALL AllocateBPD(BPD1)
  BPD1%xl = xMinLocal
  BPD1%xr = xMaxLocal
  CALL GridMaker(BPD1%xPoints,BPD1%xNumPoints,BPD1%xl,BPD1%xr,"quadratic")
  CALL MakeBasis(BPD1)
  WRITE(6,*) "First few xPoints (knots): ", BPD1%xPoints(1:5)
  BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
  CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)
  BPD1%Pot = 0d0   ! TEMPORARY: zero U+Q/(2mu) entirely, isolating the kinetic + wavefunction-
                   ! reduction-term (alphaTerm/R^2) structure alone, per direct request.

  ! TEMPORARY diagnostic: dump the actual in-memory Pot(1,1,lx,kx) at a handful of R
  ! values, to diff directly against an independent Wolfram evaluation of what the
  ! potential SHOULD be there -- removes any ambiguity about what SetAdiabaticPotential
  ! actually computed vs. what the source code appears to say it computes.
  OPEN(unit=77,file='PotDiag.dat',status='replace')
  WRITE(77,'(A)') '# R  Pot(1,1)'
  DO kx = 1,BPD1%xNumPoints-1
     DO lx = 1,LegPoints
        WRITE(77,'(1P,2E24.16)') BPD1%x(lx,kx), BPD1%Pot(1,1,lx,kx)
     ENDDO
  ENDDO
  CLOSE(77)

  ! --- kLeft fixed to the TRUE near-R=0 log-derivative (leading-order Frobenius), NOT
  !     BuildBox1CombinedBasis's own hyperrjry/Leff-based (far-field) default. Replicates
  !     that subroutine's own u1=B1+B2/constLeft combination inline with this kLeft. ---
  kLeft = 1.0d0/(2.0d0*xMinLocal)
  WRITE(6,*) "kLeft (true near-origin, leading Frobenius) = ", kLeft

  xDim2 = BPD1%xDim + 1
  ALLOCATE(xBounds2(BPD1%xNumPoints+2*BPD1%Order))
  ALLOCATE(uAux(LegPoints,BPD1%xNumPoints,xDim2),uxAux(LegPoints,BPD1%xNumPoints,xDim2))
  CALL CalcBasisFuncsBP(2,BPD1%Right,1d0,BPD1%kr,BPD1%Order,BPD1%xPoints,LegPoints,xLeg, &
       xDim2,xBounds2,BPD1%xNumPoints,0,uAux)
  CALL CalcBasisFuncsBP(2,BPD1%Right,1d0,BPD1%kr,BPD1%Order,BPD1%xPoints,LegPoints,xLeg, &
       xDim2,xBounds2,BPD1%xNumPoints,1,uxAux)

  constLeft = (MyBSpline(BPD1%Order,1,BPD1%xNumPoints,BPD1%xPoints,2,BPD1%xPoints(1)) &
       - MyBSpline(BPD1%Order,0,BPD1%xNumPoints,BPD1%xPoints,2,BPD1%xPoints(1))*kLeft) &
       / (MyBSpline(BPD1%Order,0,BPD1%xNumPoints,BPD1%xPoints,1,BPD1%xPoints(1))*kLeft &
       - MyBSpline(BPD1%Order,1,BPD1%xNumPoints,BPD1%xPoints,1,BPD1%xPoints(1)))

  ALLOCATE(u1Box1(LegPoints,BPD1%xNumPoints,NumChannels))
  ALLOCATE(ux1Box1(LegPoints,BPD1%xNumPoints,NumChannels))
  u1Box1(:,:,1) = uAux(:,:,1) + uAux(:,:,2)/constLeft
  ux1Box1(:,:,1) = uxAux(:,:,1) + uxAux(:,:,2)/constLeft

  EIG%MatrixDim = BPD1%MatrixDim
  CALL AllocateEIG(EIG)
  CALL CalcGamLamBox1(BPD1,EIG,0,NumOpenChannels,u1Box1,ux1Box1,(/kLeft/))
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

  WRITE(6,*) "NumOpen=",NumOpen," Ncc=",Ncc," evalRed=",evalRed

  ALLOCATE(copen(NumOpen),cclosed(Ncc),cfull(BPD1%MatrixDim))
  copen = evecRed(:,1)
  cclosed = -MATMUL(Xout,copen)
  DO i = 1,NumOpen
     cfull(OpenIdxOut(i)) = copen(i)
  ENDDO
  DO i = 1,Ncc
     cfull(ClosedIdxOut(i)) = cclosed(i)
  ENDDO

  Bnull%NumOpenL = 0
  Bnull%NumOpenR = 0
  CALL AllocateBox(Bnull)
  CALL InitZeroBox(Bnull)
  CALL BoxMatch(Bnull, Box1, BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
  CALL AllocateScat(SD,NumChannels)
  CALL CalcK(Box1,BPD1,SD,reducedmass,EffDim,AlphaFactor,EnergyLocal,Threshold)
  qWave = DSQRT(2d0*reducedmass*EnergyLocal-Threshold(1))
  WRITE(6,*) "q = ", qWave, " Energy = ", EnergyLocal
  WRITE(6,*) "K (from CalcK) = ", SD%K(1,1)
  WRITE(6,*) "delta = atan(K) = ", DATAN(SD%K(1,1))
  WRITE(6,*) "Box1%Zf(1,1) [u(xMax)]  = ", Box1%Zf(1,1)
  WRITE(6,*) "Box1%Zfp(1,1) [u'(xMax)] = ", Box1%Zfp(1,1)
  WRITE(6,*) "L = Zfp/Zf at xMax = ", Box1%Zfp(1,1)/Box1%Zf(1,1)
  WRITE(6,*) "xMax = ", BPD1%xr

  OPEN(unit=99,file='FEM_wavefunction_ch1_robin.dat',status='replace')
  DO kx = 1,BPD1%xNumPoints-1
     DO lx = 1,LegPoints
        BLOCK
          DOUBLE PRECISION Rval, uval
          Rval = BPD1%x(lx,kx)
          uval = 0d0
          DO ix = 1,BPD1%xDim
             IF (ix.EQ.1) THEN
                uval = uval + cfull((ix-1)*NumChannels+1)*u1Box1(lx,kx,1)
             ELSE
                uval = uval + cfull((ix-1)*NumChannels+1)*BPD1%u(lx,kx,ix)
             ENDIF
          ENDDO
          WRITE(99,*) Rval, uval
        END BLOCK
     ENDDO
  ENDDO
  CLOSE(99)
  WRITE(6,*) "Wrote FEM_wavefunction_ch1_robin.dat"

END PROGRAM PlotWavefunctionCh1Robin
