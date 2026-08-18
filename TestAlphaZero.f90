!****************************************************************************************************
! Sidesteps the wavefunction-reduction-term investigation entirely instead of ablating it: uses
! AlphaFactor=0 (the OTHER valid convention noted in RMatPropCore.f90's own GlobalVars comment --
! "u(R) = R^(AlphaFactor) Psi(R)... typical choice is either AlphaFactor=0 ... or (EffDim-1)/2"),
! for which AlphaFactor*(AlphaFactor-EffDim+2)=0 identically -- the term vanishes at the source,
! with NO code modification/ablation switch needed. This exercises the completely unmodified
! production CalcGamLam/PartitionAndEliminate/BoxMatch pipeline, just with the kinetic weight
! R^(EffDim-1-2*AlphaFactor)=R^(EffDim-1)=R^1 (D=2) instead of R^0.
!
! With U=Q=0 (set after SetAdiabaticPotential) and AlphaFactor=0, D=2, the TRUE (unreduced)
! equation is exactly the ordinary Bessel equation of order 0:
!   Psi'' + (1/R)Psi' + k^2 Psi = 0   (k^2=2*mu*E)
! whose regular solution is Psi(R)=J0(kR), giving log-derivative L(R)=Psi'/Psi=-k*J1(kR)/J0(kR)
! -- a clean, closed-form target with NO sqrt(R) prefactor (unlike the AlphaFactor=0.5 tests).
!
! Left=2/NumOpenL=0 (open/variational-elimination box 1, unchanged from the rest of this repo's
! convention) -- deliberately NOT Left=0/Dirichlet or Left=3/Robin, since Left=2 doesn't hard-code
! a specific near-origin power law and should find the correct (Neumann-like, Psi'(0)=0) regular
! solution on its own.
!****************************************************************************************************
PROGRAM TestAlphaZero
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
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim, EnergyLocal, xMinLocal, xMaxLocal
  INTEGER NumDataPoints, NumOpenChannels, xNumPoints, LegPointsLocal, i, kx, lx
  INTEGER OrderLocal

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
  AlphaFactor = 0d0   ! <<< the key change: overrides the dataset's own alpha=0.5
  EffDim = ddim
  Order = OrderLocal

  ALLOCATE(CoulombC(NumChannels))
  CoulombC = 0d0

  Threshold = 2d0*reducedmass*Threshold
  NumOpenChannels = 0
  DO i = 1,NumChannels
     IF (2d0*reducedmass*EnergyLocal.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
  ENDDO
  WRITE(6,*) "NumOpenChannels=",NumOpenChannels
  WRITE(6,*) "AlphaFactor=",AlphaFactor," EffDim=",EffDim

  BPD1%NumChannels = NumChannels
  BPD1%Order = OrderLocal
  BPD1%Left = 2    ! open/variational-elimination, matching this repo's standard box-1 convention
  BPD1%Right = 2
  BPD1%xNumPoints = xNumPoints
  BPD1%kl = 1d20
  BPD1%kr = 1d20
  CALL AllocateBPD(BPD1)
  BPD1%xl = xMinLocal
  BPD1%xr = xMaxLocal
  CALL GridMaker(BPD1%xPoints,BPD1%xNumPoints,BPD1%xl,BPD1%xr,"quadratic")
  CALL MakeBasis(BPD1)
  BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
  CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)
  ! Keep the REAL U1(R)+Q11(R)/(2mu) potential as-is (no zeroing) -- test AlphaFactor=0 as a
  ! genuine workaround for the actual physics, not just the U=Q=0 toy problem.

  EIG%MatrixDim = BPD1%MatrixDim
  CALL AllocateEIG(EIG)
  CALL CalcGamLam(BPD1,EIG,0,NumOpenChannels)   ! completely unmodified production code
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
  WRITE(6,*) "K (from CalcK) = ", SD%K(1,1)
  WRITE(6,*) "Box1%Zf(1,1) [u(xMax)]  = ", Box1%Zf(1,1)
  WRITE(6,*) "Box1%Zfp(1,1) [u'(xMax)] = ", Box1%Zfp(1,1)
  WRITE(6,*) "xMax = ", BPD1%xr, " reducedmass = ", reducedmass, " EnergyLocal = ", EnergyLocal

  OPEN(unit=99,file='FEM_wavefunction_alphazero.dat',status='replace')
  DO kx = 1,BPD1%xNumPoints-1
     DO lx = 1,LegPoints
        BLOCK
          DOUBLE PRECISION Rval, uval
          Rval = BPD1%x(lx,kx)
          uval = 0d0
          DO ix = 1,BPD1%xDim
             uval = uval + cfull((ix-1)*NumChannels+1)*BPD1%u(lx,kx,ix)
          ENDDO
          WRITE(99,*) Rval, uval
        END BLOCK
     ENDDO
  ENDDO
  CLOSE(99)
  WRITE(6,*) "Wrote FEM_wavefunction_alphazero.dat"

END PROGRAM TestAlphaZero
