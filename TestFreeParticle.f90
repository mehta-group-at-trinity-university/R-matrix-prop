!****************************************************************************************************
! Cleanest possible diagnostic: single box [xMin~0, xMax], Left=0 (Dirichlet, u(0)=0) via ordinary
! CalcGamLam (no combined-basis/kLeft machinery needed for Left=0), U=Q=0 (set after
! SetAdiabaticPotential), and the wavefunction-reduction term ablated via RMatPropCore.f90's
! ABLATE_WAVEFN_REDUCTION_TERM switch. With both the potential and that term zeroed, the equation
! is exactly u''=-k^2 u (k^2=2*mu*E), and the regular (u(0)=0) solution is EXACTLY sin(kR), giving
! log-derivative L(R) = k*cos(kR)/sin(kR) = k*cot(kR) -- a clean, closed-form target with no
! near-origin subtlety of any kind.
!****************************************************************************************************
PROGRAM TestFreeParticle
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
  DOUBLE PRECISION, PARAMETER :: TestCoeff = 0.25d0
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
  AlphaFactor = alpha
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

  BPD1%NumChannels = NumChannels
  BPD1%Order = OrderLocal
  BPD1%Left = 0    ! Dirichlet: u(xMin)=0 automatically, no combined basis/kLeft needed
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
  ! Keep the REAL U1(R)+Q11(R)/(2mu) potential as-is (no zeroing, no TestCoeff folding).
  ! ABLATE_WAVEFN_REDUCTION_TERM=.TRUE. in RMatPropCore.f90 keeps the separate
  ! -0.25/(2mu R^2) line OFF -- direct log-derivative comparison against NDSolve solving
  ! the SAME (real-potential-only, no wavefn-reduction term) equation.

  EIG%MatrixDim = BPD1%MatrixDim
  CALL AllocateEIG(EIG)
  CALL CalcGamLam(BPD1,EIG,0,NumOpenChannels)   ! ordinary CalcGamLam -- no CalcGamLamBox1 needed for Left=0

  ! TEMPORARY diagnostic: dump everything needed to independently recompute Gam0(100,100)
  ! (a generic interior diagonal element, far from any boundary) in Python, bypassing
  ! PartitionAndEliminate/BoxMatch/CalcK entirely.
  BLOCK
    INTEGER, PARAMETER :: ixTest = 100
    INTEGER :: kxT, lxT
    WRITE(6,*) "=== MATRIX ELEMENT DIAGNOSTIC ==="
    WRITE(6,*) "Gam0(ixTest,ixTest) = ", EIG%Gam0(ixTest,ixTest)
    WRITE(6,*) "Overlap(ixTest,ixTest) = ", EIG%Overlap(ixTest,ixTest)
    WRITE(6,*) "kxMin,kxMax = ", BPD1%kxMin(ixTest,ixTest), BPD1%kxMax(ixTest,ixTest)
    OPEN(unit=88,file='MatrixElementDiag.dat',status='replace')
    WRITE(88,'(A)') '# kx  lx  xPoints(kx)  xPoints(kx+1)  R  wLeg  u  ux  Pot(1,1)'
    DO kxT = BPD1%kxMin(ixTest,ixTest),BPD1%kxMax(ixTest,ixTest)
       DO lxT = 1,LegPoints
          WRITE(88,'(2I6,1P,7E20.12)') kxT, lxT, BPD1%xPoints(kxT), BPD1%xPoints(kxT+1), &
               BPD1%x(lxT,kxT), wLeg(lxT), BPD1%u(lxT,kxT,ixTest), BPD1%ux(lxT,kxT,ixTest), &
               BPD1%Pot(1,1,lxT,kxT)
       ENDDO
    ENDDO
    CLOSE(88)
    WRITE(6,*) "Wrote MatrixElementDiag.dat"
  END BLOCK

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

  OPEN(unit=99,file='FEM_wavefunction_freeparticle.dat',status='replace')
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
  WRITE(6,*) "Wrote FEM_wavefunction_freeparticle.dat"

END PROGRAM TestFreeParticle
