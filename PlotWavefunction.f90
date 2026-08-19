!****************************************************************************************************
! Standalone diagnostic: reconstructs the actual FEM (B-spline) wavefunction u(R) for a single
! energy, single box (NumBoxes=1, regularity at xMin) calculation, using PartitionAndEliminateFull's
! Xout/OpenIdxOut/ClosedIdxOut to recover the FULL spline coefficient vector (not just the box's
! boundary amplitude Zf/Zfp that CalcK uses). Written to compare directly against an independent
! RK45 integration and the near-origin Whittaker/Coulomb power series, for a single hard-coded
! DataDir/EnergyLocal (see below) -- not a general-purpose tool.
!****************************************************************************************************
PROGRAM PlotWavefunction
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
  DOUBLE PRECISION, ALLOCATABLE :: u1Box1(:,:,:), ux1Box1(:,:,:), kLeftBox1(:)
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

  DataDir = '3bodydata_delta_5ch'
  OrderLocal = 5
  xNumPoints = 40
  LegPointsLocal = 10
  xMinLocal = 0.01d0
  xMaxLocal = 20.0d0
  EnergyLocal = -0.99d0

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
  DO i = 1,NumChannels
     IF (Threshold(i).LT.0d0) THEN
        BLOCK
          DOUBLE PRECISION gamma0, B2, a2
          gamma0 = -12d0/(3d0**0.25d0*DSQRT(2d0)*Pi)
          B2 = -Threshold(i)
          a2 = 1d0/DSQRT(2d0*0.5d0*B2)
          CoulombC(i) = -gamma0/(2d0*reducedmass*a2)
        END BLOCK
     ENDIF
  ENDDO

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
  BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
  CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)

  ALLOCATE(u1Box1(LegPoints,BPD1%xNumPoints,NumChannels))
  ALLOCATE(ux1Box1(LegPoints,BPD1%xNumPoints,NumChannels))
  ALLOCATE(kLeftBox1(NumChannels))
  CALL BuildBox1CombinedBasis(BPD1,Leff,reducedmass,EffDim,AlphaFactor,EnergyLocal,Threshold, &
       u1Box1,ux1Box1,kLeftBox1,CoulombC,UInterp,QInterp)

  EIG%MatrixDim = BPD1%MatrixDim
  CALL AllocateEIG(EIG)
  CALL CalcGamLamBox1(BPD1,EIG,0,NumOpenChannels,u1Box1,ux1Box1,kLeftBox1)
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

  ! Reconstruct the physical solution's full spline coefficient vector (unique up to
  ! overall scale, since NumOpen=1 here -- only one open DOF, the box's right edge).
  ALLOCATE(copen(NumOpen),cclosed(Ncc),cfull(BPD1%MatrixDim))
  copen = evecRed(:,1)
  cclosed = -MATMUL(Xout,copen)
  DO i = 1,NumOpen
     cfull(OpenIdxOut(i)) = copen(i)
  ENDDO
  DO i = 1,Ncc
     cfull(ClosedIdxOut(i)) = cclosed(i)
  ENDDO

  ! Cross-check: compute K via the normal BoxMatch+CalcK pathway, to confirm this
  ! reconstruction is consistent with the code's own reported K-matrix.
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
  WRITE(6,*) "xMax = ", BPD1%xr, " reducedmass = ", reducedmass, " EnergyLocal = ", EnergyLocal, &
       " Threshold(1) = ", Threshold(1)

  ! Evaluate u(R) = sum_ix cfull((ix-1)*NumChannels+1)*basis_ix(R) at all quadrature
  ! points spanning the box -- channel 1 only (the only open channel here). cfull's
  ! index follows the SAME global flat convention as EIG%Gam0/Overlap/Lam everywhere
  ! else in RMatPropCore.f90 ((ix-1)*NumChannels+channel, e.g. CalcGamLam at
  ! RMatPropCore.f90:309) -- indexing it by ix alone (as an earlier version of this
  ! loop did) silently picked up a scramble of different channels' coefficients at
  ! different spatial points whenever NumChannels>1, not channel 1's own coefficients.
  OPEN(unit=99,file='FEM_wavefunction.dat',status='replace')
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
  WRITE(6,*) "Wrote FEM_wavefunction.dat"

END PROGRAM PlotWavefunction
