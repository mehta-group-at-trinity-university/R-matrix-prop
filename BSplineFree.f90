!****************************************************************************************************
! B-spline single-box free-particle test, matching the 2019-era ("Independent of AlphaFactor",
! commit 38c0e72) and DeltaScat/DeltaScat.f90's validated box-1 recipe -- NOT the newer
! BuildBox1CombinedBasis/CalcGamLamBox1 machinery (which post-dates both of those and has never
! itself been validated against an exact result):
!   - Left=2, Right=2 (two independent boundary basis functions at each edge)
!   - NumOpenL=0 at the box's own left edge (variational elimination only -- a "natural"
!     boundary condition, no explicit Robin/log-derivative matching)
!   - xStart=0.0 exactly (matching DeltaScat.f90's own BPD1%xl=xStart=0.0d0) -- safe because
!     Gauss-Legendre quadrature points are strictly interior to each element and never land
!     exactly on R=0
!   - plain CalcGamLam (not CalcGamLamBox1)
!
! Potential data comes directly from OneDimChannels at the box's own Gauss-Legendre quadrature
! points (matching SVDBound.f90/SVDRmat.f90's approach) instead of SetAdiabaticPotential's
! tabulated-interpolation pipeline, sidestepping 3bodydata_free_5ch's Uad.dat domain limit
! (only tabulated from R=0.5) entirely -- this dataset's Uad is exact (verified), but the
! *tabulation* doesn't reach the origin, unlike DeltaScat's own analytic KMSFormulas.f evaluation.
!
! CORRECTION_SIGN toggles CalcGamLam's AlphaFactor*(AlphaFactor-EffDim+2)/R^2 correction term's
! sign for this test only (RMatPropCore.f90's own CalcGamLam is left untouched) -- set to -1d0 to
! match the current repo (established this session as the empirically-working choice when paired
! with BuildBox1CombinedBasis regularity matching), +1d0 to match the 2019-validated code (and
! this session's own differential-equation derivation).
!****************************************************************************************************
PROGRAM BSplineFree
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE Scattering
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: CORRECTION_SIGN = +1d0   ! -1d0 = current repo; +1d0 = 2019/derivation

  CHARACTER(LEN=64) DataDir
  TYPE(BPData) BPD
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull, Box1
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Leff(:), Threshold(:)
  DOUBLE PRECISION, ALLOCATABLE :: Uad(:,:,:), Psi(:,:,:), S_out(:,:), Rflat(:)
  INTEGER, ALLOCATABLE :: flatIdx(:,:)

  DOUBLE PRECISION :: xMin1D, xMax1D, Shift1D, muLocal, EnergyLocal, Rmax
  INTEGER :: Order1D, Left1D, Right1D, xNumPoints1D, LegPoints1D
  INTEGER :: OrderLocal, xNumPoints, NumOpenChannels
  INTEGER :: PsiDim, HalfBandWidth1D, Ufile
  INTEGER :: i, j, kx, lx, mch, RSteps

  DataDir = '3bodydata_free_5ch'
  OrderLocal = 5
  xNumPoints = 40
  LegPoints1D = 10
  Order1D = 5
  Left1D = 1
  Right1D = 0
  xNumPoints1D = 200
  xMin1D = 0.d0
  xMax1D = 1.570796326794897d0
  Shift1D = -1.d0
  Rmax = 10.0d0
  EnergyLocal = 1.0d0

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = LegPoints1D
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  OPEN(unit=8,file=TRIM(DataDir)//'/masses.dat',status='old')
  READ(8,*)
  READ(8,*) muLocal
  CLOSE(8)
  reducedmass = muLocal

  NumChannels = 5
  ALLOCATE(Leff(NumChannels),Threshold(NumChannels))
  OPEN(unit=9,file=TRIM(DataDir)//'/FitLeff.data',status='old')
  READ(9,*)
  READ(9,*)
  READ(9,*)
  DO i = 1,NumChannels
     READ(9,*) j, Threshold(i), Leff(i)
  ENDDO
  CLOSE(9)
  Threshold = 2d0*reducedmass*Threshold

  EffDim = 2.0d0
  AlphaFactor = 0.5d0
  Order = OrderLocal

  NumOpenChannels = 0
  DO i = 1,NumChannels
     IF (2d0*reducedmass*EnergyLocal.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
  ENDDO
  WRITE(6,*) "NumChannels=",NumChannels," NumOpenChannels=",NumOpenChannels," Leff=",Leff

  !--- box: Left=2/Right=2, xl=0.0 exactly, plain CalcGamLam, NumOpenL=0 (2019/DeltaScat style) ---
  BPD%NumChannels = NumChannels
  BPD%Order = OrderLocal
  BPD%Left = 2
  BPD%Right = 2
  BPD%xNumPoints = xNumPoints
  BPD%kl = 1d20
  BPD%kr = 1d20
  CALL AllocateBPD(BPD)
  BPD%xl = 0.0d0
  BPD%xr = Rmax
  CALL GridMaker(BPD%xPoints,BPD%xNumPoints,BPD%xl,BPD%xr,"quadratic")
  CALL MakeBasis(BPD)
  WRITE(6,*) "MakeBasis done, xDim=",BPD%xDim

  !--- fill BPD%x (Gauss-Legendre quadrature points) directly, matching SetAdiabaticPotential's
  !    own R=xScale*xLeg+xScaledZero formula, but WITHOUT calling it (avoids its tabulated-
  !    interpolation domain limit) ---
  RSteps = LegPoints*(BPD%xNumPoints-1)
  ALLOCATE(Rflat(RSteps),flatIdx(LegPoints,BPD%xNumPoints-1))
  i = 0
  DO kx = 1,BPD%xNumPoints-1
     DO lx = 1,LegPoints
        i = i+1
        BPD%x(lx,kx) = 0.5d0*(BPD%xPoints(kx+1)-BPD%xPoints(kx))*xLeg(lx) &
             + 0.5d0*(BPD%xPoints(kx+1)+BPD%xPoints(kx))
        Rflat(i) = BPD%x(lx,kx)
        flatIdx(lx,kx) = i
     ENDDO
  ENDDO

  !--- Uad at every quadrature point, direct from OneDimChannels (no interpolation) ---
  PsiDim = xNumPoints1D + Order1D - 3
  IF (Left1D.EQ.2)  PsiDim = PsiDim + 1
  IF (Right1D.EQ.2) PsiDim = PsiDim + 1
  HalfBandWidth1D = Order1D
  ALLOCATE(Uad(RSteps,NumChannels,2),Psi(RSteps,PsiDim,NumChannels),S_out(HalfBandWidth1D+1,PsiDim))
  Ufile = 71
  OPEN(unit=Ufile,file=TRIM(DataDir)//'/BSplineFreeUad.dat',status='replace')
  CALL OneDimChannels(NumChannels,PsiDim,0,0,LegendreFile,LegPoints1D,Shift1D, &
       Order1D,Left1D,Right1D,reducedmass,xNumPoints1D,xMin1D,xMax1D, &
       Rflat,RSteps,Ufile,72,73,74,75,Uad,Psi,S_out,DataDir)
  CLOSE(Ufile)
  WRITE(6,*) "OneDimChannels done (",RSteps," quadrature points, direct evaluation)."

  BPD%Pot = 0d0
  BPD%Pcoup = 0d0
  DO kx = 1,BPD%xNumPoints-1
     DO lx = 1,LegPoints
        DO mch = 1,NumChannels
           BPD%Pot(mch,mch,lx,kx) = Uad(flatIdx(lx,kx),mch,1)
        ENDDO
     ENDDO
  ENDDO

  BPD%lam(1:NumChannels) = Leff(1:NumChannels)

  EIG%MatrixDim = BPD%MatrixDim
  CALL AllocateEIG(EIG)
  CALL CalcGamLamSimple(BPD,EIG,0,NumOpenChannels,CORRECTION_SIGN)
  CALL CombineGam(EIG,EnergyLocal)

  Box1%NumOpenL = 0
  Box1%NumOpenR = NumOpenChannels
  CALL AllocateBox(Box1)
  CALL InitZeroBox(Box1)

  ALLOCATE(evalRed(Box1%betaMax),evecRed(Box1%betaMax,Box1%betaMax),LamooDiag(Box1%betaMax))
  CALL PartitionAndEliminate(BPD,EIG,Box1%NumOpenL,Box1%NumOpenR,evalRed,evecRed,LamooDiag)
  WRITE(6,*) "evalRed=",evalRed

  Bnull%NumOpenL = 0
  Bnull%NumOpenR = 0
  CALL AllocateBox(Bnull)
  CALL InitZeroBox(Bnull)
  CALL BoxMatch(Bnull, Box1, BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
  CALL AllocateScat(SD,NumChannels)
  CALL CalcK(Box1,BPD,SD,reducedmass,EffDim,AlphaFactor,EnergyLocal,Threshold)
  WRITE(6,*) "CORRECTION_SIGN=",CORRECTION_SIGN
  WRITE(6,*) "Full K diagonal = ", (SD%K(i,i), i=1,NumChannels)
  WRITE(6,*) "Full K off-diag max = ", MAXVAL(ABS(SD%K)) - MAXVAL(ABS((/(SD%K(i,i),i=1,NumChannels)/)))

CONTAINS

  ! Local copy of CalcGamLam with a sign-toggle on the AlphaFactor correction term (see
  ! header). Everything else identical to RMatPropCore.f90's CalcGamLam (channel-diagonal
  ! part only -- this test has no off-diagonal Q coupling or P coupling, both exactly zero
  ! for the free particle, so those loops are omitted).
  SUBROUTINE CalcGamLamSimple(BPD,EIG,NumOpenL,NumOpenR,csign)
    USE DataStructures
    USE Quadrature
    USE GlobalVars
    IMPLICIT NONE
    TYPE(BPData), INTENT(in) :: BPD
    TYPE(GenEigVal) :: EIG
    INTEGER, INTENT(in) :: NumOpenL, NumOpenR
    DOUBLE PRECISION, INTENT(in) :: csign
    INTEGER ix,ixp,lx,kx,mch
    DOUBLE PRECISION a, ax, bx, xIntScale(BPD%xNumPoints), TempG, TempS, xScaledZero

    EIG%Gam0 = 0d0
    EIG%Overlap = 0d0
    EIG%Lam = 0d0

    DO ix = 1,BPD%xDim
       DO ixp = MAX(1,ix-BPD%Order),MIN(BPD%xDim,ix+BPD%Order)
          DO mch = 1,BPD%NumChannels
             DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                ax = BPD%xPoints(kx)
                bx = BPD%xPoints(kx+1)
                xIntScale(kx) = 0.5d0*(bx-ax)
                xScaledZero = 0.5d0*(bx+ax)
                TempG = 0.0d0
                TempS = 0.0d0
                DO lx = 1,LegPoints
                   a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                   TempG = TempG - a*(BPD%ux(lx,kx,ix)*BPD%ux(lx,kx,ixp))
                   a = a*BPD%u(lx,kx,ix)*BPD%u(lx,kx,ixp)*2d0*reducedmass
                   TempS = TempS + a
                   TempG = TempG - a*BPD%Pot(mch,mch,lx,kx)
                   TempG = TempG + csign*a*AlphaFactor*(AlphaFactor-EffDim+2d0)/(2d0*reducedmass*BPD%x(lx,kx)**2)
                ENDDO
                EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) = &
                     EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) + TempG
                EIG%Overlap((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) = &
                     EIG%Overlap((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) + TempS
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    EIG%Lam = 0.d0
    IF (BPD%Right.EQ.2) THEN
       DO mch = 1,NumOpenR
          EIG%Lam((BPD%xDim-1)*BPD%NumChannels+mch,(BPD%xDim-1)*BPD%NumChannels+mch) = &
               BPD%xr**(EffDim-1d0-2d0*AlphaFactor)
       ENDDO
    ENDIF
    DO mch = 1,NumOpenL
       EIG%Lam(mch,mch) = BPD%xl**(EffDim-1d0-2d0*AlphaFactor)
    ENDDO

  END SUBROUTINE CalcGamLamSimple

END PROGRAM BSplineFree
