!****************************************************************************************************
! Free-particle (V2Depth=0) variant of PlotWavefunction.f90: reconstructs the actual FEM
! (B-spline) wavefunction u(R) for the LOWEST channel (Leff=1) of the '3bodydata_free_5ch'
! dataset, single box (NumBoxes=1, regularity at xMinLocal via BuildBox1CombinedBasis),
! at a positive scattering energy where all 5 channels are open. Written to compare
! directly against the exact free-particle solution u(R) = sqrt(k*R)*J_1(k*R) (the
! reduced-wavefunction alpha=0.5 form of the regular hyperspherical Bessel function,
! order=Leff=1 for d=2 -- see hyperrjry/besselnew.f). Not a general-purpose tool.
!
! Differs from PlotWavefunction.f90 in three ways:
!  - DataDir/xMinLocal/xMaxLocal/EnergyLocal target the free-particle dataset instead of
!    the delta-function one (xMaxLocal=9.99, not 10.0 -- FitRangeMax=10.0 exactly, and
!    letting a quadrature point land at R=10.000000000000002 crashes the interpolator).
!  - CoulombC is forced to zero unconditionally: PlotWavefunction.f90's Threshold<0 check
!    (a delta-function-specific "genuine Coulomb-type bound-pair channel" convention, see
!    RMATPROPAdiabatic.f90's own comment on that assumption) would otherwise misfire here
!    -- this dataset's Threshold(:) are all ~1e-13 negative purely from ARPACK numerical
!    noise, not a real bound pair.
!  - With 5 open channels (not 1), PartitionAndEliminateFull's reduced eigenproblem has
!    betaMax=5 eigenvectors, not a single one -- but since P=Q=0 exactly for this test,
!    Goo/Lamoo are exactly block-diagonal by channel, so each eigenvector is (up to
!    normalization) a pure single-channel excitation. The channel-1 column is selected by
!    searching for the eigenvector with its dominant component at row 1.
!****************************************************************************************************
PROGRAM PlotWavefunctionFree
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
  INTEGER NumOpen, Ncc, ix, betaCh1

  DataDir = '3bodydata_free_5ch'
  OrderLocal = 5
  xNumPoints = 40
  LegPointsLocal = 10
  ! 3bodydata_free_5ch's Uad.dat is only tabulated from R=0.5, AND ReadAdiabaticData
  ! unconditionally drops the first tabulated point when building UInterp/PInterp/
  ! QInterp (NPts=NumDataPoints-1, x=Rdata(2:...) -- a Runge-ringing workaround tuned
  ! for the delta-function dataset's near-origin data, unconditionally applied here
  ! too) -- so the true usable interpolation domain starts at the SECOND tabulated R,
  ! 0.59596, not 0.5. xMinLocal must stay above that or SetAdiabaticPotential's
  ! UInterp/QInterp calls hit dbsval's domain-edge error. BuildBox1CombinedBasis's
  ! log-derivative matching is an EXACT Bessel-connection match at xMinLocal (not a
  ! near-origin approximation), so 0.65 is just as valid a regularity-matching point
  ! as a smaller R would be.
  xMinLocal = 0.65d0
  xMaxLocal = 9.99d0
  EnergyLocal = 1.0d0

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
  CoulombC = 0d0   ! no genuine Coulomb-type channel in the free-particle test

  Threshold = 2d0*reducedmass*Threshold

  NumOpenChannels = 0
  DO i = 1,NumChannels
     IF (2d0*reducedmass*EnergyLocal.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
  ENDDO
  WRITE(6,*) "NumChannels=",NumChannels," NumOpenChannels=",NumOpenChannels
  WRITE(6,*) "Leff=",Leff
  WRITE(6,*) "Threshold(raw, 2*mu*E units)=",Threshold

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
  WRITE(6,*) "CHECKPOINT 1: GridMaker done, xPoints(xNumPoints)=",BPD1%xPoints(BPD1%xNumPoints)
  FLUSH(6)
  CALL MakeBasis(BPD1)
  WRITE(6,*) "CHECKPOINT 2: MakeBasis done"
  FLUSH(6)
  BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
  CALL SetAdiabaticPotential(BPD1,muLocal,UInterp,PInterp,QInterp,CoulombC)
  WRITE(6,*) "CHECKPOINT 3: SetAdiabaticPotential done"
  FLUSH(6)

  ALLOCATE(u1Box1(LegPoints,BPD1%xNumPoints,NumChannels))
  ALLOCATE(ux1Box1(LegPoints,BPD1%xNumPoints,NumChannels))
  ALLOCATE(kLeftBox1(NumChannels))
  CALL BuildBox1CombinedBasis(BPD1,Leff,reducedmass,EffDim,AlphaFactor,EnergyLocal,Threshold, &
       u1Box1,ux1Box1,kLeftBox1,CoulombC,UInterp,QInterp)
  WRITE(6,*) "CHECKPOINT 4: BuildBox1CombinedBasis done"
  FLUSH(6)

  EIG%MatrixDim = BPD1%MatrixDim
  CALL AllocateEIG(EIG)
  CALL CalcGamLamBox1(BPD1,EIG,0,NumOpenChannels,u1Box1,ux1Box1,kLeftBox1)
  WRITE(6,*) "CHECKPOINT 5: CalcGamLamBox1 done"
  FLUSH(6)
  CALL CombineGam(EIG,EnergyLocal)
  WRITE(6,*) "CHECKPOINT 6: CombineGam done"
  FLUSH(6)

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
  WRITE(6,*) "CHECKPOINT 7: PartitionAndEliminateFull done"
  FLUSH(6)

  WRITE(6,*) "NumOpen=",NumOpen," Ncc=",Ncc
  WRITE(6,*) "evalRed=",evalRed

  ! P=Q=0 exactly for this test, so Goo/Lamoo are exactly block-diagonal by channel --
  ! each evecRed(:,beta) column is (up to normalization) a pure single-channel
  ! excitation. Pick the column whose dominant component sits at row 1 (channel 1).
  betaCh1 = 1
  DO i = 1,NumOpen
     IF (ABS(evecRed(1,i)).GT.ABS(evecRed(1,betaCh1))) betaCh1 = i
  ENDDO
  WRITE(6,*) "Selected beta (channel-1-pure mode) = ",betaCh1
  WRITE(6,*) "evecRed(:,betaCh1) = ",evecRed(:,betaCh1)

  ! Reconstruct the physical solution's full spline coefficient vector for this mode.
  ALLOCATE(copen(NumOpen),cclosed(Ncc),cfull(BPD1%MatrixDim))
  copen = evecRed(:,betaCh1)
  cclosed = -MATMUL(Xout,copen)
  DO i = 1,NumOpen
     cfull(OpenIdxOut(i)) = copen(i)
  ENDDO
  DO i = 1,Ncc
     cfull(ClosedIdxOut(i)) = cclosed(i)
  ENDDO

  ! Cross-check: compute K via the normal BoxMatch+CalcK pathway.
  Bnull%NumOpenL = 0
  Bnull%NumOpenR = 0
  CALL AllocateBox(Bnull)
  CALL InitZeroBox(Bnull)
  CALL BoxMatch(Bnull, Box1, BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
  CALL AllocateScat(SD,NumChannels)
  CALL CalcK(Box1,BPD1,SD,reducedmass,EffDim,AlphaFactor,EnergyLocal,Threshold)
  WRITE(6,*) "Full K diagonal = ", (SD%K(i,i), i=1,NumChannels)
  WRITE(6,*) "Box1%Zf(1,betaCh1) [u_ch1(xMax)]  = ", Box1%Zf(1,betaCh1)
  WRITE(6,*) "Box1%Zfp(1,betaCh1) [u_ch1'(xMax)] = ", Box1%Zfp(1,betaCh1)
  WRITE(6,*) "xMax = ", BPD1%xr, " reducedmass = ", reducedmass, " EnergyLocal = ", EnergyLocal

  ! Evaluate u(R) = sum_ix cfull((ix-1)*NumChannels+1)*basis_ix(R) -- channel 1 only.
  OPEN(unit=99,file='FEM_wavefunction_free.dat',status='replace')
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
  WRITE(6,*) "Wrote FEM_wavefunction_free.dat"

END PROGRAM PlotWavefunctionFree
