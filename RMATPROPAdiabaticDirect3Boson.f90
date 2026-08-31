!****************************************************************************************************
! No-interpolation analog of RMATPROPAdiabatic.f90 for the equal-mass 3-identical-boson sech^2
! problem (domain [0,pi/6], Left=Right=1): fills every box's potential via
! SechAdiabaticPotentialDirect.f90's SetAdiabaticPotentialDirect (real per-quadrature-point ARPACK
! diagonalization, SechFunctionPlugins physics) instead of AdiabaticPotential.f90's
! tabulate-then-interpolate pipeline -- eliminating the interpolant-domain-floor/near-origin-
! resolution problem entirely, exactly mirroring DeltaScat.f90's own no-interpolation design.
! Threshold/Leff/mu/alpha/ddim are still read from the existing Fit.data/FitLeff.data/masses.dat
! (large-R asymptotic fits, unaffected by near-origin potential evaluation method).
!
! Structurally identical to RMATPROPAdiabatic.f90 EXCEPT: (1) no CoulombC/ApplyCoulombC at all
! (this dataset never needs it -- see that file's own header comment on why), and (2) the last box
! -- which RMATPROPAdiabatic.f90 rebuilds via SetAdiabaticPotential every single energy, harmless
! there since interpolation lookup is cheap -- is instead computed ONCE via a dedicated BPDLast
! (mirroring BPD1's own one-time-build-then-reuse pattern), since SetAdiabaticPotentialDirect's
! ARPACK diagonalization is NOT cheap enough to redo NumEnergies times for no physics reason (the
! potential itself doesn't depend on Energy at all -- only CalcGamLam's NumOpenR argument does, and
! that's a cheap dense-linear-algebra refresh, not a re-diagonalization).
!****************************************************************************************************
PROGRAM main
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE SechAdiabaticPotentialDirect

  IMPLICIT NONE
  CHARACTER(LEN=64) DataDir
  CHARACTER(LEN=16) EnergyGridType
  INTEGER ios
  TYPE(BPData) BPD,BPD0,BPD1,BPDLast
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: xprim(:)
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Boxes(:,:,:), OverlapBoxes(:,:,:), LamBoxes(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Gam0Box1(:,:), OverlapBox1(:,:), LamBox1(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  DOUBLE PRECISION muLocal, alpha, ddim, dum
  DOUBLE PRECISION Emin, Emax, qOurs, Kelem, rawDelta, delta, prevDelta
  INTEGER NumOpenChannels, xNumPoints, LegPointsLocal, NumDataPoints
  INTEGER i, j, n, iBox, lx, kx, NumEnergies, ie, nBranch, PhaseFile, KMatFile, RecombFile
  LOGICAL havePrevDelta
  COMPLEX*16, ALLOCATABLE :: Amat(:,:), Smat(:,:)
  DOUBLE PRECISION recombProb
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0

  !----------------------------------------------------------------------
  ! Read run parameters from the same .inp format RMATPROPAdiabatic.f90 uses (minus the
  ! trailing ApplyCoulombC line, never needed here); Threshold/Leff/mu/alpha/ddim still come
  ! from DataDir's own Fit.data/FitLeff.data/masses.dat -- only Uad.dat/Pmat.dat/Qmat.dat (the
  ! near-origin-sensitive tabulated curves) are bypassed.
  !----------------------------------------------------------------------
  OPEN(unit=7,file='RMATPROPAdiabaticDirect3Boson.inp',status='old')
  READ(7,*)
  READ(7,*) DataDir
  READ(7,*)
  READ(7,*)
  READ(7,*) Order, xNumPoints, LegPointsLocal
  READ(7,*)
  READ(7,*)
  READ(7,*) NumBoxes, xStart, xEnd
  READ(7,*)
  READ(7,*)
  READ(7,*) Emin, Emax, NumEnergies
  READ(7,*,IOSTAT=ios) EnergyGridType
  IF (ios.NE.0) EnergyGridType = "linear"
  CLOSE(7)

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = LegPointsLocal
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  !--- Fit.data header: NumChannels, NumDataPoints, alpha, ddim (same format ReadAdiabaticData
  !    uses -- NumDataPoints is meaningless here, no tabulated array, but Fit.data still has
  !    the column so it must still be read positionally). ---
  OPEN(unit=21,file=TRIM(DataDir)//'/Fit.data',status='old')
  READ(21,*) NumChannels, NumDataPoints, alpha, ddim
  CLOSE(21)
  OPEN(unit=22,file=TRIM(DataDir)//'/masses.dat',status='old')
  READ(22,*)
  READ(22,*) muLocal
  CLOSE(22)
  ALLOCATE(Threshold(NumChannels),Leff(NumChannels))
  OPEN(unit=26,file=TRIM(DataDir)//'/FitLeff.data',status='old')
  READ(26,*)
  READ(26,*)
  READ(26,*)
  DO n = 1,NumChannels
     READ(26,*) j, Threshold(n), dum, Leff(n)
  ENDDO
  CLOSE(26)

  reducedmass = muLocal
  AlphaFactor = alpha
  EffDim = ddim
  Threshold = 2d0*reducedmass*Threshold
  WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," ddim = ",EffDim," mu = ",reducedmass
  WRITE(6,*) "Energy scan: Emin = ",Emin," Emax = ",Emax," NumEnergies = ",NumEnergies

  PhaseFile = 50
  OPEN(unit=PhaseFile,file='PhaseShift.dat',status='replace')
  WRITE(PhaseFile,'(A)') '# q  Energy  K  delta(rad, unwrapped, no anchor/offset applied)'
  WRITE(PhaseFile,'(A)') '@    xaxis label "q"'
  WRITE(PhaseFile,'(A)') '@    yaxis label "delta (rad)"'

  KMatFile = 51
  OPEN(unit=KMatFile,file='KMatrixFull.dat',status='replace')

  RecombFile = 52
  OPEN(unit=RecombFile,file='Recombination.dat',status='replace')
  WRITE(RecombFile,'(A)') '# Energy  NumOpenChannels  1-|S(1,1)|^2  (McGuire exact-unitarity test: should -> 0)'
  WRITE(RecombFile,'(A)') '@    xaxis label "Energy"'
  WRITE(RecombFile,'(A)') '@    yaxis label "1 - |S11|^2"'

  !---------------------------------------------------------------------
  ! ONE-TIME setup (energy-independent): box 1, interior boxes 2..NumBoxes-1, AND (unlike
  ! RMATPROPAdiabatic.f90) the last box too -- all built and their real potential evaluated via
  ! SetAdiabaticPotentialDirect exactly once here, never again inside the energy loop.
  !---------------------------------------------------------------------
  BPD1%NumChannels = NumChannels
  BPD1%Order = Order
  BPD1%Left = 2
  BPD1%Right = 2
  BPD1%xNumPoints = xNumPoints
  BPD1%kl = 1d20
  BPD1%kr = 1d20

  BPD0%NumChannels = NumChannels
  BPD0%Order = Order
  BPD0%Left = 2
  BPD0%Right = 2
  BPD0%xNumPoints = xNumPoints
  BPD0%kl = 1d20
  BPD0%kr = 1d20
  CALL AllocateBPD(BPD0)
  ALLOCATE(xprim(xNumPoints))
  CALL GridMakerLinear(BPD0%xNumPoints,0d0,1d0,xprim)
  BPD0%xPoints = xprim
  CALL MakeBasis(BPD0)
  BPD0%lam(1:NumChannels) = Leff(1:NumChannels)

  BPD = BPD0
  BPDLast = BPD0

  EIG%MatrixDim = BPD%MatrixDim
  CALL AllocateEIG(EIG)

  ! BOX BUILD ORDER IS ASCENDING-R ON PURPOSE.  SetAdiabaticPotentialDirect threads the
  ! eigenvector phase from one box into the next (OneDimChannelsGeneric's PhaseSeed), and that
  ! chain is only a continuity statement if consecutive calls are radially adjacent -- building
  ! 2..N-1 before box 1, as this once did, would seed box 1 from box N-1.  Keep box 1 first.
  WRITE(6,*) "Building box 1..."
  IF (NumBoxes.GE.2) THEN
     ALLOCATE(Gam0Box1(BPD%MatrixDim,BPD%MatrixDim))
     ALLOCATE(OverlapBox1(BPD%MatrixDim,BPD%MatrixDim))
     ALLOCATE(LamBox1(BPD%MatrixDim,BPD%MatrixDim))
     CALL AllocateBPD(BPD1)
     BPD1%xl = xStart
     BPD1%xr = xStart + (xEnd-xStart)*(1d0/DBLE(NumBoxes))**BoxSpacingPower
     CALL GridMaker(BPD1%xPoints,BPD1%xNumPoints,BPD1%xl,BPD1%xr,"quadratic")
     CALL MakeBasis(BPD1)
     BPD1%lam(1:NumChannels) = Leff(1:NumChannels)
     CALL SetAdiabaticPotentialDirect(BPD1,muLocal)
     CALL CalcGamLam(BPD1,EIG,0,NumChannels)
     Gam0Box1 = EIG%Gam0
     OverlapBox1 = EIG%Overlap
     LamBox1 = EIG%Lam
  ENDIF

  !---------------------------------------------------------------------
  ! Last box: built and its potential evaluated ONCE here via a dedicated BPDLast (NOT
  ! recomputed per energy -- see this file's own header for why). Only CalcGamLam's
  ! NumOpenR argument genuinely depends on Energy (via NumOpenChannels), and that's applied
  ! fresh each energy inside the loop below, cheaply, without re-touching BPDLast%Pot/Pcoup.
  !---------------------------------------------------------------------
  WRITE(6,*) "Building interior boxes (SetAdiabaticPotentialDirect, one ARPACK diagonalization per quadrature point)..."
  IF (NumBoxes.GE.3) THEN
     ALLOCATE(Gam0Boxes(BPD%MatrixDim,BPD%MatrixDim,2:NumBoxes-1))
     ALLOCATE(OverlapBoxes(BPD%MatrixDim,BPD%MatrixDim,2:NumBoxes-1))
     ALLOCATE(LamBoxes(BPD%MatrixDim,BPD%MatrixDim,2:NumBoxes-1))
     DO iBox = 2,NumBoxes-1
        BPD%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xPoints = BPD%xl + (BPD%xr-BPD%xl)*xprim
        DO kx = 1,BPD%xNumPoints
           DO lx = 1,LegPoints
              BPD%ux(lx,kx,1:BPD%xDim) = BPD0%ux(lx,kx,1:BPD%xDim)/(BPD%xr-BPD%xl)
           ENDDO
        ENDDO
        ! Fills Pot(m,n)=Qtilde_mn(R)/(2*mu)+delta_mn*U_m(R), Pcoup=P(R) at this box's own
        ! quadrature points via a fresh ARPACK diagonalization at each R (OneDimChannelsGeneric),
        ! not interpolation -- same target quantities as SetAdiabaticPotential, no tabulated data.
        CALL SetAdiabaticPotentialDirect(BPD,muLocal)
        ! Gam0_ij=Int[B_i'B_j' - 2*mu*Pot*B_iB_j - alpha(alpha-D+2)/(2*mu*R^2)*B_iB_j]
        ! R^(D-1-2alpha)dR (+ Pcoup coupling term); Overlap_ij=Int[B_iB_j]R^(D-1-2alpha)dR.
        CALL CalcGamLam(BPD,EIG,NumChannels,NumChannels)
        Gam0Boxes(:,:,iBox) = EIG%Gam0
        OverlapBoxes(:,:,iBox) = EIG%Overlap
        LamBoxes(:,:,iBox) = EIG%Lam
        WRITE(6,'(A,I5,A,I5)') "  interior box ",iBox," / ",NumBoxes-1
     ENDDO
  ENDIF

  WRITE(6,*) "Building last box..."
  ! BPDLast was already allocated by the "BPDLast = BPD0" deep-copy assignment above (line 148)
  ! -- same pattern as the interior-box loop's "BPD = BPD0", which likewise never calls
  ! AllocateBPD(BPD) again. A redundant CALL AllocateBPD(BPDLast) here previously crashed
  ! -fcheck=all with "Attempting to allocate already allocated variable".
  BPDLast%xl = xStart + (xEnd-xStart)*(DBLE(NumBoxes-1)/DBLE(NumBoxes))**BoxSpacingPower
  BPDLast%xr = xStart + (xEnd-xStart)*(DBLE(NumBoxes)/DBLE(NumBoxes))**BoxSpacingPower
  BPDLast%xPoints = BPDLast%xl + (BPDLast%xr-BPDLast%xl)*xprim
  DO kx = 1,BPDLast%xNumPoints
     DO lx = 1,LegPoints
        BPDLast%ux(lx,kx,1:BPDLast%xDim) = BPD0%ux(lx,kx,1:BPDLast%xDim)/(BPDLast%xr-BPDLast%xl)
     ENDDO
  ENDDO
  CALL SetAdiabaticPotentialDirect(BPDLast,muLocal)
  WRITE(6,*) "All boxes built. Starting energy scan..."

  havePrevDelta = .FALSE.

  DO ie = 1,NumEnergies
     IF (NumEnergies.GT.1) THEN
        IF (EnergyGridType.EQ."quadratic") THEN
           Energy = Emin + (Emax-Emin)*(DBLE(ie-1)/DBLE(NumEnergies-1))**2
        ELSE
           Energy = Emin + (Emax-Emin)*DBLE(ie-1)/DBLE(NumEnergies-1)
        ENDIF
     ELSE
        Energy = Emin
     ENDIF

     NumOpenChannels = 0
     DO i = 1,NumChannels
        IF (2d0*reducedmass*Energy.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
     ENDDO

     ALLOCATE(Boxes(NumBoxes))
     DO i = 1,NumBoxes
        Boxes(i)%NumOpenL = MERGE(0,NumChannels,i.EQ.1)
        IF (i.EQ.NumBoxes) THEN
           Boxes(i)%NumOpenR = NumOpenChannels
        ELSE
           Boxes(i)%NumOpenR = NumChannels
        ENDIF
        IF (i.EQ.1) THEN
           Boxes(i)%xl = xStart
        ELSE
           Boxes(i)%xl = Boxes(i-1)%xr
        ENDIF
        Boxes(i)%xr = xStart + (xEnd-xStart)*(DBLE(i)/DBLE(NumBoxes))**BoxSpacingPower
        CALL AllocateBox(Boxes(i))
        CALL InitZeroBox(Boxes(i))
     ENDDO

     Bnull%NumOpenL = 0
     Bnull%NumOpenR = 0
     CALL AllocateBox(Bnull)
     CALL InitZeroBox(Bnull)

     CALL AllocateScat(SD,NumChannels)

     !--- box 1 (only when NumBoxes.GE.2 -- BPD1/Gam0Box1 are only built in that case, see the
     !    one-time setup above. For NumBoxes==1 the sole box plays both the box-1 role
     !    (NumOpenL=0, matched against Bnull) AND the last-box role (NumOpenR=NumOpenChannels,
     !    potential from BPDLast) at once -- handled entirely in the "last box" section below,
     !    which branches on NumBoxes to match against Bnull instead of Boxes(iBox-1) in that case. ---
     IF (NumBoxes.GE.2) THEN
        EIG%Gam0 = Gam0Box1
        EIG%Overlap = OverlapBox1
        EIG%Lam = LamBox1
        ! Gam(E) = Gam0 + E*Overlap.
        CALL CombineGam(EIG,Energy)
        ALLOCATE(evalRed(Boxes(1)%betaMax),evecRed(Boxes(1)%betaMax,Boxes(1)%betaMax))
        ALLOCATE(LamooDiag(Boxes(1)%betaMax))
        ! Schur-eliminates interior/closed DOF (Omega_oo=Gam_oo-Gam_oc*Gam_cc^-1*Gam_co), then
        ! solves the small generalized eigenproblem Omega_oo*c=beta*Lam_oo*c on the boundary-
        ! retained ("open") DOF -- evalRed=beta, evecRed=c.
        CALL PartitionAndEliminate(BPD1,EIG,Boxes(1)%NumOpenL,Boxes(1)%NumOpenR,evalRed,evecRed,LamooDiag)
        ! Builds Z,Z'=-beta*Z and matches to the previous box via a generalized eigenvalue match
        ! enforcing R-matrix/log-derivative continuity (Burke thesis Eqs. 3.13-3.14).
        CALL BoxMatch(Bnull, Boxes(1), BPD1, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        DEALLOCATE(evalRed,evecRed,LamooDiag)
     ENDIF

     !--- interior boxes ---
     DO iBox = 2,NumBoxes-1
        BPD%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
        BPD%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
        EIG%Gam0 = Gam0Boxes(:,:,iBox)
        EIG%Overlap = OverlapBoxes(:,:,iBox)
        EIG%Lam = LamBoxes(:,:,iBox)
        CALL CombineGam(EIG,Energy)
        ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
        ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
        CALL PartitionAndEliminate(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR,evalRed,evecRed,LamooDiag)
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
        DEALLOCATE(evalRed,evecRed,LamooDiag)
     ENDDO

     !--- last box: CalcGamLam refreshed fresh every energy (cheap -- BPDLast%Pot/Pcoup were
     !    already filled once above, NOT re-touched here), since only NumOpenR can change.
     !    Always runs (unlike the old NumBoxes.GE.2 gate, which left Boxes(NumBoxes) never
     !    matched at all -- and therefore SD garbage -- when NumBoxes==1). When NumBoxes==1
     !    this box also plays the box-1 role (NumOpenL=0 already set on Boxes(1) above), so it
     !    is matched directly against Bnull instead of a nonexistent Boxes(0). ---
     iBox = NumBoxes
     CALL CalcGamLam(BPDLast,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
     CALL CombineGam(EIG,Energy)
     ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
     ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
     CALL PartitionAndEliminate(BPDLast,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR,evalRed,evecRed,LamooDiag)
     IF (NumBoxes.GE.2) THEN
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPDLast, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
     ELSE
        CALL BoxMatch(Bnull, Boxes(iBox), BPDLast, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
     ENDIF
     DEALLOCATE(evalRed,evecRed,LamooDiag)

     ! Rmat(i,j)=sum_beta Zf(i,beta)Zf(j,beta)/bf(beta) at xEnd, matched to asymptotic
     ! Riccati-Bessel-like reference functions s,c (order=Leff): K=(c+cp*Rmat)*(s+sp*Rmat)^-1.
     CALL CalcK(Boxes(NumBoxes),BPDLast,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)

     IF (NumOpenChannels.EQ.NumChannels .AND. NumChannels.GT.1) THEN
        WRITE(KMatFile,'(A,1P,E20.12)') '# Energy = ',Energy
        DO i = 1,NumOpenChannels
           WRITE(KMatFile,'(1P,20E20.12)') SD%K(i,1:NumOpenChannels)
        ENDDO
     ENDIF

     ALLOCATE(Amat(NumOpenChannels,NumOpenChannels),Smat(NumOpenChannels,NumOpenChannels))
     Amat = (0d0,0d0)
     Smat = (0d0,0d0)
     DO i = 1,NumOpenChannels
        Amat(i,1:NumOpenChannels) = -(0d0,1d0)*SD%K(i,1:NumOpenChannels)
        Smat(i,1:NumOpenChannels) = (0d0,1d0)*SD%K(i,1:NumOpenChannels)
        Amat(i,i) = Amat(i,i) + (1d0,0d0)
        Smat(i,i) = Smat(i,i) + (1d0,0d0)
     ENDDO
     CALL CompSqrMatInv(Amat,NumOpenChannels)
     Smat = MATMUL(Amat,Smat)
     recombProb = 1d0 - ABS(Smat(1,1))**2
     WRITE(RecombFile,'(1P,E20.12,I5,E20.12)') Energy, NumOpenChannels, recombProb
     DEALLOCATE(Amat,Smat)

     IF (NumOpenChannels.EQ.1) THEN
        Kelem = SD%K(1,1)
        rawDelta = DATAN(Kelem)
        IF (havePrevDelta) THEN
           nBranch = NINT((prevDelta-rawDelta)/Pi)
           delta = rawDelta + DBLE(nBranch)*Pi
        ELSE
           delta = rawDelta
        ENDIF
        prevDelta = delta
        havePrevDelta = .TRUE.
        qOurs = DSQRT(2d0*reducedmass*Energy-Threshold(1))
        WRITE(PhaseFile,'(1P,4E20.12)') qOurs, Energy, Kelem, delta
        WRITE(6,'(A,I5,A,1P,E14.6,A,E14.6,A,E14.6)') " ie=",ie," Energy=",Energy," q=",qOurs," delta=",delta
     ELSE
        WRITE(6,'(A,I5,A,1P,E14.6,A,I3,A)') " ie=",ie," Energy=",Energy," NumOpenChannels=",NumOpenChannels, &
             " (skipped -- phase shift only written for exactly 1 open channel)"
     ENDIF

     DO i = 1,NumBoxes
        CALL DeAllocateBox(Boxes(i))
     ENDDO
     DEALLOCATE(Boxes)
     CALL DeAllocateBox(Bnull)
     CALL DeAllocateScat(SD)
  ENDDO

  CALL DeAllocateBPD(BPD0)
  IF (NumBoxes.GE.2) THEN
     CALL DeAllocateBPD(BPD1)
     DEALLOCATE(Gam0Box1,OverlapBox1,LamBox1)
  ENDIF
  CALL DeAllocateBPD(BPDLast)
  DEALLOCATE(xprim)
  CALL DeAllocateEIG(EIG)
  IF (ALLOCATED(Gam0Boxes)) DEALLOCATE(Gam0Boxes,OverlapBoxes,LamBoxes)

  CLOSE(PhaseFile)
  CLOSE(KMatFile)
  CLOSE(RecombFile)
  WRITE(6,*)
  WRITE(6,*) "Phase shift written to PhaseShift.dat"

END PROGRAM main
