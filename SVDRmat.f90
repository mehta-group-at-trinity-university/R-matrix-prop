!****************************************************************************************************
! Single-box R-matrix in Wigner-Eisenbud form, built directly from SVDBound.f90's validated
! T+V construction -- per the user's explicit instruction to stop solving for the log-
! derivative b (PartitionAndEliminate/BoxMatch's eigenchannel approach, RMatPropCore.f90)
! and instead use the proper Wigner-Eisenbud recipe from FEMSVD/DVR-RmatProp.nb and
! "From Jia"/Code4body1D/CalcScatteringSVD.m:
!
!   1. Diagonalize H=T+V ONCE (energy-independent!) over the box's DVR grid, with the
!      LEFT boundary (R=0) Dirichlet-eliminated (regularity, matching SVDBound.f90's own
!      convention) but the RIGHT boundary (R=Rmax) RETAINED as a genuine DOF (not
!      eliminated) -- exactly the notebook's `Diagonalize[T+V, ...]` and the .m file's
!      first-sector `eig(Hmatrix(index_a,index_a))` with index_a = interior + surface DOF.
!      No Bloch-operator/Lam surface-correction term anywhere -- H is bare T+V, exactly
!      as SVDBound.f90 already builds and validates it.
!   2. Surface amplitudes gamma_{lambda,c} = EVecs(boundary-DOF-row-for-channel-c,lambda)/
!      sqrt(w_Rmax) -- the DVR basis's own Kronecker-delta property (matching the .m
!      file's `ub=transpose(uunow(end-nright+1:end,:))*xchb`).
!   3. R_cc'(E) = (1/2mu) * Sum_lambda gamma_{lambda,c}*gamma_{lambda,c'}/(EVals(lambda)-E)
!      -- reusable for EVERY energy without re-diagonalizing (matching the .m file's
!      `Rleft(...)=(ub'*diag(1./(diag(vvnow)-Eplot(kk)))*ub)/(2.d0*mu)`).
!   4. K(E) via the SAME outer hyperrjry Bessel matching CalcK uses (RMatPropCore.f90),
!      with the sqrt(Leff^2+1/4) order fix (see SVDBound.f90's header) -- reimplemented
!      directly here rather than shoehorned through CalcK's BoxData/Zf/bf interface,
!      since R(E) is now a dense matrix from a pole sum, not a Zf/bf pair.
!****************************************************************************************************
PROGRAM SVDRmat
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE SVDChannelBasis, ONLY: SVDCalcDX
  IMPLICIT NONE

  CHARACTER(LEN=64) DataDir
  CHARACTER*64 LobattoFile64
  INTEGER :: LobattoPoints, sR, sRc, NumStates
  INTEGER :: Order1D, Left1D, Right1D, xNumPoints1D, LegPoints1D
  DOUBLE PRECISION :: xMin1D, xMax1D, Shift1D, Rmin, Rmax, muLocal
  DOUBLE PRECISION :: Emin, Emax, EnergyLocal
  INTEGER :: NumEnergies

  DOUBLE PRECISION, ALLOCATABLE :: nodesR(:), weightsR(:)
  DOUBLE PRECISION, ALLOCATABLE :: dXr(:,:), TmatR(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Uad(:,:,:), Psi(:,:,:), S_out(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: O(:,:), TempPsi(:)
  DOUBLE PRECISION, ALLOCATABLE :: H(:,:), EVals(:), EVecs(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: gamma(:,:)
  INTEGER, ALLOCATABLE :: indexOf(:,:)
  DOUBLE PRECISION, EXTERNAL :: ddot

  DOUBLE PRECISION, ALLOCATABLE :: Leff(:), Threshold(:)
  DOUBLE PRECISION, ALLOCATABLE :: Rmat(:,:), s(:), c(:), sp(:), cp(:)
  DOUBLE PRECISION, ALLOCATABLE :: tmp1(:,:), tmp2(:,:), Kmat(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: k(:)
  DOUBLE PRECISION :: rhypj, rhypy, rhypjp, rhypyp, dum

  INTEGER :: PsiDim, HalfBandWidth1D, Ufile
  INTEGER :: i, j, iL, jL, nu, mu2, inu, jmu, probSize, ie
  INTEGER :: KMatFile, NumEigsKeep, NKeep, ios

  OPEN(unit=7,file='SVDRmat.inp',status='old')
  READ(7,*)
  READ(7,*) DataDir
  READ(7,*)
  READ(7,*)
  READ(7,*) LobattoPoints, NumStates, Order1D, Left1D, Right1D
  READ(7,*)
  READ(7,*)
  READ(7,*) xNumPoints1D, xMin1D, xMax1D, LegPoints1D
  READ(7,*)
  READ(7,*)
  READ(7,*) Shift1D
  READ(7,*)
  READ(7,*)
  READ(7,*) Rmax
  READ(7,*)
  READ(7,*)
  READ(7,*) Emin, Emax, NumEnergies
  NumEigsKeep = -1   ! -1 => keep all (default/back-compat if the block below is absent)
  READ(7,*,IOSTAT=ios)
  IF (ios.EQ.0) THEN
     READ(7,*,IOSTAT=ios)
     IF (ios.EQ.0) READ(7,*,IOSTAT=ios) NumEigsKeep
     IF (ios.NE.0) NumEigsKeep = -1
  ENDIF
  CLOSE(7)

  Rmin = 0.d0

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

  ! Leff/Threshold: read directly from FitLeff.data (same 4-field format
  ! ReadAdiabaticData/AdiabaticPotential.f90 uses) -- this tool needs Leff for the outer
  ! Bessel matching but, like SVDBound.f90, never touches the tabulated Uad.dat/Pmat.dat/
  ! Qmat.dat interpolation pipeline (Uad/Psi/S_out all come fresh from OneDimChannels).
  ALLOCATE(Leff(NumStates),Threshold(NumStates))
  OPEN(unit=9,file=TRIM(DataDir)//'/FitLeff.data',status='old')
  READ(9,*)
  READ(9,*)
  READ(9,*)
  DO i = 1,NumStates
     READ(9,*) j, Threshold(i), Leff(i)
  ENDDO
  CLOSE(9)
  Threshold = 2d0*reducedmass*Threshold
  WRITE(6,*) "DataDir=",TRIM(DataDir)," mu=",muLocal," NumStates=",NumStates," Leff=",Leff

  sR = LobattoPoints
  sRc = sR - 1   ! only R=0 (node 1) is Dirichlet-eliminated; R=Rmax (node sR) is RETAINED
  ALLOCATE(nodesR(sR),weightsR(sR))
  LobattoFile64 = 'GaussLobatto.dat'
  CALL GetGaussLobattoFactors(LobattoFile64,sR,nodesR,weightsR)
  nodesR = 0.5d0*(Rmin+Rmax) + 0.5d0*(Rmax-Rmin)*nodesR
  weightsR = 0.5d0*(Rmax-Rmin)*weightsR
  WRITE(6,*) "Radial grid: Rmin=",Rmin," Rmax=",Rmax," sR=",sR," sRc(retained)=",sRc

  !--- DVR derivative matrix, exactly SVDBound.f90's construction (see its own comments) ---
  ALLOCATE(dXr(sR,sRc))
  DO iL = 1,sR
     DO jL = 1,sRc
        CALL SVDCalcDX(nodesR,weightsR,jL+1,iL,sR,dXr(iL,jL))
     ENDDO
  ENDDO

  ALLOCATE(TmatR(sRc,sRc))
  DO iL = 1,sRc
     DO jL = 1,sRc
        TmatR(iL,jL) = (0.5d0/reducedmass) * SUM(weightsR(:)*dXr(:,iL)*dXr(:,jL))
     ENDDO
  ENDDO

  !--- per-R adiabatic channel data at the sRc retained nodes (interior + Rmax boundary) ---
  PsiDim = xNumPoints1D + Order1D - 3
  IF (Left1D.EQ.2)  PsiDim = PsiDim + 1
  IF (Right1D.EQ.2) PsiDim = PsiDim + 1
  HalfBandWidth1D = Order1D
  ALLOCATE(Uad(sRc,NumStates,2),Psi(sRc,PsiDim,NumStates),S_out(HalfBandWidth1D+1,PsiDim))
  Ufile = 71
  OPEN(unit=Ufile,file=TRIM(DataDir)//'/SVDRmatUad.dat',status='replace')
  CALL OneDimChannels(NumStates,PsiDim,0,0,LegendreFile,LegPoints1D,Shift1D, &
       Order1D,Left1D,Right1D,reducedmass,xNumPoints1D,xMin1D,xMax1D, &
       nodesR(2:sR),sRc,Ufile,72,73,74,75,Uad,Psi,S_out,DataDir)
  CLOSE(Ufile)
  WRITE(6,*) "OneDimChannels done."

  probSize = sRc*NumStates
  ALLOCATE(indexOf(sRc,NumStates))
  i = 1
  DO iL = 1,sRc
     DO nu = 1,NumStates
        indexOf(iL,nu) = i
        i = i+1
     ENDDO
  ENDDO

  ALLOCATE(O(probSize,probSize),TempPsi(PsiDim))
  DO iL = 1,sRc
     DO nu = 1,NumStates
        CALL dsbmv('U',PsiDim,HalfBandWidth1D,1.0d0,S_out,HalfBandWidth1D+1, &
             Psi(iL,:,nu),1,0.0d0,TempPsi,1)
        inu = indexOf(iL,nu)
        DO jL = 1,sRc
           DO mu2 = 1,NumStates
              jmu = indexOf(jL,mu2)
              O(inu,jmu) = ddot(PsiDim,TempPsi,1,Psi(jL,:,mu2),1)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE(TempPsi)
  WRITE(6,*) "Channel-overlap matrix O done."

  ALLOCATE(H(probSize,probSize))
  H = 0d0
  DO iL = 1,sRc
     DO nu = 1,NumStates
        inu = indexOf(iL,nu)
        DO jL = 1,sRc
           DO mu2 = 1,NumStates
              jmu = indexOf(jL,mu2)
              H(inu,jmu) = TmatR(iL,jL)*O(inu,jmu)
           ENDDO
        ENDDO
        H(inu,inu) = H(inu,inu) + Uad(iL,nu,1)
     ENDDO
  ENDDO
  WRITE(6,*) "Hamiltonian assembled, probSize=",probSize

  ALLOCATE(EVals(probSize),EVecs(probSize,probSize))
  WRITE(6,*) "Diagonalizing (once, energy-independent)..."
  CALL Mydsyev(H,probSize,EVals,EVecs)
  WRITE(6,*) "Done."

  !--- surface amplitudes at the Rmax boundary (retained index sRc = full node sR) ---
  ALLOCATE(gamma(NumStates,probSize))
  DO nu = 1,NumStates
     gamma(nu,:) = EVecs(indexOf(sRc,nu),:)/DSQRT(weightsR(sR))
  ENDDO

  OPEN(unit=15,file='SVDRmat_Evals.dat',status='replace')
  DO i = 1,probSize
     WRITE(15,'(I8,1P,100E20.12)') i, EVals(i), (gamma(nu,i), nu=1,NumStates)
  ENDDO
  CLOSE(15)

  ! Only the lowest fraction of a DVR spectrum is physically converged -- the highest
  ! eigenstates are quadrature/aliasing artifacts of the finite grid, not real solutions
  ! (confirmed directly: EVals(1:~280) agree to 1e-9+ between LobattoPoints=100 and 200,
  ! then diverge sharply, reaching ~90% relative error by EVals(~450) of 495 total).
  ! Including the spurious high states in the Wigner-Eisenbud pole sum corrupts R(E) for
  ! every energy, since Sum_lambda runs over the WHOLE spectrum. NumEigsKeep truncates it;
  ! -1 (or omitted from SVDRmat.inp) keeps every state (the previous, uncorrected behavior).
  IF (NumEigsKeep.LT.0 .OR. NumEigsKeep.GT.probSize) THEN
     NKeep = probSize
  ELSE
     NKeep = NumEigsKeep
  ENDIF
  WRITE(6,*) "NKeep=",NKeep," out of probSize=",probSize

  !--- energy scan: build R(E) via the Wigner-Eisenbud pole sum, then K(E) via the SAME
  !    outer hyperrjry Bessel matching CalcK uses (RMatPropCore.f90), order=sqrt(Leff^2+1/4)
  !    per SVDBound.f90's validation ---
  ALLOCATE(Rmat(NumStates,NumStates),k(NumStates))
  ALLOCATE(s(NumStates),c(NumStates),sp(NumStates),cp(NumStates))
  ALLOCATE(tmp1(NumStates,NumStates),tmp2(NumStates,NumStates),Kmat(NumStates,NumStates))

  KMatFile = 50
  OPEN(unit=KMatFile,file='SVDRmat_KMatrix.dat',status='replace')
  DO ie = 1,NumEnergies
     IF (NumEnergies.EQ.1) THEN
        EnergyLocal = Emin
     ELSE
        EnergyLocal = Emin + (Emax-Emin)*DBLE(ie-1)/DBLE(NumEnergies-1)
     ENDIF

     Rmat = 0d0
     DO i = 1,NumStates
        DO j = 1,NumStates
           Rmat(i,j) = (1d0/(2d0*reducedmass)) * &
                SUM(gamma(i,1:NKeep)*gamma(j,1:NKeep)/(EVals(1:NKeep)-EnergyLocal))
        ENDDO
     ENDDO

     DO i = 1,NumStates
        k(i) = DSQRT(2d0*reducedmass*EnergyLocal-Threshold(i))
        CALL hyperrjry(2,0.5d0,DSQRT(Leff(i)**2+0.25d0),k(i)*Rmax,rhypj,rhypy,rhypjp,rhypyp)
        s(i)  =  1d0/DSQRT(Pi*k(i))*rhypj
        c(i)  = -1d0/DSQRT(Pi*k(i))*rhypy
        sp(i) =  k(i)*1d0/DSQRT(Pi*k(i))*rhypjp
        cp(i) = -k(i)*1d0/DSQRT(Pi*k(i))*rhypyp
     ENDDO

     tmp1 = 0d0
     tmp2 = 0d0
     DO i = 1,NumStates
        DO j = 1,NumStates
           tmp1(i,j) = cp(i)*Rmat(i,j)
           tmp2(i,j) = sp(i)*Rmat(i,j)
        ENDDO
        tmp1(i,i) = tmp1(i,i) + c(i)
        tmp2(i,i) = tmp2(i,i) + s(i)
     ENDDO
     CALL SqrMatInv(tmp2,NumStates)
     Kmat = MATMUL(tmp1,tmp2)

     WRITE(KMatFile,'(A,1P,E20.12)') '# Energy = ',EnergyLocal
     DO i = 1,NumStates
        WRITE(KMatFile,'(1P,100E20.12)') (Kmat(i,j), j=1,NumStates)
     ENDDO
  ENDDO
  CLOSE(KMatFile)
  WRITE(6,*) "Wrote SVDRmat_KMatrix.dat"

END PROGRAM SVDRmat
