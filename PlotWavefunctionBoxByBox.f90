!****************************************************************************************************
! Box-by-box wavefunction reconstruction, single channel, comparing the tabulated/adiabatic
! B-spline box chain against the SVD/direct box chain.
!
! FIXED METHOD (the naive "use Boxes(i)%Zf(1,1) directly as a boundary VALUE" approach in this
! file's earlier version was wrong: BoxMatch always unit-renormalizes Zf/Zfp at every box --
! |Zf|=1 exactly for a single channel -- so it carries no cross-box AMPLITUDE information at all,
! only the log-derivative ratio bf=-Zfp/Zf is physically meaningful there. Using Zf as if it were
! u(R) discarded genuine relative-scale information and produced the wild interior "kinks").
!
! Correct approach: forward sweep, matching VALUE *and* DERIVATIVE at every interface (not two
! edge VALUES). Seed box 1 with an arbitrary overall scale (its own regular solution -- no
! external constraint at the origin). For each subsequent box, BoxMatch's own (pre-renormalization)
! internal quantities
!   Z(i,beta)  =  x^(0.5*RPowExp) * evecRed(i,beta) / Nbeta(beta)     [RPowExp=EffDim-1-2*AlphaFactor]
!   ZP(i,beta) = -evalRed(beta) * Z(i,beta)
! (RMatPropCore.f90's BoxMatch, reproduced here since BoxMatch itself only returns the final
! renormalized Zf/Zfp, discarding this) let us solve, given KNOWN u(xL) and u'(xL) entering the
! box, the linear system  sum_beta Zleft(beta)*c(beta) = u(xL),  sum_beta ZPleft(beta)*c(beta) =
! u'(xL)  for the coefficients c(beta) in the box's own eigenbasis. copen = sum_beta
! evecRed(:,beta)*c(beta) is then the CORRECT open-DOF vector for Xout back-substitution, and
! u(xR)=sum_beta Zright(beta)*c(beta), u'(xR)=sum_beta ZPright(beta)*c(beta) propagate to the next
! box. This is exactly the psi_m -> f_m - g_m*K_mn asymptotic-matching idea applied at every
! interior interface instead of only the outermost one: one value + one derivative, not two values.
!****************************************************************************************************
PROGRAM PlotWavefunctionBoxByBox
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE AdiabaticPotential
  USE AdiabaticInterfaces
  USE SechFunctionPlugins
  USE GenericSVDChannelBasis, ONLY: BuildSVDBoxGeneric
  USE SVDChannelBasis, ONLY: GetGaussRadauFactors
  IMPLICIT NONE

  CHARACTER(LEN=64) DataDir
  TYPE(BPData) BPD0, BPD
  TYPE(GenEigVal) EIG
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(BoxData) Bnull
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:), CoulombC(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION muLocal, alpha, ddim, EnergyLocal, xMinLocal, xEndLocal, xStartLocal
  INTEGER NumDataPoints, xNumPoints, LegPointsLocal, OrderLocal
  INTEGER, PARAMETER :: NBoxesLocal = 10
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  INTEGER iBox, ix, lx, kx

  TYPE BoxRecon
     DOUBLE PRECISION, ALLOCATABLE :: Xout(:,:), cfull(:)
     INTEGER, ALLOCATABLE :: OpenIdxOut(:), ClosedIdxOut(:)
     INTEGER NumOpen, Ncc, MatrixDim
     DOUBLE PRECISION xl, xr
  END TYPE BoxRecon

  ! ==================================================================
  ! Path A: tabulated/adiabatic, B-spline box chain
  ! ==================================================================
  BLOCK
    TYPE(BoxRecon), ALLOCATABLE :: Recon(:)
    DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
    DOUBLE PRECISION, ALLOCATABLE :: copen(:), cclosed(:), cbeta(:)
    DOUBLE PRECISION, ALLOCATABLE :: Nbeta(:), Zleft(:), ZPleft(:), Zright(:), ZPright(:)
    DOUBLE PRECISION uL, uLp, uR, uRp, RPowExp, det
    INTEGER NOL, NOR, beta, k

    DataDir = '3bodydata_sech2_3boson_ch1'
    OrderLocal = 5
    xNumPoints = 40
    LegPointsLocal = 10
    xStartLocal = 0.0d0
    xEndLocal = 40.0d0
    EnergyLocal = 0.5d0

    Pi = dacos(-1d0)
    LegendreFile = 'Legendre.dat'
    LegPoints = LegPointsLocal
    ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
    CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

    CALL ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alpha,ddim, &
         Threshold,Leff,UInterp,PInterp,QInterp)
    reducedmass = muLocal
    AlphaFactor = 0d0
    EffDim = ddim
    RPowExp = EffDim - 1d0 - 2d0*AlphaFactor
    ALLOCATE(CoulombC(NumChannels))
    CoulombC = 0d0

    ALLOCATE(Boxes(NBoxesLocal),Recon(NBoxesLocal))
    Bnull%NumOpenL = 0
    Bnull%NumOpenR = 0
    CALL AllocateBox(Bnull)
    CALL InitZeroBox(Bnull)

    uL = 0d0
    uLp = 0d0
    DO iBox = 1,NBoxesLocal
       BPD%NumChannels = NumChannels
       BPD%Order = OrderLocal
       BPD%Left = 2
       BPD%Right = 2
       BPD%xNumPoints = xNumPoints
       BPD%kl = 1d20
       BPD%kr = 1d20
       CALL AllocateBPD(BPD)
       BPD%xl = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox-1)/DBLE(NBoxesLocal))**BoxSpacingPower
       BPD%xr = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox)/DBLE(NBoxesLocal))**BoxSpacingPower
       CALL GridMaker(BPD%xPoints,BPD%xNumPoints,BPD%xl,BPD%xr,"linear")
       CALL MakeBasis(BPD)
       BPD%lam(1:NumChannels) = Leff(1:NumChannels)
       CALL SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp,CoulombC)

       EIG%MatrixDim = BPD%MatrixDim
       CALL AllocateEIG(EIG)
       Boxes(iBox)%NumOpenL = MERGE(0,1,iBox.EQ.1)
       Boxes(iBox)%NumOpenR = 1
       NOL = Boxes(iBox)%NumOpenL
       NOR = Boxes(iBox)%NumOpenR
       CALL CalcGamLam(BPD,EIG,NOL,NOR)
       CALL CombineGam(EIG,EnergyLocal)
       CALL AllocateBox(Boxes(iBox))
       CALL InitZeroBox(Boxes(iBox))

       Recon(iBox)%NumOpen = Boxes(iBox)%betaMax
       Recon(iBox)%Ncc = BPD%MatrixDim - Recon(iBox)%NumOpen
       Recon(iBox)%MatrixDim = BPD%MatrixDim
       Recon(iBox)%xl = BPD%xl
       Recon(iBox)%xr = BPD%xr
       ALLOCATE(evalRed(Recon(iBox)%NumOpen),evecRed(Recon(iBox)%NumOpen,Recon(iBox)%NumOpen))
       ALLOCATE(LamooDiag(Recon(iBox)%NumOpen))
       ALLOCATE(Recon(iBox)%Xout(Recon(iBox)%Ncc,Recon(iBox)%NumOpen))
       ALLOCATE(Recon(iBox)%OpenIdxOut(Recon(iBox)%NumOpen),Recon(iBox)%ClosedIdxOut(Recon(iBox)%Ncc))
       CALL PartitionAndEliminateFull(BPD,EIG,NOL,NOR, &
            evalRed,evecRed,LamooDiag,Recon(iBox)%Xout,Recon(iBox)%OpenIdxOut,Recon(iBox)%ClosedIdxOut)

       IF (iBox.EQ.1) THEN
          CALL BoxMatch(Bnull, Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
       ELSE
          CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
       ENDIF

       ! Reproduce BoxMatch's own (pre-renormalization) Z/ZP so we can solve for the physically
       ! correct eigenbasis coefficients c(beta) given a KNOWN (value,derivative) at xL, instead
       ! of assuming Zf directly gives a boundary value (see this file's header).
       ALLOCATE(Nbeta(Recon(iBox)%NumOpen))
       ALLOCATE(Zleft(Recon(iBox)%NumOpen),ZPleft(Recon(iBox)%NumOpen))
       ALLOCATE(Zright(Recon(iBox)%NumOpen),ZPright(Recon(iBox)%NumOpen))
       DO beta = 1,Recon(iBox)%NumOpen
          Nbeta(beta) = 0d0
          DO k = 1,Recon(iBox)%NumOpen
             Nbeta(beta) = Nbeta(beta) + evecRed(k,beta)*LamooDiag(k)*evecRed(k,beta)
          ENDDO
          Nbeta(beta) = DSQRT(Nbeta(beta))
          ! NOTE: raw evecRed, NOT BoxMatch's own x^(0.5*RPowExp)/Nbeta-scaled convention --
          ! that scaling exists only to make BoxMatch's OWN Zf/Zfp output unit-normalized; it
          ! is not the basis copen/cfull are expressed in. Using it here (an earlier attempt)
          ! silently reintroduced the same box-to-box scale mismatch this rewrite is fixing.
          IF (NOL.GT.0) THEN
             Zleft(beta) = evecRed(1,beta)
             ! Burke's b is defined w.r.t. the OUTWARD normal, which flips sign between a
             ! box's left/right edges -- BoxMatch's own ZP=-b*Z formula (applied identically
             ! to both edges) is therefore already in outward-normal convention, equal to
             ! +du/dR at the right edge but -du/dR at the left edge. Flip sign here so
             ! ZPleft represents the ordinary coordinate derivative du/dR, matching uLp's
             ! own convention (an ordinary derivative carried over from the previous box).
             ZPleft(beta) = +evalRed(beta)*Zleft(beta)
          ENDIF
          Zright(beta) = evecRed(NOL+1,beta)
          ZPright(beta) = -evalRed(beta)*Zright(beta)
       ENDDO

       ALLOCATE(cbeta(Recon(iBox)%NumOpen))
       IF (iBox.EQ.1) THEN
          ! No external constraint at the origin (Bnull) -- box 1 has a single mode
          ! (betaMax=NumOpenR=1); pick an arbitrary overall scale, c(1)=1.
          cbeta(1) = 1d0
       ELSE
          ! Solve [Zleft(1) Zleft(2); ZPleft(1) ZPleft(2)] * [c1;c2] = [uL;uLp] (2x2, single
          ! channel -- NumOpenL=NumOpenR=1 gives exactly 2 equations for betaMax=2 unknowns).
          det = Zleft(1)*ZPleft(2) - Zleft(2)*ZPleft(1)
          cbeta(1) = (uL*ZPleft(2) - uLp*Zleft(2))/det
          cbeta(2) = (Zleft(1)*uLp - ZPleft(1)*uL)/det
       ENDIF

       ALLOCATE(copen(Recon(iBox)%NumOpen),cclosed(Recon(iBox)%Ncc),Recon(iBox)%cfull(Recon(iBox)%MatrixDim))
       copen = MATMUL(evecRed,cbeta)
       uR = SUM(Zright*cbeta)
       uRp = SUM(ZPright*cbeta)
       cclosed = -MATMUL(Recon(iBox)%Xout,copen)
       DO ix = 1,Recon(iBox)%NumOpen
          Recon(iBox)%cfull(Recon(iBox)%OpenIdxOut(ix)) = copen(ix)
       ENDDO
       DO ix = 1,Recon(iBox)%Ncc
          Recon(iBox)%cfull(Recon(iBox)%ClosedIdxOut(ix)) = cclosed(ix)
       ENDDO
       DEALLOCATE(copen,cclosed,cbeta,Nbeta,Zleft,ZPleft,Zright,ZPright)
       uL = uR
       uLp = uRp

       OPEN(unit=93,file='wavefunction_adiabatic_boxbybox.dat',status='unknown', &
            position=MERGE('REWIND','APPEND',iBox.EQ.1))
       DO kx = 1,BPD%xNumPoints-1
          DO lx = 1,LegPoints
             BLOCK
               DOUBLE PRECISION Rval, uval
               Rval = BPD%x(lx,kx)
               uval = 0d0
               DO ix = 1,BPD%xDim
                  uval = uval + Recon(iBox)%cfull((ix-1)*NumChannels+1)*BPD%u(lx,kx,ix)
               ENDDO
               WRITE(93,*) Rval, uval, iBox
             END BLOCK
          ENDDO
       ENDDO
       CLOSE(93)

       DEALLOCATE(evalRed,evecRed,LamooDiag)
       CALL DeAllocateEIG(EIG)
       CALL DeAllocateBPD(BPD)
    ENDDO
    WRITE(6,*) "=== Path A done: wrote wavefunction_adiabatic_boxbybox.dat ==="
    WRITE(6,*) "Final box Zf(1,1) [unit-normalized] = ",Boxes(NBoxesLocal)%Zf(1,1)
    ! uL here equals copen at the boundary index directly, which for the B-spline "Left=2"
    ! combined boundary basis function IS u(xEnd) with no further conversion -- unlike Path
    ! B below, where the analogous raw coefficient still needs /sqrt(weight). Don't compare
    ! this number directly against Path B's uL print -- see wavefunction_*_boxbybox.dat
    ! (the actual reconstructed files, both correctly converted) for a meaningful comparison.
    WRITE(6,*) "Final box u(xEnd) [physical, this chain's scale] = ",uL

    ! ------------------------------------------------------------------
    ! Absolute physical normalization: psi -> f - g*K at xEnd (user-requested check).
    ! Derived directly from CalcK's own formula (K=(c+cp*Rmat)/(s+sp*Rmat), Rmat=Zf^2/bf)
    ! to find which combination of s,c it corresponds to: Psi=Rmat*Psi' with ansatz Psi=s-c*K
    ! gives K=(s-Rmat*sp)/(c-Rmat*cp) -- does NOT match. Psi=-Rmat*Psi' with Psi=c-s*K gives
    ! K=(c+Rmat*cp)/(s+Rmat*sp) -- MATCHES CalcK exactly. So in the user's own f_m-g_m*K_mn
    ! notation (with legacy's f=this file's c [cosine-like/irregular], g=s [sine-like/regular],
    ! per CalcK's own header comment), u_true(xEnd) = c - s*K, u_true'(xEnd) = cp - sp*K.
    BLOCK
      DOUBLE PRECISION :: kwave, s0,c0,sp0,cp0, rj,ry,rjp,ryp, RmatFinal, Kmat, uTrue, uTrueP, rescale
      kwave = DSQRT(2d0*muLocal*EnergyLocal - 2d0*muLocal*Threshold(1))
      CALL hyperrjry(INT(EffDim),AlphaFactor,Leff(1),kwave*xEndLocal,rj,ry,rjp,ryp)
      s0 = 1d0/DSQRT(Pi*kwave)*rj
      c0 = -1d0/DSQRT(Pi*kwave)*ry
      sp0 = kwave*1d0/DSQRT(Pi*kwave)*rjp
      cp0 = -kwave*1d0/DSQRT(Pi*kwave)*ryp
      RmatFinal = Boxes(NBoxesLocal)%Zf(1,1)**2/Boxes(NBoxesLocal)%bf(1)
      Kmat = (c0+cp0*RmatFinal)/(s0+sp0*RmatFinal)
      uTrue = c0 - s0*Kmat
      uTrueP = cp0 - sp0*Kmat
      rescale = uTrue/uL
      WRITE(6,*) "K (from f-gK matching)          = ",Kmat
      WRITE(6,*) "u_true(xEnd) [asymptotic, f-g*K] = ",uTrue
      WRITE(6,*) "rescale factor (uTrue/uL)        = ",rescale
    END BLOCK
  END BLOCK

  ! ==================================================================
  ! Path B: SVD/direct, single box chain, NumStates=1, GridType='Radau'/'Lobatto'
  ! ==================================================================
  BLOCK
    TYPE(BoxRecon), ALLOCATABLE :: ReconB(:)
    TYPE(BoxData), ALLOCATABLE :: BoxesB(:)
    TYPE(BoxData) BnullB
    TYPE(BPData) BPDSVD
    TYPE(GenEigVal) EIGSVD
    DOUBLE PRECISION, ALLOCATABLE :: evalRedB(:), evecRedB(:,:), LamooDiagB(:)
    DOUBLE PRECISION, ALLOCATABLE :: copenB(:), cclosedB(:), cbetaB(:)
    DOUBLE PRECISION, ALLOCATABLE :: NbetaB(:), ZleftB(:), ZPleftB(:), ZrightB(:), ZPrightB(:)
    DOUBLE PRECISION, ALLOCATABLE :: ODiagL(:), ODiagR(:), nodesR(:), weightsR(:)
    DOUBLE PRECISION :: wsq1Out, wsqLOut, ShiftLocal, uL, uLp, uR, uRp, RPowExp, det
    DOUBLE PRECISION :: lastNodeWeight
    CHARACTER*64 :: LegendreFile64, datadirLoc
    CHARACTER(LEN=16) :: GType
    INTEGER, PARAMETER :: LSVD = 24
    INTEGER :: j, NOL, NOR, beta, k

    LegendreFile64 = 'Legendre.dat'
    datadirLoc = '.'
    ShiftLocal = -5.0d0
    RPowExp = EffDim - 1d0 - 2d0*AlphaFactor
    ALLOCATE(ODiagL(1),ODiagR(1))

    ALLOCATE(BoxesB(NBoxesLocal),ReconB(NBoxesLocal))
    BnullB%NumOpenL = 0
    BnullB%NumOpenR = 0
    CALL AllocateBox(BnullB)
    CALL InitZeroBox(BnullB)

    uL = 0d0
    uLp = 0d0
    DO iBox = 1,NBoxesLocal
       BLOCK
         DOUBLE PRECISION :: a1,a2
         a1 = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox-1)/DBLE(NBoxesLocal))**BoxSpacingPower
         a2 = xStartLocal + (xEndLocal-xStartLocal)*(DBLE(iBox)/DBLE(NBoxesLocal))**BoxSpacingPower
         GType = MERGE('Radau   ','Lobatto ',iBox.EQ.1)

         IF (ALLOCATED(BPDSVD%lam)) DEALLOCATE(BPDSVD%lam)
         CALL BuildSVDBoxGeneric(a1,a2,LSVD,1,5,1,1, &
              200,0d0,dacos(-1d0)/6.0d0,LegendreFile64,10,ShiftLocal,datadirLoc, &
              .FALSE.,SechPotential,SechGridMaker,SechRobinCoeff, &
              MERGE(0,1,iBox.EQ.1),1,TRIM(GType),BPDSVD,EIGSVD,ODiagL,ODiagR,wsq1Out,wsqLOut)

         CALL CombineGam(EIGSVD,EnergyLocal)

         BoxesB(iBox)%NumOpenL = MERGE(0,1,iBox.EQ.1)
         BoxesB(iBox)%NumOpenR = 1
         NOL = BoxesB(iBox)%NumOpenL
         NOR = BoxesB(iBox)%NumOpenR
         CALL AllocateBox(BoxesB(iBox))
         CALL InitZeroBox(BoxesB(iBox))

         ReconB(iBox)%NumOpen = BoxesB(iBox)%betaMax
         ReconB(iBox)%Ncc = BPDSVD%MatrixDim - ReconB(iBox)%NumOpen
         ReconB(iBox)%MatrixDim = BPDSVD%MatrixDim
         ReconB(iBox)%xl = a1
         ReconB(iBox)%xr = a2
         ALLOCATE(evalRedB(ReconB(iBox)%NumOpen),evecRedB(ReconB(iBox)%NumOpen,ReconB(iBox)%NumOpen))
         ALLOCATE(LamooDiagB(ReconB(iBox)%NumOpen))
         ALLOCATE(ReconB(iBox)%Xout(ReconB(iBox)%Ncc,ReconB(iBox)%NumOpen))
         ALLOCATE(ReconB(iBox)%OpenIdxOut(ReconB(iBox)%NumOpen),ReconB(iBox)%ClosedIdxOut(ReconB(iBox)%Ncc))
         CALL PartitionAndEliminateFull(BPDSVD,EIGSVD,NOL,NOR, &
              evalRedB,evecRedB,LamooDiagB,ReconB(iBox)%Xout,ReconB(iBox)%OpenIdxOut,ReconB(iBox)%ClosedIdxOut)

         IF (iBox.EQ.1) THEN
            CALL BoxMatch(BnullB, BoxesB(iBox), BPDSVD, evalRedB, evecRedB, LamooDiagB, EffDim, AlphaFactor)
         ELSE
            CALL BoxMatch(BoxesB(iBox-1), BoxesB(iBox), BPDSVD, evalRedB, evecRedB, LamooDiagB, EffDim, AlphaFactor)
         ENDIF

         ALLOCATE(NbetaB(ReconB(iBox)%NumOpen))
         ALLOCATE(ZleftB(ReconB(iBox)%NumOpen),ZPleftB(ReconB(iBox)%NumOpen))
         ALLOCATE(ZrightB(ReconB(iBox)%NumOpen),ZPrightB(ReconB(iBox)%NumOpen))
         DO beta = 1,ReconB(iBox)%NumOpen
            NbetaB(beta) = 0d0
            DO k = 1,ReconB(iBox)%NumOpen
               NbetaB(beta) = NbetaB(beta) + evecRedB(k,beta)*LamooDiagB(k)*evecRedB(k,beta)
            ENDDO
            NbetaB(beta) = DSQRT(NbetaB(beta))
            ! Raw evecRedB, not BoxMatch's own scaled convention -- see Path A's identical note.
            IF (NOL.GT.0) THEN
               ZleftB(beta) = evecRedB(1,beta)
               ! Same outward-normal sign flip as Path A -- see that block's own comment.
               ZPleftB(beta) = +evalRedB(beta)*ZleftB(beta)
            ENDIF
            ZrightB(beta) = evecRedB(NOL+1,beta)
            ZPrightB(beta) = -evalRedB(beta)*ZrightB(beta)
         ENDDO

         ALLOCATE(cbetaB(ReconB(iBox)%NumOpen))
         IF (iBox.EQ.1) THEN
            cbetaB(1) = 1d0
         ELSE
            det = ZleftB(1)*ZPleftB(2) - ZleftB(2)*ZPleftB(1)
            cbetaB(1) = (uL*ZPleftB(2) - uLp*ZleftB(2))/det
            cbetaB(2) = (ZleftB(1)*uLp - ZPleftB(1)*uL)/det
         ENDIF

         ALLOCATE(copenB(ReconB(iBox)%NumOpen),cclosedB(ReconB(iBox)%Ncc), &
              ReconB(iBox)%cfull(ReconB(iBox)%MatrixDim))
         copenB = MATMUL(evecRedB,cbetaB)
         uR = SUM(ZrightB*cbetaB)
         uRp = SUM(ZPrightB*cbetaB)
         cclosedB = -MATMUL(ReconB(iBox)%Xout,copenB)
         DO j = 1,ReconB(iBox)%NumOpen
            ReconB(iBox)%cfull(ReconB(iBox)%OpenIdxOut(j)) = copenB(j)
         ENDDO
         DO j = 1,ReconB(iBox)%Ncc
            ReconB(iBox)%cfull(ReconB(iBox)%ClosedIdxOut(j)) = cclosedB(j)
         ENDDO
         DEALLOCATE(copenB,cclosedB,cbetaB,NbetaB,ZleftB,ZPleftB,ZrightB,ZPrightB)
         uL = uR
         uLp = uRp

         BLOCK
           DOUBLE PRECISION, ALLOCATABLE :: nodesQ(:), weightsQ(:)
           CHARACTER*64 :: LobattoFile64
           INTEGER :: LQuad, QOffset
           LQuad = LSVD
           QOffset = 0
           ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
           IF (iBox.EQ.1) THEN
              CALL GetGaussRadauFactors(LSVD,nodesQ,weightsQ)
           ELSE
              LobattoFile64 = 'GaussLobatto.dat'
              CALL GetGaussLobattoFactors(LobattoFile64,LSVD,nodesQ,weightsQ)
           ENDIF
           nodesQ = 0.5d0*(a1+a2) + 0.5d0*(a2-a1)*nodesQ
           weightsQ = 0.5d0*(a2-a1)*weightsQ

           OPEN(unit=94,file='wavefunction_svd_boxbybox.dat',status='unknown', &
                position=MERGE('REWIND','APPEND',iBox.EQ.1))
           DO j = 1,LSVD
              WRITE(94,*) nodesQ(QOffset+j), ReconB(iBox)%cfull(j)/DSQRT(weightsQ(QOffset+j)), iBox
           ENDDO
           CLOSE(94)
           IF (iBox.EQ.NBoxesLocal) lastNodeWeight = weightsQ(QOffset+LSVD)
           DEALLOCATE(nodesQ,weightsQ)
         END BLOCK

         DEALLOCATE(evalRedB,evecRedB,LamooDiagB)
         CALL DeAllocateEIG(EIGSVD)
       END BLOCK
    ENDDO
    WRITE(6,*) "=== Path B done: wrote wavefunction_svd_boxbybox.dat ==="
    WRITE(6,*) "Final box Zf(1,1) [unit-normalized] = ",BoxesB(NBoxesLocal)%Zf(1,1)
    ! uL here is the RAW DVR coefficient (copen at the boundary index), NOT yet divided by
    ! sqrt(weight) at that node -- unlike Path A's analogous print, this is NOT u(xEnd)
    ! directly (the actual reconstructed file, wavefunction_svd_boxbybox.dat, correctly
    ! applies that conversion at every node via cfull/sqrt(weight); this raw value does not
    ! and should not be compared against Path A's print).
    WRITE(6,*) "Final box coefficient (raw, NOT yet /sqrt(weight)) = ",uL

    ! Absolute physical normalization, same derivation as Path A's own block (see its
    ! comment for the full f-g*K derivation) -- using the properly DVR-converted
    ! u(xEnd)=uL/sqrt(lastNodeWeight), not the raw coefficient uL itself.
    BLOCK
      DOUBLE PRECISION :: kwave, s0,c0,sp0,cp0, rj,ry,rjp,ryp, RmatFinal, Kmat, uTrue, uChainB, rescale
      kwave = DSQRT(2d0*muLocal*EnergyLocal - 2d0*muLocal*Threshold(1))
      CALL hyperrjry(INT(EffDim),AlphaFactor,Leff(1),kwave*xEndLocal,rj,ry,rjp,ryp)
      s0 = 1d0/DSQRT(Pi*kwave)*rj
      c0 = -1d0/DSQRT(Pi*kwave)*ry
      sp0 = kwave*1d0/DSQRT(Pi*kwave)*rjp
      cp0 = -kwave*1d0/DSQRT(Pi*kwave)*ryp
      RmatFinal = BoxesB(NBoxesLocal)%Zf(1,1)**2/BoxesB(NBoxesLocal)%bf(1)
      Kmat = (c0+cp0*RmatFinal)/(s0+sp0*RmatFinal)
      uTrue = c0 - s0*Kmat
      uChainB = uL/DSQRT(lastNodeWeight)
      rescale = uTrue/uChainB
      WRITE(6,*) "K (from f-gK matching)          = ",Kmat
      WRITE(6,*) "u_chain(xEnd) [DVR-converted]    = ",uChainB
      WRITE(6,*) "u_true(xEnd) [asymptotic, f-g*K] = ",uTrue
      WRITE(6,*) "rescale factor (uTrue/uChainB)   = ",rescale
    END BLOCK
  END BLOCK

END PROGRAM PlotWavefunctionBoxByBox
