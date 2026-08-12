!****************************************************************************************************
! Third variant of SVDChannelBasis.f90's BuildSVDBox, specific to the equal-mass three-boson
! delta-function problem: gets its per-R channel data (Uad, cross-R channel overlap) from a
! GENUINE numerical diagonalization -- unlike AnalyticSVDChannelBasis.f90's BuildSVDBoxAnalytic,
! which bypasses diagonalization entirely via KMSFormulas.f's closed-form SolveQ/EvalOverlapAt --
! but WITHOUT approximating the delta function by a smooth potential well the way
! SVDChannelBasis.f90's own BuildSVDBox does via OneDimChannels' Poschl-Teller well (which can
! never exactly reproduce a genuine discontinuous-derivative interaction, however narrow).
!
! Instead: between phi=0 and phi=pi/6 the delta-function problem has NO potential at all -- its
! entire effect is the jump condition at phi=pi/6, which reduces (for the wedge-restricted,
! even-parity solution) to a boundary log-derivative Phi'(pi/6)/Phi(pi/6) = sqrt(2*mu)*R, exactly
! (verified numerically to machine precision against KMSFormulas.f's SolveQ, channel-independent
! as expected since the jump condition is a property of the interaction, not the channel -- see
! CheckLogDeriv.f90 in this session's scratch tests). This subroutine builds a B-spline basis over
! phi in [0,pi/6] via the shared library's CalcBasisFuncsBP (Left=1 Neumann at phi=0, matching the
! even-parity Phi'(0)=0 condition; Right=3 Robin/mixed at phi=pi/6, kr=sqrt(2*mu)*R, exact), then
! diagonalizes the FREE (V=0) angular kinetic operator -- WITH the surface/Bloch term the Robin
! condition requires (T_ij -= kr*u_i(pi/6)*u_j(pi/6); the boundary-combined basis function is
! exactly 1 at its own edge by the standard clamped-B-spline normalization, so this is a single
! scalar correction to T's last diagonal entry -- omitting it silently gives the WRONG, pure
! Neumann-Neumann free-particle spectrum instead; verified this exact failure mode directly before
! finding the fix). Confirmed to reproduce KMSFormulas.f's exact spectrum (eigenvalues -q1^2,
! kappa_2^2,...,kappa_5^2 simultaneously, from ONE diagonalization) to machine precision for the
! lowest few channels, degrading to ~1e-8 for the highest -- ordinary B-spline basis-set
! truncation, not a sign of anything wrong.
!
! Since the angular DOMAIN (and B-spline knot placement) never moves with R -- only kr does -- the
! usual "shared angular basis, R-varying coefficients" projection trick (BuildSVDBox's own S_out
! machinery) does not directly apply here (the BASIS FUNCTIONS themselves, not just the
! coefficients, depend on R through kr). Instead, each R's diagonalized eigenfunction is evaluated
! POINTWISE at the shared angular Gauss-Legendre quadrature grid (valid regardless of which
! R-dependent basis produced it), and the cross-R channel overlap is a plain quadrature-weighted
! dot product of those pointwise values -- no basis-projection matrix needed at all.
!****************************************************************************************************
MODULE RobinSVDChannelBasis
  IMPLICIT NONE
CONTAINS
  !****************************************************************************************************
  ! Signature mirrors AnalyticSVDChannelBasis.f90's BuildSVDBoxAnalytic, with two extra inputs
  ! (OrderPhi, xNumPointsPhi) controlling the ANGULAR B-spline basis resolution -- independent of
  ! L/GridType, which describe the RADIAL DVR grid as usual.
  SUBROUTINE BuildSVDBoxRobin(a1,a2,L,NumStates,NumOpenL,NumOpenR,GridType,OrderPhi,xNumPointsPhi, &
       BPD,EIG,ODiagL,ODiagR,wsq1Out,wsqLOut)
    USE DataStructures
    USE GlobalVars
    USE Quadrature
    USE SVDChannelBasis, ONLY: GetGaussRadauFactors, SVDCalcDX
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: a1,a2
    INTEGER, INTENT(in) :: L,NumStates,NumOpenL,NumOpenR,OrderPhi,xNumPointsPhi
    CHARACTER(LEN=*), INTENT(in) :: GridType
    TYPE(BPData) :: BPD
    TYPE(GenEigVal) :: EIG
    DOUBLE PRECISION, INTENT(out) :: ODiagL(NumStates), ODiagR(NumStates)
    DOUBLE PRECISION, INTENT(out) :: wsq1Out, wsqLOut

    TYPE(BPData) :: BPDPhi
    DOUBLE PRECISION, ALLOCATABLE :: nodesR(:), weightsR(:)
    DOUBLE PRECISION, ALLOCATABLE :: nodesQ(:), weightsQ(:), wpowQ(:)
    DOUBLE PRECISION, ALLOCATABLE :: Uad(:,:), Psi(:,:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: O(:,:), dXr(:,:), TmatBulk(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Tphi(:,:), Sphi(:,:), evalPhi(:), evecPhi(:,:), work(:)
    INTEGER :: iL,jL,nu,mu2,inu,jmu2,LQuad,QOffset,lwork,info,kxp,lxp
    DOUBLE PRECISION :: RPowExp, RpowIL, wsq1, wsqL, Rval, kRight, wPhi
    CHARACTER*64 :: LobattoFile64

    RPowExp = EffDim - 1d0 - 2d0*AlphaFactor

    !--- geometry: identical to BuildSVDBox/BuildSVDBoxAnalytic ---
    BPD%xl = a1
    BPD%xr = a2
    BPD%NumChannels = NumStates
    BPD%xDim = L
    BPD%MatrixDim = L*NumStates
    BPD%Order = L-1
    ALLOCATE(BPD%lam(NumStates))
    BPD%lam = 0d0

    EIG%MatrixDim = BPD%MatrixDim
    CALL AllocateEIG(EIG)

    !--- radial DVR grid over [a1,a2]: identical logic to BuildSVDBox/BuildSVDBoxAnalytic ---
    ALLOCATE(nodesR(L),weightsR(L))
    IF (GridType.EQ.'Radau') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBoxRobin: GridType='Radau' requires NumOpenL=0 (no left-boundary DOF exists)"
          STOP
       ENDIF
       LQuad = L
       QOffset = 0
       ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
       CALL GetGaussRadauFactors(L,nodesQ,weightsQ)
    ELSEIF (GridType.EQ.'LobattoDirichlet') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBoxRobin: GridType='LobattoDirichlet' requires NumOpenL=0 (a1 DOF is eliminated)"
          STOP
       ENDIF
       LQuad = L+1
       QOffset = 1
       ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
       LobattoFile64 = 'GaussLobatto.dat'
       CALL GetGaussLobattoFactors(LobattoFile64,LQuad,nodesQ,weightsQ)
    ELSE
       LQuad = L
       QOffset = 0
       ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
       LobattoFile64 = 'GaussLobatto.dat'
       CALL GetGaussLobattoFactors(LobattoFile64,L,nodesQ,weightsQ)
    ENDIF
    nodesQ = 0.5d0*(a1+a2) + 0.5d0*(a2-a1)*nodesQ
    weightsQ = 0.5d0*(a2-a1)*weightsQ
    nodesR = nodesQ(QOffset+1:QOffset+L)
    weightsR = weightsQ(QOffset+1:QOffset+L)
    ALLOCATE(wpowQ(LQuad))
    wpowQ = weightsQ*nodesQ**RPowExp

    !--- angular B-spline basis over phi in [0,pi/6]: geometry is FIXED (only kr varies with R
    !    below), so build it once here (AllocateBPD + xPoints), then just update BPDPhi%kr and
    !    call MakeBasis fresh inside the R-loop -- cheap, no reallocation needed each time. ---
    BPDPhi%NumChannels = 1
    BPDPhi%Order = OrderPhi
    BPDPhi%Left = 1
    BPDPhi%Right = 3
    BPDPhi%xNumPoints = xNumPointsPhi
    BPDPhi%kl = 1d20
    CALL AllocateBPD(BPDPhi)
    BLOCK
      INTEGER :: iphi
      DO iphi = 1,xNumPointsPhi
         BPDPhi%xPoints(iphi) = (Pi/6.0d0)*DBLE(iphi-1)/DBLE(xNumPointsPhi-1)
      ENDDO
    END BLOCK

    ALLOCATE(Tphi(BPDPhi%xDim,BPDPhi%xDim),Sphi(BPDPhi%xDim,BPDPhi%xDim))
    ALLOCATE(evalPhi(BPDPhi%xDim),evecPhi(BPDPhi%xDim,BPDPhi%xDim))
    lwork = MAX(1,3*BPDPhi%xDim-1)
    ALLOCATE(work(lwork))

    ALLOCATE(Uad(L,NumStates),Psi(L,LegPoints,xNumPointsPhi-1,NumStates))

    DO iL = 1,L
       Rval = nodesR(iL)
       kRight = DSQRT(2d0*reducedmass)*Rval
       BPDPhi%kr = kRight
       CALL MakeBasis(BPDPhi)

       !--- bulk kinetic + overlap matrices, plus the Robin surface/Bloch term (see this
       !    file's own header for the derivation) -- WITHOUT it the Right boundary condition
       !    is silently dropped and the spectrum comes out as the wrong, pure Neumann-Neumann
       !    free-particle one instead. ---
       Tphi = 0d0
       Sphi = 0d0
       DO jL = 1,BPDPhi%xDim
          DO mu2 = 1,BPDPhi%xDim
             DO kxp = 1,BPDPhi%xNumPoints-1
                DO lxp = 1,LegPoints
                   wPhi = 0.5d0*(BPDPhi%xPoints(kxp+1)-BPDPhi%xPoints(kxp))*wLeg(lxp)
                   Tphi(jL,mu2) = Tphi(jL,mu2) + wPhi*BPDPhi%ux(lxp,kxp,jL)*BPDPhi%ux(lxp,kxp,mu2)
                   Sphi(jL,mu2) = Sphi(jL,mu2) + wPhi*BPDPhi%u(lxp,kxp,jL)*BPDPhi%u(lxp,kxp,mu2)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       Tphi(BPDPhi%xDim,BPDPhi%xDim) = Tphi(BPDPhi%xDim,BPDPhi%xDim) - kRight

       CALL DSYGV(1,'V','U',BPDPhi%xDim,Tphi,BPDPhi%xDim,Sphi,BPDPhi%xDim,evalPhi,work,lwork,info)
       IF (info.NE.0) THEN
          WRITE(6,*) "BuildSVDBoxRobin: DSYGV failed, info=",info," at R=",Rval
          STOP
       ENDIF
       evecPhi = Tphi   ! DSYGV overwrites the first matrix argument with the eigenvectors

       ! DSYGV returns each eigenvector up to an arbitrary overall sign, chosen
       ! independently at every R -- left uncorrected, the cross-R overlap between
       ! adjacent DVR nodes can pick up a spurious sign flip wherever LAPACK happens to
       ! choose opposite signs at neighboring R's (confirmed directly: this showed up as
       ! Gam0 entries matching BuildSVDBoxAnalytic in MAGNITUDE but with a flipped sign,
       ! not a small basis-set-truncation-level difference).
       !
       ! First fix attempt: enforce Phi(0)>0 always, using Phi(0)=evecPhi(1,nu) directly
       ! (the Left=1 Neumann basis's first function is exactly 1 at phi=0 by clamped-knot
       ! normalization). This is safe for the continuum (cos) channels nu>1, since
       ! Phi(0)=B*cos(0)=B is never suppressed there. But channel 1 (nu=1, the bound-pair
       ! cosh branch) broke it: Phi(0)=B*cosh(0)=B while Phi(pi/6)=B*cosh(q1*pi/6) grows
       ! exponentially with q1(R) -- so at large R (large q1, e.g. R>~65 in the 20-box/
       ! xEnd=100 delta-function run), Phi(0) is exponentially SUPPRESSED relative to the
       ! eigenvector's overall O(1) normalization and decays into roundoff noise
       ! (confirmed directly: evecPhi(1,1) ~3e-5 vs ~2 for other channels at R~65-70,
       ! sign essentially random node-to-node) -- causing exactly this spurious sign-flip
       ! failure mode, just relocated to channel 1 specifically and only at large R (this
       ! is what caused RMATPROPHybridSVDRobin.x's box-14-onward blowup: max|Onum-Oana|=2,
       ! i.e. a full sign flip, isolated via DiagBox14.f90 in this session's scratch tests).
       !
       ! Correct fix: use whichever boundary point is NEVER small for that channel's
       ! branch. cos(kappa*phi) channels (nu>1): Phi(0)=B is robust (Phi(pi/6) is NOT --
       ! it genuinely crosses zero whenever kappa(R)*pi/6 sweeps through pi/2+n*pi, which
       ! is what broke the very first fix attempt, pinning ALL channels via Phi(pi/6)).
       ! cosh(q*phi) channel (nu=1): Phi(pi/6)=B*cosh(q*pi/6) is robust instead (cosh is
       ! monotonically growing and never crosses zero, and is the numerically DOMINANT
       ! endpoint value, never suppressed into roundoff).
       DO nu = 1,NumStates
          IF (nu.EQ.1) THEN
             IF (evecPhi(BPDPhi%xDim,nu).LT.0d0) evecPhi(:,nu) = -evecPhi(:,nu)
          ELSE
             IF (evecPhi(1,nu).LT.0d0) evecPhi(:,nu) = -evecPhi(:,nu)
          ENDIF
       ENDDO

       DO nu = 1,NumStates
          Uad(iL,nu) = evalPhi(nu)/(2d0*reducedmass*Rval**2)
          DO kxp = 1,xNumPointsPhi-1
             DO lxp = 1,LegPoints
                Psi(iL,lxp,kxp,nu) = SUM(evecPhi(1:BPDPhi%xDim,nu)*BPDPhi%u(lxp,kxp,1:BPDPhi%xDim))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !--- channel-overlap matrix O(l,nu;l',nu') for ALL DVR-point pairs: plain quadrature-
    !    weighted dot product of the pointwise-evaluated eigenfunctions (valid regardless of
    !    which R-dependent basis produced them -- see this file's own header). ---
    ALLOCATE(O(L*NumStates,L*NumStates))
    O = 0d0
    DO iL = 1,L
       DO nu = 1,NumStates
          inu = (iL-1)*NumStates+nu
          DO jL = 1,L
             DO mu2 = 1,NumStates
                jmu2 = (jL-1)*NumStates+mu2
                DO kxp = 1,xNumPointsPhi-1
                   DO lxp = 1,LegPoints
                      wPhi = 0.5d0*(BPDPhi%xPoints(kxp+1)-BPDPhi%xPoints(kxp))*wLeg(lxp)
                      O(inu,jmu2) = O(inu,jmu2) + wPhi*Psi(iL,lxp,kxp,nu)*Psi(jL,lxp,kxp,mu2)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !--- DVR derivative matrix / bulk kinetic matrix: IDENTICAL to BuildSVDBox/
    !    BuildSVDBoxAnalytic -- properties of the radial DVR grid alone. ---
    ALLOCATE(dXr(LQuad,L))
    DO iL = 1,LQuad
       DO jL = 1,L
          CALL SVDCalcDX(nodesQ,weightsQ,jL+QOffset,iL,LQuad,dXr(iL,jL))
       ENDDO
    ENDDO
    ALLOCATE(TmatBulk(L,L))
    DO iL = 1,L
       DO jL = 1,L
          TmatBulk(iL,jL) = SUM(wpowQ(:)*dXr(:,iL)*dXr(:,jL))
       ENDDO
    ENDDO

    !--- assemble EIG%Gam0/Overlap/Lam -- IDENTICAL formulas to BuildSVDBox/BuildSVDBoxAnalytic,
    !    fed the numerically-diagonalized O/Uad. ---
    EIG%Gam0 = 0d0
    EIG%Overlap = 0d0
    EIG%Lam = 0d0
    DO iL = 1,L
       DO nu = 1,NumStates
          inu = (iL-1)*NumStates+nu
          DO jL = 1,L
             DO mu2 = 1,NumStates
                jmu2 = (jL-1)*NumStates+mu2
                EIG%Gam0(inu,jmu2) = EIG%Gam0(inu,jmu2) - TmatBulk(iL,jL)*O(inu,jmu2)
             ENDDO
          ENDDO
          RpowIL = nodesR(iL)**RPowExp
          EIG%Overlap(inu,inu) = 2d0*reducedmass*RpowIL*O(inu,inu)
          EIG%Gam0(inu,inu) = EIG%Gam0(inu,inu) - 2d0*reducedmass*RpowIL*Uad(iL,nu)
       ENDDO
    ENDDO

    wsq1 = 1d0/weightsR(1)
    DO nu = 1,NumOpenL
       inu = nu
       EIG%Lam(inu,inu) = wsq1*(a1**RPowExp)*O(inu,inu)
    ENDDO
    wsqL = 1d0/weightsR(L)
    DO nu = 1,NumOpenR
       inu = (L-1)*NumStates+nu
       EIG%Lam(inu,inu) = wsqL*(a2**RPowExp)*O(inu,inu)
    ENDDO

    DO nu = 1,NumStates
       ODiagL(nu) = O(nu,nu)
       ODiagR(nu) = O((L-1)*NumStates+nu,(L-1)*NumStates+nu)
    ENDDO
    wsq1Out = wsq1
    wsqLOut = wsqL

    CALL DeAllocateBPD(BPDPhi)
    DEALLOCATE(nodesR,weightsR,nodesQ,weightsQ,wpowQ,Uad,Psi,O,dXr,TmatBulk)
    DEALLOCATE(Tphi,Sphi,evalPhi,evecPhi,work)
  END SUBROUTINE BuildSVDBoxRobin
END MODULE RobinSVDChannelBasis
