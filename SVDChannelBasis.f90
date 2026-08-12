!****************************************************************************************************
! Builds one SVD (slow-variable-discretization) R-matrix box: a hyperradial sector
! [a1,a2] covered by a single Gauss-Lobatto DVR grid, combined with per-R adiabatic
! channel functions from Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90
! (OneDimChannels, CouplingFlag=0 -- SVD needs no d/dR coupling, that's the whole
! point of the method). Plays the same role MakeBasis+SetAdiabaticPotential+CalcGamLam
! play for a B-spline box (RMatPropCore.f90): on return, BPD's geometry/Order and
! EIG%Gam0/Overlap/Lam are ready for CombineGam/PartitionAndEliminate/BoxMatch,
! unchanged. See Suno, PRA 109, 042814 (2024), Eqs. 12-20 for the underlying physics;
! single-box template adapted from ~/Documents/GitHub/4BodySVD/4BodySVD.f90 (GetGaussLobattoFactors
! grid, CalcDX derivative matrix, CalcOMatrix channel overlap).
!
! Design notes (see plan for full derivation):
!  - PartitionAndEliminate only ever reads the OPEN-block DIAGONAL of EIG%Lam
!    (RMatPropCore.f90:867: Lamoo(ic,ic)=EIG%Lam(OpenIdx(ic),OpenIdx(ic))) -- so Lam
!    must be diagonal, exactly mirroring CalcGamLam's own boundary convention
!    (Lam(boundary dof,itself) = [basis fn value at its own boundary node]^2 * R^power),
!    NOT Suno's literal surface-term formula (which would populate full boundary rows).
!    For the normalized DVR basis pi_l(R)=L_l(R)/sqrt(w_l), pi_l(R_l)=1/sqrt(w_l), so
!    the l=1/l=L boundary Lam entries pick up a 1/w factor that has no B-spline analog
!    (B-splines are normalized to exactly 1 at their own boundary point instead).
!  - CalcGamLam bakes the hyperradial measure R**(EffDim-1-2*AlphaFactor) into its
!    Legendre-quadrature weight and reuses it for both the kinetic and potential/overlap
!    terms; matched here via wpow(k)=weightsR(k)*nodesR(k)**RPowExp for the kinetic
!    term, and nodesR(l)**RPowExp directly for the (DVR-diagonal-exact) overlap/
!    potential terms -- the plain quadrature weight cancels there because the DVR
!    orthogonality property collapses the l/=l' angular-overlap "plain" integral to a
!    Kronecker delta at exactly weight R_l**RPowExp (the w_l's cancel against the
!    1/sqrt(w_l) basis normalization). Only the kinetic term stays dense in (l,l').
!
! ANGULAR-GRID CONSISTENCY (found via a free-particle V=0 diagnostic, where the exact
! answer -- K=0, zero cross-channel coupling -- is known in closed form): this box's
! channel-overlap matrix O(l,nu;l',nu') is computed via a SINGLE angular banded overlap
! matrix S_out, captured once from OneDimChannels, then reused via dsbmv for every
! (l,l') pair -- correct only if every R point in the box's DVR grid shares the SAME
! angular basis/knot placement. adiabaticSolver1D.f90's stock OneDimChannels calls
! GridMakerHHL fresh at every R, which for R > xRswitch = 10*PotRange/Pi switches to a
! grid whose knot placement itself depends continuously on R (deltax=PotRange/R) --
! silently mismatching the basis actually used at earlier R points in that regime.
! Confirmed empirically: with the stock solver, an SVD box straddling or lying entirely
! above xRswitch showed self-overlap deviating from 1 by up to 20% and spurious cross-
! channel O values up to 0.46 (should be exactly 0 for non-interacting channels).
!
! Fixed at the source rather than by constraining callers, and WITHOUT giving up
! GridMakerHHL's coalescence-adjacent refinement (needed for accuracy with a real,
! short-ranged two-body potential): R-matrix-prop's local copy of OneDimChannels
! (adiabaticSolver1D.f90, used only by BuildSVDBox) now calls GridMakerHHL exactly
! ONCE per call -- not once per R the way stock OneDimChannels does -- anchored at
! R(RSteps), the box's own largest/rightmost R (BuildSVDBox always passes its L DVR
! nodes in ascending order, so this is the box's right edge a2). GridMakerHHL's
! refined-region width narrows monotonically with R (deltax=PotRange/R), so anchoring
! on the box's largest R gives a grid at least as finely resolved as any smaller R in
! the same box would independently request. Every R point in a box therefore shares
! one grid -- whichever GridMakerHHL would have chosen (uniform or refined) for that
! box's own R(RSteps) -- with no dependence on PotRange/xRswitch beyond ordinary grid-
! quality tradeoffs, and no constraint on where a box's [a1,a2] can sit.
! (The upstream Adiabatic-Scattering-BoundStates copy is untouched -- its own callers,
! e.g. dataset generation via DriveAdiabaticSolver1D.x, evaluate each R independently
! and never need cross-R-point consistency, so the fully R-dependent grid is fine to
! leave as-is there.) Verified via the free-particle diagnostic with boxes spanning
! both sides of xRswitch (so some anchor in the uniform branch, others in the refined
! branch): self-overlap/cross-channel-O are machine-precision correct in every box,
! and the resulting K-matrix matches a plain-uniform-grid-everywhere version to
! floating-point noise (~1e-13 relative) -- confirming the physics is correctly
! independent of which angular basis GridMakerHHL happens to choose.
!****************************************************************************************************
MODULE SVDChannelBasis
  IMPLICIT NONE
CONTAINS
  !****************************************************************************************************
  ! Gauss-Radau nodes/weights (Legendre weight, w(x)=1) on [-1,1], Points total nodes,
  ! ONE FIXED at x0=+1 (right endpoint), the remaining Points-1 strictly interior --
  ! used for the leftmost SVD box so no DVR node sits exactly at the true R=0 origin
  ! (GridMakerHHL inside OneDimChannels has a r0/R term that diverges there).
  !
  ! Golub (1973)-style modified-Jacobi-matrix eigenvalue method, derived here from the
  ! standard orthonormal Legendre 3-term recurrence (alpha_k=0, beta_k=k/sqrt(4k^2-1)):
  ! replace ONLY the last diagonal entry of the standard N-point Jacobi matrix with a
  ! value alpha_tilde chosen so that x0 becomes an exact eigenvalue (root of the
  ! degree-N orthogonal polynomial). Solving beta_N*q_N(x0)=(x0-alpha_tilde)*q_{N-1}(x0)
  ! - beta_{N-1}*q_{N-2}(x0) = 0 for alpha_tilde gives alpha_tilde = x0 -
  ! beta_{N-1}*q_{N-2}(x0)/q_{N-1}(x0), where q_{N-1}(x0),q_{N-2}(x0) come from the
  ! UNMODIFIED recurrence (they don't depend on the entry being modified). Verified
  ! independently in Python against exact-degree-(2N-2) polynomial integration (the
  ! defining property of N-point Gauss-Radau) to machine precision for N=3..24 before
  ! porting here.
  SUBROUTINE GetGaussRadauFactors(Points,x,w)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: Points
    DOUBLE PRECISION, INTENT(out) :: x(Points), w(Points)
    INTEGER :: k, info
    DOUBLE PRECISION :: x0, betakm1, qkm2, qkm1, qk, bk, bkp1
    DOUBLE PRECISION :: D(Points), E(Points-1), Z(Points,Points), WORK(MAX(1,2*Points-2))

    x0 = 1.0d0

    ! q_{-1}=0, q_0=1, then build q_1(x0),...,q_{N-1}(x0) via the UNMODIFIED recurrence
    ! beta_{k+1} q_{k+1} = (x0-alpha_k) q_k - beta_k q_{k-1}, alpha_k=0 always.
    qkm2 = 0d0
    qkm1 = 1d0
    DO k = 0,Points-2
       bk = MERGE(0d0, DBLE(k)/DSQRT(4d0*DBLE(k)**2-1d0), k.EQ.0)
       bkp1 = DBLE(k+1)/DSQRT(4d0*DBLE(k+1)**2-1d0)
       qk = (x0*qkm1 - bk*qkm2)/bkp1
       qkm2 = qkm1
       qkm1 = qk
    ENDDO
    ! qkm1 = q_{N-1}(x0), qkm2 = q_{N-2}(x0) (N=Points)

    D = 0d0
    DO k = 1,Points-1
       E(k) = DBLE(k)/DSQRT(4d0*DBLE(k)**2-1d0)
    ENDDO
    betakm1 = DBLE(Points-1)/DSQRT(4d0*DBLE(Points-1)**2-1d0)
    D(Points) = x0 - betakm1*qkm2/qkm1

    CALL DSTEV('V',Points,D,E,Z,Points,WORK,info)
    IF (info.NE.0) THEN
       WRITE(6,*) "GetGaussRadauFactors: DSTEV failed, info=",info
       STOP
    ENDIF
    x = D
    w = 2d0*Z(1,:)**2
  END SUBROUTINE GetGaussRadauFactors
  !****************************************************************************************************
  SUBROUTINE BuildSVDBox(a1,a2,L,NumStates,Order1D,Left1D,Right1D, &
       xNumPoints1D,xMin1D,xMax1D,LegendreFile1D,LegPoints1D,Shift,datadir, &
       NumOpenL,NumOpenR,GridType,BPD,EIG,ODiagL,ODiagR,wsq1Out,wsqLOut)
    USE DataStructures
    USE GlobalVars
    USE Quadrature, ONLY: GetGaussLobattoFactors
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: a1,a2,xMin1D,xMax1D,Shift
    INTEGER, INTENT(in) :: L,NumStates,Order1D,Left1D,Right1D
    INTEGER, INTENT(in) :: xNumPoints1D,LegPoints1D,NumOpenL,NumOpenR
    CHARACTER(LEN=*), INTENT(in) :: LegendreFile1D, datadir
    ! 'Lobatto' (both a1,a2 are DVR nodes -- ordinary interior/last SVD box), 'Radau'
    ! (only a2 is a DVR node, a1 excluded so OneDimChannels never sees it -- but leaves
    ! the box UNCONSTRAINED/open at a1, not Dirichlet -- see BuildSVDBox's grid-setup
    ! comment; correct only when the true regularity condition at a1 actually IS open,
    ! e.g. AlphaFactor=0), or 'LobattoDirichlet' (a1 IS a formal Lobatto node but its DOF
    ! is eliminated -- true Dirichlet u(a1)=0, still never evaluating OneDimChannels at
    ! a1 -- correct when AlphaFactor=(EffDim-1)/2, per delta_function_validation.md).
    ! Caller must pass NumOpenL=0 for 'Radau' or 'LobattoDirichlet' -- neither has a
    ! left-boundary DOF to represent as open.
    CHARACTER(LEN=*), INTENT(in) :: GridType
    TYPE(BPData) :: BPD
    TYPE(GenEigVal) :: EIG
    ! Boundary O-diagonal (O(1,nu;1,nu) and O(L,nu;L,nu)) and the wsq1=1/w(1),
    ! wsqL=1/w(L) normalization factors -- everything CalcLastBoxLam needs to rebuild
    ! JUST the Lam boundary term for a NEW NumOpenR (e.g. a "last box" whose open-
    ! channel count changes every energy) WITHOUT redoing the expensive ARPACK
    ! per-R diagonalization above. Gam0/Overlap never depend on NumOpenL/NumOpenR at
    ! all (see header note), so they stay valid/cacheable regardless.
    DOUBLE PRECISION, INTENT(out) :: ODiagL(NumStates), ODiagR(NumStates)
    DOUBLE PRECISION, INTENT(out) :: wsq1Out, wsqLOut

    DOUBLE PRECISION, ALLOCATABLE :: nodesR(:), weightsR(:)
    DOUBLE PRECISION, ALLOCATABLE :: nodesQ(:), weightsQ(:), wpowQ(:)
    DOUBLE PRECISION, ALLOCATABLE :: Uad(:,:,:), Psi(:,:,:), S_out(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: O(:,:), dXr(:,:), TmatBulk(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: TempPsi(:)
    INTEGER :: PsiDim, HalfBandWidth1D
    INTEGER :: iL,jL,nu,mu2,inu,jmu2,Ufile,LQuad,QOffset
    DOUBLE PRECISION :: RPowExp, RpowIL, wsq1, wsqL
    DOUBLE PRECISION, EXTERNAL :: ddot
    ! GetGaussLobattoFactors' File dummy is CHARACTER*64 (fixed-length, old-style) --
    ! passing a shorter literal directly leaves the dummy's tail unpadded/undefined
    ! (confirmed: caused a "cannot open file ''" crash). Assigning into a like-kind
    ! CHARACTER*64 local blank-pads it correctly before the call.
    CHARACTER*64 :: LobattoFile64

    RPowExp = EffDim - 1d0 - 2d0*AlphaFactor

    !--- geometry: BPD descriptor for a dense (Order=xDim-1) SVD box.
    !    Only xl/xr/NumChannels/xDim/MatrixDim/Order/lam are ever read downstream
    !    (PartitionAndEliminate, BoxMatch -- verified they never touch BPD%u/ux/Pot/
    !    Pcoup/x/xPoints/Left/Right), so those are the only fields set here; the
    !    B-spline-only arrays AllocateBPD would otherwise allocate are skipped.
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

    !--- DVR grid over [a1,a2]: 'Lobatto' (both edges are nodes, ordinary interior/last
    !    box), 'Radau' (only a2 is a node -- avoids ever evaluating OneDimChannels at
    !    R=a1, but does NOT impose Dirichlet there: the retained basis is left
    !    unconstrained/extrapolated at a1, which is the WRONG regularity condition when
    !    AlphaFactor=(EffDim-1)/2 -- confirmed via the free-particle diagnostic and the
    !    Adiabatic-Scattering-BoundStates delta_function_validation.md memory: at that
    !    alpha, u~sqrt(R)*Psi requires strict Dirichlet u(0)=0, not an open boundary), or
    !    'LobattoDirichlet' (true Dirichlet at a1: build an L+1-point Lobatto grid that
    !    DOES include a1 as a formal node, then eliminate that node's DOF entirely --
    !    matching the validated B-spline Left=0 treatment for alpha=(EffDim-1)/2).
    !    LQuad/QOffset let the same dXr/TmatBulk code below work for all three cases: the
    !    quadrature SUM always runs over the FULL LQuad-point grid (an eliminated node's
    !    weight still contributes to the bulk kinetic integral), while the retained BASIS
    !    (nodesR, size L) is always just the last L of those LQuad points -- so
    !    OneDimChannels/Uad is still never evaluated at a1 for 'Radau' or
    !    'LobattoDirichlet' either. ---
    ALLOCATE(nodesR(L),weightsR(L))
    IF (GridType.EQ.'Radau') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBox: GridType='Radau' requires NumOpenL=0 (no left-boundary DOF exists)"
          STOP
       ENDIF
       LQuad = L
       QOffset = 0
       ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
       CALL GetGaussRadauFactors(L,nodesQ,weightsQ)
    ELSEIF (GridType.EQ.'LobattoDirichlet') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBox: GridType='LobattoDirichlet' requires NumOpenL=0 (a1 DOF is eliminated)"
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

    !--- per-R adiabatic channel data (CouplingFlag=0: no d/dR coupling needed) ---
    PsiDim = xNumPoints1D + Order1D - 3
    IF (Left1D.EQ.2)  PsiDim = PsiDim + 1
    IF (Right1D.EQ.2) PsiDim = PsiDim + 1
    HalfBandWidth1D = Order1D
    ALLOCATE(Uad(L,NumStates,2),Psi(L,PsiDim,NumStates),S_out(HalfBandWidth1D+1,PsiDim))
    Ufile = 71
    OPEN(unit=Ufile,file=TRIM(datadir)//'/SVDBoxUad.dat',status='replace')
    CALL OneDimChannels(NumStates,PsiDim,0,0,LegendreFile1D,LegPoints1D,Shift, &
         Order1D,Left1D,Right1D,reducedmass,xNumPoints1D,xMin1D,xMax1D, &
         nodesR,L,Ufile,72,73,74,75,Uad,Psi,S_out,datadir)
    CLOSE(Ufile)

    !--- channel-overlap matrix O(l,nu;l',nu') for ALL DVR-point pairs in the box
    !    (Suno Eq. 20; adapted from 4BodySVD.f90:398-415 CalcOMatrix, but not
    !    restricted to l,l' -- a box here IS the whole sector, unlike FEMSVD's
    !    finer global sector chain) ---
    ALLOCATE(O(L*NumStates,L*NumStates))
    ALLOCATE(TempPsi(PsiDim))
    DO iL = 1,L
       DO nu = 1,NumStates
          CALL dsbmv('U',PsiDim,HalfBandWidth1D,1.0d0,S_out,HalfBandWidth1D+1, &
               Psi(iL,:,nu),1,0.0d0,TempPsi,1)
          inu = (iL-1)*NumStates+nu
          DO jL = 1,L
             DO mu2 = 1,NumStates
                jmu2 = (jL-1)*NumStates+mu2
                O(inu,jmu2) = ddot(PsiDim,TempPsi,1,Psi(jL,:,mu2),1)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(TempPsi)

    !--- DVR derivative matrix: basis functions are Lagrange polynomials through ALL
    !    LQuad quadrature nodes (including any eliminated a1 node for 'LobattoDirichlet'
    !    -- dropping a DOF does not remove that node from the Lagrange interpolant the
    !    OTHER basis functions are built from), so dXr must be evaluated against the FULL
    !    grid even though only the retained L basis indices (jL, mapped to full index
    !    jL+QOffset) survive as independent DOF. For 'Radau'/'Lobatto' (QOffset=0,
    !    LQuad=L) this is identical to the previous all-L-points construction; an R-matrix
    !    box must retain the boundary points as independent DOF there (unlike
    !    4BodySVD.f90's CalcDX, which restricts to bound-state interior points only). ---
    ALLOCATE(dXr(LQuad,L))
    DO iL = 1,LQuad
       DO jL = 1,L
          CALL SVDCalcDX(nodesQ,weightsQ,jL+QOffset,iL,LQuad,dXr(iL,jL))
       ENDDO
    ENDDO

    !--- bulk DVR kinetic matrix, R-power-weighted quadrature (Suno Eq. 18's
    !    integral term; CalcGamLam analog, RMatPropCore.f90:301, no explicit 2*mu
    !    factor -- that convention is applied uniformly below). Quadrature sum runs over
    !    ALL LQuad points -- an eliminated a1 node still contributes its weight to the
    !    bulk integral even though it's not an independent DOF. ---
    ALLOCATE(TmatBulk(L,L))
    DO iL = 1,L
       DO jL = 1,L
          TmatBulk(iL,jL) = SUM(wpowQ(:)*dXr(:,iL)*dXr(:,jL))
       ENDDO
    ENDDO

    !--- assemble EIG%Gam0/Overlap/Lam ---
    EIG%Gam0 = 0d0
    EIG%Overlap = 0d0
    EIG%Lam = 0d0
    DO iL = 1,L
       DO nu = 1,NumStates
          inu = (iL-1)*NumStates+nu
          ! kinetic term: dense over ALL (jL,mu2) -- Suno Eq. 18 first term
          DO jL = 1,L
             DO mu2 = 1,NumStates
                jmu2 = (jL-1)*NumStates+mu2
                EIG%Gam0(inu,jmu2) = EIG%Gam0(inu,jmu2) - TmatBulk(iL,jL)*O(inu,jmu2)
             ENDDO
          ENDDO
          ! overlap + adiabatic-potential terms: DVR-diagonal exact (l=l',nu=nu' only)
          RpowIL = nodesR(iL)**RPowExp
          EIG%Overlap(inu,inu) = 2d0*reducedmass*RpowIL*O(inu,inu)
          EIG%Gam0(inu,inu) = EIG%Gam0(inu,inu) - 2d0*reducedmass*RpowIL*Uad(iL,nu,1)
       ENDDO
    ENDDO

    !--- boundary Lam: diagonal only, matching CalcGamLam's own convention (see header
    !    note) -- pi_l(R_l)=1/sqrt(w_l) for the normalized DVR basis, so the boundary
    !    self-value-squared factor is 1/w instead of B-spline's implicit 1 ---
    wsq1 = 1d0/weightsR(1)
    DO nu = 1,NumOpenL
       inu = nu   ! (1-1)*NumStates+nu
       EIG%Lam(inu,inu) = wsq1*(a1**RPowExp)*O(inu,inu)
    ENDDO
    wsqL = 1d0/weightsR(L)
    DO nu = 1,NumOpenR
       inu = (L-1)*NumStates+nu
       EIG%Lam(inu,inu) = wsqL*(a2**RPowExp)*O(inu,inu)
    ENDDO

    ! stash what CalcLastBoxLam needs to redo this box's Lam for a different NumOpenR
    DO nu = 1,NumStates
       ODiagL(nu) = O(nu,nu)
       ODiagR(nu) = O((L-1)*NumStates+nu,(L-1)*NumStates+nu)
    ENDDO
    wsq1Out = wsq1
    wsqLOut = wsqL

    DEALLOCATE(nodesR,weightsR,nodesQ,weightsQ,wpowQ,Uad,Psi,S_out,O,dXr,TmatBulk)
  END SUBROUTINE BuildSVDBox
  !****************************************************************************************************
  ! Cheaply rebuilds JUST the Lam boundary term for an SVD box acting as the LAST box
  ! in the chain (NumOpenR=NumOpenChannels(Energy), which changes across an energy
  ! scan) -- using the ODiagL/ODiagR/wsq1/wsqL stashed by the original BuildSVDBox
  ! call, without repeating the expensive per-R ARPACK diagonalization. EIG%Gam0/
  ! Overlap are untouched (never depended on NumOpenL/NumOpenR to begin with).
  SUBROUTINE RebuildSVDBoxLam(EIG,ODiagL,ODiagR,wsq1,wsqL,a1,a2,NumStates,L,NumOpenL,NumOpenR)
    USE DataStructures
    USE GlobalVars
    IMPLICIT NONE
    TYPE(GenEigVal) :: EIG
    DOUBLE PRECISION, INTENT(in) :: ODiagL(NumStates), ODiagR(NumStates), wsq1, wsqL, a1, a2
    INTEGER, INTENT(in) :: NumStates, L, NumOpenL, NumOpenR
    INTEGER :: nu, inu
    DOUBLE PRECISION :: RPowExp

    RPowExp = EffDim - 1d0 - 2d0*AlphaFactor
    EIG%Lam = 0d0
    DO nu = 1,NumOpenL
       inu = nu
       EIG%Lam(inu,inu) = wsq1*(a1**RPowExp)*ODiagL(nu)
    ENDDO
    DO nu = 1,NumOpenR
       inu = (L-1)*NumStates+nu
       EIG%Lam(inu,inu) = wsqL*(a2**RPowExp)*ODiagR(nu)
    ENDDO
  END SUBROUTINE RebuildSVDBoxLam
  !****************************************************************************************************
  ! Deallocates only what BuildSVDBox allocated on BPD (lam) -- NOT DeAllocateBPD,
  ! which assumes the B-spline-only arrays (u,ux,Pot,Pcoup,x,xPoints,xBounds,kxMin,
  ! kxMax) were allocated too, which they weren't for an SVD box.
  SUBROUTINE DeAllocateSVDBPD(BPD)
    USE DataStructures
    IMPLICIT NONE
    TYPE(BPData) :: BPD
    DEALLOCATE(BPD%lam)
  END SUBROUTINE DeAllocateSVDBPD
  !****************************************************************************************************
  ! Derivative of the i-th Gauss-Lobatto Lagrange DVR basis function, evaluated at
  ! node l, both ranging over the FULL set of sz nodes (no interior-only restriction --
  ! see header note). Adapted from ~/Documents/GitHub/4BodySVD/4BodySVD.f90:284-306
  ! (CalcDX), which computes the same quantity but only calls it for i restricted to
  ! interior nodes (its bound-state Dirichlet elimination). Includes the standard DVR
  ! normalization factor 1/sqrt(w_i) (pi_i(R)=L_i(R)/sqrt(w_i)).
  SUBROUTINE SVDCalcDX(nodes,weights,i,l,sz,res)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: i,l,sz
    DOUBLE PRECISION, INTENT(in) :: nodes(sz),weights(sz)
    DOUBLE PRECISION, INTENT(out) :: res
    INTEGER j,k
    DOUBLE PRECISION holder(sz),val

    res = 0d0
    DO k = 1,sz
       IF (k.NE.i) THEN
          holder = 0d0
          DO j = 1,sz
             IF (j.NE.i) holder(j) = (nodes(l)-nodes(j))/(nodes(i)-nodes(j))
          ENDDO
          holder(i) = 1d0
          holder(k) = 1d0
          val = PRODUCT(holder)/(nodes(i)-nodes(k))
          res = res + val
       ENDIF
    ENDDO
    res = res*(weights(i)**(-0.5d0))
  END SUBROUTINE SVDCalcDX
  !****************************************************************************************************
END MODULE SVDChannelBasis
