!****************************************************************************************************
! Analytic drop-in replacement for SVDChannelBasis.f90's BuildSVDBox, specific to the equal-mass
! three-boson delta-function problem of Mehta, Esry & Greene, PRA 76, 022711 (2007) -- the same
! system DeltaScat/KMSFormulas.f already solves exactly. Where BuildSVDBox gets its per-R channel
! data (Uad, and the channel-function overlap needed to relate different R's) by numerically
! diagonalizing OneDimChannels' B-spline Hamiltonian at every DVR node and projecting through a
! shared angular B-spline overlap matrix, this file evaluates q_nu(R)/kappa_nu(R) directly via
! KMSFormulas.f's SolveQ (exact, warm-started Newton solve of the transcendental equation) and gets
! the cross-R channel overlap directly via KMSFormulas.f's EvalOverlapAt (exact closed form, see
! that subroutine's own header) -- no B-spline basis, no ARPACK diagonalization, no basis-set error
! anywhere. Everything else (the DVR radial grid, derivative matrix, kinetic bulk matrix, and the
! final Gam0/Overlap/Lam assembly) is IDENTICAL to BuildSVDBox -- only the source of the per-R
! channel data (Uad, O) differs, so this box should be interchangeable with a BuildSVDBox-built one
! anywhere downstream (CombineGam/PartitionAndEliminate/BoxMatch, RebuildSVDBoxLam).
!
! Purpose: eliminates ALL sources of numerical/basis-set error from the hybrid SVD/adiabatic
! propagator's box construction for this specific benchmark problem, so a hybrid-propagator run can
! be compared directly against DeltaScat/the exact Mehta-Shepard result the same way
! RMATPROPAdiabatic_seeded.f90 isolated near-origin interpolation error from the B-spline pipeline.
!****************************************************************************************************
MODULE AnalyticSVDChannelBasis
  IMPLICIT NONE
CONTAINS
  !****************************************************************************************************
  ! Signature mirrors BuildSVDBox with every B-spline-only argument (Order1D/Left1D/Right1D/
  ! xNumPoints1D/xMin1D/xMax1D/LegendreFile1D/LegPoints1D/Shift/datadir) dropped -- none of them
  ! are needed since there is no OneDimChannels call here at all. reducedmass/AlphaFactor/EffDim
  ! are still read from GlobalVars, exactly as BuildSVDBox does, and must already be set to this
  ! problem's fixed values (see AnalyticDeltaPotential.f90's DeltaMu/DeltaAlpha/DeltaDdim) by the
  ! caller before this is invoked.
  SUBROUTINE BuildSVDBoxAnalytic(a1,a2,L,NumStates,NumOpenL,NumOpenR,GridType, &
       BPD,EIG,ODiagL,ODiagR,wsq1Out,wsqLOut)
    USE DataStructures
    USE GlobalVars
    USE SVDChannelBasis, ONLY: GetGaussRadauFactors, SVDCalcDX
    USE Quadrature, ONLY: GetGaussLobattoFactors
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: a1,a2
    INTEGER, INTENT(in) :: L,NumStates,NumOpenL,NumOpenR
    ! Same three GridType conventions as BuildSVDBox -- see that subroutine's own header comment.
    CHARACTER(LEN=*), INTENT(in) :: GridType
    TYPE(BPData) :: BPD
    TYPE(GenEigVal) :: EIG
    DOUBLE PRECISION, INTENT(out) :: ODiagL(NumStates), ODiagR(NumStates)
    DOUBLE PRECISION, INTENT(out) :: wsq1Out, wsqLOut

    DOUBLE PRECISION, ALLOCATABLE :: nodesR(:), weightsR(:)
    DOUBLE PRECISION, ALLOCATABLE :: nodesQ(:), weightsQ(:), wpowQ(:)
    DOUBLE PRECISION, ALLOCATABLE :: qAll(:,:), Uad(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: O(:,:), dXr(:,:), TmatBulk(:,:)
    INTEGER :: iL,jL,nu,mu2,inu,jmu2,LQuad,QOffset,nEscalate,niter
    DOUBLE PRECISION :: RPowExp, RpowIL, wsq1, wsqL, qguess, Rval
    DOUBLE PRECISION, EXTERNAL :: EvalU
    CHARACTER*64 :: LobattoFile64

    RPowExp = EffDim - 1d0 - 2d0*AlphaFactor

    !--- geometry: identical to BuildSVDBox ---
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

    !--- DVR grid over [a1,a2]: identical logic/three GridType cases to BuildSVDBox ---
    ALLOCATE(nodesR(L),weightsR(L))
    IF (GridType.EQ.'Radau') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBoxAnalytic: GridType='Radau' requires NumOpenL=0 (no left-boundary DOF exists)"
          STOP
       ENDIF
       LQuad = L
       QOffset = 0
       ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
       CALL GetGaussRadauFactors(L,nodesQ,weightsQ)
    ELSEIF (GridType.EQ.'LobattoDirichlet') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBoxAnalytic: GridType='LobattoDirichlet' requires NumOpenL=0 (a1 DOF is eliminated)"
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

    !--- exact per-R channel data: q_nu(R_l) via SolveQ (KMSFormulas.f), warm-started in
    !    DESCENDING R order per channel -- matches AnalyticDeltaPotential.f90's own convention,
    !    needed because channel 1's q grows ~R*sqrt(2*mu) at large R (a fixed small-R guess fails
    !    to converge there), while channels >1's kappa approaches a fixed constant (6*nu-9) --
    !    see AnalyticDeltaPotential.f90's qguess initialization for the same two branches. ---
    ALLOCATE(qAll(L,NumStates),Uad(L,NumStates))
    DO nu = 1,NumStates
       IF (nu.EQ.1) THEN
          qguess = nodesR(L)*DSQRT(2d0*reducedmass)
       ELSE
          qguess = 6.0d0*DBLE(nu)-9.0d0+1d-8
       ENDIF
       DO iL = L,1,-1
          Rval = nodesR(iL)
          CALL SolveQ(nu-1,Rval,reducedmass,qguess,qAll(iL,nu),nEscalate,niter)
          qguess = qAll(iL,nu)
          Uad(iL,nu) = EvalU(nu,Rval,reducedmass,qAll(iL,nu))
       ENDDO
    ENDDO

    !--- channel-overlap matrix O(l,nu;l',nu') for ALL DVR-point pairs -- EXACT closed form
    !    (EvalOverlapAt), no B-spline projection/ARPACK diagonalization anywhere. ---
    ALLOCATE(O(L*NumStates,L*NumStates))
    DO iL = 1,L
       DO nu = 1,NumStates
          inu = (iL-1)*NumStates+nu
          DO jL = 1,L
             DO mu2 = 1,NumStates
                jmu2 = (jL-1)*NumStates+mu2
                CALL EvalOverlapAt(nu,qAll(iL,nu),mu2,qAll(jL,mu2),O(inu,jmu2))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !--- DVR derivative matrix / bulk kinetic matrix: IDENTICAL to BuildSVDBox -- these are
    !    properties of the radial DVR grid alone, independent of the angular/channel basis. ---
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

    !--- assemble EIG%Gam0/Overlap/Lam -- IDENTICAL formulas to BuildSVDBox, just fed the
    !    analytically-computed O/Uad instead of the numerically-computed ones. ---
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

    !--- boundary Lam: diagonal only, identical convention to BuildSVDBox ---
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

    ! stash what RebuildSVDBoxLam needs to redo this box's Lam for a different NumOpenR --
    ! RebuildSVDBoxLam itself is reused UNCHANGED from SVDChannelBasis (it only consumes these
    ! stashed values plus geometry, never the numeric-vs-analytic origin of O/Uad).
    DO nu = 1,NumStates
       ODiagL(nu) = O(nu,nu)
       ODiagR(nu) = O((L-1)*NumStates+nu,(L-1)*NumStates+nu)
    ENDDO
    wsq1Out = wsq1
    wsqLOut = wsqL

    DEALLOCATE(nodesR,weightsR,nodesQ,weightsQ,wpowQ,qAll,Uad,O,dXr,TmatBulk)
  END SUBROUTINE BuildSVDBoxAnalytic
END MODULE AnalyticSVDChannelBasis
