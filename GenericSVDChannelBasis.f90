!****************************************************************************************************
! Generalized drop-in analog of SVDChannelBasis.f90's BuildSVDBox: builds one SVD
! (slow-variable-discretization) R-matrix box exactly the same way (same DVR grid/
! GridType conventions, same channel-overlap/Gam0/Overlap/Lam assembly), but gets its
! per-R adiabatic channel data from AdiabaticSolverGeneric.f90's OneDimChannelsGeneric
! instead of the hardcoded adiabaticSolver1D.f90 OneDimChannels -- so any system with a
! PotentialProc/GridMakerProcI/RobinCoeffProc plugin (lib/AdiabaticInterfaces.f90) can
! build SVD boxes without needing its own hand-rolled per-R diagonalization the way
! RobinSVDChannelBasis.f90 did for the delta-function benchmark.
!
! Grid-consistency requirement is IDENTICAL to BuildSVDBox's own (see that file's header):
! every R point in a box must share the same angular basis, or the shared channel-
! overlap-matrix trick below (dense O over the whole box, single banded S_out projected
! via dsbmv) is wrong. BuildSVDBox enforced this by anchoring its local OneDimChannels'
! GridMakerHHL call at the box's own rightmost R (R(L)) instead of calling fresh every R.
! GridRebuildEveryR here is the same idea, generalized: pass .FALSE. (the safe default
! for any box-building use, matching BuildSVDBox's own hardcoded choice) and this
! subroutine anchors OneDimChannelsGeneric's grid at nodesR(L) automatically -- the
! caller's GridMakerProc is still free to be R-independent entirely (as it is for the
! atom-ion problem, where the angular domain never moves and only a Robin/BC parameter
! would vary with R, if the problem needed one at all). Pass .TRUE. only if a system
! genuinely needs a different angular grid at every DVR node within the SAME box (rare --
! breaks the shared-S_out-projection trick below, would need a per-node O computation
! like RobinSVDChannelBasis.f90's pointwise-evaluation approach instead).
!****************************************************************************************************
MODULE GenericSVDChannelBasis
  IMPLICIT NONE
CONTAINS
  SUBROUTINE BuildSVDBoxGeneric(a1,a2,L,NumStates,Order1D,Left1D,Right1D, &
       xNumPoints1D,xMin1D,xMax1D,LegendreFile1D,LegPoints1D,Shift,datadir, &
       GridRebuildEveryR,MyPotential,MyGridMaker,MyRobinCoeff, &
       NumOpenL,NumOpenR,GridType,BPD,EIG,ODiagL,ODiagR,wsq1Out,wsqLOut)
    USE DataStructures
    USE GlobalVars
    USE Quadrature, ONLY: GetGaussLobattoFactors
    USE AdiabaticInterfaces
    USE AdiabaticSolverGeneric, ONLY: OneDimChannelsGeneric
    USE SVDChannelBasis, ONLY: GetGaussRadauFactors, SVDCalcDX
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: a1,a2,xMin1D,xMax1D,Shift
    INTEGER, INTENT(in) :: L,NumStates,Order1D,Left1D,Right1D
    INTEGER, INTENT(in) :: xNumPoints1D,LegPoints1D,NumOpenL,NumOpenR
    CHARACTER(LEN=*), INTENT(in) :: LegendreFile1D, datadir
    LOGICAL, INTENT(in) :: GridRebuildEveryR
    PROCEDURE(PotentialProc) :: MyPotential
    PROCEDURE(GridMakerProcI) :: MyGridMaker
    PROCEDURE(RobinCoeffProc) :: MyRobinCoeff
    ! Same three GridType conventions as BuildSVDBox -- see that subroutine's own header.
    CHARACTER(LEN=*), INTENT(in) :: GridType
    TYPE(BPData) :: BPD
    TYPE(GenEigVal) :: EIG
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

    !--- DVR grid over [a1,a2]: identical to BuildSVDBox ---
    ALLOCATE(nodesR(L),weightsR(L))
    IF (GridType.EQ.'Radau') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBoxGeneric: GridType='Radau' requires NumOpenL=0 (no left-boundary DOF exists)"
          STOP
       ENDIF
       LQuad = L
       QOffset = 0
       ALLOCATE(nodesQ(LQuad),weightsQ(LQuad))
       CALL GetGaussRadauFactors(L,nodesQ,weightsQ)
    ELSEIF (GridType.EQ.'LobattoDirichlet') THEN
       IF (NumOpenL.NE.0) THEN
          WRITE(6,*) "BuildSVDBoxGeneric: GridType='LobattoDirichlet' requires NumOpenL=0 (a1 DOF is eliminated)"
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

    !--- per-R adiabatic channel data via the pluggable engine (CouplingFlag=0: no d/dR
    !    coupling needed -- SVD needs none, that's the whole point of the method).
    !    GridRebuildEveryR=.FALSE. anchors the angular grid at nodesR(L), the box's own
    !    rightmost R -- same consistency requirement/convention as BuildSVDBox's own
    !    hardcoded GridMakerHHL-anchoring fix (see this file's header). ---
    PsiDim = xNumPoints1D + Order1D - 3
    IF (Left1D.EQ.2)  PsiDim = PsiDim + 1
    IF (Right1D.EQ.2) PsiDim = PsiDim + 1
    HalfBandWidth1D = Order1D
    ALLOCATE(Uad(L,NumStates,2),Psi(L,PsiDim,NumStates),S_out(HalfBandWidth1D+1,PsiDim))
    Ufile = 71
    OPEN(unit=Ufile,file=TRIM(datadir)//'/SVDBoxUad.dat',status='replace')
    CALL OneDimChannelsGeneric(NumStates,PsiDim,0,0,LegendreFile1D,LegPoints1D,Shift, &
         Order1D,Left1D,Right1D,reducedmass,xNumPoints1D,xMin1D,xMax1D, &
         nodesR,L,GridRebuildEveryR,nodesR(L), &
         MyPotential,MyGridMaker,MyRobinCoeff, &
         Ufile,72,73,74,75,Uad,Psi,S_out,datadir)
    CLOSE(Ufile)

    !--- channel-overlap matrix O(l,nu;l',nu') for ALL DVR-point pairs in the box --
    !    identical to BuildSVDBox: a box IS the whole sector here, dense over all L. ---
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

    !--- DVR derivative matrix / bulk kinetic matrix: identical to BuildSVDBox ---
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

    !--- assemble EIG%Gam0/Overlap/Lam -- identical formulas to BuildSVDBox ---
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
          EIG%Gam0(inu,inu) = EIG%Gam0(inu,inu) - 2d0*reducedmass*RpowIL*Uad(iL,nu,1)
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

    DEALLOCATE(nodesR,weightsR,nodesQ,weightsQ,wpowQ,Uad,Psi,S_out,O,dXr,TmatBulk)
  END SUBROUTINE BuildSVDBoxGeneric
END MODULE GenericSVDChannelBasis
