!****************************************************************************************************
! Generalized drop-in analog of Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90's
! OneDimChannels, parameterized via AdiabaticInterfaces.f90's PotentialProc/GridMakerProcI/
! RobinCoeffProc instead of a hardcoded `use potential`+sumpairwisepot(DR)+GridMakerHHL+CalcBasisFuncs
! pipeline. Built as a NEW, separate subroutine (OneDimChannelsGeneric) rather than editing
! OneDimChannels in place -- that subroutine is not a module procedure, so adding optional/pluggable
! arguments to it directly would require every existing caller across every repo that links it
! (DriveAdiabaticSolver1D.f90, SVDChannelBasis.f90, PlotWavefunctionSVD.f90, SVDBound.f90, SVDRmat.f90,
! ...) to pick up an explicit interface (i.e. a USE statement) at the same time -- out of scope for a
! first validation of this design. This file lets the design be proven out (see the delta-function/
! Robin-BC test this was built for) with zero risk to anything already depending on OneDimChannels.
!
! Both subroutines here ARE module procedures (unlike OneDimChannels/CalcHamiltonian themselves,
! which are bare external subroutines) so that every caller gets real compile-time argument-list
! checking (types/ranks/count, and which concrete procedure gets passed to which PROCEDURE(...)
! dummy) via `USE AdiabaticSolverGeneric`, rather than relying on hand-matched implicit-interface
! external linkage the way the rest of this codebase's F77-style subroutines do. Worth doing now,
! while there's only one caller (this file's own delta-function test) -- much cheaper than
! retrofitting after other callers exist.
!
! Everything NOT specific to potential/grid/BC construction (ARPACK shift-invert setup via MyDsband,
! phase-fixing via FixPhase, Hellmann-Feynman coupling via CalcPHF, CalcEigenErrors, CalcOverlap) is
! reused UNCHANGED from adiabaticSolver1D.f90/matrix_stuff.f90 -- both already linked into this
! repo's build, so those external subroutines are called directly, not reimplemented here (calling a
! plain external subroutine from inside a module procedure needs no special treatment).
!
! Two behavioral generalizations beyond the three AdiabaticInterfaces axes:
!   - Grid rebuild cadence: OneDimChannels always rebuilds the angular grid every R (calling
!     GridMakerHHL fresh each iteration). Real use cases need other cadences too (R-matrix-prop's
!     own SVD-box fork needs "once, anchored at the box's rightmost R"; this file's own
!     delta-function test needs "once, ever" since its grid never depends on R at all) -- controlled
!     here by GridRebuildEveryR (logical) + RAnchor (used only when .FALSE., before the R-loop).
!   - Basis-function rebuild is NOT purely gated on "did the grid change" the way canonical's
!     CalcNewBasisFunc flag assumes: a Robin BC's basis functions depend on kLeft/kRight, which can
!     be genuinely R-dependent (this file's own test: kRight=sqrt(2*mu)*R) even when the geometric
!     grid itself is completely fixed. So whenever Left=3 .OR. Right=3, u/uxx are rebuilt fresh every
!     R regardless of GridRebuildEveryR -- exactly matching RobinSVDChannelBasis.f90's own
!     unconditional per-R MakeBasis call, which this file's test validates against.
!****************************************************************************************************
MODULE AdiabaticSolverGeneric
  IMPLICIT NONE
CONTAINS
SUBROUTINE OneDimChannelsGeneric(NumStates,PsiDim,PsiFlag,CouplingFlag,LegendreFile,&
     LegPoints,Shift,Order,Left,Right,mu,xNumPoints,xMin,xMax,&
     R,RSteps,GridRebuildEveryR,RAnchor,&
     MyPotential,MyGridMaker,MyRobinCoeff,&
     Ufile,Pfile,Qfile,dPfile,VQfile,Uad,Psi,S_out,datadir)
  USE AdiabaticInterfaces
  USE quadrature, ONLY: GetGaussFactors
  IMPLICIT NONE
  INTEGER Ufile,Pfile,Qfile,dPfile,VQfile,PsiDim
  INTEGER LegPoints,xNumPoints
  INTEGER NumStates,PsiFlag,Order,Left,Right
  INTEGER RSteps,CouplingFlag,CalcNewBasisFunc
  LOGICAL, INTENT(in) :: GridRebuildEveryR
  DOUBLE PRECISION, INTENT(in) :: RAnchor
  PROCEDURE(PotentialProc) :: MyPotential
  PROCEDURE(GridMakerProcI) :: MyGridMaker
  PROCEDURE(RobinCoeffProc) :: MyRobinCoeff
  DOUBLE PRECISION Shift
  DOUBLE PRECISION xMin,xMax
  DOUBLE PRECISION R(RSteps)
  DOUBLE PRECISION, ALLOCATABLE :: xPoints(:)
  DOUBLE PRECISION :: Psi(RSteps,PsiDim,NumStates),Uad(RSteps,NumStates,2)
  LOGICAL, ALLOCATABLE :: Select(:)

  INTEGER iparam(11),ncv,info
  INTEGER i,j,iR,nQ
  INTEGER LeadDim,MatrixDim,HalfBandWidth
  INTEGER xDim
  INTEGER, ALLOCATABLE :: iwork(:)
  INTEGER, ALLOCATABLE :: xBounds(:)
  DOUBLE PRECISION Tol
  DOUBLE PRECISION mu, kLeft, kRight

  DOUBLE PRECISION, ALLOCATABLE :: LUFac(:,:),workl(:)
  DOUBLE PRECISION, ALLOCATABLE :: workd(:),Residuals(:)
  DOUBLE PRECISION, ALLOCATABLE :: xLeg(:),wLeg(:)
  DOUBLE PRECISION, ALLOCATABLE :: u(:,:,:),uxx(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: S(:,:),H(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: lPsi(:,:),mPsi(:,:),Energies(:,:)

  DOUBLE PRECISION, ALLOCATABLE :: P(:,:),Q(:,:),dP(:,:),prevP(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: V_band(:,:), dV_band(:,:)

  CHARACTER*64 LegendreFile
  CHARACTER(LEN=64), INTENT(in) :: datadir
  CHARACTER(LEN=128) :: PsiFilename
  DOUBLE PRECISION :: S_out(Order+1, *)

  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  xDim = xNumPoints+Order-3
  IF (Left .EQ. 2) xDim = xDim + 1
  IF (Right .EQ. 2) xDim = xDim + 1

  MatrixDim = xDim
  HalfBandWidth = Order
  LeadDim = 3*HalfBandWidth+1

  ALLOCATE(xPoints(xNumPoints))
  ALLOCATE(xBounds(xNumPoints+2*Order))
  ALLOCATE(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim))
  ALLOCATE(S(HalfBandWidth+1,MatrixDim),H(HalfBandWidth+1,MatrixDim))
  ALLOCATE(P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates))

  ncv = MIN(MatrixDim,10*NumStates)
  LeadDim = 3*HalfBandWidth+1
  ALLOCATE(iwork(MatrixDim))
  ALLOCATE(Select(ncv))
  ALLOCATE(LUFac(LeadDim,MatrixDim))
  ALLOCATE(workl(ncv*ncv+8*ncv))
  ALLOCATE(workd(3*MatrixDim))
  ALLOCATE(lPsi(MatrixDim,ncv),mPsi(MatrixDim,ncv))
  ALLOCATE(Residuals(MatrixDim))
  ALLOCATE(Energies(ncv,2))
  ALLOCATE(V_band(HalfBandWidth+1,MatrixDim), dV_band(HalfBandWidth+1,MatrixDim))
  ALLOCATE(prevP(NumStates,NumStates))
  info = 0
  iR = 1
  CalcNewBasisFunc = 1
  Tol = 1e-20
  OPEN(unit=16, file=TRIM(datadir)//'/QuickVQmat.dat', status='replace')

  ! Grid built once, before the loop, when not adaptive per-R (see this file's header).
  IF (.NOT. GridRebuildEveryR) THEN
     CALL MyGridMaker(RAnchor,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
  ENDIF

  DO iR = 1,RSteps
     IF (GridRebuildEveryR) THEN
        CALL MyGridMaker(R(iR),xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
     ENDIF
     CALL MyRobinCoeff(R(iR),kLeft,kRight)
     ! Robin basis functions depend on kLeft/kRight, not just grid geometry -- always rebuild
     ! when a Robin BC is in play, regardless of GridRebuildEveryR (see this file's header).
     IF (CalcNewBasisFunc.EQ.1 .OR. Left.EQ.3 .OR. Right.EQ.3) THEN
        CALL CalcBasisFuncsBP(Left,Right,kLeft,kRight,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,u)
        CALL CalcBasisFuncsBP(Left,Right,kLeft,kRight,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,uxx)
     ENDIF

     CALL CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,xBounds,HalfBandWidth,S)

     CALL CalcHamiltonianGeneric(R(iR),mu,Order,xPoints,LegPoints,xLeg,wLeg,xDim,&
          xNumPoints,u,uxx,xBounds,HalfBandWidth,H,V_band,dV_band,MyPotential)

     CALL MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,&
          H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,NumStates,&
          Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,&
          ncv*ncv+8*ncv,iwork,info)
     IF (iR .GT. 1) CALL FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,lPsi,mPsi)
     lPsi = mPsi
     CALL CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,&
          HalfBandWidth,NumStates,mPsi,Energies,ncv)
     IF (info.EQ.0 .AND. DABS(Energies(1,1)).GT.1d-3) THEN
        IF (Energies(1,1).LT.0d0) THEN
           Shift = 1.2d0*Energies(1,1)
        ELSE
           Shift = 0.8d0*Energies(1,1)
        ENDIF
     ENDIF
     WRITE(Ufile,20) R(iR),(Energies(i,1), i = 1,MIN(NumStates,iparam(5)))
     WRITE(6,20) R(iR),(Energies(i,1), i = 1,MIN(NumStates,iparam(5)))

     IF (CouplingFlag .NE. 0) THEN
        nQ = MIN(NumStates,iparam(5))
        CALL CalcPHF(nQ,MatrixDim,HalfBandWidth,R(iR),mu,&
             Energies,ncv,mPsi,S,V_band,dV_band,1.d-6,P)
        CALL dgemm('N','N',NumStates,NumStates,NumStates,&
             -1.d0,P,NumStates,P,NumStates,0.d0,Q,NumStates)
        IF (iR .GT. 1) THEN
           dP = (P - prevP) / (R(iR) - R(iR-1))
        ELSE
           dP = 0.d0
        ENDIF
        prevP = P

        WRITE(Pfile,*) R(iR)
        WRITE(Qfile,*) R(iR)
        WRITE(dPfile,*) R(iR)
        WRITE(VQfile,*) R(iR)
        DO i = 1,nQ
           WRITE(Pfile,20)  (P(i,j), j = 1,nQ)
           WRITE(Qfile,20)  (Q(i,j), j = 1,nQ)
           WRITE(dPfile,20) (dP(i,j), j = 1,nQ)
           WRITE(VQfile,20) (MERGE(Energies(i,1), 0.d0, i==j) &
                + (0.5d0/mu)*Q(i,j), j=1,nQ)
        ENDDO
        WRITE(16,20) R(iR), &
             ((Energies(i,1) - (0.5d0/mu)*Q(i,i)), i=1,nQ)
     ENDIF
     IF (PsiFlag .NE. 0) THEN
        OPEN(unit=97, file=TRIM(datadir)//'/PhiGrid.dat', status='replace')
        DO i = 1,xNumPoints
           WRITE(97,*) xPoints(i)
        ENDDO
        CLOSE(unit=97)
        WRITE(PsiFilename,'(a,i4.4,a)') TRIM(datadir)//'/Psi_R', iR, '.dat'
        OPEN(unit=999+iR, file=TRIM(PsiFilename), status='replace')
        DO i = 1,MatrixDim
           WRITE(999+iR,20) (mPsi(i,j), j = 1,NumStates)
        ENDDO
        CLOSE(unit=999+iR)
     ENDIF

     DO i=1,NumStates
        Uad(iR,i,1)=Energies(i,1)
        Uad(iR,i,2)=Energies(i,2)
     ENDDO
     DO i=1,PsiDim
        DO j=1,NumStates
           Psi(iR,i,j)=mPsi(i,j)
        ENDDO
     ENDDO
  ENDDO
  CLOSE(unit=16)

  S_out(1:HalfBandWidth+1, 1:MatrixDim) = S
  DEALLOCATE(S,H)
  DEALLOCATE(Energies)
  DEALLOCATE(iwork)
  DEALLOCATE(Select)
  DEALLOCATE(LUFac)
  DEALLOCATE(workl)
  DEALLOCATE(workd)
  DEALLOCATE(lPsi,mPsi)
  DEALLOCATE(Residuals)
  DEALLOCATE(P,Q,dP,prevP)
  DEALLOCATE(V_band,dV_band)
  DEALLOCATE(xPoints)
  DEALLOCATE(xLeg,wLeg)
  DEALLOCATE(xBounds)
  DEALLOCATE(u,uxx)

20 FORMAT(1P,100E16.8)
END SUBROUTINE OneDimChannelsGeneric
!****************************************************************************************************
! Generalized drop-in analog of adiabaticSolver1D.f90's CalcHamiltonian: identical kinetic-matrix
! (<u|u''>) and potential-matrix-element assembly, but potential VALUES come from a caller-supplied
! PotentialProc instead of a hardcoded `use potential`+sumpairwisepot(DR)+pairwise-distance pipeline.
! Dropping that pipeline also drops the r12/r13/r23/cosx/cosxp/cosxm machinery entirely -- a
! concrete PotentialProc implementation owns whatever geometry conversion IT needs internally,
! given just the (R,phi) quadrature points.
!
! NOTE on the <u|u''> convention: this does NOT need an explicit Robin/Bloch surface-term
! correction the way a <u'|u'> formulation (e.g. RobinSVDChannelBasis.f90's own Tphi) does.
! Integration by parts gives INT(u*u'') = [u*u']_boundary - INT(u'*u'). As long as the basis
! functions themselves already satisfy the boundary condition pointwise (CalcBasisFuncsBP's
! Right=3 combined function is built so u'(edge)=kRight*u(edge) exactly), the <u|u''> integral
! already has that boundary information baked in via the basis functions' own derivative values --
! no separate correction term needed. Verified numerically (not just by this derivation) against
! RobinSVDChannelBasis.f90/KMSFormulas.f's SolveQ in this file's own validation test.
!****************************************************************************************************
SUBROUTINE CalcHamiltonianGeneric(R,mu,Order,xPoints,LegPoints,xLeg,wLeg,xDim,&
     xNumPoints,u,uxx,xBounds,HalfBandWidth,H,V_band,dV_band,MyPotential)
  USE AdiabaticInterfaces
  IMPLICIT NONE
  INTEGER Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
  DOUBLE PRECISION R,mu
  DOUBLE PRECISION xPoints(*),xLeg(*),wLeg(*)
  DOUBLE PRECISION H(HalfBandWidth+1,xDim)
  DOUBLE PRECISION V_band(HalfBandWidth+1,xDim), dV_band(HalfBandWidth+1,xDim)
  DOUBLE PRECISION u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)
  PROCEDURE(PotentialProc) :: MyPotential

  INTEGER ix,ixp,kx,lx
  INTEGER Row,NewRow,Col
  INTEGER, ALLOCATABLE :: kxMin(:,:),kxMax(:,:)
  DOUBLE PRECISION a,m
  DOUBLE PRECISION xTempT,xTempV,ax,bx,xScaledZero
  DOUBLE PRECISION, ALLOCATABLE :: xIntScale(:),xT(:,:),xV(:,:),xdV(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: x0(:,:),Pot(:,:),dPot(:,:)

  ALLOCATE(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim),xdV(xDim,xDim))
  ALLOCATE(x0(LegPoints,xNumPoints))
  ALLOCATE(kxMin(xDim,xDim),kxMax(xDim,xDim))
  ALLOCATE(Pot(LegPoints,xNumPoints), dPot(LegPoints,xNumPoints))

  m = -1.0d0/(2.0d0*mu*R*R)

  DO kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     DO lx = 1,LegPoints
        x0(lx,kx) = xIntScale(kx)*xLeg(lx)+xScaledZero
     ENDDO
  ENDDO

  DO ix = 1,xDim
     DO ixp = 1,xDim
        kxMin(ixp,ix) = MAX(xBounds(ix),xBounds(ixp))
        kxMax(ixp,ix) = MIN(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
     ENDDO
  ENDDO

  CALL MyPotential(R,LegPoints,xNumPoints,x0,Pot,dPot)

  DO ix = 1,xDim
     DO ixp = MAX(1,ix-Order),MIN(xDim,ix+Order)
        xT(ix,ixp) = 0.0d0
        xV(ix,ixp) = 0.0d0
        xdV(ix,ixp) = 0.0d0
        DO kx = kxMin(ixp,ix),kxMax(ixp,ix)
           xTempT = 0.0d0
           xTempV = 0.0d0
           DO lx = 1,LegPoints
              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
              xTempT = xTempT + a*uxx(lx,kx,ixp)
              xTempV = xTempV + a*(Pot(lx,kx))*u(lx,kx,ixp)
              xdV(ix,ixp) = xdV(ix,ixp) + a*(dPot(lx,kx))*u(lx,kx,ixp)
           ENDDO
           xT(ix,ixp) = xT(ix,ixp) + xTempT
           xV(ix,ixp) = xV(ix,ixp) + xTempV
        ENDDO
     ENDDO
  ENDDO

  H = 0.0d0
  V_band = 0.0d0
  dV_band = 0.0d0
  DO ix = 1,xDim
     Row=ix
     DO ixp = MAX(1,ix-Order),MIN(xDim,ix+Order)
        Col = ixp
        IF (Col .GE. Row) THEN
           NewRow = HalfBandWidth+1+Row-Col
           H(NewRow,Col) = (m*xT(ix,ixp)+xV(ix,ixp))
           V_band(NewRow,Col)  = xV(ix,ixp)
           dV_band(NewRow,Col) = xdV(ix,ixp)
        ENDIF
     ENDDO
  ENDDO

  DEALLOCATE(Pot,dPot,x0)
  DEALLOCATE(xIntScale,xT,xV,xdV)
  DEALLOCATE(kxMin,kxMax)
END SUBROUTINE CalcHamiltonianGeneric
END MODULE AdiabaticSolverGeneric
