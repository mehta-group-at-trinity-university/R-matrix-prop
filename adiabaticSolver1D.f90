Subroutine OneDimChannels(NumStates,PsiDim,PsiFlag,CouplingFlag,LegendreFile,&
     LegPoints,Shift,Order,Left,Right,mu,xNumPoints,xMin,xMax,&
     R,RSteps,&
     !RFirst,RLast,&
     !rhogrid,phigrid,krho,kphi,rhoknot,phiknot,Nrho,Nphi,Ubcoef,&
     Ufile,Pfile,Qfile,dPfile,VQfile,Uad,Psi,S_out,datadir)
  use quadrature, only: GetGaussFactors
  use potential

  implicit none
!  integer Nrho, Nphi,krho,kphi
  integer Ufile,Pfile,Qfile,dPfile,VQfile,iphi,PsiDim
!  double precision rhogrid(Nrho),phigrid(Nphi),rhoknot(Nrho+krho),phiknot(Nphi+kphi),Ubcoef(Nrho*Nphi)
  integer LegPoints,xNumPoints
  integer NumStates,PsiFlag,Order,Left,Right
  integer RSteps,CouplingFlag,CalcNewBasisFunc
  double precision mass,Shift,xtest
  DOUBLE PRECISION RFirst,RLast,XFirst,XLast,StepX
  double precision xMin,xMax
  !  double precision, allocatable :: R(:)
  double precision R(RSteps)
  double precision, allocatable :: xPoints(:)
  double precision :: Psi(RSteps,PsiDim,NumStates),Uad(RSteps,NumStates,2)
  logical, allocatable :: Select(:)

  integer iparam(11),ncv,info
  integer i,j,k,iR,nQ
  integer LeadDim,MatrixDim,HalfBandWidth
  integer xDim
  integer, allocatable :: iwork(:)
  integer, allocatable :: xBounds(:)
  double precision Tol
  double precision TotalMemory
  double precision mu, mu12, r0diatom, dDiatom, etaOVERpi, Pi

  double precision, allocatable :: LUFac(:,:),workl(:)
  double precision, allocatable :: workd(:),Residuals(:)
  double precision, allocatable :: xLeg(:),wLeg(:)
  double precision, allocatable :: u(:,:,:),uxx(:,:,:)
  double precision, allocatable :: S(:,:),H(:,:)
  double precision, allocatable :: lPsi(:,:),mPsi(:,:),Energies(:,:)

  double precision, allocatable :: P(:,:),Q(:,:),dP(:,:),prevP(:,:)
  double precision, allocatable :: V_band(:,:), dV_band(:,:)

  common/MassInfo/mu12,r0diatom,dDiatom

  character*64 LegendreFile
  character(len=64), intent(in) :: datadir
  character(len=128) :: PsiFilename
  double precision :: S_out(Order+1, *)   ! returned B-spline overlap (same at all R)
!  write(6,*) LegendreFile

  !     SET UP A LOG-GRID IN THE HYPERRADIUS
  !XFirst = dlog10(RFirst)
  !XLast = dlog10(RLast)
  !StepX=(XLast-XFirst)/(RSteps-1.d0)

  !allocate(R(RSteps))
  !do i = 1,RSteps
     !     read(5,*) R(i)
     !     R(i)= (XFirst+(i-1)*StepX)**3
    ! R(i)= 10.d0**(XFirst+(i-1)*StepX)
  !enddo

  !      if (mod(xNumPoints,2) .ne. 0) then
  !         write(6,*) 'xNumPoints not divisible by 2'
  !         xNumPoints = (xNumPoints/2)*2
  !         write(6,*) '   truncated to ',xNumPoints
  !      endif

  ! test the interpolated function
!  do i = 1, RSteps
!     write(3,*) '#', R(i)
!     do iphi = 1, Nphi
!        write(3,*) R(i), phigrid(iphi), dbs2vl(R(i),phigrid(iphi),krho,kphi,rhoknot,phiknot,Nrho,Nphi,Ubcoef)
!     enddo
!     write(3,*) ' '
!  enddo


  !     mu=0.5
!  mu=mass/dsqrt(3.0d0)
  !      mu=2.d0*mu12
  !      mu=mass
  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  xDim = xNumPoints+Order-3
  if (Left .eq. 2) xDim = xDim + 1
  if (Right .eq. 2) xDim = xDim + 1

  MatrixDim = xDim
  HalfBandWidth = Order
  LeadDim = 3*HalfBandWidth+1


  TotalMemory = 2.0d0*(HalfBandWidth+1)*MatrixDim ! S, H
  TotalMemory = TotalMemory + 2.0d0*LegPoints*(Order+2)*xDim ! x splines
  TotalMemory = TotalMemory + 2.0d0*NumStates*NumStates ! P and Q matrices
  TotalMemory = TotalMemory + LeadDim*MatrixDim ! LUFac
  TotalMemory = TotalMemory + 4.0d0*NumStates*MatrixDim ! channel functions
  TotalMemory = TotalMemory + 4*xDim*xDim ! (CalcHamiltonian)
  TotalMemory = TotalMemory + LegPoints**2*xNumPoints ! (CalcHamiltonian)
  TotalMemory = 8.0d0*TotalMemory/(1024.0d0*1024.0d0)

!  write(6,*)
!  write(6,*) 'MatrixDim ',MatrixDim
!  write(6,*) 'HalfBandWidth ',HalfBandWidth
!  write(6,*) 'Approximate peak memory usage (in Mb) ',TotalMemory
!  write(6,*)

  allocate(xPoints(xNumPoints))
  allocate(xBounds(xNumPoints+2*Order))
  allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim))
  allocate(S(HalfBandWidth+1,MatrixDim),H(HalfBandWidth+1,MatrixDim))
  allocate(P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates))

  ! Widened from 3*NumStates -- see the same fix in
  ! Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90 for the reasoning.
  ncv = MIN(MatrixDim,10*NumStates)
  LeadDim = 3*HalfBandWidth+1
  allocate(iwork(MatrixDim))
  allocate(Select(ncv))
  allocate(LUFac(LeadDim,MatrixDim))
  allocate(workl(ncv*ncv+8*ncv))
  allocate(workd(3*MatrixDim))
  allocate(lPsi(MatrixDim,ncv),mPsi(MatrixDim,ncv))
  allocate(Residuals(MatrixDim))
  allocate(Energies(ncv,2))
  allocate(V_band(HalfBandWidth+1,MatrixDim), dV_band(HalfBandWidth+1,MatrixDim))
  allocate(prevP(NumStates,NumStates))
  info = 0
  iR=1
  CalcNewBasisFunc=1
  Tol=1e-20
  ! Shift was previously hardcoded to -20.d0 here unconditionally, silently discarding
  ! whatever the caller passed in. For a genuinely repulsive problem (all eigenvalues
  ! positive -- e.g. V2Depth=0, a free-particle test case) that mismatch makes
  ! ARPACK's shift-invert Lanczos fail to build its factorization (dsaupd/dsband
  ! info=-9999) at essentially every R. The caller's own Shift is now the starting
  ! point instead; the adaptive per-R update below is unchanged. Kept in sync with the
  ! same fix in Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90.
  open(unit=16, file=trim(datadir)//'/QuickVQmat.dat', status='replace')

  ! Angular grid/basis built ONCE here, outside the R-loop -- NOT once per R the way
  ! stock OneDimChannels does it. This local copy is only ever called by
  ! SVDChannelBasis.f90's BuildSVDBox with the L DVR points of a SINGLE R-matrix box,
  ! and needs the SAME angular basis at every one of them (it captures S_out once and
  ! reuses it via dsbmv for every cross-R-point channel overlap -- Suno Eq. 20).
  ! GridMakerHHL's own grid is R-dependent above xRswitch=10*PotRange/Pi (its refined-
  ! near-coalescence branch has deltax=PotRange/R, shrinking continuously with R), so
  ! calling it fresh per R (as stock OneDimChannels does) breaks that requirement --
  ! confirmed empirically via a free-particle diagnostic (self-overlap off by up to
  ! 20%, spurious cross-channel coupling up to 0.46, for a box reaching above
  ! xRswitch). Still want the coalescence-adjacent refinement for a REAL potential,
  ! though (that's the resolution GridMakerHHL exists to provide) -- so rather than
  ! falling back to a plain uniform grid unconditionally, build the grid ONCE per box
  ! via a single GridMakerHHL call at R(RSteps), the box's own largest/rightmost R
  ! (BuildSVDBox always passes its L DVR nodes in ascending order, so this is the
  ! box's right edge a2). GridMakerHHL's refined region narrows monotonically with R
  ! (deltax=PotRange/R), so anchoring on the box's largest R gives the grid that is at
  ! least as finely resolved as any smaller R in the same box would independently
  ! request -- erring toward more resolution rather than less for the rest of the box.
  call GridMakerHHL(R(RSteps),PotRange,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,u)
  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,uxx)
  ! S depends only on the (now fixed) basis, not on R -- computed once, matching what
  ! the S_out output argument's own comment already documented as the assumption
  ! ("returned B-spline overlap (same at all R)").
  call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,xBounds,HalfBandWidth,S)

  do iR = 1,RSteps
     call CalcHamiltonian(R(iR),mu,Order,xPoints,LegPoints,xLeg,wLeg,xDim,&
          xNumPoints,u,uxx,xBounds,HalfBandWidth,H,V_band,dV_band)

     call MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,&
          H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,NumStates,&
          Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,&
          ncv*ncv+8*ncv,iwork,info)
     if (iR .gt. 1) call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,lPsi,mPsi)
     lPsi = mPsi   ! save as phase reference for next R-step
     call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,&
          HalfBandWidth,NumStates,mPsi,Energies,ncv)
     ! Track the ground state's scale on WHICHEVER side of zero it's on, staying
     ! safely on the same side (never crossing zero) so ARPACK's shift-invert keeps
     ! targeting the lowest NumStates eigenvalues. Previously only the negative
     ! branch was handled, leaving Shift stuck at its (possibly badly mismatched)
     ! starting value for any purely-repulsive channel (e.g. V2Depth=0, a free-
     ! particle test) -- that caused dsaupd/dsband to fail (info=-9999) at most R.
     ! Kept in sync with the same fix in Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90.
     ! Floored at 1d-3 in magnitude -- see the same fix in
     ! Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90 for the reasoning
     ! (a near-zero ground state otherwise chases Shift onto the eigenvalue itself,
     ! making (H-Shift*S) numerically singular for shift-invert).
     IF (info.EQ.0 .AND. DABS(Energies(1,1)).GT.1d-3) THEN
        IF (Energies(1,1).LT.0d0) THEN
           Shift = 1.2d0*Energies(1,1)
        ELSE
           Shift = 0.8d0*Energies(1,1)
        ENDIF
     ENDIF
     write(Ufile,20) R(iR),(Energies(i,1), i = 1,min(NumStates,iparam(5)))
     write(6,20) R(iR),(Energies(i,1), i = 1,min(NumStates,iparam(5)))

     if (CouplingFlag .ne. 0) then
        nQ = min(NumStates,iparam(5))
        call CalcPHF(nQ,MatrixDim,HalfBandWidth,R(iR),mu,&
             Energies,ncv,mPsi,S,V_band,dV_band,1.d-6,P)
        call dgemm('N','N',NumStates,NumStates,NumStates,&
             -1.d0,P,NumStates,P,NumStates,0.d0,Q,NumStates)
        if (iR .gt. 1) then
           dP = (P - prevP) / (R(iR) - R(iR-1))
        else
           dP = 0.d0
        endif
        prevP = P

        write(Pfile,*) R(iR)
        write(Qfile,*) R(iR)
        write(dPfile,*) R(iR)
        write(VQfile,*) R(iR)
        do i = 1,nQ
           write(Pfile,20)  (P(i,j), j = 1,nQ)
           write(Qfile,20)  (Q(i,j), j = 1,nQ)
           write(dPfile,20) (dP(i,j), j = 1,nQ)
!           write(VQfile,20) (merge(Energies(i,1) - 0.25d0/(2d0*mu*R(iR)**2), 0.d0, i==j) &
!                + (0.5d0/mu)*Q(i,j), j=1,nQ)
           write(VQfile,20) (merge(Energies(i,1), 0.d0, i==j) &
                + (0.5d0/mu)*Q(i,j), j=1,nQ)
        enddo
        write(16,20) R(iR), &
!             ((Energies(i,1) - 0.25d0/(2d0*mu*R(iR)**2) - (0.5d0/mu)*Q(i,i)), i=1,nQ)
             ((Energies(i,1) - (0.5d0/mu)*Q(i,i)), i=1,nQ)
     endif
     !         write(400,20) R(iR),R(iR)**3.0d0*(Energies(2,1)-Q(2,2))
     if (PsiFlag .ne. 0) then
        open(unit=97, file=trim(datadir)//'/PhiGrid.dat', status='replace')
        do i = 1,xNumPoints
           write(97,*) xPoints(i)
        enddo
        close(unit=97)
        write(PsiFilename,'(a,i4.4,a)') trim(datadir)//'/Psi_R', iR, '.dat'
        open(unit=999+iR, file=trim(PsiFilename), status='replace')
        do i = 1,MatrixDim
           write(999+iR,20) (mPsi(i,j), j = 1,NumStates)
        enddo
        close(unit=999+iR)
     endif

     do i=1,NumStates
        Uad(iR,i,1)=Energies(i,1)
        Uad(iR,i,2)=Energies(i,2)
     end do
     do i=1,PsiDim
        do j=1,NumStates
           Psi(iR,i,j)=mPsi(i,j)
        end do
     end do

     
  enddo
  close(unit=16)

  S_out(1:HalfBandWidth+1, 1:MatrixDim) = S
  deallocate(S,H)
  deallocate(Energies)
  deallocate(iwork)
  deallocate(Select)
  deallocate(LUFac)
  deallocate(workl)
  deallocate(workd)
  deallocate(lPsi,mPsi)
  deallocate(Residuals)
  deallocate(P,Q,dP,prevP)
  deallocate(V_band,dV_band)
  deallocate(xPoints)
  deallocate(xLeg,wLeg)
  deallocate(xBounds)
  deallocate(u,uxx)

  !deallocate(R)

10 format(1P,100e25.15)
19 format(1P,100e16.8)
20 format(1P,100e22.12)
1002 format(a64)


end subroutine OneDimChannels
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,&
     xNumPoints,u,xBounds,HalfBandWidth,S)
  implicit none
  integer Order,LegPoints,xDim,xNumPoints,xBounds(xNumPoints+2*Order),HalfBandWidth
  double precision xPoints(*),xLeg(*),wLeg(*)
  double precision S(HalfBandWidth+1,xDim)
  double precision u(LegPoints,xNumPoints,xDim)

  integer ix,ixp,kx,lx
  integer i1,i1p
  integer Row,NewRow,Col
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  double precision a,b,m
  double precision xTempS
  double precision ax,bx
  double precision, allocatable :: xIntScale(:),xS(:,:)

  allocate(xIntScale(xNumPoints),xS(xDim,xDim))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))

  S = 0.0d0

  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
  enddo

  !      do ix=1,xNumPoints+2*Order
  !         print*, ix, xBounds(ix)
  !      enddo

  do ix = 1,xDim
     do ixp = 1,xDim
        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
     enddo
  enddo

  do ix = 1,xDim
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        xS(ixp,ix) = 0.0d0
        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
           xTempS = 0.0d0
           do lx = 1,LegPoints
              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
              b = a*u(lx,kx,ixp)
              xTempS = xTempS + b
           enddo
           xS(ixp,ix) = xS(ixp,ix) + xTempS
        enddo
     enddo
  enddo

  do ix = 1,xDim
     Row=ix
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        Col = ixp
        if (Col .ge. Row) then
           NewRow = HalfBandWidth+1+Row-Col
           S(NewRow,Col) = xS(ixp,ix)
           !               write(26,*) ix,ixp,S(NewRow,Col)
        endif
     enddo
  enddo

  deallocate(xIntScale,xS)
  deallocate(kxMin,kxMax)

  return
end subroutine CalcOverlap
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcHamiltonian(R,mu,Order,xPoints,LegPoints,xLeg,wLeg,xDim,&
     xNumPoints,u,uxx,xBounds,HalfBandWidth,H,V_band,dV_band)
     !rhogrid,phigrid,krho,kphi,rhoknot,phiknot,Nrho,Nphi,Ubcoef

  use potential
  implicit none
  !integer Nrho, Nphi,krho,kphi
  !double precision rhogrid(Nrho),phigrid(Nphi),rhoknot(Nrho+krho),phiknot(Nphi+kphi),Ubcoef(Nrho*Nphi)
  integer Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
  double precision R,mu
  double precision xPoints(*),xLeg(*),wLeg(*)
  double precision H(HalfBandWidth+1,xDim)
  double precision V_band(HalfBandWidth+1,xDim), dV_band(HalfBandWidth+1,xDim)
  double precision u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)

  integer ix,ixp,kx,lx
  integer i1,i1p
  integer Row,NewRow,Col
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  double precision a,b,m,Pi
  double precision Rall,r12,r12a,r23,r23a,r23b,r23c,r13,r13a,r13b,r13c,r14,r24,r34
  double precision u1,sys_ss_pot,V12,V23,V31
  double precision VInt,VTempInt,potvalue, xTempV
  !     double precision TempPot,VInt,VTempInt
  double precision x,ax,bx,xScaledZero,xTempT,xTempS,xInt
  double precision, allocatable :: Pot(:,:), dPot(:,:)
  double precision, allocatable :: xIntScale(:),xT(:,:),xV(:,:),xdV(:,:)
  double precision, allocatable :: cosx0(:,:),x0(:,:),cosxp(:,:),cosxm(:,:),cosx(:,:),sinx(:,:)

  double precision mu12,r0diatom,dDiatom
  !      common/MassInfo/mu12,r0diatom,dDiatom

  allocate(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim),xdV(xDim,xDim))
  allocate(cosx0(LegPoints,xNumPoints),x0(LegPoints,xNumPoints),cosxp(LegPoints,xNumPoints),cosxm(LegPoints,xNumPoints))
  allocate(cosx(LegPoints,xNumPoints))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
  allocate(Pot(LegPoints,xNumPoints), dPot(LegPoints,xNumPoints))

  Pi = 3.1415926535897932385d0

  m = -1.0d0/(2.0d0*mu*R*R)

  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     do lx = 1,LegPoints
        x0(lx,kx) = xIntScale(kx)*xLeg(lx)+xScaledZero
        x = x0(lx,kx)
!        write(6,*) 'x0(',lx,',',kx,') = ', x0(lx,kx)
!        write(2,*)  x0(lx,kx)
        cosx(lx,kx) = dcos(x)
        cosxp(lx,kx) = dcos(x+Pi/3.0d0)
        cosxm(lx,kx) = dcos(x-Pi/3.0d0)
     enddo
  enddo

  do ix = 1,xDim
     do ixp = 1,xDim
        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
     enddo
  enddo

  do kx = 1,xNumPoints-1
     do lx = 1,LegPoints
        !     comment this out when the surface is a pure function of R and phi (i.e. x)
        r12 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosx(lx,kx))
        r13 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxm(lx,kx))
        r23 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxp(lx,kx))
        call sumpairwisepot(r12, r13, r23, potvalue)
        !potvalue = dbs2vl(R,x0(lx,kx),krho,kphi,rhoknot,phiknot,Nrho,Nphi,Ubcoef)
        Pot(lx,kx) = potvalue
        call sumpairwisepotDR(r12, r13, r23, R, potvalue)
        dPot(lx,kx) = potvalue
        !            write(24,*) kx, lx, Pot(lx,kx)
     enddo
  enddo

  do ix = 1,xDim
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        xT(ix,ixp) = 0.0d0
        xV(ix,ixp) = 0.0d0
        xdV(ix,ixp) = 0.0d0
        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
           xTempT = 0.0d0
           xTempV = 0.0d0
           do lx = 1,LegPoints
              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
              xTempT = xTempT + a*uxx(lx,kx,ixp)
              xTempV = xTempV + a*(Pot(lx,kx))*u(lx,kx,ixp)
              xdV(ix,ixp) = xdV(ix,ixp) + a*(dPot(lx,kx))*u(lx,kx,ixp)
           enddo
           xT(ix,ixp) = xT(ix,ixp) + xTempT
           xV(ix,ixp) = xV(ix,ixp) + xTempV
        enddo
     enddo
  enddo

  H = 0.0d0
  V_band = 0.0d0
  dV_band = 0.0d0
  do ix = 1,xDim
     Row=ix
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        Col = ixp
        if (Col .ge. Row) then
           NewRow = HalfBandWidth+1+Row-Col
           H(NewRow,Col) = (m*xT(ix,ixp)+xV(ix,ixp))
           V_band(NewRow,Col)  = xV(ix,ixp)
           dV_band(NewRow,Col) = xdV(ix,ixp)
        endif
     enddo
  enddo



  deallocate(Pot, dPot)
  deallocate(xIntScale,xT,xV,xdV)
  deallocate(cosx0,x0,cosxp,cosxm,cosx)
  deallocate(kxMin,kxMax)


  return
end subroutine CalcHamiltonian

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcPMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim
  double precision RDelt
  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
  double precision S(HalfBandWidth+1,MatrixDim)
  double precision P(NumStates,NumStates)

  integer i,j,k
  double precision a,ddot
  double precision, allocatable :: TempPsi1(:),TempPsi2(:)

  allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))

  a = 0.5d0/RDelt

  do j = 1,NumStates
     do k = 1,MatrixDim
        TempPsi1(k) = rPsi(k,j)-lPsi(k,j)
     enddo
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
     do i = 1,NumStates
        P(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  return
end subroutine CalcPMatrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcQMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,Q)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim
  double precision RDelt
  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
  double precision S(HalfBandWidth+1,MatrixDim)
  double precision Q(NumStates,NumStates)

  integer i,j,k
  double precision a,ddot
  double precision, allocatable :: TempPsi1(:),TempPsi2(:)

  allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))

  a = 1.0d0/(RDelt**2)

  do j = 1,NumStates
     do k = 1,MatrixDim
        TempPsi1(k) = lPsi(k,j)+rPsi(k,j)-2.0d0*mPsi(k,j)
     enddo
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
     do i = 1,NumStates
        Q(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  return
end subroutine CalcQMatrix
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!      
subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,mPsi,rPsi)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim,ncv
  double precision S(HalfBandWidth+1,MatrixDim),Psi(MatrixDim,ncv)
  double precision mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv)

  integer i,j
  double precision Phase,ddot
  double precision, allocatable :: TempPsi(:)

  allocate(TempPsi(MatrixDim))

  do i = 1,NumStates
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rPsi(1,i),1,0.0d0,TempPsi,1)
     Phase = ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
     if (Phase .lt. 0.0d0) then
        do j = 1,MatrixDim
           rPsi(j,i) = -rPsi(j,i)
        enddo
     endif
  enddo

  deallocate(TempPsi)

  return
end subroutine FixPhase
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMaker222(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
  implicit none
  integer xNumPoints,yNumPoints
  double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

  integer i,j,k
  double precision Pi
  double precision r0New
  double precision xRswitch,yRswitch
  double precision xDelt,x0,x1,x2
  double precision yDelt,y0,y1,y2

  Pi = 3.1415926535897932385d0

  r0New = 3.0d0*r0

  xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0))))

  if (R .gt. xRswitch) then
     x0 = xMin
     x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
     x2 = xMax
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints/2)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x0
        k = k + 1
     enddo
     xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x1
        k = k + 1
     enddo
  else
     x0 = xMin
     x1 = xMax
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints-1)
     do i = 1,xNumPoints
        xPoints(k) = (i-1)*xDelt + x0
        k = k + 1
     enddo
  endif

  yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))

  if (R .gt. yRswitch) then
     y0 = yMin
     y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
     y2 = yMax
     k = 1
     yDelt = (y1-y0)/dfloat(yNumPoints/2)
     do i = 1,yNumPoints/2
        yPoints(k) = (i-1)*yDelt + y0
        k = k + 1
     enddo
     yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
     do i = 1,yNumPoints/2
        yPoints(k) = (i-1)*yDelt + y1
        k = k + 1
     enddo
  else
     y0 = yMin
     y1 = yMax
     k = 1
     yDelt = (y1-y0)/dfloat(yNumPoints-1)
     do i = 1,yNumPoints
        yPoints(k) = (i-1)*yDelt + y0
        k = k + 1
     enddo
  endif

  return
end subroutine GridMaker222

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine GridMaker111(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
  implicit none
  integer xNumPoints,yNumPoints
  double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)


  integer i,j,k
  double precision Pi
  double precision r0New
  double precision xRswitch,yRswitch
  double precision xDelt,x0,x1,x2
  double precision yDelt,y0,y1,y2

  Pi = 3.1415926535897932385d0


  r0New = 3.0d0*r0

  xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0)) ))


  if (R .gt. xRswitch) then
     x0 = xMin
     x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
     x2 = xMax
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints/2)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x0
        k = k + 1
     enddo
     xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x1
        k = k + 1
     enddo
  else
     x0 = xMin
     x1 = xMax
     k = 1
     xDelt = (x1-x0)/dfloat(xNumPoints-1)
     do i = 1,xNumPoints
        xPoints(k) = (i-1)*xDelt + x0
        k = k + 1
     enddo
  endif

  yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))

  if (R .gt. yRswitch) then
     y0 = yMin
     y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
     y2 = yMax
     k = 1
     yDelt = (y1-y0)/dfloat(yNumPoints/2)
     do i = 1,yNumPoints/2
        yPoints(k) = (i-1)*yDelt + y0
        k = k + 1
     enddo
     yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
     do i = 1,yNumPoints/2
        yPoints(k) = (i-1)*yDelt + y1
        k = k + 1
     enddo
  else
     y0 = yMin
     y1 = yMax
     k = 1
     yDelt = (y1-y0)/dfloat(yNumPoints-1)
     do i = 1,yNumPoints
        yPoints(k) = (i-1)*yDelt + y0
        k = k + 1
     enddo
  endif

  return
end subroutine GridMaker111

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! NOTE: upstream adiabaticSolver1D.f90's GridMaker(mu,R,r0,...) subroutine was removed
! from this local copy -- it is dead code (never called by OneDimChannels, which uses
! GridMakerHHL instead; confirmed only a commented-out reference remains upstream) and
! its name collides at link time with RMatPropCore.f90's own, actively-used
! GridMaker(grid,numpts,E1,E2,scale). Removing this unreachable duplicate avoids the
! collision without touching either upstream file's real functionality.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerHHL(R, r0, xNumPoints, xMin, xMax, xPoints, CalcNewBasisFunc)
  ! Non-uniform phi grid that concentrates knots near the two coalescence angles:
  !   phi23 = pi/6  (r23=0, interior, for equal-mass 1D 3-body coordinates)
  !   xMax  = pi/2  (r12=0, right boundary)
  ! For small R (R <= xRswitch) a uniform grid is used instead.
  ! The angular refinement half-width is deltax = r0/R, so it narrows with R.
  implicit none
  integer, intent(in)    :: xNumPoints, CalcNewBasisFunc
  double precision, intent(in)  :: R, r0, xMin, xMax
  double precision, intent(out) :: xPoints(xNumPoints)

  double precision :: Pi, phi23, deltax, xRswitch
  double precision :: x0, x1, x2, x3, x4, xDelt
  integer :: i, k

  Pi     = 3.1415926535897932385d0
  phi23  = Pi/6.0d0                    ! r23=0 coalescence angle (equal-mass coords)
  deltax = r0/R                        ! angular half-width of refined region
  xRswitch = 10.0d0*r0/Pi             ! refine for R > xRswitch (~3.2 for r0=1)

  if (R .gt. xRswitch) then
     x0 = xMin
     x1 = phi23 - deltax
     x2 = phi23 + deltax
     x3 = xMax  - deltax
     x4 = xMax
     ! Guard against deltax pushing segments to degenerate or reversed widths
     if (x1 .le. x0) x1 = x0 + (phi23 - x0)*0.5d0
     if (x2 .ge. x3) x2 = x3 - (x3 - phi23)*0.5d0
     k = 1
     ! Segment 1: [xMin, phi23-deltax] — coarse
     xDelt = (x1-x0)/dble(xNumPoints/4)
     do i = 1, xNumPoints/4
        xPoints(k) = x0 + (i-1)*xDelt
        k = k + 1
     enddo
     ! Segment 2: [phi23-deltax, phi23+deltax] — refined near r23=0
     xDelt = (x2-x1)/dble(xNumPoints/4)
     do i = 1, xNumPoints/4
        xPoints(k) = x1 + (i-1)*xDelt
        k = k + 1
     enddo
     ! Segment 3: [phi23+deltax, xMax-deltax] — coarse
     xDelt = (x3-x2)/dble(xNumPoints/4)
     do i = 1, xNumPoints/4
        xPoints(k) = x2 + (i-1)*xDelt
        k = k + 1
     enddo
     ! Segment 4: [xMax-deltax, xMax] — refined near r12=0
     xDelt = (x4-x3)/dble(xNumPoints/4-1)
     do i = 1, xNumPoints/4
        xPoints(k) = x3 + (i-1)*xDelt
        k = k + 1
     enddo
  else
     ! Uniform grid for small R where coalescence singularities are less severe
     x0 = xMin
     xDelt = (xMax-x0)/dble(xNumPoints-1)
     do i = 1, xNumPoints
        xPoints(i) = x0 + (i-1)*xDelt
     enddo
  endif

  return
end subroutine GridMakerHHL
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcCoupling(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P,Q,dP)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDim
  double precision RDelt
  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
  double precision S(HalfBandWidth+1,MatrixDim),testorth
  double precision P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates)

  integer i,j,k
  double precision aP,aQ,ddot
  double precision, allocatable :: lDiffPsi(:),rDiffPsi(:),TempPsi(:),TempPsiB(:),rSumPsi(:)
  double precision, allocatable :: TempmPsi(:)

  allocate(lDiffPsi(MatrixDim),rDiffPsi(MatrixDim),TempPsi(MatrixDim),&
       TempPsiB(MatrixDim),rSumPsi(MatrixDim))
  allocate(TempmPsi(MatrixDim))

      aP = 0.5d0/RDelt
      aQ = aP*aP

      do j = 1,NumStates
       do k = 1,MatrixDim
        rDiffPsi(k) = rPsi(k,j)-lPsi(k,j)
	rSumPsi(k)  = lPsi(k,j)+mPsi(k,j)+rPsi(k,j)
       enddo
       call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)
       call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rSumPsi,1,0.0d0,TempPsiB,1)
      do i = 1,NumStates
        P(i,j) = aP*ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
	dP(i,j)= ddot(MatrixDim,mPsi(1,i),1,TempPsiB,1)
        do k = 1,MatrixDim
         lDiffPsi(k) = rPsi(k,i)-lPsi(k,i)
        enddo
        Q(i,j) = -aQ*ddot(MatrixDim,lDiffPsi,1,TempPsi,1)
       enddo
      enddo

	do j=1,NumStates
	 do i=j,NumStates
	  dP(i,j)=2.d0*aQ*(dP(i,j)-dP(j,i))
	  dP(j,i)=-dP(i,j)
	 enddo
	enddo


  deallocate(lDiffPsi,rDiffPsi,TempPsi,rSumPsi,TempPsiB,TempmPsi)

  return
end subroutine CalcCoupling

subroutine CalcPHF(NumStates, MatrixDim, HalfBandWidth, R, mu, &
                   Energies, ncv, Psi, S, V_band, dV_band, degTol, P_out)
  ! Compute P_{ij} = [(2/R)*<Phi_i|V|Phi_j> + <Phi_i|dV/dR|Phi_j>] / (E_j - E_i)
  ! using the Hellmann-Feynman theorem on the scaled Hamiltonian h~ = -T + (2muR^2)V.
  ! Degenerate pairs (|E_j - E_i| < degTol) and diagonal elements are set to 0.
  implicit none
  integer, intent(in)    :: NumStates, MatrixDim, HalfBandWidth, ncv
  double precision, intent(in) :: R, mu, degTol
  double precision, intent(in) :: Energies(ncv,2)
  double precision, intent(in) :: Psi(MatrixDim,ncv)
  double precision, intent(in) :: S(HalfBandWidth+1,MatrixDim)
  double precision, intent(in) :: V_band(HalfBandWidth+1,MatrixDim)
  double precision, intent(in) :: dV_band(HalfBandWidth+1,MatrixDim)
  double precision, intent(out) :: P_out(NumStates,NumStates)

  integer :: i, j
  double precision :: Vij, dVij, dE, ddot
  double precision, allocatable :: tmpV(:), tmpdV(:)

  allocate(tmpV(MatrixDim), tmpdV(MatrixDim))

  P_out = 0.d0
  do j = 1, NumStates
     call dsbmv('U', MatrixDim, HalfBandWidth, 1.d0, V_band,  HalfBandWidth+1, Psi(1,j), 1, 0.d0, tmpV,  1)
     call dsbmv('U', MatrixDim, HalfBandWidth, 1.d0, dV_band, HalfBandWidth+1, Psi(1,j), 1, 0.d0, tmpdV, 1)
     do i = 1, NumStates
        if (i .eq. j) cycle
        dE = Energies(j,1) - Energies(i,1)
        if (dabs(dE) .lt. degTol) cycle
        Vij  = ddot(MatrixDim, Psi(1,i), 1, tmpV,  1)
        dVij = ddot(MatrixDim, Psi(1,i), 1, tmpdV, 1)
        P_out(i,j) = ((2.d0/R)*Vij + dVij) / dE
     enddo
  enddo

  deallocate(tmpV, tmpdV)
end subroutine CalcPHF



      











  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$  (NumStates,PsiFlag,Shift,&
!!$     Order,Left,Right,alpha,mu,NumMasses,&
!!$     xNumPoints,xMin,xMax,&
!!$     RSteps,RDerivDelt,RFirst,RLast,&
!!$     R,Uad,Psi,eDim,psiDim,S,sDim,r0)
!!$  use quadrature
!!$  implicit none
!!$
!!$  integer xNumPoints,NumMasses
!!$  integer NumStates,PsiFlag,Order,Left,Right
!!$  integer RSteps,CouplingFlag,CalcNewBasisFunc
!!$  double precision alpha,mass,Shift,m1,m2,m3
!!$  double precision RLeft,RRight,RDerivDelt
!!$  double precision RFirst,RLast,XFirst,XLast,StepX
!!$  double precision xMin,xMax,R(RSteps)
!!$  double precision, allocatable :: xPoints(:)
!!$
!!$  logical, allocatable :: Select(:)
!!$
!!$  integer iparam(11),ncv,info
!!$  integer i,j,k,iR
!!$  integer LeadDim,MatrixDim,HalfBandWidth
!!$  integer xDim,eDim,psiDim,sDim
!!$  integer, allocatable :: iwork(:)
!!$  integer, allocatable :: xBounds(:)
!!$  double precision Tol
!!$  double precision TotalMemory
!!$  double precision r0,massarray(NumMasses)
!!$  double precision mu, mu12,mu123,r0diatom, dDiatom, etaOVERpi, Pi, phi12, phi23, mgamma
!!$
!!$  double precision, allocatable :: LUFac(:,:),workl(:)
!!$  double precision, allocatable :: workd(:),Residuals(:)
!!$  double precision, allocatable :: u(:,:,:),uxx(:,:,:),H(:,:)
!!$  double precision, allocatable :: S(:,:)
!!$  double precision, allocatable :: lPsi(:,:),mPsi(:,:),rPsi(:,:),Energies(:,:)
!!$  double precision, allocatable :: P(:,:),Q(:,:),dP(:,:)
!!$  double precision :: Psi(RSteps,psiDim,eDim),Uad(RSteps,eDim,2)
!!$
!!$  common/MassInfo/mu12,r0diatom,dDiatom
!!$
!!$!  character*64 LegendreFile
!!$
!!$  Pi=dacos(-1.d0)
!!$  write(6,*) 'Pi=',Pi
!!$
!!$!  m1 = massarray(1)
!!$!  m2 = massarray(2)
!!$!  m3 = massarray(3)
!!$!  mu = 3.d0**(-0.5d0)
!!$  !  mu12=m1*m2/(m1+m2)
!!$  !  mu123=(m1+m2)*m3/(m1+m2+m3)
!!$  !  mu=dsqrt(mu12*mu123)
!!$  Pi=dacos(-1.d0)
!!$  write(6,*) 'Pi=',Pi, 'mu = ', mu
!!$!  mgamma = mu/m1
!!$!  phi12=Pi/2
!!$!  phi23=datan(mgamma)
!!$  r0 = 1.d0
!!$!  stop
!!$
!!$
!!$  !      if (mod(xNumPoints,2) .ne. 0) then
!!$  !         write(6,*) 'xNumPoints not divisible by 2'
!!$  !         xNumPoints = (xNumPoints/2)*2
!!$  !         write(6,*) '   truncated to ',xNumPoints
!!$  !      endif
!!$
!!$  RLeft = 0.0d0
!!$  allocate(xLeg(LegPoints),wLeg(LegPoints))
!!$
!!$  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)
!!$
!!$  xDim = xNumPoints+Order-3
!!$  if (Left .eq. 2) xDim = xDim + 1
!!$  if (Right .eq. 2) xDim = xDim + 1
!!$
!!$  MatrixDim = xDim
!!$  HalfBandWidth = Order
!!$  LeadDim = 3*HalfBandWidth+1
!!$
!!$
!!$  TotalMemory = 2.0d0*(HalfBandWidth+1)*MatrixDim ! S, H
!!$  TotalMemory = TotalMemory + 2.0d0*LegPoints*(Order+2)*xDim ! x splines
!!$  TotalMemory = TotalMemory + 2.0d0*NumStates*NumStates ! P and Q matrices
!!$  TotalMemory = TotalMemory + LeadDim*MatrixDim ! LUFac
!!$  TotalMemory = TotalMemory + 4.0d0*NumStates*MatrixDim ! channel functions
!!$  TotalMemory = TotalMemory + 4*xDim*xDim ! (CalcHamiltonian)
!!$  TotalMemory = TotalMemory + LegPoints**2*xNumPoints ! (CalcHamiltonian)
!!$  TotalMemory = 8.0d0*TotalMemory/(1024.0d0*1024.0d0)
!!$
!!$  write(6,*)
!!$  write(6,*) 'MatrixDim ',MatrixDim
!!$  write(6,*) 'HalfBandWidth ',HalfBandWidth
!!$  write(6,*) 'Approximate peak memory usage (in Mb) ',TotalMemory
!!$  write(6,*)
!!$  write(6,*) "xNumPoints = ", xNumPoints
!!$  allocate(xPoints(xNumPoints))
!!$  allocate(xBounds(xNumPoints+2*Order))
!!$  allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim))
!!$  allocate(S(HalfBandWidth+1,MatrixDim))
!!$  allocate(H(HalfBandWidth+1,MatrixDim))
!!$  allocate(P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates))
!!$
!!$  ncv = 2*NumStates
!!$  LeadDim = 3*HalfBandWidth+1
!!$  allocate(iwork(MatrixDim))
!!$  allocate(Select(ncv))
!!$  allocate(LUFac(LeadDim,MatrixDim))
!!$  allocate(workl(ncv*ncv+8*ncv))
!!$  allocate(workd(3*MatrixDim))
!!$  allocate(lPsi(MatrixDim,ncv),mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv))
!!$  allocate(Residuals(MatrixDim))
!!$  allocate(Energies(ncv,2))
!!$  info = 0
!!$  iR=1
!!$
!!$  Tol=1e-20
!!$  print*, 'calling GridMaker'
!!$  call GridMaker(mu,R(1),2.0d0, xNumPoints,xMin,xMax,xPoints)
!!$  !call GridMakerHHL(mu,mu12,mu123,phi23,R(iR),r0,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
!!$  
!!$  print*, 'done... Calculating Basis functions'
!!$  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,u)
!!$  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,uxx)
!!$  !     must move this block inside the loop if the grid is adaptive
!!$  print*, 'done... Calculating overlap matrix'
!!$  call CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,xBounds,HalfBandWidth,S)
!!$!
!!$  print*, "Done.  Now entering the radial loop"
!!$  
!!$!  write(6,*) "array sizes", size(S,1), size(S,2), size(H,1), size(H,2)
!!$!  stop
!!$  do iR = 1,RSteps     
!!$!     write(6,*) "About to calculate Hamiltonian.  iR = ", iR
!!$     call CalcHamiltonian(alpha,R(iR),mu,Order,xPoints,&
!!$          LegPoints,xLeg,wLeg,xDim,&
!!$          xNumPoints,u,uxx,xBounds,&
!!$          HalfBandWidth,H)
!!$ !    write(6,*) "Done.  Now diagonalizing..."
!!$     call MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,&
!!$          H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,NumStates,&
!!$          Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,&
!!$          ncv*ncv+8*ncv,iwork,info)
!!$ !    write(6,*) "Done.  Fixing phase..."
!!$     if (iR .ne. 1) call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,Psi(iR-1,:,:),mPsi)
!!$     
!!$     call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,HalfBandWidth,NumStates,mPsi,Energies,ncv)
!!$ !    write(6,*) "Done."
!!$!     IF(R(iR).GT. 2.2d0) Shift = 0.93d0*Energies(1,1)
!!$     !         do i = 1,min(NumStates,iparam(5))
!!$     !     Energies(i,1) = Energies(i,1) + 1.875d0/(mu*R(iR)*R(iR))
!!$     !         enddo
!!$     write(200,20) R(iR),(Energies(i,1), i = 1,min(NumStates,iparam(5)))
!!$     write(6,20) R(iR),(Energies(i,1), i = 1,min(NumStates,iparam(5)))
!!$!     write(6,*)
!!$  !   write(6,*) R(iR)
!!$
!!$     do i=1,eDim
!!$        Uad(iR,i,1)=Energies(i,1)
!!$        Uad(iR,i,2)=Energies(i,2)
!!$     end do
!!$     do i=1,MatrixDim
!!$        do j=1,eDim
!!$           Psi(iR,i,j)=mPsi(i,j)
!!$        end do
!!$     end do
!!$     
!!$     !         write(400,20) R(iR),R(iR)**3.0d0*(Energies(2,1)-Q(2,2))
!!$     if (PsiFlag .ne. 0) then
!!$        do i = 1,xNumPoints
!!$           write(97,*) xPoints(i)
!!$        enddo
!!$        do i = 1,MatrixDim
!!$           write(999+iR,20) (mPsi(i,j), j = 1,NumStates)
!!$        enddo
!!$        close(unit=999+iR)
!!$     endif
!!$
!!$  enddo
!!$
!!$  deallocate(S)
!!$  deallocate(H)
!!$  deallocate(Energies)
!!$  deallocate(iwork)
!!$  deallocate(Select)
!!$  deallocate(LUFac)
!!$  deallocate(workl)
!!$  deallocate(workd)
!!$  deallocate(lPsi,mPsi,rPsi)
!!$  deallocate(Residuals)
!!$  deallocate(P,Q,dP)
!!$  deallocate(xPoints)
!!$  deallocate(xLeg,wLeg)
!!$  deallocate(xBounds)
!!$  deallocate(u,uxx)
!!$
!!$!  deallocate(R)
!!$
!!$10 format(1P,100e25.15)
!!$20 format(1P,100e22.12)
!!$1002 format(a64)
!!$
!!$  stop
!!$end Subroutine adiabaticSolver1D
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine CalcOverlap(Order,xPoints,LegPoints,xLeg,wLeg,xDim,xNumPoints,u,xBounds,HalfBandWidth,S)
!!$  implicit none
!!$  integer Order,LegPoints,xDim,xNumPoints,xBounds(xNumPoints+2*Order),HalfBandWidth
!!$  double precision xPoints(xNumPoints),xLeg(LegPoints),wLeg(LegPoints)
!!$  double precision S(HalfBandWidth+1,xDim)
!!$  double precision u(LegPoints,xNumPoints,xDim)
!!$
!!$  integer ix,ixp,kx,lx
!!$  integer i1,i1p
!!$  integer Row,NewRow,Col
!!$  integer, allocatable :: kxMin(:,:),kxMax(:,:)
!!$  double precision a,b,m
!!$  double precision xTempS
!!$  double precision ax,bx
!!$  double precision, allocatable :: xIntScale(:),xS(:,:)
!!$
!!$  allocate(xIntScale(xNumPoints),xS(xDim,xDim))
!!$  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
!!$
!!$  S = 0.0d0
!!$
!!$  do kx = 1,xNumPoints-1
!!$     ax = xPoints(kx)
!!$     bx = xPoints(kx+1)
!!$     xIntScale(kx) = 0.5d0*(bx-ax)
!!$  enddo
!!$
!!$  !      do ix=1,xNumPoints+2*Order
!!$  !         print*, ix, xBounds(ix)
!!$  !      enddo
!!$
!!$  do ix = 1,xDim
!!$     do ixp = 1,xDim
!!$        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
!!$        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
!!$     enddo
!!$  enddo
!!$
!!$  do ix = 1,xDim
!!$     do ixp = max(1,ix-Order),min(xDim,ix+Order)
!!$        xS(ixp,ix) = 0.0d0
!!$        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
!!$           xTempS = 0.0d0
!!$           do lx = 1,LegPoints
!!$              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
!!$              b = a*u(lx,kx,ixp)
!!$              xTempS = xTempS + b
!!$           enddo
!!$           xS(ixp,ix) = xS(ixp,ix) + xTempS
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  do ix = 1,xDim
!!$     Row=ix
!!$     do ixp = max(1,ix-Order),min(xDim,ix+Order)
!!$        Col = ixp
!!$        if (Col .ge. Row) then
!!$           NewRow = HalfBandWidth+1+Row-Col
!!$           S(NewRow,Col) = xS(ixp,ix)
!!$           write(6,*) ix,ixp,S(NewRow,Col)
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$  deallocate(xIntScale,xS)
!!$  deallocate(kxMin,kxMax)
!!$
!!$  return
!!$end subroutine CalcOverlap
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine CalcHamiltonian(alpha,R,mu,Order,xPoints,LegPoints,xLeg,wLeg,xDim,&
!!$     xNumPoints,u,uxx,xBounds,HalfBandWidth,H)
!!$  implicit none
!!$  integer Order,LegPoints,xDim,xNumPoints,xBounds(*),HalfBandWidth
!!$  double precision alpha,R,mu
!!$  double precision xPoints(xNumPoints),xLeg(LegPoints),wLeg(LegPoints)
!!$  double precision H(HalfBandWidth+1,xDim)
!!$  double precision u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)
!!$
!!$  integer ix,ixp,kx,lx
!!$  integer i1,i1p
!!$  integer Row,NewRow,Col
!!$  integer, allocatable :: kxMin(:,:),kxMax(:,:)
!!$  double precision a,b,m,Pi
!!$  double precision Rall,r12,r12a,r23,r23a,r23b,r23c,r13,r13a,r13b,r13c,r14,r24,r34
!!$  double precision u1,sys_ss_pot,V12,V23,V31
!!$  double precision VInt,VTempInt,potvalue, xTempV
!!$  
!!$  double precision x,ax,bx,xScaledZero,xTempT,xTempS,xInt
!!$  double precision, allocatable :: Pot(:,:)
!!$  double precision, allocatable :: xIntScale(:),xT(:,:),xV(:,:)
!!$  double precision, allocatable :: cosx0(:,:),cosxp(:,:),cosxm(:,:),cosx(:,:),sinx(:,:)
!!$
!!$  double precision mu12,r0diatom,dDiatom
!!$  !      common/MassInfo/mu12,r0diatom,dDiatom
!!$
!!$  allocate(xIntScale(xNumPoints),xT(xDim,xDim),xV(xDim,xDim))
!!$  allocate(cosx0(LegPoints,xNumPoints),cosxp(LegPoints,xNumPoints),cosxm(LegPoints,xNumPoints))
!!$  allocate(cosx(LegPoints,xNumPoints))
!!$  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
!!$  allocate(Pot(LegPoints,xNumPoints))
!!$
!!$  Pi = 3.1415926535897932385d0
!!$
!!$  m = -1.0d0/(2.0d0*mu*R*R)
!!$
!!$  do kx = 1,xNumPoints-1
!!$
!!$     ax = xPoints(kx)
!!$     x = xPoints(kx+1)
!!$     xIntScale(kx) = 0.5d0*(bx-ax)
!!$     xScaledZero = 0.5d0*(bx+ax)
!!$
!!$     do lx = 1,LegPoints
!!$        x = xIntScale(kx)*xLeg(lx)+xScaledZero
!!$        cosx(lx,kx) = dcos(x)
!!$        cosxp(lx,kx) = dcos(x+Pi/3.0d0)
!!$        cosxm(lx,kx) = dcos(x-Pi/3.0d0)
!!$     enddo
!!$  enddo
!!$
!!$  do ix = 1,xDim
!!$     do ixp = 1,xDim
!!$        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
!!$        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
!!$     enddo
!!$  enddo
!!$
!!$  do kx = 1,xNumPoints-1
!!$     do lx = 1,LegPoints
!!$!        r12 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosx(lx,kx))
!!$!        r13 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxm(lx,kx))
!!$!        r23 = dsqrt(2.d0/dsqrt(3.0d0))*R*dabs(cosxp(lx,kx))
!!$        !        call  sumpairwisepot(r12, r13, r23, V2Depth, PotRange, potvalue)
!!$        potvalue = 0.5d0*mu*(R*cosx(lx,kx))**2
!!$        Pot(lx,kx) = 0d0!alpha*potvalue
!!$     enddo
!!$  enddo
!!$
!!$  xT = 0.0d0
!!$  xV = 0.0d0
!!$  
!!$  do ix = 1,xDim
!!$     do ixp = max(1,ix-Order),min(xDim,ix+Order)
!!$        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
!!$           xTempT = 0.0d0
!!$           xTempV = 0.0d0
!!$           do lx = 1,LegPoints
!!$              a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
!!$              xTempT = xTempT + a*uxx(lx,kx,ixp)
!!$              xTempV = xTempV + a*(Pot(lx,kx))*u(lx,kx,ixp)
!!$           enddo
!!$           xT(ix,ixp) = xT(ix,ixp) + xTempT
!!$           xV(ix,ixp) = xV(ix,ixp) + xTempV
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  H = 0.0d0      
!!$  do ix = 1,xDim
!!$     Row=ix
!!$     do ixp = max(1,ix-Order),min(xDim,ix+Order)
!!$        Col = ixp
!!$        if (Col .ge. Row) then
!!$           NewRow = HalfBandWidth+1+Row-Col
!!$           H(NewRow,Col) = (m*xT(ix,ixp)+xV(ix,ixp))
!!$!           write(25,*) ix,ixp,H(NewRow,Col)
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$  deallocate(Pot)
!!$  deallocate(xIntScale,xT,xV)
!!$  deallocate(cosx0,cosxp,cosxm,cosx)
!!$  deallocate(kxMin,kxMax)
!!$!  stop
!!$
!!$  return
!!$end subroutine CalcHamiltonian
!!$
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine CalcPMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P)
!!$  implicit none
!!$  integer NumStates,HalfBandWidth,MatrixDim
!!$  double precision RDelt
!!$  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
!!$  double precision S(HalfBandWidth+1,MatrixDim)
!!$  double precision P(NumStates,NumStates)
!!$
!!$  integer i,j,k
!!$  double precision a,ddot
!!$  double precision, allocatable :: TempPsi1(:),TempPsi2(:)
!!$
!!$  allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))
!!$
!!$  a = 0.5d0/RDelt
!!$
!!$  do j = 1,NumStates
!!$     do k = 1,MatrixDim
!!$        TempPsi1(k) = rPsi(k,j)-lPsi(k,j)
!!$     enddo
!!$     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
!!$     do i = 1,NumStates
!!$        P(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
!!$     enddo
!!$  enddo
!!$
!!$  deallocate(TempPsi1,TempPsi2)
!!$
!!$  return
!!$end subroutine CalcPMatrix
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine CalcQMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,Q)
!!$  implicit none
!!$  integer NumStates,HalfBandWidth,MatrixDim
!!$  double precision RDelt
!!$  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
!!$  double precision S(HalfBandWidth+1,MatrixDim)
!!$  double precision Q(NumStates,NumStates)
!!$
!!$  integer i,j,k
!!$  double precision a,ddot
!!$  double precision, allocatable :: TempPsi1(:),TempPsi2(:)
!!$
!!$  allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))
!!$
!!$  a = 1.0d0/(RDelt**2)
!!$
!!$  do j = 1,NumStates
!!$     do k = 1,MatrixDim
!!$        TempPsi1(k) = lPsi(k,j)+rPsi(k,j)-2.0d0*mPsi(k,j)
!!$     enddo
!!$     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
!!$     do i = 1,NumStates
!!$        Q(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
!!$     enddo
!!$  enddo
!!$
!!$  deallocate(TempPsi1,TempPsi2)
!!$
!!$  return
!!$end subroutine CalcQMatrix
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!!$subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,mPsi,rPsi)
!!$  implicit none
!!$  integer NumStates,HalfBandWidth,MatrixDim,ncv
!!$  double precision S(HalfBandWidth+1,MatrixDim),Psi(MatrixDim,ncv)
!!$  double precision mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv)
!!$
!!$  integer i,j
!!$  double precision Phase,ddot
!!$  double precision, allocatable :: TempPsi(:)
!!$
!!$  allocate(TempPsi(MatrixDim))
!!$
!!$  do i = 1,NumStates
!!$     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rPsi(1,i),1,0.0d0,TempPsi,1)
!!$     Phase = ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
!!$     if (Phase .lt. 0.0d0) then
!!$        do j = 1,MatrixDim
!!$           rPsi(j,i) = -rPsi(j,i)
!!$        enddo
!!$     endif
!!$  enddo
!!$
!!$  deallocate(TempPsi)
!!$
!!$  return
!!$end subroutine FixPhase
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine GridMakerHHL(mu,mu12,mu123,phi23,R,r0,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
!!$  implicit none
!!$  integer xNumPoints,CalcNewBasisFunc
!!$  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
!!$  double precision mu12,mu123,phi23
!!$  integer i,j,k,OPGRID
!!$  double precision Pi
!!$  double precision r0New
!!$  double precision xRswitch
!!$  double precision xDelt,x0,x1,x2,x3,x4,deltax
!!$
!!$
!!$  Pi = 3.1415926535897932385d0
!!$  r0New=r0*2.0d0
!!$  deltax = r0/R
!!$
!!$  xRswitch = 20.0d0*r0New/Pi
!!$  OPGRID=1
!!$  write (6,*) 'xMin = ', xMin, 'xMax = ', xMax
!!$  if((OPGRID.eq.1).and.(R.gt.xRswitch)) then
!!$     print*, 'R>xRswitch!! using modified grid!!'
!!$     x0 = xMin
!!$     x1 = phi23 - deltax  
!!$     x2 = phi23 + deltax
!!$     x3 = xMax-deltax
!!$     x4 = xMax
!!$     k = 1
!!$     xDelt = (x1-x0)/dfloat(xNumPoints/4)
!!$     do i = 1,xNumPoints/4
!!$        xPoints(k) = (i-1)*xDelt + x0
!!$        !     write(6,*) k, xPoints(k)
!!$        k = k + 1
!!$     enddo
!!$     xDelt = (x2-x1)/dfloat(xNumPoints/4)
!!$     do i = 1,xNumPoints/4
!!$        xPoints(k) = (i-1)*xDelt + x1
!!$        !     write(6,*) k, xPoints(k)
!!$        k = k + 1
!!$     enddo
!!$     xDelt = (x3-x2)/dfloat(xNumPoints/4)
!!$     do i = 1,xNumPoints/4
!!$        xPoints(k) = (i-1)*xDelt + x2
!!$        !     write(6,*) k, xPoints(k)
!!$        k = k + 1
!!$     enddo
!!$     xDelt = (x4-x3)/dfloat(xNumPoints/4-1)
!!$     do i = 1, xNumPoints/4
!!$        xPoints(k) = (i-1)*xDelt + x3
!!$        !     write(6,*) k, xPoints(k)
!!$        k = k + 1
!!$     enddo
!!$
!!$     !     FOR SMALL R, USE A LINEAR GRID
!!$  else
!!$     k = 1
!!$     xDelt = (xMax-xMin)/dfloat(xNumPoints-1)
!!$     do i = 1,xNumPoints
!!$        xPoints(k) = (i-1)*xDelt + x0
!!$        k = k + 1
!!$     enddo
!!$  endif
!!$
!!$  !     Smooth Grid 
!!$
!!$  write(20,*) 1, xPoints(1)
!!$  do i = 2, xNumPoints-1
!!$     xPoints(i)=(xPoints(i-1)+2.d0*xPoints(i)+xPoints(i+1))/4.d0
!!$     write(20,*) i, xPoints(i)
!!$  enddo
!!$  write(20,*) xNumPoints, xPoints(xNumPoints)
!!$  write(20,*) ' ' 
!!$  !      write(96,15) (xPoints(k),k=1,xNumPoints)
!!$
!!$15 format(6(1x,1pd12.5))
!!$
!!$
!!$
!!$  return
!!$end subroutine GridMakerHHL
!!$!      cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine GridMaker222(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
!!$  implicit none
!!$  integer xNumPoints,yNumPoints
!!$  double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)
!!$
!!$  integer i,j,k
!!$  double precision Pi
!!$  double precision r0New
!!$  double precision xRswitch,yRswitch
!!$  double precision xDelt,x0,x1,x2
!!$  double precision yDelt,y0,y1,y2
!!$
!!$  Pi = 3.1415926535897932385d0
!!$
!!$  r0New = 3.0d0*r0
!!$
!!$  xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0))))
!!$
!!$  if (R .gt. xRswitch) then
!!$     x0 = xMin
!!$     x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
!!$     x2 = xMax
!!$     k = 1
!!$     xDelt = (x1-x0)/dfloat(xNumPoints/2)
!!$     do i = 1,xNumPoints/2
!!$        xPoints(k) = (i-1)*xDelt + x0
!!$        k = k + 1
!!$     enddo
!!$     xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
!!$     do i = 1,xNumPoints/2
!!$        xPoints(k) = (i-1)*xDelt + x1
!!$        k = k + 1
!!$     enddo
!!$  else
!!$     x0 = xMin
!!$     x1 = xMax
!!$     k = 1
!!$     xDelt = (x1-x0)/dfloat(xNumPoints-1)
!!$     do i = 1,xNumPoints
!!$        xPoints(k) = (i-1)*xDelt + x0
!!$        k = k + 1
!!$     enddo
!!$  endif
!!$
!!$  yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))
!!$
!!$  if (R .gt. yRswitch) then
!!$     y0 = yMin
!!$     y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
!!$     y2 = yMax
!!$     k = 1
!!$     yDelt = (y1-y0)/dfloat(yNumPoints/2)
!!$     do i = 1,yNumPoints/2
!!$        yPoints(k) = (i-1)*yDelt + y0
!!$        k = k + 1
!!$     enddo
!!$     yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
!!$     do i = 1,yNumPoints/2
!!$        yPoints(k) = (i-1)*yDelt + y1
!!$        k = k + 1
!!$     enddo
!!$  else
!!$     y0 = yMin
!!$     y1 = yMax
!!$     k = 1
!!$     yDelt = (y1-y0)/dfloat(yNumPoints-1)
!!$     do i = 1,yNumPoints
!!$        yPoints(k) = (i-1)*yDelt + y0
!!$        k = k + 1
!!$     enddo
!!$  endif
!!$
!!$  return
!!$end subroutine GridMaker222
!!$
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!!$subroutine GridMaker111(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
!!$  implicit none
!!$  integer xNumPoints,yNumPoints
!!$  double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)
!!$
!!$
!!$  integer i,j,k
!!$  double precision Pi
!!$  double precision r0New
!!$  double precision xRswitch,yRswitch
!!$  double precision xDelt,x0,x1,x2
!!$  double precision yDelt,y0,y1,y2
!!$
!!$  Pi = 3.1415926535897932385d0
!!$
!!$
!!$  r0New = 3.0d0*r0
!!$
!!$  xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0)) ))
!!$
!!$
!!$  if (R .gt. xRswitch) then
!!$     x0 = xMin
!!$     x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
!!$     x2 = xMax
!!$     k = 1
!!$     xDelt = (x1-x0)/dfloat(xNumPoints/2)
!!$     do i = 1,xNumPoints/2
!!$        xPoints(k) = (i-1)*xDelt + x0
!!$        k = k + 1
!!$     enddo
!!$     xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
!!$     do i = 1,xNumPoints/2
!!$        xPoints(k) = (i-1)*xDelt + x1
!!$        k = k + 1
!!$     enddo
!!$  else
!!$     x0 = xMin
!!$     x1 = xMax
!!$     k = 1
!!$     xDelt = (x1-x0)/dfloat(xNumPoints-1)
!!$     do i = 1,xNumPoints
!!$        xPoints(k) = (i-1)*xDelt + x0
!!$        k = k + 1
!!$     enddo
!!$  endif
!!$
!!$  yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))
!!$
!!$  if (R .gt. yRswitch) then
!!$     y0 = yMin
!!$     y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
!!$     y2 = yMax
!!$     k = 1
!!$     yDelt = (y1-y0)/dfloat(yNumPoints/2)
!!$     do i = 1,yNumPoints/2
!!$        yPoints(k) = (i-1)*yDelt + y0
!!$        k = k + 1
!!$     enddo
!!$     yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
!!$     do i = 1,yNumPoints/2
!!$        yPoints(k) = (i-1)*yDelt + y1
!!$        k = k + 1
!!$     enddo
!!$  else
!!$     y0 = yMin
!!$     y1 = yMax
!!$     k = 1
!!$     yDelt = (y1-y0)/dfloat(yNumPoints-1)
!!$     do i = 1,yNumPoints
!!$        yPoints(k) = (i-1)*yDelt + y0
!!$        k = k + 1
!!$     enddo
!!$  endif
!!$
!!$  return
!!$end subroutine GridMaker111
!!$
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine GridMaker(mu,R,r0,xNumPoints,xMin,xMax,xPoints)
!!$  implicit none
!!$  integer xNumPoints
!!$  double precision mu,R,r0,xMin,xMax,xPoints(xNumPoints)
!!$  integer i,j,k,OPGRID
!!$  double precision Pi
!!$  double precision r0New
!!$  double precision xRswitch
!!$  double precision xDelt,x0,x1,x2
!!$
!!$  write(6,*) "In gridmaker..."
!!$  Pi = 3.1415926535897932385d0
!!$  write(6,*) xMin, xMax
!!$  x0 = xMin
!!$  x1 = xMax
!!$  write(96,*) 'x0,x1=',x0,x1
!!$  
!!$!  r0New=10.0d0*r0           !/2.0d0
!!$  !r0New=Pi/12.0*R
!!$ ! xRswitch = 12.0d0*r0New/Pi
!!$ ! OPGRID=0
!!$  !write(6,*) "entering conditional"
!!$  !if((OPGRID.eq.1).and.(R.gt.xRswitch)) then
!!$  !   print*, 'R>xRswitch!! using modified grid!!'
!!$  !   x0 = xMin
!!$  !   x1 = xMax - r0New/R  
!!$  !   x2 = xMax
!!$  !   k = 1
!!$  !   xDelt = (x1-x0)/dfloat(xNumPoints/2)
!!$  !   do i = 1,xNumPoints/2
!!$  !      xPoints(k) = (i-1)*xDelt + x0
!!$  !      !            print*, k, xPoints(k), xDelt
!!$  !      k = k + 1
!!$  !   enddo
!!$  !   xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
!!$  !   do i = 1,xNumPoints/2
!!$  !      xPoints(k) = (i-1)*xDelt + x1
!!$  !      !            print*, k, xPoints(k), xDelt
!!$  !      k = k + 1
!!$  !   enddo
!!$  !else
!!$  !   x0 = xMin
!!$  !   x1 = xMax
!!$  k = 1
!!$  xDelt = (x1-x0)/dble(xNumPoints-1)
!!$  write(6,*) "entering loop"
!!$  do i = 1,xNumPoints
!!$     xPoints(k) = (i-1)*xDelt + x0
!!$     print*, k, xPoints(k), xDelt
!!$     k = k + 1
!!$  enddo
!!$ ! endif
!!$
!!$  !      write(96,15) (xPoints(k),k=1,xNumPoints)
!!$15 format(6(1x,1pd12.5))
!!$
!!$
!!$
!!$  return
!!$end subroutine GridMaker
!!$!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine CalcCoupling(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P,Q,dP)
!!$  implicit none
!!$  integer NumStates,HalfBandWidth,MatrixDim
!!$  double precision RDelt
!!$  double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
!!$  double precision S(HalfBandWidth+1,MatrixDim),testorth
!!$  double precision P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates)
!!$
!!$  integer i,j,k
!!$  double precision aP,aQ,ddot
!!$  double precision, allocatable :: lDiffPsi(:),rDiffPsi(:),TempPsi(:),TempPsiB(:),rSumPsi(:)
!!$  double precision, allocatable :: TempmPsi(:)
!!$
!!$  allocate(lDiffPsi(MatrixDim),rDiffPsi(MatrixDim),TempPsi(MatrixDim),&
!!$       TempPsiB(MatrixDim),rSumPsi(MatrixDim))
!!$  allocate(TempmPsi(MatrixDim))
!!$
!!$  aP = 0.5d0/RDelt
!!$  aQ = aP*aP
!!$
!!$  do j = 1,NumStates
!!$     do k = 1,MatrixDim
!!$        rDiffPsi(k) = rPsi(k,j)-lPsi(k,j)
!!$        rSumPsi(k)  = lPsi(k,j)+mPsi(k,j)+rPsi(k,j)
!!$        !            rSumPsi(k)  = lPsi(k,j)-2.0d0*mPsi(k,j)+rPsi(k,j)
!!$        !            rSumPsi(k)  = lPsi(k,j)+rPsi(k,j)
!!$     enddo
!!$     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)
!!$     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rSumPsi,1,0.0d0,TempPsiB,1)
!!$     !         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,mPsi(1,j),1,0.0d0,TempmPsi,1)
!!$
!!$     do i = 1,NumStates
!!$
!!$        !            testorth=ddot(MatrixDim,mPsi(1,i),1,TempmPsi,1)
!!$        !            write(309,*) i,j, '   testorth=',testorth
!!$
!!$        P(i,j) = aP*ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
!!$        dP(i,j)= ddot(MatrixDim,mPsi(1,i),1,TempPsiB,1)
!!$
!!$        do k = 1,MatrixDim
!!$           lDiffPsi(k) = rPsi(k,i)-lPsi(k,i)
!!$        enddo
!!$        Q(i,j) = -aQ*ddot(MatrixDim,lDiffPsi,1,TempPsi,1)
!!$     enddo
!!$  enddo
!!$
!!$  do j=1,NumStates
!!$     do i=j,NumStates
!!$        dP(i,j)=2.d0*aQ*(dP(i,j)-dP(j,i))
!!$        dP(j,i)=-dP(i,j)
!!$     enddo
!!$  enddo
!!$
!!$  deallocate(lDiffPsi,rDiffPsi,TempPsi,rSumPsi,TempPsiB,TempmPsi)
!!$
!!$  return
!!$end subroutine CalcCoupling
!!$
!!$
!!$
!!$
!!$      
