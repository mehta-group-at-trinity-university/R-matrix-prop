!****************************************************************************************************
MODULE GlobalVars
  IMPLICIT NONE
  INTEGER NumParticles, NumChannels,  NumAllChan, Order
  INTEGER StartBC, EndBC, xTotNumPoints, NumBoxes,lmax
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION AlphaFactor ! This is the parameter that appears in the reduced wavefunction u(R) = R^(AlphaFactor) Psi(R)
  ! Typical choice is either AlphaFactor = 0 (reduced wavefunction = wavefunction), or AlphaFactor = (EffDim - 1)/2 (eliminates 1st derivative terms from KE)
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION reducedmass, xStart, xEnd, energy,kStart,kEnd
  DOUBLE PRECISION SpatialDim, EffDim
  DOUBLE PRECISION, ALLOCATABLE :: mass(:)
  CHARACTER*64 InputFile
  COMPLEX*16 II
  PARAMETER(II=(0.0d0,1.0d0))
  !****************************************************************************************************
CONTAINS
  SUBROUTINE ReadGlobal
    USE Quadrature
    IMPLICIT NONE
    INTEGER n
    ! Be sure to match the format of the input file to the format of the read statements below
    OPEN(unit=7,file=InputFile(1:INDEX(InputFile,' ')-1),action='read')
    READ(7,*)
    READ(7,*) NumParticles, NumChannels, SpatialDim, NumAllChan, Order
    ALLOCATE(mass(NumParticles))
    READ(7,*)
    READ(7,*)
    READ(7,*) (mass(n), n=1,NumParticles)
    READ(7,*)
    READ(7,*)
    READ(7,*) xStart, xEnd, xTotNumPoints, LegPoints
    READ(7,*)
    READ(7,*)
    READ(7,*) NumBoxes, StartBC, EndBC, kStart, kEnd
    READ(7,*)
    READ(7,*)
    READ(7,*) Energy

    lmax = 2*NumChannels
    CLOSE(unit=7)
    EffDim = NumParticles*SpatialDim - SpatialDim
    !AlphaFactor = 0d0
    AlphaFactor = (EffDim-1d0)/2d0

    IF (NumParticles.EQ.2) THEN
       reducedmass = mass(1)*mass(2)/(mass(1)+mass(2))
    ELSE
       WRITE(6,*) "Reduced mass not set. Must set reduced mass"
       STOP
    END IF


  END SUBROUTINE ReadGlobal
  !****************************************************************************************************
END MODULE GlobalVars
!****************************************************************************************************
!****************************************************************************************************
MODULE DataStructures
  IMPLICIT NONE
  !****************************************************************************************************
  TYPE BPData
     !This data type contains the array of basis functions and potential matrix
     INTEGER xNumPoints, xDim, Left, Right, Order, NumChannels, MatrixDim
     DOUBLE PRECISION, ALLOCATABLE :: u(:,:,:),ux(:,:,:),xPoints(:),x(:,:)
     DOUBLE PRECISION, ALLOCATABLE :: Pot(:,:,:,:),lam(:)
     DOUBLE PRECISION kl, kr, xl,xr
     INTEGER, ALLOCATABLE :: xBounds(:), kxMin(:,:), kxMax(:,:)
  END TYPE BPData
  !****************************************************************************************************
  TYPE GenEigVal
     INTEGER MatrixDim
     DOUBLE PRECISION, ALLOCATABLE :: Gam(:,:), Lam(:,:), evec(:,:), eval(:)
  END TYPE GenEigVal
  !****************************************************************************************************
  TYPE BoxData
     INTEGER NumOpenL, NumOpenR, betaMax
     DOUBLE PRECISION xL, xR
     DOUBLE PRECISION, ALLOCATABLE :: Z(:,:), ZP(:,:), NBeta(:), Norm(:,:)
     DOUBLE PRECISION, ALLOCATABLE :: b(:), bf(:), Zf(:,:), Zfp(:,:)
     !double precision, allocatable :: RF(:,:),K(:,:)

  END TYPE BoxData
  TYPE ScatData
     DOUBLE PRECISION, ALLOCATABLE :: K(:,:), R(:,:), sigma(:,:)
     complex*16, ALLOCATABLE :: f(:,:), S(:,:), T(:,:)

  END TYPE ScatData

CONTAINS
  !****************************************************************************************************
  SUBROUTINE AllocateScat(SD,N)
    IMPLICIT NONE
    TYPE(ScatData) SD
    INTEGER N
    ALLOCATE(SD%K(N,N),SD%R(N,N),SD%sigma(N,N))
    ALLOCATE(SD%f(N,N),SD%S(N,N),SD%T(N,N))
    SD%K=0d0
    SD%R=0d0
    SD%f=0d0
    SD%sigma=0d0
    SD%S=(0d0,0d0)
    SD%T=(0d0,0d0)
  END SUBROUTINE AllocateScat
  !****************************************************************************************************
  SUBROUTINE DeAllocateScat(SD)
    IMPLICIT NONE
    TYPE(ScatData) SD

    DEALLOCATE(SD%K,SD%R,SD%f,SD%sigma,SD%S,SD%T)

  END SUBROUTINE DeAllocateScat
  !****************************************************************************************************
  SUBROUTINE AllocateEIG(EIG)
    IMPLICIT NONE
    TYPE(GenEigVal) :: EIG
    ALLOCATE(EIG%Gam(EIG%MatrixDim,EIG%MatrixDim),EIG%Lam(EIG%MatrixDim,EIG%MatrixDim))
    ALLOCATE(EIG%eval(EIG%MatrixDim),EIG%evec(EIG%MatrixDim,EIG%MatrixDim))
    EIG%Gam=0d0
    EIG%Lam=0d0
    EIG%eval=0d0
    EIG%Evec=0d0

  END SUBROUTINE AllocateEIG
  !****************************************************************************************************
  SUBROUTINE DeAllocateEIG(EIG)
    IMPLICIT NONE
    TYPE(GenEigVal) :: EIG
    DEALLOCATE(EIG%Gam,EIG%Lam,EIG%eval,EIG%evec)

  END SUBROUTINE DeAllocateEIG
  !****************************************************************************************************
  SUBROUTINE AllocateBox(B)
    IMPLICIT NONE
    TYPE(BoxData) :: B

    B%betaMax = B%NumOpenR + B%NumOpenL

    ALLOCATE(B%Z(B%betaMax,B%betaMax),B%ZP(B%betaMax,B%betaMax))
    ALLOCATE(B%b(B%betaMax))
    ALLOCATE(B%NBeta(B%betaMax))
    ALLOCATE(B%Norm(B%betaMax,B%betaMax))
    ALLOCATE(B%Zf(B%NumOpenR,B%NumOpenR),B%Zfp(B%NumOpenR,B%NumOpenR))
    ALLOCATE(B%bf(B%NumOpenR))

  END SUBROUTINE AllocateBox
  !****************************************************************************************************
  SUBROUTINE InitZeroBox(B)
    IMPLICIT NONE
    TYPE(BoxData) :: B
    B%Z=0d0
    B%ZP=0d0
    B%b=0d0
    B%NBeta=0d0
    B%Norm=0d0
    B%Zf=0d0
    B%Zfp=0d0
    B%bf=0d0

  END SUBROUTINE InitZeroBox
  !****************************************************************************************************
  SUBROUTINE DeAllocateBox(B)
    IMPLICIT NONE
    TYPE(BoxData) :: B

    DEALLOCATE(B%Z,B%ZP)
    DEALLOCATE(B%b)
    DEALLOCATE(B%NBeta,B%Norm,B%Zf,B%bf,B%Zfp)

  END SUBROUTINE DeAllocateBox
  !****************************************************************************************************
  SUBROUTINE AllocateBPD(BPD)
    USE Quadrature
    IMPLICIT NONE
    TYPE(BPData) BPD
    INTEGER xDimMin
    !Determine the basis dimension for these boundary conditions
    xDimMin=BPD%xNumPoints+BPD%Order-3
    BPD%xDim=xDimMin
    IF (BPD%Left .EQ. 2) BPD%xDim = BPD%xDim + 1
    IF (BPD%Right .EQ. 2) BPD%xDim = BPD%xDim + 1
    BPD%MatrixDim = BPD%NumChannels*BPD%xDim

    ALLOCATE(BPD%xPoints(BPD%xNumPoints))
    ALLOCATE(BPD%xBounds(BPD%xNumPoints + 2*BPD%Order))   ! allocate xBounds before calling CalcBasisFuncs
    ALLOCATE(BPD%x(LegPoints,BPD%xNumPoints-1))

    ALLOCATE(BPD%u(LegPoints,BPD%xNumPoints,BPD%xDim),BPD%ux(LegPoints,BPD%xNumPoints,BPD%xDim)) ! allocate memory for basis functions
    ALLOCATE(BPD%kxMin(BPD%xDim,BPD%xDim),BPD%kxMax(BPD%xDim,BPD%xDim)) !
    ALLOCATE(BPD%Pot(BPD%NumChannels,BPD%NumChannels,LegPoints,BPD%xNumPoints-1))

    ALLOCATE(BPD%lam(BPD%NumChannels))

    BPD%lam=0.d0

  END SUBROUTINE AllocateBPD
  !****************************************************************************************************
  SUBROUTINE DeAllocateBPD(BPD)
    IMPLICIT NONE
    TYPE(BPData) BPD
    DEALLOCATE(BPD%xPoints,BPD%xBounds,BPD%x,BPD%u,BPD%ux,BPD%kxMin,BPD%kxMax,BPD%Pot,BPD%lam)
  END SUBROUTINE DeAllocateBPD
  !****************************************************************************************************
  SUBROUTINE MakeBasis(BPD)
    USE Quadrature
    IMPLICIT NONE

    TYPE(BPData) BPD
    INTEGER ix, ixp

    CALL CalcBasisFuncsBP(BPD%Left,BPD%Right,BPD%kl,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         BPD%xDim,BPD%xBounds,BPD%xNumPoints,0,BPD%u)
    CALL CalcBasisFuncsBP(BPD%Left,BPD%Right,BPD%kl,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         BPD%xDim,BPD%xBounds,BPD%xNumPoints,1,BPD%ux)

    ! Determine the bounds for integration of matrix elements
    DO ix = 1,BPD%xDim
       DO ixp = 1,BPD%xDim
          BPD%kxMin(ixp,ix) = MAX(BPD%xBounds(ix),BPD%xBounds(ixp))
          BPD%kxMax(ixp,ix) = MIN(BPD%xBounds(ix+BPD%Order+1),BPD%xBounds(ixp+BPD%Order+1))-1 !
       ENDDO
    ENDDO

  END SUBROUTINE MakeBasis
  !****************************************************************************************************
  SUBROUTINE checkpot(BPD,file)
    USE Quadrature
    IMPLICIT NONE
    INTEGER file,mch,nch,lx,kx
    TYPE(BPData) BPD

    DO mch = 1,BPD%NumChannels
       DO nch = 1,mch
          DO kx = 1, BPD%xNumPoints-1
             DO lx = 1,LegPoints
                WRITE(file,*) BPD%x(lx,kx), BPD%Pot(mch,nch,lx,kx)
             ENDDO
          ENDDO
          WRITE(file,*)
       ENDDO
       WRITE(file,*)
    ENDDO

  END SUBROUTINE checkpot

END MODULE DataStructures
!****************************************************************************************************
SUBROUTINE CalcGamLam(BPD,EIG)
  USE DataStructures
  USE Quadrature
  USE GlobalVars
  IMPLICIT NONE
  TYPE(BPData), INTENT(in) :: BPD
  TYPE(GenEigVal) :: EIG
  INTEGER ix,ixp,lx,kx,mch,nch
  DOUBLE PRECISION a, ax, bx, xIntScale(BPD%xNumPoints), TempG, xScaledZero
  EIG%Gam=0d0
  EIG%Lam=0d0
  EIG%eval=0d0
  EIG%evec=0d0
  !Calculate the Gamma Matrix elements
  DO ix = 1,BPD%xDim
     DO ixp = MAX(1,ix-BPD%Order),MIN(BPD%xDim,ix+BPD%Order)
        DO mch = 1,BPD%NumChannels
           ! Do the channel-diagonal part first:
           DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
              ax = BPD%xPoints(kx)
              bx = BPD%xPoints(kx+1)
              xIntScale(kx) = 0.5d0*(bx-ax)
              xScaledZero = 0.5d0*(bx+ax)
              TempG = 0.0d0
              DO lx = 1,LegPoints
                 a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                 ! The KE matrix elements
                 TempG = TempG - a*(BPD%ux(lx,kx,ix)*BPD%ux(lx,kx,ixp))
                 ! The diagonal "overlap"*energy and additional term from reducing the wavefunction
                 TempG = TempG + a*BPD%u(lx,kx,ix)*2*reducedmass*(Energy - (BPD%Pot(mch,mch,lx,kx) - &
                      AlphaFactor*(AlphaFactor-EffDim+2)/(2d0*reducedmass*BPD%x(lx,kx)**2)))*BPD%u(lx,kx,ixp)
              ENDDO
              EIG%Gam((mch-1)*BPD%xDim+ix,(mch-1)*BPD%xDim+ixp) = EIG%Gam((mch-1)*BPD%xDim+ix,(mch-1)*BPD%xDim+ixp) + TempG ! place values into Gamma0
           ENDDO

           ! Now do the off-diagonal parts
           DO nch = 1, mch-1
              DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                 ax = BPD%xPoints(kx)
                 bx = BPD%xPoints(kx+1)
                 xIntScale(kx) = 0.5d0*(bx-ax)
                 xScaledZero = 0.5d0*(bx+ax)
                 TempG = 0.0d0
                 DO lx = 1,LegPoints
                    a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                    ! The potential matrix elements calculated here
                    TempG = TempG - a*BPD%u(lx,kx,ix)*2.0d0*reducedmass*BPD%Pot(mch,nch,lx,kx)*BPD%u(lx,kx,ixp)
                 ENDDO
                 EIG%Gam((mch-1)*BPD%xDim+ix,(nch-1)*BPD%xDim+ixp) = EIG%Gam((mch-1)*BPD%xDim +ix, (nch-1)*BPD%xDim +ixp) + TempG ! place values into EIG%Gam
              ENDDO
              EIG%Gam((nch-1)*BPD%xDim+ix,(mch-1)*BPD%xDim+ixp) = EIG%Gam((mch-1)*BPD%xDim+ix,(nch-1)*BPD%xDim+ixp)  ! fill in symmetric element mch <--> nch
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !Calculate the Lambda matrix elements
  EIG%Lam = 0.d0

  IF(BPD%Right.EQ.2) THEN
     DO mch = 1, BPD%NumChannels
        EIG%Lam( (mch-1)*BPD%xDim + BPD%xDim, (mch-1)*BPD%xDim + BPD%xDim ) = BPD%xr**(EffDim-1d0-2d0*AlphaFactor)
     ENDDO
  ENDIF
  IF(BPD%Left.EQ.2) THEN
     DO mch = 1, NumChannels
        EIG%Lam((mch-1)*BPD%xDim + 1,(mch-1)*BPD%xDim + 1) = BPD%xl**(EffDim-1d0-2d0*AlphaFactor)
     ENDDO
  ENDIF

END SUBROUTINE CalcGamLam
!****************************************************************************************************
MODULE Scattering
  USE DataStructures
CONTAINS
  SUBROUTINE CalcK(B,BPD,SD,mu,d,alpha,EE,Eth)
    IMPLICIT NONE
    TYPE(BoxData), INTENT(in) :: B
    TYPE(BPData), INTENT(in) :: BPD
    TYPE(ScatData) :: SD
    DOUBLE PRECISION mu, EE,rm,d,alpha
    DOUBLE PRECISION, ALLOCATABLE :: s(:),c(:),sp(:),cp(:)
    DOUBLE PRECISION, ALLOCATABLE :: Imat(:,:),Jmat(:,:)
    DOUBLE PRECISION rhypj,rhypy,rhypjp,rhypyp,Pi
    DOUBLE PRECISION k(BPD%NumChannels),Eth(BPD%NumChannels)
    complex*16, allocatable :: tmp(:,:),Identity(:,:)
    complex*16  II
    INTEGER i,j,no,nw,nc,beta

    II=(0d0,1d0)
    Pi=dacos(-1d0)
    rm=BPD%xr
    no=0
    nw=0
    nc=0

    DO i = 1,BPD%NumChannels
       IF (EE.GE.Eth(i)) THEN
          k(i) = dsqrt(2d0*mu*(EE-Eth(i))) ! k is real
          no=no+1
       ELSE
          k(i) = dsqrt(2d0*mu*(Eth(i)-EE)) ! k->kappa, kappa is real
          IF( (k(i)*rm).LT.10d0) nw = nw+1
          IF( (k(i)*rm).GE.10d0) nc = nc+1
       ENDIF
    ENDDO
    !      write(6,*) "no = ", no
    IF((no+nw+nc).NE.BPD%NumChannels) THEN
       WRITE(6,*) "Channel miscount in calcK"
       STOP
    ENDIF

    ALLOCATE(s(no),c(no),Imat(no,no),Jmat(no,no),tmp(no,no))
    ALLOCATE(sp(no),cp(no))
    allocate(Identity(no,no))
    Identity = 0d0;

    DO i = 1,no
      Identity(i,i) = 1d0
      CALL hyperrjry(INT(d),alpha,BPD%lam(i),k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
      s(i) = rhypj/dsqrt(Pi*k(i))
      c(i) = -rhypy/dsqrt(Pi*k(i))
      sp(i) = dsqrt(k(i)/Pi)*rhypjp
      cp(i) = -dsqrt(k(i)/Pi)*rhypyp
    ENDDO
    Imat=0d0
    Jmat=0d0
    DO i=1,no
       DO beta=1,no
          Imat(i,beta) = (B%Zf(i,beta)*cp(i) - B%Zfp(i,beta)*c(i))/(s(i)*cp(i)-c(i)*sp(i))
          Jmat(i,beta) = (B%Zf(i,beta)*sp(i) - B%Zfp(i,beta)*s(i))/(c(i)*sp(i)-s(i)*cp(i))
       ENDDO
    ENDDO

    CALL SqrMatInv(Imat, no)
    SD%K = MATMUL(Jmat,Imat)

    SD%R=0.0d0
    DO i=1,B%NumOpenR
       DO j=1,B%NumOpenR
          DO beta = 1, B%NumOpenR
             SD%R(i,j) = SD%R(i,j) - B%Zf(i,beta)*B%Zf(j,beta)/B%bf(beta) !
          ENDDO
       ENDDO
    ENDDO

    !call sqrmatinv(B%Zfp,B%NumOpenR)  !This gives the same result as the code segment  executed above
    !SD%R = matmul(B%Zf,B%Zfp)

    !tmp = Identity - II*SD%K
    !write(6,*) "Identity matrix:"
    !call printmatrix(Identity,no,no,6)
    !write(6,*) "K matrix:"
    !call printmatrix(SD%K,no,no,6)
    !write(6,*) "real part of 1 - II*K"
    !call printmatrix(realpart(tmp),no,no,6)
    !write(6,*) "Imaginary part of 1 - II*K"
    !call printmatrix(aimag(tmp),no,no,6)
    !SD%S = Identity + II*SD%K
    !call CompSqrMatInv(tmp,no)
    !SD%S = MATMUL(SD%S,tmp)
    !SD%T = -II*0.5d0*(SD%S-Identity)

    DEALLOCATE(s,c,Imat,Jmat,sp,cp)
  END SUBROUTINE CalcK

END MODULE Scattering

!****************************************************************************************************
SUBROUTINE makeEgrid(Egrid,NumE,E1,E2,scale)
  DOUBLE PRECISION Egrid(NumE)
  DOUBLE PRECISION E1,E2,LE1,LE2,DE,DLE
  INTEGER NumE, iE
  CHARACTER(LEN=*), INTENT(IN) :: scale
  !--------------------------------------------
  ! Linear grid:
  !--------------------------------------------
  IF(scale.EQ."linear") THEN
     DE=(E2-E1)/DBLE(NumE-1)
     DO iE=1,NumE
        Egrid(iE) = E1 + (iE-1)*DE
     ENDDO
  ENDIF
  !--------------------------------------------
  ! Log grid:
  !--------------------------------------------
  IF(scale.EQ."log") THEN
     LE1=dlog(E1)
     LE2=dlog(E2)
     LDE=(LE2-LE1)/DBLE(NumE-1)
     DO iE=1,NumE
        Egrid(iE) = dexp(LE1 + (iE-1)*LDE)
     ENDDO
  ENDIF

END SUBROUTINE makeEgrid
!****************************************************************************************************
PROGRAM main
  USE DataStructures
  USE GlobalVars
  USE BalujaParameters
  USE Quadrature
  USE MorsePotential
  USE scattering
  USE DipoleDipole

  IMPLICIT NONE
  TYPE(BPData) BPD,BPD1,BPD0
  TYPE(GenEigVal) EIG,EIG1
  TYPE(BoxData) BA, BB, Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  TYPE(Morse) :: M
  TYPE(DPData) :: DP
  DOUBLE PRECISION, ALLOCATABLE :: evec(:,:), eval(:)!, temp0(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Egrid(:), xprim(:)

  DOUBLE PRECISION xDelt
  INTEGER NumE, iE, beta, i, iBox,lx,kx

  !----------------------------------------------------------------------
  ! Read information from the input file
  !----------------------------------------------------------------------
  InputFile = 'RMATPROP.inp'
  CALL ReadGlobal()
  !--------------------------------------------------
  ! Read in the quadrature data stored in Quadrature module.
  !--------------------------------------------------
  LegendreFile = 'Legendre.dat'
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)
  !---------------------------------------------------------------------
  ! allocate the data for the Boxes
  !---------------------------------------------------------------------
  xDelt=(xEnd-xStart)/DBLE(NumBoxes)
  ALLOCATE(Boxes(NumBoxes))
  Boxes(1)%NumOpenL = 0
  Boxes(1)%NumOpenR = NumChannels
  Boxes(1)%betaMax = Boxes(1)%NumOpenL + Boxes(1)%NumOpenR
  Boxes(1)%xl=xStart
  Boxes(1)%xr=xDelt
  CALL AllocateBox(Boxes(1))
  CALL InitZeroBox(Boxes(1))
  WRITE(6,*) "Box",1,Boxes(1)%xl,Boxes(1)%xr
  CALL printmatrix(Boxes(1)%Zf,Boxes(1)%NumOpenR,Boxes(1)%NumOpenR,6)
  DO i = 2,NumBoxes
     Boxes(i)%NumOpenL = NumChannels
     Boxes(i)%NumOpenR = NumChannels
     Boxes(i)%betaMax = Boxes(1)%NumOpenL + Boxes(1)%NumOpenR
     Boxes(i)%xl=Boxes(i-1)%xr
     Boxes(i)%xr=Boxes(i-1)%xr+xDelt
     CALL AllocateBox(Boxes(i))
     CALL InitZeroBox(Boxes(i))
     WRITE(6,*) "Box",i,Boxes(i)%xl,Boxes(i)%xr
     CALL printmatrix(Boxes(i)%Zf,Boxes(i)%NumOpenR,Boxes(i)%NumOpenR,6)
  ENDDO

  DP%lmax=lmax
  CALL AllocateDP(DP)

  CALL MakeDipoleDipoleCouplingMatrix(DP)
!  WRITE(6,*) "cllp for m=0:"
!  CALL printmatrix(DP%cllp(0:DP%lmax,0:DP%lmax,0),DP%lmax,DP%lmax,6)
!  STOP

  !-------------------------------------------------------------------
  ! Intitializes some BPD1 variables to the input values.
  ! This is the basis set used at the left edge (r=0) in Boxes(1)
  !-------------------------------------------------------------------
  BPD1%NumChannels = NumChannels
  BPD1%Order = Order
  BPD1%Left = 0
  BPD1%Right = 2
  BPD1%xNumPoints = xTotNumPoints
  BPD1%kl = kStart ! only relevant for Left = 3. This is the normal log-derivative at BPD%xl
  BPD1%kr = kEnd   ! only relevant for Right = 3. This is the normal log-derivative at BPD%xr

  !---------------------------------------------------------------------
  ! Intitializes some BPD0 variables to the input values.
  ! This is the "primitive" basis set defined over range [0,1] (see below.)
  !---------------------------------------------------------------------
  BPD0%NumChannels = NumChannels
  BPD0%Order = Order
  BPD0%Left = 2
  BPD0%Right = 2
  BPD0%xNumPoints = xTotNumPoints
  BPD0%kl = kStart ! only relevant for Left = 3. This is the normal log-derivative at BPD%xl
  BPD0%kr = kEnd ! only relevant for Right = 3. This is the normal log-derivative at BPD%xr

  !---------------------------------------------------------------------
  ! Copy over the properties of the primitive set to the workhorse set BPD that will
  ! actually be used for propagation.
  !---------------------------------------------------------------------
  BPD=BPD0
  ! Allocate memory for the basis sets and potential arrays
  CALL AllocateBPD(BPD1)
  CALL AllocateBPD(BPD)
  CALL AllocateBPD(BPD0)

  !---------------------------------------------------------------------
  ! Make and fill the arrays for the basis used in Boxes(1)
  !---------------------------------------------------------------------
  BPD1%xl=Boxes(1)%xl
  BPD1%xr=Boxes(1)%xr
  CALL GridMakerLinear(BPD1%xNumPoints,BPD1%xl,BPD1%xr,BPD1%xPoints)
  CALL Makebasis(BPD1)

  !---------------------------------------------------------------------
  ! Make and fill the arrays for the primitive set
  !---------------------------------------------------------------------
  ALLOCATE(xprim(xTotNumPoints))
  CALL GridMakerLinear(BPD0%xNumPoints,0d0,1d0,xprim)
  BPD0%xPoints=xprim
  CALL MakeBasis(BPD0)
  !Copy over the primitive set to the workhose set
  BPD=BPD0
  ! Initialize the potential data
  CALL InitMorse(M)

  ! make the energy grid
  NumE=2000
  ALLOCATE(Egrid(NumE))
  CALL makeEgrid(Egrid,NumE,M%Eth(2)+0.01d0,M%Eth(3)-0.01d0,"linear")

  ! Fill the arrays for the potential matrix for the first box
  CALL SetMorsePotential(BPD1,M)

  ! Allocate memory for the eigenvalue problem
  EIG1%MatrixDim=BPD1%MatrixDim
  CALL AllocateEIG(EIG1)
  WRITE(6,*) "EIG1%MatrixDim = ",EIG1%MatrixDim

  EIG%MatrixDim=BPD%MatrixDim
  CALL AllocateEIG(EIG)
  WRITE(6,*) "EIG%MatrixDim = ",EIG%MatrixDim

  ! Allocate memory for the scattering data
  CALL AllocateScat(SD,Boxes(NumBoxes)%NumOpenR)

  DO iE = 1, NumE
     Energy = Egrid(iE)
     iBox = 1
     CALL CalcGamLam(BPD1,EIG1)
     CALL Mydggev(EIG1%MatrixDim,EIG1%Gam,EIG1%MatrixDim,EIG1%Lam,EIG1%MatrixDim,EIG1%eval,EIG1%evec)
     CALL BoxMatch(Bnull, Boxes(iBox), BPD1, EIG1, EffDim, AlphaFactor)
     CALL CalcK(Boxes(iBox),BPD1,SD,reducedmass,EffDim,AlphaFactor,Egrid(iE),M%Eth)
     DO iBox = 2, NumBoxes
        !BPD=BPD0  ! recall the stored BPD0 into BPD in case some things are being rewritten
        BPD%xl=Boxes(iBox)%xl
        BPD%xr=Boxes(iBox)%xr
        BPD%xPoints = BPD%xl + (BPD%xr-BPD%xl)*xprim
        DO kx = 1,BPD%xNumPoints
           DO lx = 1,LegPoints
              BPD%ux(lx,kx,1:BPD%xDim) = BPD0%ux(lx,kx,1:BPD%xDim)/(BPD%xr-BPD%xl)
           ENDDO
        ENDDO
        CALL SetMorsePotential(BPD,M)  ! would like to change this so we don't have to do it for each energy.
        CALL CalcGamLam(BPD,EIG)
        CALL Mydggev(EIG%MatrixDim,EIG%Gam,EIG%MatrixDim,EIG%Lam,EIG%MatrixDim,EIG%eval,EIG%evec)
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, EIG, EffDim, AlphaFactor)
        SD%K=0d0
        CALL CalcK(Boxes(iBox),BPD,SD,reducedmass,EffDim,AlphaFactor,Egrid(iE),M%Eth)
     ENDDO
     !WRITE(6,*) Energy, SD%sigma(1,1), SD%sigma(1,2), SD%sigma(2,1), SD%sigma(2,2)
     WRITE(6,*)  "K-matrix:"
     call printmatrix(SD%K,2,2,6)
     WRITE(10,*) Energy, SD%K
  ENDDO

20 FORMAT(1P,100e14.8)
END PROGRAM main
!****************************************************************************************************
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)

  DO j = 1,nr
     WRITE(file,20) (M(j,k), k = 1,nc)
  ENDDO

20 FORMAT(1P,100D20.12)
30 FORMAT(100F12.6)
END SUBROUTINE printmatrix
!****************************************************************************************************
SUBROUTINE GridMakerLinear(xNumPoints,x1,x2,xPoints)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: xNumPoints
  DOUBLE PRECISION, INTENT(in) :: x1,x2
  DOUBLE PRECISION, INTENT(out) :: xPoints(xNumPoints)
  INTEGER i
  DOUBLE PRECISION xDelt
  xDelt = (x2-x1)/DBLE(xNumPoints-1)
  DO i = 1,xNumPoints
     xPoints(i) = (i-1)*xDelt + x1 ! Simple linear grid
  ENDDO
END SUBROUTINE GridMakerLinear
!****************************************************************************************************
SUBROUTINE BoxMatch(BA, BB, BPD, EIG, dim, alphafact)
  USE DataStructures
  IMPLICIT NONE
  TYPE(BPData), INTENT(in) :: BPD
  TYPE(GenEigVal), INTENT(in) :: EIG
  TYPE(BoxData), INTENT(in) :: BA
  TYPE(BoxData) BB

  DOUBLE PRECISION, ALLOCATABLE :: tempnorm(:,:),  temp0(:,:), ZT(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: BigZA(:,:), BigZB(:,:),Dvec(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Dval(:)
  INTEGER, ALLOCATABLE :: Dikeep(:)
  DOUBLE PRECISION NewNorm, dim, alphafact, absd
  INTEGER i,j,k,beta,betaprime,nch,mch!,betamax
  INTEGER ikeep(BB%betaMax)

  j=1
  ikeep = 0
  DO i = 1, BPD%MatrixDim
     IF(ABS(EIG%eval(i)).GE.1e-12) THEN
        BB%b(j) = EIG%eval(i)
        ikeep(j)=i
        j = j+1
     ENDIF
  ENDDO

  ALLOCATE(tempnorm(BPD%MatrixDim,BB%betaMax))
  BB%Norm=0d0
  tempnorm=0d0
  BB%Nbeta=0d0
  DO beta=1,BB%betaMax
     DO betaprime=1,BB%betaMax
        BB%Norm(beta,betaprime)=0d0
        DO i=1,BPD%MatrixDim
           tempnorm(i,betaprime)=0.0d0
           DO j=1,BPD%MatrixDim
              tempnorm(i,betaprime) = tempnorm(i,betaprime) + EIG%Lam(i,j)*EIG%evec(j,ikeep(betaprime)) !
           ENDDO
           BB%Norm(beta,betaprime) = BB%Norm(beta,betaprime) + EIG%evec(i,ikeep(beta))*tempnorm(i,betaprime) !
        ENDDO
     ENDDO
     BB%Nbeta(beta) = dsqrt(BB%Norm(beta,beta))
  ENDDO

  DO beta = 1,BB%betaMax
     DO i = 1,BB%NumOpenL
        BB%Z(i,beta) = BPD%xl**(0.5d0*(dim-1d0-2d0*alphafact))*&
             EIG%evec((i-1)*BPD%xDim + 1, ikeep(beta))/BB%Nbeta(beta)!
        BB%ZP(i,beta) = -BB%b(beta)*BB%Z(i,beta)
     ENDDO
     DO i = 1, BB%NumOpenR
        BB%Z(i+BB%NumOpenL,beta) = BPD%xr**(0.5d0*(dim-1d0-2d0*alphafact))*&
             EIG%evec((i-1)*BPD%xDim + BPD%xDim,ikeep(beta))/BB%Nbeta(beta) !
        BB%ZP(i+BB%NumOpenL,beta) = -BB%b(beta)*BB%Z(i+BB%NumOpenL,beta)
     ENDDO
  ENDDO

  ALLOCATE(ZT(BB%betaMax,BB%betaMax))
  ALLOCATE(temp0(BB%betaMax,BB%betaMax))
  ZT = BB%Z
  temp0 = 0d0
  temp0 = MATMUL(TRANSPOSE(ZT),BB%Z)

  ALLOCATE(Dikeep(BB%betaMax + BB%NumOpenL))
  ALLOCATE(BigZA(BB%NumOpenL + BB%betaMax,BB%NumOpenL+BB%betaMax),BigZB(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax)) !
  ALLOCATE(Dvec(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax))
  ALLOCATE(Dval(BB%NumOpenL+BB%betaMax))

  Dvec=0d0
  BigZA=0d0
  BigZB=0d0

  DO i=1,BB%NumOpenL
     DO beta=1,BB%NumOpenL
        BigZA(i,beta)=BA%Zf(i,beta)
        BigZA(i+BB%NumOpenL,beta)=BA%Zfp(i,beta)
     ENDDO
  ENDDO

  DO i=1,BB%NumOpenL
     DO beta=1,BB%betaMax
        BigZA(i, BB%NumOpenL+beta) = -BB%Z(i,beta)
        BigZA(i+BB%NumOpenL,beta+BB%NumOpenL) = BB%ZP(i,beta)
     ENDDO
  ENDDO

  DO i=1,BB%NumOpenR
     DO beta=1,BB%betaMax
        BigZA(i+2*BB%NumOpenL, BB%NumOpenL+beta)=-BB%ZP(i+BB%NumOpenL,beta)
        BigZB(i+2*BB%NumOpenL, BB%NumOpenL+beta)=BB%Z(i+BB%NumOpenL,beta)
     ENDDO
  ENDDO

  CALL Mydggev(BB%NumOpenL+BB%betaMax,BigZA,BB%NumOpenL+BB%betaMax,BigZB,BB%NumOpenL+BB%betaMax,Dval,Dvec) !
  j=1
  DO i = 1,BB%NumOpenL+BB%betaMax
     absd=ABS(Dval(i))
     IF((absd.GE.1d-12).AND.(absd.LT.1d12)) THEN
        Dikeep(j)=i
        BB%bf(j)=Dval(i)
        j=j+1
     ENDIF
  ENDDO

  DO beta=1,BB%NumOpenR
     !for each beta, construct the final state with constant log-derivative at the right boundary
     DO i=1, BB%NumOpenR
        BB%Zf(i,beta) = 0.0d0
        DO betaprime = 1,BB%betaMax
           BB%Zf(i,beta) =  BB%Zf(i,beta) + BB%Z(BB%NumOpenL+i,betaprime)*Dvec(BB%NumOpenL+betaprime,Dikeep(beta)) !
        ENDDO
     ENDDO
     ! calculate the normalization of the final state
     NewNorm=0.0d0
     DO i=1,BB%NumOpenR
        NewNorm=NewNorm+BB%Zf(i,beta)**2
     ENDDO
     NewNorm=dsqrt(NewNorm)
     ! normalize and set the derivative Zfp
     DO i=1,BB%NumOpenR
        BB%Zf(i,beta) = BB%Zf(i,beta)/NewNorm
        BB%Zfp(i,beta) = -BB%Zf(i,beta)*BB%bf(beta)
     ENDDO
  ENDDO

  DEALLOCATE(Dikeep,Dvec,Dval,BigZA,BigZB)
  DEALLOCATE(temp0)
  DEALLOCATE(tempnorm)
  DEALLOCATE(ZT)
END SUBROUTINE BoxMatch
!****************************************************************************************************
SUBROUTINE checkbessel(xmin,xmax,npts,lam,d,file)
  IMPLICIT NONE
  DOUBLE PRECISION xmin, xmax
  DOUBLE PRECISION lam,hypj,hypy,hypjp,hypyp,hypi,hypk,hypip,hypkp
  DOUBLE PRECISION, ALLOCATABLE :: x(:)
  INTEGER d,file
  INTEGER n, npts

  ALLOCATE(x(npts))
  DO n = 1, npts
     x(n) = xmin + (n-1)*(xmax-xmin)/DBLE(npts-1)
  ENDDO
  DO n=1,npts
     !call hyperjy(d,lam,x(n),hypj,hypy,hypjp,hypyp)
     CALL hyperrjry(d,DBLE(0.5d0*(d-1d0)),lam,x(n),hypj,hypy,hypjp,hypyp)
     !call hyperik(d,lam,x(n),hypj,hypy,hypjp,hypyp)
     !call hyperrirk(d,dble(0.5d0*(d-1d0)),lam,x(n),hypj,hypy,hypjp,hypyp)
     WRITE(file,20) x(n), hypj, hypy, hypjp, hypyp
  ENDDO
20 FORMAT(1P,100D12.4)
END SUBROUTINE checkbessel
!****************************************************************************************************
SUBROUTINE SetZeroPotential(BPD)
  USE DataStructures
  USE Quadrature
  IMPLICIT NONE
  TYPE(BPData) BPD
  INTEGER kx,lx,mch,nch
  DOUBLE PRECISION ax,bx,xScaledZero,pot(3,3)
  DOUBLE PRECISION xScale(BPD%xNumPoints)

  DO kx = 1,BPD%xNumPoints-1
     ax = BPD%xPoints(kx)
     bx = BPD%xPoints(kx+1)
     xScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     DO lx = 1,LegPoints
        BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
        BPD%Pot(:,:,lx,kx) = 0d0
     ENDDO
  ENDDO
END SUBROUTINE SetZeroPotential
!****************************************************************************************************
SUBROUTINE testCG
  IMPLICIT NONE
  INTEGER J1D,J2D,JD,M1D,M2D,MD
  DOUBLE PRECISION cg
  DOUBLE PRECISION, EXTERNAL :: CLEBSH

  CALL SETUPAM

  J1D=1
  J2D=1
  JD=2
  M1D=1
  M2D=1
  MD=2

  CG = CLEBSH(J1D,J2D,JD,M1D,M2D,MD)
  WRITE(6,*) J1D, J2D, JD, M1D, M2D, MD, CG
END SUBROUTINE testCG
!****************************************************************************************************
