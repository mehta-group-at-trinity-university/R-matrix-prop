MODULE Quadrature
  IMPLICIT NONE
  INTEGER LegPoints
  DOUBLE PRECISION, ALLOCATABLE :: xLeg(:),wLeg(:)
  CHARACTER*64 LegendreFile

CONTAINS
  SUBROUTINE GetGaussFactors(File,Points,x,w)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c     This subroutine retrieves the points and weights for
    !c      Gaussian Quadrature from a given file
    !c
    !c     Variables:
    !c      File		name of file containing points and
    !c			 weights
    !c      Points		number of points to be used for
    !c			quadrature
    !c      x		array of nodes
    !c      w		array of weights
    !c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    IMPLICIT NONE
    INTEGER Points
    DOUBLE PRECISION x(Points),w(Points)
    CHARACTER*64 File
    INTEGER i,tempPoints

    OPEN(unit=7,file=File(1:INDEX(File,' ')-1))
    DO i = 1,18
       READ(7,*)
    ENDDO
    READ(7,*) tempPoints
    DO WHILE (tempPoints .NE. Points)
       DO i = 1,tempPoints
          READ(7,*)
       ENDDO
       READ(7,*) tempPoints
    ENDDO
    DO i = 1,Points
       READ(7,*) x(i),w(i)
    ENDDO
    CLOSE(unit=7)
    RETURN
  END SUBROUTINE GetGaussFactors

END MODULE Quadrature
  !****************************************************************************************************
MODULE GlobalVars
  IMPLICIT NONE
  INTEGER NumParticles, NumChannels,  NumAllChan, Order
  INTEGER StartBC, EndBC, xTotNumPoints, NumBoxes
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
  TYPE BPData !This data type contains the array of basis functions and potential matrix
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
     DOUBLE PRECISION, ALLOCATABLE :: K(:,:), R(:,:), f(:,:), sigma(:,:)
     COMPLEX*8, ALLOCATABLE :: S(:,:), T(:,:)

  END TYPE ScatData

CONTAINS
  !****************************************************************************************************
  SUBROUTINE AllocateScat(SD,N)
    IMPLICIT NONE
    TYPE(ScatData) SD
    INTEGER N
    ALLOCATE(SD%K(N,N),SD%R(N,N),SD%f(N,N),SD%sigma(N,N))
    ALLOCATE(SD%S(N,N),SD%T(N,N))
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
  SUBROUTINE Makebasis(BPD)
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

  END SUBROUTINE Makebasis
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
      DOUBLE PRECISION rhypj,rhypy,rhypjp,rhypyp
      DOUBLE PRECISION k(BPD%NumChannels),Eth(BPD%NumChannels)
      INTEGER i,j,no,nw,nc,beta

      rm=BPD%xr
      no=0
      nw=0
      nc=0

      DO i = 1,BPD%NumChannels
         IF (EE.GE.Eth(i)) THEN
            k(i) = dsqrt(2d0*mu*(EE)-Eth(i)) ! k is real
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

      ALLOCATE(s(no),c(no),Imat(no,no),Jmat(no,no))
      ALLOCATE(sp(no),cp(no))
      DO i = 1,no
         !        write(6,*) d, k(i), rm, k(i)*rm,BPD%lam(i)
         CALL hyperrjry(INT(d),alpha,BPD%lam(i),k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
         s(i) = dsqrt(mu)*rhypj  ! the factor of sqrt(mu) is for energy normalization
         c(i) = -dsqrt(mu)*rhypy ! the factor of sqrt(mu) is for energy normalization
         sp(i) = k(i)*dsqrt(mu)*rhypjp
         cp(i) = -k(i)*dsqrt(mu)*rhypyp
      ENDDO
      Imat=0d0
      Jmat=0d0
      DO i=1,no
         DO beta=1,no
            Imat(i,beta) = (B%Zf(i,beta)*cp(i) - B%Zfp(i,beta)*c(i))/(s(i)*cp(i)-c(i)*sp(i))
            Jmat(i,beta) = (B%Zf(i,beta)*sp(i) - B%Zfp(i,beta)*s(i))/(c(i)*sp(i)-s(i)*cp(i))
         ENDDO
      ENDDO
!!$      write(6,*) "Imat:"
!!$      call printmatrix(Imat,no,no,6)
!!$      write(6,*) "Jmat:"
!!$      call printmatrix(Jmat,no,no,6)
      CALL SqrMatInv(Imat, no)
      SD%K = MATMUL(Jmat,Imat)
      !WRITE(6,*) "B%NumOpenR = ", B%NumOpenR
      DO i=1,B%NumOpenR
         DO j=1,B%NumOpenR
            SD%R(i,j)=0.0d0
            DO beta = 1, B%NumOpenR
               SD%R(i,j) = SD%R(i,j) - B%Zf(i,beta)*B%Zf(j,beta)/B%bf(beta) !
            ENDDO
         ENDDO
      ENDDO
      !call sqrmatinv(BB%Zfp,BB%NumOpenR)  This gives the same result as the code segment  executed above
      !SD%R = matmul(BB%Zf,BB%Zfp)

      DEALLOCATE(s,c,Imat,Jmat,sp,cp)
    END SUBROUTINE CalcK

  END MODULE Scattering
  !****************************************************************************************************
  MODULE BalujaParameters
  IMPLICIT NONE
    DOUBLE PRECISION BalujaMCmat(3,3,2), BalujaEth(3)              ! this is the Baluja et al Multipole coupling matrix and thresholds
    DOUBLE PRECISION RMatBaluja1(3,3)                                   ! this is the Baluja et al R-matrix at r=5 for their test case
    DOUBLE PRECISION RMatBaluja2(3,3)                                   ! this is the Baluja et al R-matrix at r=15 for their test case
    DOUBLE PRECISION Znet                                               ! this is the charge seen by the electron
    DOUBLE PRECISION lmom(3)                                            ! partial wave for each channel
  CONTAINS
    !----------------------------------------------------------------------------------------------------
    ! This subroutine sets the multipole coupling elements for the test case described in Baluja et al CPC (1982) paper
    !----------------------------------------------------------------------------------------------------
    SUBROUTINE SetMultipoleCoup()
      IMPLICIT NONE

      lmom(1)=0d0
      lmom(2)=2d0
      lmom(3)=1d0

      Znet = 1d0

      !I beleive the thresholds here are really threshold values of 2*mu*Eth
      BalujaEth(1)=0.0d0
      BalujaEth(2)=0.0d0
      BalujaEth(3)=2.1689316d0

      ! These are the aij's
      BalujaMCmat(1,1,1)=0.0d0
      BalujaMCmat(1,2,1)=0.0d0
      BalujaMCmat(1,3,1)=-0.5682977d0

      BalujaMCmat(2,1,1)=0.0d0
      BalujaMCmat(2,2,1)=0.0d0
      BalujaMCmat(2,3,1)=-0.8036944d0

      BalujaMCmat(3,1,1)=-0.5682977d0
      BalujaMCmat(3,2,1)=-0.8036944d0
      BalujaMCmat(3,3,1)=0.0d0

      BalujaMCmat(1,1,2)=0.0d0
      BalujaMCmat(1,2,2)=-0.5554260d0
      BalujaMCmat(1,3,2)=0.0d0

      BalujaMCmat(2,1,2)=-0.5554260d0
      BalujaMCmat(2,2,2)=-0.3927455d0
      BalujaMCmat(2,3,2)=0.0d0

      BalujaMCmat(3,1,2)=0.0d0
      BalujaMCmat(3,2,2)=0.0d0
      BalujaMCmat(3,3,2)=0.0d0

      !THe starting R-matrix
      RMatBaluja1(1,1) = 0.145599d0
      RMatBaluja1(1,2) = 0.00559130d0
      RMatBaluja1(1,3) = 0.000514900d0
      RMatBaluja1(2,1) = 0.00559130d0
      RMatBaluja1(2,2) = -0.00811540d0
      RMatBaluja1(2,3) = -0.01735370d0
      RMatBaluja1(3,1) = 0.000514900d0
      RMatBaluja1(3,2) = -0.01735370d0
      RMatBaluja1(3,3) = 0.00829450d0

      !The ending R-matrix
      RMatBaluja2(1,1) = -0.392776d0
      RMatBaluja2(1,2) = -0.0134222d0
      RMatBaluja2(1,3) = -0.00283486d0
      RMatBaluja2(2,1) = -0.0134222d0
      RMatBaluja2(2,2) = 0.0211938d0
      RMatBaluja2(2,3) = 0.01156956d0
      RMatBaluja2(3,1) = -0.00283486d0
      RMatBaluja2(3,2) = 0.01156956d0
      RMatBaluja2(3,3) = 0.145266d0


    END SUBROUTINE SetMultipoleCoup
    !****************************************************************************************************
    SUBROUTINE SetBalujaPotential(BPD)
      USE DataStructures
      USE Quadrature
      IMPLICIT NONE
      TYPE(BPData) BPD
      INTEGER kx,lx,mch,nch
      DOUBLE PRECISION ax,bx,xScaledZero,pot(3,3)
      DOUBLE PRECISION xScale(BPD%xNumPoints)

      CALL SetMultipoleCoup()

      DO kx = 1,BPD%xNumPoints-1
         ax = BPD%xPoints(kx)
         bx = BPD%xPoints(kx+1)
         xScale(kx) = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         DO lx = 1,LegPoints
            BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
!            call POTLMX(3,BPD%x(lx,kx),pot)
            !BPD%Pot(:,:,lx,kx) = pot
            DO mch=1,BPD%NumChannels
               BPD%Pot(mch,mch,lx,kx) = BalujaPotential(mch,mch,BPD%x(lx,kx))
               DO nch=1,mch-1
                  BPD%Pot(mch,nch,lx,kx) = BalujaPotential(mch,nch,BPD%x(lx,kx))
                  BPD%Pot(nch,mch,lx,kx) = BPD%Pot(mch,nch,lx,kx) ! Potential is symmetric
               ENDDO
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE SetBalujaPotential
    !****************************************************************************************************
    DOUBLE PRECISION FUNCTION BalujaPotential(mch,nch,R)
      !----------------------------------------------------------------------------------------------------
      ! This subroutine returns the Baluja et al potential
      !----------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mch, nch
      INTEGER Lam
      DOUBLE PRECISION, INTENT(in) :: R
      DOUBLE PRECISION reducedmass

      reducedmass = 1.d0
      BalujaPotential = 0d0
      DO Lam = 1,2
         BalujaPotential = BalujaPotential + BalujaMCmat(mch,nch,Lam)*R**(-Lam-1d0)/(2.d0*reducedmass)
      ENDDO
      IF(mch.EQ.nch) THEN
         BalujaPotential = BalujaPotential + &
              lmom(mch)*(lmom(mch)+1d0)/(2d0*reducedmass*R**2) - Znet/R + BalujaEth(mch)/(2.d0*reducedmass) !
      END IF

      RETURN

    END FUNCTION BalujaPotential
    !****************************************************************************************************
  END MODULE BalujaParameters
  !****************************************************************************************************
      MODULE MorsePotential
      USE DataStructures
      USE Quadrature
      IMPLICIT NONE
      TYPE Morse
         DOUBLE PRECISION D(3),V(3,3),rc(3,3),a(3,3)
         DOUBLE PRECISION Eth(3)
         INTEGER NumMorseOpen,l
      END TYPE Morse

    CONTAINS
      SUBROUTINE InitMorse(M)
        TYPE(Morse) :: M
        M%D(1)=1.d0
        M%D(2)=4.d0
        M%D(3)=6.d0
        M%V=0d0
        M%V(1,2)=0.2d0
        M%V(2,1)=M%V(1,2)
        M%V(2,3)=0.2d0
        M%V(3,2)=M%V(2,3)
        M%a=1d0
        M%a(1,2)=0.5d0
        M%a(2,1)=M%a(1,2)
        M%a(2,3)=0.5d0
        M%a(3,2)=M%a(2,3)
        M%rc=1d0                !initialize coupling ranges to 1
        M%rc(1,2) = 0.5d0       !set any coupling ranges different from 1d0
        M%rc(2,1) = M%rc(1,2)
        M%l=0
        M%Eth(1) = 0d0
        M%Eth(2) = 3d0
        M%Eth(3) = 5d0

      END SUBROUTINE InitMorse

      FUNCTION Morse1(a,D,re,r)
        DOUBLE PRECISION Morse1, a, D, re, r
        Morse1 = D*(1d0 - dexp(-(r - re)/a))**2 - D
        RETURN
      END FUNCTION Morse1

      SUBROUTINE SetMorsePotential(BPD,M)

        IMPLICIT NONE
        TYPE(BPData) BPD
        TYPE(Morse), INTENT(in) :: M
        INTEGER kx,lx,mch,nch,NChan
        DOUBLE PRECISION ax,bx,xScaledZero,pot(3,3)
        DOUBLE PRECISION xScale(BPD%xNumPoints)
        Nchan=3
        BPD%Pot(:,:,:,:) = 0d0

!        write(6,*) "V:"
!        call printmatrix(M%V,3,3,6)
!        write(6,*) "D:"
!        call printmatrix(M%D,3,1,6)

        DO kx = 1,BPD%xNumPoints-1
           ax = BPD%xPoints(kx)
           bx = BPD%xPoints(kx+1)
           xScale(kx) = 0.5d0*(bx-ax)
           xScaledZero = 0.5d0*(bx+ax)
           DO lx = 1,LegPoints
              BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
              DO mch = 1,NChan
                 BPD%Pot(mch,mch,lx,kx) = Morse1(M%a(mch,mch),M%D(mch),M%rc(mch,mch),BPD%x(lx,kx)) + M%Eth(mch)
                 DO nch = 1, mch-1
                    BPD%Pot(mch,nch,lx,kx) = M%V(mch,nch)*dexp(-(BPD%x(lx,kx)-M%rc(mch,nch))**2/M%a(mch,nch)**2)
                    BPD%Pot(nch,mch,lx,kx) = BPD%Pot(mch,nch,lx,kx)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO

      END SUBROUTINE SetMorsePotential
    END MODULE MorsePotential
    !****************************************************************************************************
    SUBROUTINE makeEgrid(Egrid,NumE,E1,E2)
      DOUBLE PRECISION Egrid(NumE)
      DOUBLE PRECISION E1,E2
      INTEGER NumE, iE

      DO iE=1,NumE
         Egrid(iE) = E1 + (iE-1)*(E2-E1)/(NumE-1)
      ENDDO
    END SUBROUTINE makeEgrid
!****************************************************************************************************
PROGRAM main
  USE DataStructures
  USE GlobalVars
  USE BalujaParameters
  USE Quadrature
  USE MorsePotential
  USE scattering

  IMPLICIT NONE
  TYPE(BPData) BPD,BPD0
  TYPE(GenEigVal) EIG,EIG0
  TYPE(BoxData) BA, BB, Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  TYPE(Morse) :: M
  DOUBLE PRECISION, ALLOCATABLE :: evec(:,:), eval(:)!, temp0(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Egrid(:)
  DOUBLE PRECISION xDelt
  INTEGER NumE, iE, beta, i, iBox

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
  Boxes(1)%NumOpenR = 3
  Boxes(1)%betaMax = Boxes(1)%NumOpenL + Boxes(1)%NumOpenR
  Boxes(1)%xl=xStart
  Boxes(1)%xr=xDelt
  CALL AllocateBox(Boxes(1))
  CALL InitZeroBox(Boxes(1))
  WRITE(6,*) "Box",1,Boxes(1)%xl,Boxes(1)%xr
  CALL printmatrix(Boxes(1)%Zf,Boxes(1)%NumOpenR,Boxes(1)%NumOpenR,6)
  DO i = 2,NumBoxes
    Boxes(i)%NumOpenL = 3
    Boxes(i)%NumOpenR = 3
    Boxes(i)%betaMax = Boxes(1)%NumOpenL + Boxes(1)%NumOpenR
    Boxes(i)%xl=Boxes(i-1)%xr
    Boxes(i)%xr=Boxes(i-1)%xr+xDelt
    CALL AllocateBox(Boxes(i))
    CALL InitZeroBox(Boxes(i))
    WRITE(6,*) "Box",i,Boxes(i)%xl,Boxes(i)%xr
    CALL printmatrix(Boxes(i)%Zf,Boxes(i)%NumOpenR,Boxes(i)%NumOpenR,6)
  ENDDO

    ! Intitializes some BPD0 variables to the input values.  This is the basis set used at the left edge (r=0)
  BPD0%NumChannels = NumChannels
  BPD0%Order = Order
  BPD0%Left = 0
  BPD0%Right = 2
  BPD0%xNumPoints = xTotNumPoints
  BPD0%kl = kStart ! only relevant for Left = 3. This is the normal log-derivative at BPD%xl
  BPD0%kr = kEnd   ! only relevant for Right = 3. This is the normal log-derivative at BPD%xr

    ! Intitializes some BPD variables to the input values.
  BPD%NumChannels = NumChannels
  BPD%Order = Order
  BPD%Left = 2
  BPD%Right = 2
  BPD%xNumPoints = xTotNumPoints
  BPD%kl = kStart ! only relevant for Left = 3. This is the normal log-derivative at BPD%xl
  BPD%kr = kEnd ! only relevant for Right = 3. This is the normal log-derivative at BPD%xr


  CALL AllocateBPD(BPD0)
  CALL AllocateBPD(BPD)

  BPD0%xl=Boxes(1)%xl
  BPD0%xr=Boxes(1)%xr
  CALL GridMakerLinear(BPD0%xNumPoints,BPD0%xl,BPD0%xr,BPD0%xPoints)
  CALL Makebasis(BPD0)

  CALL InitMorse(M)

  !call checkpot(BPD0,100)
  !call checkpot(BPD,101)
  NumE=2000
  ALLOCATE(Egrid(NumE))
  CALL makeEgrid(Egrid,NumE,M%Eth(2)+0.01d0,M%Eth(3)-0.01d0)

  !----------------------------------------------------------------------------------------------------
  ! Comment/uncomment the next line if you want to print the basis to file fort.300
  !----------------------------------------------------------------------------------------------------
  !call CheckBasisBP(BPD0%xDim,BPD0%xNumPoints,BPD0%xBounds,BPD0%Order,300,LegPoints,BPD0%u,BPD0%x)
  !CALL CheckBasisBP(BPD%xDim,BPD%xNumPoints,BPD%xBounds,BPD%Order,301,LegPoints,BPD%u,BPD%x)

  CALL SetMorsePotential(BPD0,M)


  EIG0%MatrixDim=BPD0%MatrixDim
  CALL AllocateEIG(EIG0)
  WRITE(6,*) "EIG0%MatrixDim = ",EIG0%MatrixDim

  EIG%MatrixDim=BPD%MatrixDim
  CALL AllocateEIG(EIG)
  WRITE(6,*) "EIG%MatrixDim = ",EIG%MatrixDim

  CALL AllocateScat(SD,Boxes(NumBoxes)%NumOpenR)

  DO iE = 1, NumE
    Energy = Egrid(iE)
    DO iBox = 1, NumBoxes

      IF(iBox.eq.1) THEN

        CALL CalcGamLam(BPD0,EIG0)
        CALL Mydggev(EIG0%MatrixDim,EIG0%Gam,EIG0%MatrixDim,EIG0%Lam,EIG0%MatrixDim,EIG0%eval,EIG0%evec)
        CALL BoxMatch(Bnull, Boxes(iBox), BPD0, EIG0, EffDim, AlphaFactor)
        CALL CalcK(Boxes(iBox),BPD0,SD,reducedmass,EffDim,AlphaFactor,Egrid(iE),M%Eth)

      ELSE
        BPD%xl=Boxes(iBox)%xl
        BPD%xr=Boxes(iBox)%xr
        CALL GridMakerLinear(BPD%xNumPoints,BPD%xl,BPD%xr,BPD%xPoints)
        CALL Makebasis(BPD)
        CALL SetMorsePotential(BPD,M)

        CALL CalcGamLam(BPD,EIG)
        CALL Mydggev(EIG%MatrixDim,EIG%Gam,EIG%MatrixDim,EIG%Lam,EIG%MatrixDim,EIG%eval,EIG%evec)
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, EIG, EffDim, AlphaFactor)
        SD%K=0d0
        CALL CalcK(Boxes(iBox),BPD,SD,reducedmass,EffDim,AlphaFactor,Egrid(iE),M%Eth)



      ENDIF
    ENDDO
    write(6,*) "K-matrix:", SD%K
    WRITE(10,*) Energy, SD%K
  ENDDO

  CALL DeAllocateBPD(BPD)

  !----------------------------------------------------------------------------------------------------
  ! be sure to change the values of NumOpenL and NumOpenR above
  ! use this part for the box-matching algorithm
!   BPD%Left = 2
!   BPD%Right = 2
! !!$
!
!   BA%NumOpenL = 0
!   BA%NumOpenR = 3
!   BA%betaMax = BA%NumOpenL+BA%NumOpenR
!
!   BB%NumOpenL = 3
!   BB%NumOpenR = 3
!   energy = 2.5d0
!   BB%betaMax = 6
!   xstart = 5d0
!   xend = 5d0+10d0
!
!   reducedmass = 1d0
!
!   EIG%MatrixDim=BPD%MatrixDim
!   CALL Boxes(ieroBox(BA)
!   CALL InitZeroBox(BB)
!   CALL AllocateScat(SD,BA%NumOpenR)
!   CALL AllocateBPD(BPD)
!   BPD%xl=xStart
!   BPD%xr=xEnd
!   CALL GridMakerLinear(BPD%xNumPoints,BPD%xl,BPD%xr,BPD%xPoints)
!   CALL Makebasis(BPD)
!   CALL SetBalujaPotential(BPD)
!   CALL checkpot(BPD,101)
!   CALL AllocateEIG(EIG)
!
!   WRITE(6,*) "EIG%MatrixDim = ",EIG%MatrixDim
!
!   CALL MYDSYEV(xStart*RMatBaluja1,BA%betaMax,BA%b,BA%Z)
!   BA%Zf=BA%Z
!   BA%b = -1.d0/BA%b
!   BA%bf = BA%b
!   DO beta = 1,BA%betaMax
!      DO i = 1,BA%NumOpenR
!         BA%ZfP(i,beta) = -BA%bf(beta)*BA%Zf(i,beta)
!      ENDDO
!   ENDDO
!
!   WRITE(6,*) "calculating gamma for baluja problem"
!   CALL CalcGamLam(BPD,EIG)
!   CALL Mydggev(EIG%MatrixDim,EIG%Gam,EIG%MatrixDim,EIG%Lam,EIG%MatrixDim,EIG%eval,EIG%evec)
!   CALL BoxMatch(BA, BB, BPD, EIG, EffDim,AlphaFactor)
!   CALL CalcK(BB,BPD,SD,reducedmass,EffDim,AlphaFactor,Energy,BalujaEth)
!
!   WRITE(6,*) "The initial R-Matrix is:"
!   CALL printmatrix(RMatBaluja1,3,3,6)
!   WRITE(6,*) "Final R-matrix:"
!
!   CALL printmatrix(SD%R/xEnd,BB%NumOpenR,BB%NumOpenR,6)
! !!$  !call printmatrix(BB%RF,BB%NumOpenR,BB%NumOpenR,6)
!
!   WRITE(6,*) "The Exact R-Matrix From Baluja et al at R=15.0 is:"
!   CALL printmatrix(RMatBaluja2,3,3,6)
!

!  call checkbessel(0.0001d0,10d0,100,1d0,3,100)


20 FORMAT(1P,100e14.8)
END PROGRAM main
!****************************************************************************************************
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
  DOUBLE PRECISION NewNorm, dim, alphafact
  INTEGER i,j,k,beta,betaprime,nch,mch!,betamax
  INTEGER ikeep(BB%betaMax)


!    call printmatrix(EIG%eval,BPD%MatrixDim,1,6)
    j=1
    ikeep = 0
    DO i = 1, BPD%MatrixDim
       IF(ABS(EIG%eval(i)).GE.1e-12) THEN
          !write(6,"(A,T20,I5,T30,E12.6)") 'eval',i, EIG%eval(i)
          BB%b(j) = EIG%eval(i)
          ikeep(j)=i
          j = j+1
       ENDIF
    ENDDO
    !BB%betaMax=j-1

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
       !write(6,*) 'norm(beta,beta)=',BB%Norm(beta,beta),'Nbeta(',beta,') = ',BB%Nbeta(beta) !
    ENDDO

    !write(6,*) "Norm matrix:"
    !call printmatrix(BB%Norm,BB%betaMax,BB%betaMax,6)

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
    !write(6,*) "Z:"
!    call printmatrix(BB%Z,BB%betaMax,BB%betaMax,6)
!    call printmatrix(BB%ZP,BB%betaMax,BB%betaMax,6)

    ALLOCATE(temp0(BB%betaMax,BB%betaMax))
    ZT = BB%Z
    temp0 = 0d0
    temp0 = MATMUL(TRANSPOSE(ZT),BB%Z)
    !call dgemm('T','N',BB%betaMax,BB%betaMax,BB%betaMax,1.0d0,ZT,BB%betaMax,BB%Z,BB%betaMax,0.0d0,temp0,BB%betaMax) !

    !WRITE(6,*) "Check norm:"
    !CALL printmatrix(temp0,BB%betaMax,BB%betaMax,6)

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
          !write(6,*) "using BA%Z(), BA%ZP= = ",BA%Zf(i,beta), -BA%Zfp(i,beta),BA%Zf(i,beta)*BA%bf(beta)
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

    !print*, 'BigZA : '
    !call printmatrix(BigZA, BB%NumOpenL+BB%betaMax, BB%NumOpenL+BB%betaMax,6)
    !print*, 'BigZB : '
    !call printmatrix(BigZB, BB%NumOpenL+BB%betaMax, BB%NumOpenL+BB%betaMax,6)

    CALL Mydggev(BB%NumOpenL+BB%betaMax,BigZA,BB%NumOpenL+BB%betaMax,BigZB,BB%NumOpenL+BB%betaMax,Dval,Dvec) !
    !print*, 'Eigenvalues of box-matching eigenproblem : '
    !call printmatrix(Dval,BB%NumOpenL+BB%betaMax,1,6)

!!$    print*, 'Dvec : '
!!$    do j = 1,BB%betaMax
!!$       write(6,*) (Dvec(j,k), k = 1,BB%betaMax)
!!$    enddo
    j=1
    DO i = 1,BB%NumOpenL+BB%betaMax
       IF(ABS(Dval(i)).GE.1e-12) THEN
          Dikeep(j)=i
          BB%bf(j)=Dval(i)
          !write(6,*) i, j, BB%bf(j), (Dvec(k,i), k=BB%NumOpenL+1,BB%NumOpenL+BB%betaMax)
          j=j+1
       ENDIF
    ENDDO
    !print*, 'Final Number of Channels : ',j-1, "should be ",BB%NumOpenR
    DO beta=1,BB%NumOpenR
       !for each beta, construct the final state with constant log-derivative at the right boundary
       !write(6,*) "beta value being calculated now is: ", beta, BB%bf(beta)
       DO i=1, BB%NumOpenR
          BB%Zf(i,beta) = 0.0d0
          DO betaprime = 1,BB%betaMax
             BB%Zf(i,beta) =  BB%Zf(i,beta) + BB%Z(BB%NumOpenL+i,betaprime)*Dvec(BB%NumOpenL+betaprime,Dikeep(beta)) !
             !write(6,*) "adding...",BB%NumOpenL+i, betaprime, BB%Z(BB%NumOpenL+i,betaprime),"*",&
             !     Dvec(BB%NumOpenL+betaprime,Dikeep(beta))
          ENDDO
          !write(6,*) i,beta,"Zf = ", BB%Zf(i,beta)
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
    !PRINT*, 'Dvec : '
    !CALL printmatrix(Dvec,BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax,6)
    !print*, 'Zf'
    !call printmatrix(BB%Zf,BB%NumOpenR,BB%NumOpenR,6)

    DEALLOCATE(Dikeep,Dvec,Dval,BigZA,BigZB)
    DEALLOCATE(temp0)
    DEALLOCATE(tempnorm)
    DEALLOCATE(ZT)
  END SUBROUTINE BoxMatch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !use this routine to test if the potential is right.  Looks ok!
  SUBROUTINE POTLMX(NCHAN,R,V)
    USE BalujaParameters
    IMPLICIT NONE
    INTEGER I,J,NCHAN,NMX,LAMAX,K,ION,LCHL(3),EL2
    DOUBLE PRECISION V(NCHAN,NCHAN),R,VP,RR
    LAMAX=2
    ION=1
    NMX=3
    NCHAN=3
!C
!C      SETS UP EXAMPLE POTENTIAL FOR USE IN RPROP
!C

!C
    DO I=1,NCHAN
       EL2 = lmom(I)*(lmom(I)+1)
       DO J=1,NCHAN
          VP = 0.
          RR = 1.D0/R
          IF(I.EQ.J) VP=-2.D0*DFLOAT(ION)*RR+EL2*RR*RR + BalujaEth(I)
          DO K=1,LAMAX
             VP = VP+BalujaMCmat(I,J,K)*RR**(K+1)
          ENDDO
          V(I,J) = 0.5d0*VP
       ENDDO
    ENDDO

  END SUBROUTINE POTLMX
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
20  FORMAT(1P,100D12.4)
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

    !****************************************************************************************************
