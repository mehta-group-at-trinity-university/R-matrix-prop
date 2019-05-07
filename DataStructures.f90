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
     INTEGER NumOpenL, NumOpenR, betaMax, bfmax
     DOUBLE PRECISION xL, xR
     DOUBLE PRECISION, ALLOCATABLE :: Z(:,:), ZP(:,:), NBeta(:), Norm(:,:)
     DOUBLE PRECISION, ALLOCATABLE :: b(:), bf(:), Zf(:,:), Zfp(:,:)
     !double precision, allocatable :: RF(:,:),K(:,:)

  END TYPE BoxData
  !****************************************************************************************************
  TYPE ScatData

    DOUBLE PRECISION, ALLOCATABLE :: K(:,:), R(:,:), sigma(:,:), sigmatot(:)
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
    ALLOCATE(SD%sigmatot(-N:N))
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

    DEALLOCATE(SD%K,SD%R,SD%f,SD%sigma,SD%S,SD%T,SD%sigmatot)

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
