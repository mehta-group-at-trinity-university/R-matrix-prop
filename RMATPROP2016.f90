module Quadrature
  implicit none
  integer LegPoints
  double precision, allocatable :: xLeg(:),wLeg(:)
  character*64 LegendreFile

  contains 
    subroutine GetGaussFactors(File,Points,x,w)
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
      implicit none
      integer Points
      double precision x(Points),w(Points)
      character*64 File
      integer i,tempPoints

      open(unit=7,file=File(1:index(File,' ')-1))
      do i = 1,18
         read(7,*)
      enddo
      read(7,*) tempPoints
      do while (tempPoints .ne. Points)
         do i = 1,tempPoints
            read(7,*)
         enddo
         read(7,*) tempPoints
      enddo
      do i = 1,Points
         read(7,*) x(i),w(i)
      enddo
      close(unit=7)
      return
    End subroutine GetGaussFactors

  end module Quadrature
  !****************************************************************************************************
module GlobalVars
  implicit none
  integer NumParticles, NumChannels,  NumAllChan, Order
  integer NumSectors, StartBC, EndBC, xTotNumPoints
  !----------------------------------------------------------------------------------------------------
  double precision AlphaFactor ! This is the parameter that appears in the reduced wavefunction u(R) = R^(AlphaFactor) Psi(R)
  ! Typical choice is either AlphaFactor = 0 (reduced wavefunction = wavefunction), or AlphaFactor = (EffDim - 1)/2 (eliminates 1st derivative terms from KE)
  !----------------------------------------------------------------------------------------------------
  double precision reducedmass, xStart, xEnd, energy,kStart,kEnd
  double precision SpatialDim, EffDim
  double precision, allocatable :: mass(:)
  character*64 InputFile
  complex*16 II
  parameter(II=(0.0d0,1.0d0))
!****************************************************************************************************
contains
  subroutine ReadGlobal
    use Quadrature
    implicit none
    integer n
    ! Be sure to match the format of the input file to the format of the read statements below
    open(unit=7,file=InputFile(1:index(InputFile,' ')-1),action='read')
    read(7,*)
    read(7,*) NumParticles, NumChannels, SpatialDim, NumAllChan, Order
    allocate(mass(NumParticles))
    read(7,*)
    read(7,*)
    read(7,*) (mass(n), n=1,NumParticles)
    read(7,*)
    read(7,*)
    read(7,*) xStart, xEnd, xTotNumPoints, LegPoints
    read(7,*)
    read(7,*)
    read(7,*) NumSectors, StartBC, EndBC, kStart, kEnd
    read(7,*)
    read(7,*)
    read(7,*) Energy
    
    close(unit=7)
    EffDim = NumParticles*SpatialDim - SpatialDim
    !AlphaFactor = 0d0
    AlphaFactor = (EffDim-1d0)/2d0
    
    If (NumParticles.eq.2) then
       reducedmass = mass(1)*mass(2)/(mass(1)+mass(2))
    else
       write(6,*) "Reduced mass not set. Must set reduced mass"
       stop
    end If

    
  end subroutine ReadGlobal
  !****************************************************************************************************
end module GlobalVars
!****************************************************************************************************
!****************************************************************************************************
module DataStructures
  implicit none
  !****************************************************************************************************
  type BPData !This data type contains the array of basis functions and potential matrix
     integer xNumPoints, xDim, Left, Right, Order, NumChannels, MatrixDim
     double precision, allocatable :: u(:,:,:),ux(:,:,:),xPoints(:),x(:,:)
     double precision, allocatable :: Pot(:,:,:,:),lam(:)
     double precision kl, kr, xl,xr
     integer, allocatable :: xBounds(:), kxMin(:,:), kxMax(:,:)
  end type BPData
  !****************************************************************************************************
  type GenEigVal
     integer MatrixDim
     double precision, allocatable :: Gam(:,:), Lam(:,:), evec(:,:), eval(:)
  end type GenEigVal
  !****************************************************************************************************
  type BoxData
     integer NumOpenL, NumOpenR, betaMax
     double precision xL, xR
     double precision, allocatable :: Z(:,:), ZP(:,:), NBeta(:), Norm(:,:)
     double precision, allocatable :: b(:), bf(:), Zf(:,:), Zfp(:,:)
     !double precision, allocatable :: RF(:,:),K(:,:)
     
  end type BoxData
  type ScatData
     double precision, allocatable :: K(:,:), R(:,:), f(:,:), sigma(:,:)
     complex*8, allocatable :: S(:,:), T(:,:)
     
  end type ScatData
  
contains
  !****************************************************************************************************
  subroutine AllocateScat(SD,N)
    implicit none
    type(ScatData) SD
    integer N
    allocate(SD%K(N,N),SD%R(N,N),SD%f(N,N),SD%sigma(N,N))
    allocate(SD%S(N,N),SD%T(N,N))
    SD%K=0d0
    SD%R=0d0
    SD%f=0d0
    SD%sigma=0d0
    SD%S=(0d0,0d0)
    SD%T=(0d0,0d0)
  end subroutine AllocateScat
  !****************************************************************************************************
  subroutine DeAllocateScat(SD)
    implicit none
    type(ScatData) SD

    deallocate(SD%K,SD%R,SD%f,SD%sigma,SD%S,SD%T)
    
  end subroutine DeAllocateScat
  !****************************************************************************************************
  subroutine AllocateEIG(EIG)
    implicit none
    type(GenEigVal) :: EIG
    allocate(EIG%Gam(EIG%MatrixDim,EIG%MatrixDim),EIG%Lam(EIG%MatrixDim,EIG%MatrixDim))
    allocate(EIG%eval(EIG%MatrixDim),EIG%evec(EIG%MatrixDim,EIG%MatrixDim))
    EIG%Gam=0d0
    EIG%Lam=0d0
    EIG%eval=0d0
    EIG%Evec=0d0
    
  end subroutine AllocateEIG
  !****************************************************************************************************
  subroutine DeAllocateEIG(EIG)
    implicit none
    type(GenEigVal) :: EIG
    deallocate(EIG%Gam,EIG%Lam,EIG%eval,EIG%evec)
    
  end subroutine DeAllocateEIG
  !****************************************************************************************************
  subroutine AllocateBox(B)
    implicit none
    type(BoxData) :: B
    
    B%betaMax = B%NumOpenR + B%NumOpenL
    
    allocate(B%Z(B%betaMax,B%betaMax),B%ZP(B%betaMax,B%betaMax))
    allocate(B%b(B%betaMax))
    allocate(B%NBeta(B%betaMax))
    allocate(B%Norm(B%betaMax,B%betaMax))
    allocate(B%Zf(B%NumOpenR,B%NumOpenR),B%Zfp(B%NumOpenR,B%NumOpenR))
    allocate(B%bf(B%NumOpenR))

  end subroutine AllocateBox
  !****************************************************************************************************
  subroutine InitZeroBox(B)
    implicit none
    type(BoxData) :: B
    B%Z=0d0
    B%ZP=0d0
    B%b=0d0
    B%NBeta=0d0
    B%Norm=0d0
    B%Zf=0d0
    B%Zfp=0d0
    B%bf=0d0

  end subroutine InitZeroBox
  !****************************************************************************************************
  subroutine DeAllocateBox(B)
    implicit none
    type(BoxData) :: B

    deallocate(B%Z,B%ZP)
    deallocate(B%b)
    deallocate(B%NBeta,B%Norm,B%Zf,B%bf,B%Zfp)

    
  end subroutine DeAllocateBox
  !****************************************************************************************************
  subroutine AllocateBPD(BPD)
    use Quadrature
    implicit none
    type(BPData) BPD
    integer xDimMin
    !Determine the basis dimension for these boundary conditions
    xDimMin=BPD%xNumPoints+BPD%Order-3
    BPD%xDim=xDimMin
    if (BPD%Left .eq. 2) BPD%xDim = BPD%xDim + 1
    if (BPD%Right .eq. 2) BPD%xDim = BPD%xDim + 1
    BPD%MatrixDim = BPD%NumChannels*BPD%xDim
    
    allocate(BPD%xPoints(BPD%xNumPoints))
    allocate(BPD%xBounds(BPD%xNumPoints + 2*BPD%Order))   ! allocate xBounds before calling CalcBasisFuncs
    allocate(BPD%x(LegPoints,BPD%xNumPoints-1))
    
    allocate(BPD%u(LegPoints,BPD%xNumPoints,BPD%xDim),BPD%ux(LegPoints,BPD%xNumPoints,BPD%xDim)) ! allocate memory for basis functions
    allocate(BPD%kxMin(BPD%xDim,BPD%xDim),BPD%kxMax(BPD%xDim,BPD%xDim)) !
    allocate(BPD%Pot(BPD%NumChannels,BPD%NumChannels,LegPoints,BPD%xNumPoints-1))

    allocate(BPD%lam(BPD%NumChannels))

    BPD%lam=0.d0
    
  end subroutine AllocateBPD
  !****************************************************************************************************
  subroutine DeAllocateBPD(BPD)
    implicit none
    type(BPData) BPD
    deallocate(BPD%xPoints,BPD%xBounds,BPD%x,BPD%u,BPD%ux,BPD%kxMin,BPD%kxMax,BPD%Pot,BPD%lam)
  end subroutine DeAllocateBPD
  !****************************************************************************************************
  subroutine Makebasis(BPD)
    use Quadrature
    implicit none
    
    type(BPData) BPD
    integer ix, ixp
       
    call CalcBasisFuncsBP(BPD%Left,BPD%Right,BPD%kl,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         BPD%xDim,BPD%xBounds,BPD%xNumPoints,0,BPD%u)
    call CalcBasisFuncsBP(BPD%Left,BPD%Right,BPD%kl,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         BPD%xDim,BPD%xBounds,BPD%xNumPoints,1,BPD%ux)
    
    ! Determine the bounds for integration of matrix elements
    do ix = 1,BPD%xDim
       do ixp = 1,BPD%xDim
          BPD%kxMin(ixp,ix) = max(BPD%xBounds(ix),BPD%xBounds(ixp))
          BPD%kxMax(ixp,ix) = min(BPD%xBounds(ix+BPD%Order+1),BPD%xBounds(ixp+BPD%Order+1))-1 !
       enddo
    enddo
    
  end subroutine Makebasis
  !****************************************************************************************************
  subroutine checkpot(BPD,file)
    use Quadrature
    implicit none
    integer file,mch,nch,lx,kx
    type(BPData) BPD

    do mch = 1,BPD%NumChannels
       do nch = 1,mch
          do kx = 1, BPD%xNumPoints-1
             do lx = 1,LegPoints
                write(file,*) BPD%x(lx,kx), BPD%Pot(mch,nch,lx,kx)
             enddo
          enddo
          write(file,*)
       enddo
       write(file,*)
    enddo
    
  end subroutine checkpot
  
end module DataStructures
!****************************************************************************************************
  subroutine CalcGamLam(BPD,EIG)
    use DataStructures
    use Quadrature
    use GlobalVars
    implicit none
    type(BPData), intent(in) :: BPD
    type(GenEigVal) :: EIG
    integer ix,ixp,lx,kx,mch,nch
    double precision a, ax, bx, xIntScale(BPD%xNumPoints), TempG, xScaledZero
    EIG%Gam=0d0
    EIG%Lam=0d0
    EIG%eval=0d0
    EIG%evec=0d0
    !Calculate the Gamma Matrix elements
    do ix = 1,BPD%xDim
       do ixp = max(1,ix-BPD%Order),min(BPD%xDim,ix+BPD%Order)
          do mch = 1,BPD%NumChannels
             ! Do the channel-diagonal part first:
             do kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                ax = BPD%xPoints(kx)
                bx = BPD%xPoints(kx+1)
                xIntScale(kx) = 0.5d0*(bx-ax)
                xScaledZero = 0.5d0*(bx+ax)
                TempG = 0.0d0
                do lx = 1,LegPoints
                   a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                   ! The KE matrix elements
                   TempG = TempG - a*(BPD%ux(lx,kx,ix)*BPD%ux(lx,kx,ixp))
                   ! The diagonal "overlap"*energy and additional term from reducing the wavefunction
                   TempG = TempG + a*BPD%u(lx,kx,ix)*2*reducedmass*(Energy - (BPD%Pot(mch,mch,lx,kx) - &
                        AlphaFactor*(AlphaFactor-EffDim+2)/(2d0*reducedmass*BPD%x(lx,kx)**2)))*BPD%u(lx,kx,ixp)
                enddo
                EIG%Gam((mch-1)*BPD%xDim+ix,(mch-1)*BPD%xDim+ixp) = EIG%Gam((mch-1)*BPD%xDim+ix,(mch-1)*BPD%xDim+ixp) + TempG ! place values into Gamma0
             enddo

             ! Now do the off-diagonal parts
             do nch = 1, mch-1
                do kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                   ax = BPD%xPoints(kx)
                   bx = BPD%xPoints(kx+1)
                   xIntScale(kx) = 0.5d0*(bx-ax)
                   xScaledZero = 0.5d0*(bx+ax)
                   TempG = 0.0d0
                   do lx = 1,LegPoints
                      a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                      ! The potential matrix elements calculated here
                      TempG = TempG - a*BPD%u(lx,kx,ix)*2.0d0*reducedmass*BPD%Pot(mch,nch,lx,kx)*BPD%u(lx,kx,ixp)
                   enddo
                   EIG%Gam((mch-1)*BPD%xDim+ix,(nch-1)*BPD%xDim+ixp) = EIG%Gam((mch-1)*BPD%xDim +ix, (nch-1)*BPD%xDim +ixp) + TempG ! place values into EIG%Gam
                enddo
                EIG%Gam((nch-1)*BPD%xDim+ix,(mch-1)*BPD%xDim+ixp) = EIG%Gam((mch-1)*BPD%xDim+ix,(nch-1)*BPD%xDim+ixp)  ! fill in symmetric element mch <--> nch
             enddo
          enddo
       enddo
    enddo

    !Calculate the Lambda matrix elements
    EIG%Lam = 0.d0

    if(BPD%Right.eq.2) then
       do mch = 1, BPD%NumChannels
          EIG%Lam( (mch-1)*BPD%xDim + BPD%xDim, (mch-1)*BPD%xDim + BPD%xDim ) = BPD%xr**(EffDim-1d0-2d0*AlphaFactor)
       enddo
    endif
    if(BPD%Left.eq.2) then
       do mch = 1, NumChannels
          EIG%Lam((mch-1)*BPD%xDim + 1,(mch-1)*BPD%xDim + 1) = BPD%xl**(EffDim-1d0-2d0*AlphaFactor)
       enddo
    endif

  end subroutine CalcGamLam
  !****************************************************************************************************
  module Scattering
    use DataStructures
  contains
    subroutine CalcK(B,BPD,SD,mu,d,alpha,EE,Eth)
      implicit none
      type(BoxData), intent(in) :: B
      type(BPData), intent(in) :: BPD
      type(ScatData) :: SD
      double precision mu, EE,rm,d,alpha
      double precision, allocatable :: s(:),c(:),sp(:),cp(:)
      double precision, allocatable :: Imat(:,:),Jmat(:,:)
      double precision rhypj,rhypy,rhypjp,rhypyp      
      double precision k(BPD%NumChannels),Eth(BPD%NumChannels)
      integer i,j,no,nw,nc,beta

      rm=BPD%xr
      no=0
      nw=0
      nc=0

      do i = 1,BPD%NumChannels
         if (EE.ge.Eth(i)) then
            k(i) = dsqrt(2d0*mu*(EE)-Eth(i)) ! k is real
            no=no+1
         else
            k(i) = dsqrt(2d0*mu*(Eth(i)-EE)) ! k->kappa, kappa is real
            if( (k(i)*rm).lt.10d0) nw = nw+1
            if( (k(i)*rm).ge.10d0) nc = nc+1
         endif
      enddo
!      write(6,*) "no = ", no
      if((no+nw+nc).ne.BPD%NumChannels) then
         write(6,*) "Channel miscount in calcK"
         stop
      endif

      allocate(s(no),c(no),Imat(no,no),Jmat(no,no))
      allocate(sp(no),cp(no))
      do i = 1,no
         !        write(6,*) d, k(i), rm, k(i)*rm,BPD%lam(i)
         call hyperrjry(int(d),alpha,BPD%lam(i),k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
         s(i) = dsqrt(mu)*rhypj  ! the factor of sqrt(mu) is for energy normalization
         c(i) = -dsqrt(mu)*rhypy ! the factor of sqrt(mu) is for energy normalization
         sp(i) = k(i)*dsqrt(mu)*rhypjp
         cp(i) = -k(i)*dsqrt(mu)*rhypyp

      enddo
      Imat=0d0
      Jmat=0d0
      do i=1,no
         do beta=1,no
            Imat(i,beta) = (B%Zf(i,beta)*cp(i) - B%Zfp(i,beta)*c(i))/(s(i)*cp(i)-c(i)*sp(i))
            Jmat(i,beta) = (B%Zf(i,beta)*sp(i) - B%Zfp(i,beta)*s(i))/(c(i)*sp(i)-s(i)*cp(i))
         enddo
      enddo
!!$      write(6,*) "Imat:"
!!$      call printmatrix(Imat,no,no,6)
!!$      write(6,*) "Jmat:"
!!$      call printmatrix(Jmat,no,no,6)
      call SqrMatInv(Imat, no)
      SD%K = matmul(Jmat,Imat)
      write(6,*) "B%NumOpenR = ", B%NumOpenR
      do i=1,B%NumOpenR
         do j=1,B%NumOpenR
            SD%R(i,j)=0.0d0
            do beta = 1, B%NumOpenR
               SD%R(i,j) = SD%R(i,j) - B%Zf(i,beta)*B%Zf(j,beta)/B%bf(beta) ! 
            enddo
         enddo
      enddo
      !call sqrmatinv(BB%Zfp,BB%NumOpenR)  This gives the same result as the code segment  executed above
      !SD%R = matmul(BB%Zf,BB%Zfp)

      deallocate(s,c,Imat,Jmat,sp,cp)
    end subroutine CalcK
    
  end module Scattering
  !****************************************************************************************************
  module BalujaParameters
  implicit none
    double precision BalujaMCmat(3,3,2), BalujaEth(3)              ! this is the Baluja et al Multipole coupling matrix and thresholds
    double precision RMatBaluja1(3,3)                                   ! this is the Baluja et al R-matrix at r=5 for their test case
    double precision RMatBaluja2(3,3)                                   ! this is the Baluja et al R-matrix at r=15 for their test case
    double precision Znet                                               ! this is the charge seen by the electron
    double precision lmom(3)                                            ! partial wave for each channel
  contains
    !----------------------------------------------------------------------------------------------------
    ! This subroutine sets the multipole coupling elements for the test case described in Baluja et al CPC (1982) paper
    !----------------------------------------------------------------------------------------------------
    subroutine SetMultipoleCoup()
      implicit none

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
      
      
    end subroutine SetMultipoleCoup
    !****************************************************************************************************
    subroutine SetBalujaPotential(BPD)
      use DataStructures
      use Quadrature
      implicit none
      type(BPData) BPD
      integer kx,lx,mch,nch
      double precision ax,bx,xScaledZero,pot(3,3)
      double precision xScale(BPD%xNumPoints)

      call SetMultipoleCoup()
      
      do kx = 1,BPD%xNumPoints-1
         ax = BPD%xPoints(kx)
         bx = BPD%xPoints(kx+1)
         xScale(kx) = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do lx = 1,LegPoints
            BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
!            call POTLMX(3,BPD%x(lx,kx),pot)
            !BPD%Pot(:,:,lx,kx) = pot
            do mch=1,BPD%NumChannels
               BPD%Pot(mch,mch,lx,kx) = BalujaPotential(mch,mch,BPD%x(lx,kx))
               do nch=1,mch-1
                  BPD%Pot(mch,nch,lx,kx) = BalujaPotential(mch,nch,BPD%x(lx,kx))
                  BPD%Pot(nch,mch,lx,kx) = BPD%Pot(mch,nch,lx,kx) ! Potential is symmetric
               enddo
            enddo
         enddo
      enddo
    end subroutine SetBalujaPotential
    !****************************************************************************************************
    double precision function BalujaPotential(mch,nch,R)
      !----------------------------------------------------------------------------------------------------
      ! This subroutine returns the Baluja et al potential
      !----------------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: mch, nch
      integer Lam
      double precision, intent(in) :: R
      double precision reducedmass

      reducedmass = 1.d0
      BalujaPotential = 0d0
      do Lam = 1,2
         BalujaPotential = BalujaPotential + BalujaMCmat(mch,nch,Lam)*R**(-Lam-1d0)/(2.d0*reducedmass)
      enddo
      if(mch.eq.nch) then
         BalujaPotential = BalujaPotential + &
              lmom(mch)*(lmom(mch)+1d0)/(2d0*reducedmass*R**2) - Znet/R + BalujaEth(mch)/(2.d0*reducedmass) !
      end if
      
      return
      
    end function BalujaPotential
    !****************************************************************************************************
  end module BalujaParameters
  !****************************************************************************************************
      module MorsePotential
      use DataStructures
      use Quadrature
      implicit none
      type Morse
         double precision D(3),V(3,3),rc(3,3),a(3,3)
         double precision Eth(3)
         integer NumMorseOpen,l
      end type Morse

    contains
      subroutine InitMorse(M)
        type(Morse) :: M
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
        
      end subroutine InitMorse
      
      function Morse1(a,D,re,r)
        double precision Morse1, a, D, re, r
        Morse1 = D*(1d0 - dexp(-(r - re)/a))**2 - D
        return
      end function Morse1
      
      subroutine SetMorsePotential(BPD,M)
        
        implicit none
        type(BPData) BPD
        type(Morse), intent(in) :: M
        integer kx,lx,mch,nch,NChan
        double precision ax,bx,xScaledZero,pot(3,3)
        double precision xScale(BPD%xNumPoints)
        Nchan=3
        BPD%Pot(:,:,:,:) = 0d0

!        write(6,*) "V:"
!        call printmatrix(M%V,3,3,6)
!        write(6,*) "D:"
!        call printmatrix(M%D,3,1,6)
        
        do kx = 1,BPD%xNumPoints-1
           ax = BPD%xPoints(kx)
           bx = BPD%xPoints(kx+1)
           xScale(kx) = 0.5d0*(bx-ax)
           xScaledZero = 0.5d0*(bx+ax)
           do lx = 1,LegPoints
              BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
              do mch = 1,NChan
                 BPD%Pot(mch,mch,lx,kx) = Morse1(M%a(mch,mch),M%D(mch),M%rc(mch,mch),BPD%x(lx,kx)) + M%Eth(mch)
                 do nch = 1, mch-1
                    BPD%Pot(mch,nch,lx,kx) = M%V(mch,nch)*dexp(-(BPD%x(lx,kx)-M%rc(mch,nch))**2/M%a(mch,nch)**2)
                    BPD%Pot(nch,mch,lx,kx) = BPD%Pot(mch,nch,lx,kx)
                 enddo
              enddo
           enddo
        enddo

      end subroutine SetMorsePotential
    end module MorsePotential
    !****************************************************************************************************
    subroutine makeEgrid(Egrid,NumE,E1,E2)
      double precision Egrid(NumE)
      double precision E1,E2
      integer NumE, iE

      do iE=1,NumE
         Egrid(iE) = E1 + (iE-1)*(E2-E1)/(NumE-1)
      enddo
    end subroutine makeEgrid
!****************************************************************************************************
program main
  use DataStructures
  use GlobalVars
  use BalujaParameters
  use Quadrature
  use MorsePotential
  use scattering

  implicit none
  type(BPData) BPD,BPD0
  type(GenEigVal) EIG,EIG0
  type(BoxData) BA, BB, Bnull
  type(ScatData) SD
  type(Morse) :: M
  double precision, allocatable :: evec(:,:), eval(:)!, temp0(:,:)
  double precision, allocatable :: Egrid(:)
  integer NumE, iE, beta, i
  
  !----------------------------------------------------------------------
  ! Read information from the input file
  !----------------------------------------------------------------------
  InputFile = 'RMATPROP.inp'
  call ReadGlobal()
  !--------------------------------------------------
  ! Read in the quadrature data stored in Quadrature module.
  !--------------------------------------------------
  LegendreFile = 'Legendre.dat'
  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

    ! Intitializes some BPD0 variables to the input values.  This is the basis set used at the left edge (r=0)
  BPD0%NumChannels = NumChannels
  BPD0%Order = Order
  BPD0%Left = StartBC
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

  
  call AllocateBPD(BPD0)
  BPD0%xl=xStart
  BPD0%xr=xEnd/2d0
  call GridMakerLinear(BPD0%xNumPoints,BPD0%xl,BPD0%xr,BPD0%xPoints)
  call Makebasis(BPD0)

  call AllocateBPD(BPD)
  !BPD%xl=xStart!BPD0%xr
  BPD%xl=BPD0%xr
  BPD%xr=xEnd
  call GridMakerLinear(BPD%xNumPoints,BPD%xl,BPD%xr,BPD%xPoints)
  call Makebasis(BPD)

!!$  write(6,"(A,T25,2I3)") 'Left BC...', BPD%Left
!!$  write(6,"(A,T25,2I3)") 'Right BC...', BPD%Right
!!$  write(6,"(A,T25,2I3)") 'xNumPoints....', BPD%xNumPoints
!!$  write(6,"(A,T25,2I3)") 'Order....', BPD%Order
!!$  write(6,"(A,T25,F3.1)") 'EffDim...', EffDim
!!$  write(6,"(A,T25,F3.1)") 'Energy...', Energy
!!$  write(6,"(A,T20,F8.4)") 'Reduced mass...', reducedmass
!!$  write(6,"(A,T20,100F8.4)") 'xPoints...', BPD%xPoints
!!$  write(6,"(A,T25,2I3)") 'xDim....', BPD%xDim
!!$  write(6,"(A,T20,I8)") 'MatrixDim...', BPD%MatrixDim      


  !--------------------------------------------------
  ! Initialize the Baluja et al data and print it to the screen
  !--------------------------------------------------
  !call SetBalujaPotential(BPD)
  !call SetZeroPotential(BPD)
  call InitMorse(M)

  call SetMorsePotential(BPD0,M)
  call SetMorsePotential(BPD,M)
  !call checkpot(BPD0,100)
  !call checkpot(BPD,101)
  NumE=1000
  allocate(Egrid(NumE))
  call makeEgrid(Egrid,NumE,M%Eth(1)+0.01d0,M%Eth(2)-0.01d0)
  
  
  !----------------------------------------------------------------------------------------------------
  ! Comment/uncomment the next line if you want to print the basis to file fort.300
  !----------------------------------------------------------------------------------------------------
  !call CheckBasisBP(BPD0%xDim,BPD0%xNumPoints,BPD0%xBounds,BPD0%Order,300,LegPoints,BPD0%u,BPD0%x)
  call CheckBasisBP(BPD%xDim,BPD%xNumPoints,BPD%xBounds,BPD%Order,301,LegPoints,BPD%u,BPD%x)

!  write(6,*) '--------------------------------'
!  write(6,*) 'The starting R-Matrix at r=5.d0:'
!  call printmatrix(RMatBaluja1,3,3,6)
!  write(6,*) '--------------------------------'

  Bnull%NumOpenL = 0
  Bnull%NumOpenR = 0
  Bnull%betaMax = 0
  
  BA%NumOpenL = 0
  BA%NumOpenR = 3
  BA%betaMax = BA%NumOpenL+BA%NumOpenR
  
  BB%NumOpenL = 3
  BB%NumOpenR = 3
  BB%betaMax = BB%NumOpenL+BB%NumOpenR

  call AllocateBox(Bnull)
  call AllocateBox(BA)
  call AllocateBox(BB)

  
  EIG0%MatrixDim=BPD0%MatrixDim
  call AllocateEIG(EIG0)
  write(6,*) "EIG0%MatrixDim = ",EIG0%MatrixDim

  EIG%MatrixDim=BPD%MatrixDim
  call AllocateEIG(EIG)
  write(6,*) "EIG%MatrixDim = ",EIG%MatrixDim

  do iE = 1, NumE

     Energy = Egrid(iE)
     call InitZeroBox(BA)
     call InitZeroBox(Bnull)
     call AllocateScat(SD,BA%NumOpenR)
     call CalcGamLam(BPD0,EIG0)
     call Mydggev(EIG0%MatrixDim,EIG0%Gam,EIG0%MatrixDim,EIG0%Lam,EIG0%MatrixDim,EIG0%eval,EIG0%evec)
     call BoxMatch(Bnull, BA, BPD0, EIG0, EffDim, AlphaFactor)
     call CalcK(BA,BPD0,SD,reducedmass,EffDim,AlphaFactor,Egrid(iE),M%Eth)
     write(9,*) Energy, SD%K(1,1)
     call DeAllocateScat(SD)
     
     call AllocateScat(SD,BB%NumOpenR)
     call InitZeroBox(BB)
     call CalcGamLam(BPD,EIG)
     call Mydggev(EIG%MatrixDim,EIG%Gam,EIG%MatrixDim,EIG%Lam,EIG%MatrixDim,EIG%eval,EIG%evec)
     call BoxMatch(BA, BB, BPD, EIG, EffDim, AlphaFactor)
     SD%K=0d0
     call CalcK(BB,BPD,SD,reducedmass,EffDim,AlphaFactor,Egrid(iE),M%Eth)
     write(10,*) Energy, SD%K(1,1)
     write(6,*) Energy, SD%K(1,1)
     call DeAllocateScat(SD)
     
  enddo

  call DeAllocateBPD(BPD)
  call DeAllocateEIG(EIG)
  !----------------------------------------------------------------------------------------------------
  ! be sure to change the values of NumOpenL and NumOpenR above
  ! use this part for the box-matching algorithm
  BPD%Left = 2
  BPD%Right = 2
!!$

  BA%NumOpenL = 0
  BA%NumOpenR = 3
  BA%betaMax = BA%NumOpenL+BA%NumOpenR

  BB%NumOpenL = 3
  BB%NumOpenR = 3
  energy = 2.5d0
  BB%betaMax = 6
  xstart = 5d0
  xend = 5d0+10d0

  reducedmass = 1d0
  
  EIG%MatrixDim=BPD%MatrixDim
  call InitZeroBox(BA)
  call InitZeroBox(BB)
  call AllocateScat(SD,BA%NumOpenR)
  call AllocateBPD(BPD)
  BPD%xl=xStart
  BPD%xr=xEnd
  call GridMakerLinear(BPD%xNumPoints,BPD%xl,BPD%xr,BPD%xPoints)
  call Makebasis(BPD)
  call SetBalujaPotential(BPD)
  call checkpot(BPD,101)
  call AllocateEIG(EIG)
  
  write(6,*) "EIG%MatrixDim = ",EIG%MatrixDim

  call MYDSYEV(xStart*RMatBaluja1,BA%betaMax,BA%b,BA%Z)
  BA%Zf=BA%Z
  BA%b = -1.d0/BA%b
  BA%bf = BA%b
  do beta = 1,BA%betaMax
     do i = 1,BA%NumOpenR
        BA%ZfP(i,beta) = -BA%bf(beta)*BA%Zf(i,beta)
     enddo
  enddo

  write(6,*) "calculating gamma for baluja problem"
  call CalcGamLam(BPD,EIG)
  call Mydggev(EIG%MatrixDim,EIG%Gam,EIG%MatrixDim,EIG%Lam,EIG%MatrixDim,EIG%eval,EIG%evec)
  call BoxMatch(BA, BB, BPD, EIG, EffDim,AlphaFactor)
  call CalcK(BB,BPD,SD,reducedmass,EffDim,AlphaFactor,Energy,BalujaEth)

  write(6,*) "The initial R-Matrix is:"
  call printmatrix(RMatBaluja1,3,3,6)
  write(6,*) "Final R-matrix:"

  call printmatrix(SD%R/xEnd,BB%NumOpenR,BB%NumOpenR,6)
!!$  !call printmatrix(BB%RF,BB%NumOpenR,BB%NumOpenR,6)

  write(6,*) "The Exact R-Matrix From Baluja et al at R=15.0 is:"
  call printmatrix(RMatBaluja2,3,3,6)


!  call checkbessel(0.0001d0,10d0,100,1d0,3,100)
  

20 format(1P,100e14.8)
end program main
!****************************************************************************************************
!****************************************************************************************************
subroutine printmatrix(M,nr,nc,file)
  implicit none
  integer nr,nc,file,j,k
  double precision M(nr,nc)
  
  do j = 1,nr
     write(file,20) (M(j,k), k = 1,nc)
  enddo

20 format(1P,100D20.12)
30 format(100F12.6)
end subroutine printmatrix
!****************************************************************************************************
subroutine GridMakerLinear(xNumPoints,x1,x2,xPoints)
  implicit none
  integer, intent(in) :: xNumPoints
  double precision, intent(in) :: x1,x2
  double precision, intent(out) :: xPoints(xNumPoints)
  integer i
  double precision xDelt
  xDelt = (x2-x1)/dble(xNumPoints-1)
  do i = 1,xNumPoints
     xPoints(i) = (i-1)*xDelt + x1 ! Simple linear grid
  enddo
end subroutine GridMakerLinear
!****************************************************************************************************
subroutine BoxMatch(BA, BB, BPD, EIG, dim, alphafact)
  use DataStructures
  implicit none
  type(BPData), intent(in) :: BPD
  type(GenEigVal), intent(in) :: EIG
  type(BoxData), intent(in) :: BA
  type(BoxData) BB
  
  double precision, allocatable :: tempnorm(:,:),  temp0(:,:), ZT(:,:)
  double precision, allocatable :: BigZA(:,:), BigZB(:,:),Dvec(:,:)
  double precision, allocatable :: Dval(:)
  integer, allocatable :: Dikeep(:)
  double precision NewNorm, dim, alphafact
  integer i,j,k,beta,betaprime,nch,mch!,betamax
  integer ikeep(BB%betaMax)
  

!    call printmatrix(EIG%eval,BPD%MatrixDim,1,6)
    j=1
    ikeep = 0
    do i = 1, BPD%MatrixDim
       if(abs(EIG%eval(i)).ge.1e-12) then
          !write(6,"(A,T20,I5,T30,E12.6)") 'eval',i, EIG%eval(i)
          BB%b(j) = EIG%eval(i)
          ikeep(j)=i
          j = j+1
       endif
    enddo
    !BB%betaMax=j-1

    allocate(tempnorm(BPD%MatrixDim,BB%betaMax))
    BB%Norm=0d0
    tempnorm=0d0
    BB%Nbeta=0d0
    do beta=1,BB%betaMax
       do betaprime=1,BB%betaMax
          BB%Norm(beta,betaprime)=0d0
          do i=1,BPD%MatrixDim
             tempnorm(i,betaprime)=0.0d0
             do j=1,BPD%MatrixDim
                tempnorm(i,betaprime) = tempnorm(i,betaprime) + EIG%Lam(i,j)*EIG%evec(j,ikeep(betaprime)) ! 
             enddo
             BB%Norm(beta,betaprime) = BB%Norm(beta,betaprime) + EIG%evec(i,ikeep(beta))*tempnorm(i,betaprime) !  
          enddo
       enddo
       BB%Nbeta(beta) = dsqrt(BB%Norm(beta,beta))
       !write(6,*) 'norm(beta,beta)=',BB%Norm(beta,beta),'Nbeta(',beta,') = ',BB%Nbeta(beta) ! 
    enddo

    !write(6,*) "Norm matrix:"
    !call printmatrix(BB%Norm,BB%betaMax,BB%betaMax,6)

    do beta = 1,BB%betaMax
       do i = 1,BB%NumOpenL
          BB%Z(i,beta) = BPD%xl**(0.5d0*(dim-1d0-2d0*alphafact))*&
               EIG%evec((i-1)*BPD%xDim + 1, ikeep(beta))/BB%Nbeta(beta)! 
          BB%ZP(i,beta) = -BB%b(beta)*BB%Z(i,beta)
       enddo
       do i = 1, BB%NumOpenR
          BB%Z(i+BB%NumOpenL,beta) = BPD%xr**(0.5d0*(dim-1d0-2d0*alphafact))*&
               EIG%evec((i-1)*BPD%xDim + BPD%xDim,ikeep(beta))/BB%Nbeta(beta) ! 
          BB%ZP(i+BB%NumOpenL,beta) = -BB%b(beta)*BB%Z(i+BB%NumOpenL,beta)
       enddo
    enddo
    
    allocate(ZT(BB%betaMax,BB%betaMax))
    !write(6,*) "Z:"
!    call printmatrix(BB%Z,BB%betaMax,BB%betaMax,6)
!    call printmatrix(BB%ZP,BB%betaMax,BB%betaMax,6)

    allocate(temp0(BB%betaMax,BB%betaMax))
    ZT = BB%Z
    temp0 = 0d0
    temp0 = matmul(transpose(ZT),BB%Z)
    !call dgemm('T','N',BB%betaMax,BB%betaMax,BB%betaMax,1.0d0,ZT,BB%betaMax,BB%Z,BB%betaMax,0.0d0,temp0,BB%betaMax) !

    write(6,*) "Check norm:"
    call printmatrix(temp0,BB%betaMax,BB%betaMax,6)

    allocate(Dikeep(BB%betaMax + BB%NumOpenL))
    allocate(BigZA(BB%NumOpenL + BB%betaMax,BB%NumOpenL+BB%betaMax),BigZB(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax)) ! 
    allocate(Dvec(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax))
    allocate(Dval(BB%NumOpenL+BB%betaMax))

    Dvec=0d0
    BigZA=0d0
    BigZB=0d0
    
    do i=1,BB%NumOpenL
       do beta=1,BB%NumOpenL
          BigZA(i,beta)=BA%Zf(i,beta)
          BigZA(i+BB%NumOpenL,beta)=BA%Zfp(i,beta)
          !write(6,*) "using BA%Z(), BA%ZP= = ",BA%Zf(i,beta), -BA%Zfp(i,beta),BA%Zf(i,beta)*BA%bf(beta)
       enddo
    enddo
      
    do i=1,BB%NumOpenL
       do beta=1,BB%betaMax
          BigZA(i, BB%NumOpenL+beta) = -BB%Z(i,beta)
          BigZA(i+BB%NumOpenL,beta+BB%NumOpenL) = BB%ZP(i,beta)
       enddo
    enddo
    
    do i=1,BB%NumOpenR
       do beta=1,BB%betaMax
          BigZA(i+2*BB%NumOpenL, BB%NumOpenL+beta)=-BB%ZP(i+BB%NumOpenL,beta)
          BigZB(i+2*BB%NumOpenL, BB%NumOpenL+beta)=BB%Z(i+BB%NumOpenL,beta)
       enddo
    enddo
    
    !print*, 'BigZA : '
    !call printmatrix(BigZA, BB%NumOpenL+BB%betaMax, BB%NumOpenL+BB%betaMax,6)
    !print*, 'BigZB : '
    !call printmatrix(BigZB, BB%NumOpenL+BB%betaMax, BB%NumOpenL+BB%betaMax,6)

    call Mydggev(BB%NumOpenL+BB%betaMax,BigZA,BB%NumOpenL+BB%betaMax,BigZB,BB%NumOpenL+BB%betaMax,Dval,Dvec) !
    !print*, 'Eigenvalues of box-matching eigenproblem : '
    !call printmatrix(Dval,BB%NumOpenL+BB%betaMax,1,6)
    
!!$    print*, 'Dvec : '
!!$    do j = 1,BB%betaMax
!!$       write(6,*) (Dvec(j,k), k = 1,BB%betaMax)
!!$    enddo
    j=1
    do i = 1,BB%NumOpenL+BB%betaMax
       if(abs(Dval(i)).ge.1e-12) then
          Dikeep(j)=i
          BB%bf(j)=Dval(i)
          !write(6,*) i, j, BB%bf(j), (Dvec(k,i), k=BB%NumOpenL+1,BB%NumOpenL+BB%betaMax)
          j=j+1
       endif
    enddo
    !print*, 'Final Number of Channels : ',j-1, "should be ",BB%NumOpenR
    do beta=1,BB%NumOpenR
       !for each beta, construct the final state with constant log-derivative at the right boundary
       !write(6,*) "beta value being calculated now is: ", beta, BB%bf(beta)
       do i=1, BB%NumOpenR
          BB%Zf(i,beta) = 0.0d0
          do betaprime = 1,BB%betaMax
             BB%Zf(i,beta) =  BB%Zf(i,beta) + BB%Z(BB%NumOpenL+i,betaprime)*Dvec(BB%NumOpenL+betaprime,Dikeep(beta)) !
             !write(6,*) "adding...",BB%NumOpenL+i, betaprime, BB%Z(BB%NumOpenL+i,betaprime),"*",&
             !     Dvec(BB%NumOpenL+betaprime,Dikeep(beta))
          enddo
          !write(6,*) i,beta,"Zf = ", BB%Zf(i,beta)
       enddo
       ! calculate the normalization of the final state
       NewNorm=0.0d0
       do i=1,BB%NumOpenR
          NewNorm=NewNorm+BB%Zf(i,beta)**2
       enddo
       NewNorm=dsqrt(NewNorm)
       ! normalize and set the derivative Zfp
       do i=1,BB%NumOpenR
          BB%Zf(i,beta) = BB%Zf(i,beta)/NewNorm
          BB%Zfp(i,beta) = -BB%Zf(i,beta)*BB%bf(beta)
       enddo
    enddo
    print*, 'Dvec : '
    call printmatrix(Dvec,BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax,6)
    !print*, 'Zf'
    !call printmatrix(BB%Zf,BB%NumOpenR,BB%NumOpenR,6)

    deallocate(Dikeep,Dvec,Dval,BigZA,BigZB)         
    deallocate(temp0)
    deallocate(tempnorm)
    deallocate(ZT)
  end subroutine BoxMatch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !use this routine to test if the potential is right.  Looks ok!
  SUBROUTINE POTLMX(NCHAN,R,V)
    use BalujaParameters
    IMPLICIT none
    integer I,J,NCHAN,NMX,LAMAX,K,ION,LCHL(3),EL2
    double precision V(NCHAN,NCHAN),R,VP,RR
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
          enddo
          V(I,J) = 0.5d0*VP
       enddo
    enddo

  END SUBROUTINE POTLMX
  !****************************************************************************************************
  subroutine checkbessel(xmin,xmax,npts,lam,d,file)
    implicit none
    double precision xmin, xmax
    double precision lam,hypj,hypy,hypjp,hypyp,hypi,hypk,hypip,hypkp
    double precision, allocatable :: x(:)
    integer d,file
    integer n, npts

    allocate(x(npts))
    do n = 1, npts
       x(n) = xmin + (n-1)*(xmax-xmin)/dble(npts-1)
    enddo
    do n=1,npts
       !call hyperjy(d,lam,x(n),hypj,hypy,hypjp,hypyp)
       call hyperrjry(d,dble(0.5d0*(d-1d0)),lam,x(n),hypj,hypy,hypjp,hypyp)
       !call hyperik(d,lam,x(n),hypj,hypy,hypjp,hypyp)
       !call hyperrirk(d,dble(0.5d0*(d-1d0)),lam,x(n),hypj,hypy,hypjp,hypyp)
       write(file,20) x(n), hypj, hypy, hypjp, hypyp
    enddo
20  format(1P,100D12.4)
  end subroutine checkbessel
  !****************************************************************************************************
    subroutine SetZeroPotential(BPD)
      use DataStructures
      use Quadrature
      implicit none
      type(BPData) BPD
      integer kx,lx,mch,nch
      double precision ax,bx,xScaledZero,pot(3,3)
      double precision xScale(BPD%xNumPoints)

      do kx = 1,BPD%xNumPoints-1
         ax = BPD%xPoints(kx)
         bx = BPD%xPoints(kx+1)
         xScale(kx) = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do lx = 1,LegPoints
            BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
            BPD%Pot(:,:,lx,kx) = 0d0
         enddo
      enddo
    end subroutine SetZeroPotential
    !****************************************************************************************************

    !****************************************************************************************************
