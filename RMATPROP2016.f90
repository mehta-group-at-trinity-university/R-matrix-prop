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
    end subroutine GetGaussFactors

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
!    use DataStructures
    !    use GlobalVars
    use Quadrature
    implicit none
    integer n
    ! Be sure to match the format of the input file to the format of the read statements below
    open(unit=7,file=InputFile(1:index(InputFile,' ')-1),action='read')
    read(7,*)
    read(7,*) NumParticles, NumChannels, SpatialDim, NumAllChan
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
    AlphaFactor = 0d0!(EffDim-1d0)/2d0
    
    If (NumParticles.eq.2) then
       reducedmass = mass(1)*mass(2)/(mass(1)+mass(2))
    else
       write(6,*) "Reduced mass not set. Must set reduced mass"
       stop
    end If
    Order = 5
    
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
     double precision, allocatable :: Pot(:,:,:,:)
     double precision kl, kr
     integer, allocatable :: xBounds(:), kxMin(:,:), kxMax(:,:)
  end type BPData
  !****************************************************************************************************
  type GenEigVal
     integer MatrixDim
     double precision, allocatable :: Gam(:,:), Lam(:,:), evec(:,:), eval(:)
  end type GenEigVal
  !****************************************************************************************************
  type BoxData
     integer NumOpenL, NumOpenR, NumOpen, betaMax
     double precision xL, xR
     double precision, allocatable :: Z(:,:), ZP(:,:), NBeta(:), Norm(:,:)
     double precision, allocatable :: b(:), bf(:), Zf(:,:)
     double precision, allocatable :: RF(:,:)!, T(:,:), S(:,:), K(:,:)
     double precision, allocatable :: R11(:,:), R12(:,:), R21(:,:), R22(:,:)
  end type BoxData
end module DataStructures
!****************************************************************************************************
!  module MatrixElements
!contains

  subroutine CalcGamLam(BPD,EIG)
    use DataStructures
    use Quadrature
    use GlobalVars
    implicit none
    type(BPData) :: BPD
    type(GenEigVal) :: EIG
    integer ix,ixp,lx,kx,mch,nch
    double precision a, ax, bx, xIntScale(BPD%xNumPoints), TempG, xScaledZero, x1, x2

    
    x1=BPD%xPoints(1)
    x2=BPD%xPoints(BPD%xNumPoints)

    allocate(EIG%Gam(BPD%MatrixDim,BPD%MatrixDim))
    allocate(EIG%Lam(BPD%MatrixDim,BPD%MatrixDim))
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
                ! Ultimately, these matrix elements should be calculated elsewhere since they won't change from sector to sector.
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
          EIG%Lam( (mch-1)*BPD%xDim + BPD%xDim, (mch-1)*BPD%xDim + BPD%xDim ) = x2**(EffDim-1d0-2d0*AlphaFactor)
       enddo
    endif
    if(BPD%Left.eq.2) then
       do mch = 1, NumChannels
          EIG%Lam((mch-1)*BPD%xDim + 1,(mch-1)*BPD%xDim + 1) = -x1**(EffDim-1d0-2d0*AlphaFactor)
       enddo
    endif

  end subroutine CalcGamLam
  !****************************************************************************************************
  !end module MatrixElements
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
      
      RMatBaluja1(1,1) = 0.145599d0
      RMatBaluja1(1,2) = 0.00559130d0
      RMatBaluja1(1,3) = 0.000514900d0

      RMatBaluja1(2,1) = 0.00559130d0
      RMatBaluja1(2,2) = -0.00811540d0
      RMatBaluja1(2,3) = -0.01735370d0

      RMatBaluja1(3,1) = 0.000514900d0
      RMatBaluja1(3,2) = -0.01735370d0
      RMatBaluja1(3,3) = 0.00829450d0


      
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
            call POTLMX(3,BPD%x(lx,kx),pot)
            BPD%Pot(:,:,lx,kx) = pot
!!$            
!!$            do mch=1,BPD%NumChannels
!!$               !BPD%Pot(mch,mch,lx,kx) = BalujaPotential(mch,mch,BPD%x(lx,kx))
!!$               BPD%Pot(mch,mch,lx,kx) = pot(mch,mch)
!!$               do nch=1,mch-1
!!$                  !BPD%Pot(mch,nch,lx,kx) = BalujaPotential(mch,nch,BPD%x(lx,kx))
!!$                  BPD%Pot(mch,nch,lx,kx) = pot(mch,nch)!BalujaPotential(mch,nch,BPD%x(lx,kx))
!!$                  BPD%Pot(nch,mch,lx,kx) = BPD%Pot(mch,nch,lx,kx) ! Potential is symmetric
!!$               enddo
!!$            enddo
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
         BalujaPotential = BalujaPotential - BalujaMCmat(mch,nch,Lam)*R**(-Lam-1d0)/(2.d0*reducedmass)
      enddo
      if(mch.eq.nch) then
         BalujaPotential = BalujaPotential + &
              lmom(mch)*(lmom(mch)+1)/(2d0*reducedmass*R**2) - Znet/R + BalujaEth(mch)/(2.d0*reducedmass) !
      end if
      
      return
      
    end function BalujaPotential
    !****************************************************************************************************
  end module BalujaParameters
  !****************************************************************************************************

!****************************************************************************************************
!****************************************************************************************************
program main
  use DataStructures
  use GlobalVars
  use BalujaParameters
  use Quadrature
!  use MatrixElements
  implicit none
  type(BPData) BPD
  type(GenEigVal) EIG
  type(BoxData) BA, BB
  double precision, allocatable :: evec(:,:), eval(:)!, temp0(:,:)

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

  ! Intitializes some BPD variables to the input values.
  BPD%NumChannels = NumChannels
  BPD%Order = Order
  BPD%Left = StartBC
  BPD%Right = EndBC
  BPD%xNumPoints=xTotNumPoints
  BPD%kl = kStart ! only relevant for Left = 3. This is the normal log-derivative at x1
  BPD%kr = kEnd ! only relevant for Right = 3. This is the normal log-derivative at x2
  
  allocate(BPD%xPoints(BPD%xNumPoints))
  allocate(BPD%xBounds(BPD%xNumPoints + 2*BPD%Order))   ! allocate xBounds before calling CalcBasisFuncs
  allocate(BPD%x(LegPoints,BPD%xNumPoints-1))

  call GridMakerLinear(BPD%xNumPoints,xStart,xEnd,BPD%xPoints)
  call MakeBPData(BPD)

  write(6,"(A,T25,2I3)") 'Left BC...', BPD%Left
  write(6,"(A,T25,2I3)") 'Right BC...', BPD%Right
  write(6,"(A,T25,2I3)") 'xNumPoints....', BPD%xNumPoints
  write(6,"(A,T25,2I3)") 'Order....', BPD%Order
  write(6,"(A,T25,F3.1)") 'EffDim...', EffDim
  write(6,"(A,T25,F3.1)") 'Energy...', Energy
  write(6,"(A,T20,F8.4)") 'Reduced mass...', reducedmass
  write(6,"(A,T20,100F8.4)") 'xPoints...', BPD%xPoints
  write(6,"(A,T25,2I3)") 'xDim....', BPD%xDim
  write(6,"(A,T20,I8)") 'MatrixDim...', BPD%MatrixDim      


  !--------------------------------------------------
  ! Initialize the Baluja et al data and print it to the screen
  !--------------------------------------------------
  !call SetBalujaPotential(BPD)
  call SetZeroPotential(BPD)
  !----------------------------------------------------------------------------------------------------
  ! Comment/uncomment the next line if you want to print the basis to file fort.300
  !----------------------------------------------------------------------------------------------------
  call CheckBasisBP(BPD%xDim,BPD%xNumPoints,BPD%xBounds,BPD%Order,300,LegPoints,BPD%u,BPD%x)

  write(6,*) '--------------------------------'
  write(6,*) 'The starting R-Matrix at r=5.d0:'
  call printmatrix(RMatBaluja1,3,3,6)
  write(6,*) '--------------------------------'

  EIG%MatrixDim=BPD%MatrixDim
  !allocate(temp0(BPD%MatrixDim,BPD%MatrixDim))
  allocate(EIG%evec(EIG%MatrixDim,EIG%MatrixDim),EIG%eval(EIG%MatrixDim))
  call CalcGamLam(BPD,EIG)
  call Mydggev(EIG%MatrixDim,EIG%Gam,EIG%MatrixDim,EIG%Lam,EIG%MatrixDim,EIG%eval,EIG%evec)

  BA%NumOpenL = 0
  BA%NumOpenR = 3
  BA%NumOpen = 3
  BA%betaMax = BA%NumOpen
  
  BB%NumOpenL = 0
  BB%NumOpenR = 3
  BB%NumOpen = 3
  BB%betaMax = 3


  !BA%RF = -xStart*RMatBaluja1
  allocate(BA%RF(BB%NumOpenR,BA%NumOpenR))
  allocate(BA%b(BA%betaMax),BA%bf(BA%NumOpenR))
  allocate(BA%Z(BA%betaMax,BA%betaMax),BA%Zf(BA%NumOpenR,BA%NumOpenR))
  allocate(BB%RF(BB%NumOpenR,BB%NumOpenR),BB%Z(BB%NumOpen,BB%NumOpen),BB%ZP(BB%NumOpen,BB%NumOpen))
  allocate(BB%Zf(BB%NumOpenR,BB%NumOpenR),BB%bf(BB%NumOpenR),BB%b(BB%NumOpen))

  !----------------------------------------------------------------------------------------------------
  ! use this part to test the free particle solutions
  ! be sure to change the values of NumOpenL and NumOpenR above
  call BoxMatch(BA, BB, BPD, EIG)
  write(6,*) "Final R-matrix:"
  call printmatrix(BB%RF,BB%NumOpenR,BB%NumOpenR,6)
  
  !----------------------------------------------------------------------------------------------------
  ! be sure to change the values of NumOpenL and NumOpenR above
!!$  !use this part for baluja algorithm
!!$  BA%RF = RMatBaluja1  
!!$  call RProp(BA, BB, BPD, EIG)
!!$  write(6,*) "Final R-matrix:"
!!$  call printmatrix(BB%RF,BB%NumOpenR,BB%NumOpenR,6)



  !----------------------------------------------------------------------------------------------------
  ! be sure to change the values of NumOpenL and NumOpenR above
  ! use this part for the box-matching algorithm
!!$  BA%RF = xStart*RMatBaluja1
!!$  !BA%RF = RMatBaluja1  
!!$  call MYDSYEV(BA%RF,BA%betaMax,BA%b,BA%Z)
!!$  BA%Zf=BA%Z
!!$
!!$  BA%b = -1.d0/BA%b
!!$  BA%bf = BA%b
!!$  call BoxMatch(BA, BB, BPD, EIG)
!!$  write(6,*) "Final R-matrix:"
!!$  call printmatrix(BB%RF/xEnd,BB%NumOpenR,BB%NumOpenR,6)
!!$  !call printmatrix(BB%RF,BB%NumOpenR,BB%NumOpenR,6)
!!$
!!$  write(6,*) "The Exact R-Matrix From Baluja et al at R=15.0 is:"
!!$  call printmatrix(RMatBaluja2,3,3,6)
!!$  if(NumSectors.eq.1) then
!!$     call RMATPROP(Left,Right,xStart,xEnd)
!!$  end if


  call checkbessel(0.0001d0,10d0,100,1d0,3,100)
  

  

20 format(1P,100e14.8)
end program main
!****************************************************************************************************
!****************************************************************************************************
subroutine printmatrix(M,nr,nc,file)
  implicit none
  integer nr,nc,file,j,k
  double precision M(nr,nc)
  
  do j = 1,nr
     write(file,*) (M(j,k), k = 1,nc)
  enddo

20 format(1P,100D12.4)
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
subroutine MakeBPData(BPD)
  use Quadrature
  use DataStructures
  implicit none
  
  type(BPData) BPD
  integer xDimMin, ix, ixp
  
  !Determine the basis dimension for these boundary conditions
  xDimMin=BPD%xNumPoints+BPD%Order-3
  BPD%xDim=xDimMin
  if (BPD%Left .eq. 2) BPD%xDim = BPD%xDim + 1
  if (BPD%Right .eq. 2) BPD%xDim = BPD%xDim + 1
  BPD%MatrixDim = BPD%NumChannels*BPD%xDim
  
  allocate(BPD%u(LegPoints,BPD%xNumPoints,BPD%xDim),BPD%ux(LegPoints,BPD%xNumPoints,BPD%xDim)) ! allocate memory for basis functions
  allocate(BPD%kxMin(BPD%xDim,BPD%xDim),BPD%kxMax(BPD%xDim,BPD%xDim)) !
  allocate(BPD%Pot(BPD%NumChannels,BPD%NumChannels,LegPoints,BPD%xNumPoints-1))
  
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
  
end subroutine MakeBPData
!****************************************************************************************************
  subroutine RProp(BA, BB, BPD, EIG)
    use DataStructures
    use GlobalVars
    implicit none
    type(BPData), intent(in) :: BPD
    type(GenEigVal), intent(in) :: EIG
    type(BoxData), intent(in) :: BA
    type(BoxData) BB
    double precision, allocatable :: tempnorm(:,:), temp1(:,:), temp2(:,:), temp0(:,:), ZT(:,:),tempR(:,:)
    double precision x1, x2
    integer i,j,beta,betaprime,nch,mch!,betamax
    integer, allocatable :: ikeep(:)
    
    allocate(ikeep(BPD%MatrixDim))
!    call printmatrix(EIG%eval,BPD%MatrixDim,1,6)
    j=1
    ikeep = 0
    do i = 1, BPD%MatrixDim
       if(abs(EIG%eval(i)).ge.1e-12) then
          write(6,"(A,T20,I5,T30,E12.6)") 'eval',i, EIG%eval(i)
          ikeep(j)=i
          j = j+1
       endif
    enddo
    BB%betaMax=j-1
    allocate(BB%NBeta(BB%betaMax),BB%Norm(BB%betaMax,BB%betaMax))
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

    write(6,*) "Norm matrix:"
    call printmatrix(BB%Norm,BB%betaMax,BB%betaMax,6)

    x1=BPD%xPoints(1)
    x2=BPD%xPoints(BPD%xNumPoints)

    do beta = 1,BB%betaMax
       do i = 1,BB%NumOpenL
          BB%Z(i,beta) = EIG%evec((i-1)*BPD%xDim + 1, ikeep(beta))/BB%Nbeta(beta)! 
          BB%ZP(i,beta) = EIG%eval(ikeep(beta))*BB%Z(i,beta)
       enddo
       do i = 1, BB%NumOpenR
          BB%Z(i+BB%NumOpenL,beta) = EIG%evec((i-1)*BPD%xDim + BPD%xDim,ikeep(beta))/BB%Nbeta(beta) ! 
          BB%ZP(i+BB%NumOpenL,beta) = -EIG%eval(ikeep(beta))*BB%Z(i+BB%NumOpenL,beta)
       enddo
    enddo
    
    allocate(ZT(BB%betaMax,BB%betaMax))
    write(6,*) "Z:"
    call printmatrix(BB%Z,BB%betaMax,BB%betaMax,6)
!    call printmatrix(BB%ZP,BB%betaMax,BB%betaMax,6)

    allocate(temp0(BB%betaMax,BB%betaMax))
    ZT = BB%Z
    temp0 = 0d0
    call dgemm('T','N',BB%betaMax,BB%betaMax,BB%betaMax,1.0d0,ZT,BB%betaMax,BB%Z,BB%betaMax,0.0d0,temp0,BB%betaMax) !

    write(6,*) "Check norm:"
    call printmatrix(temp0,BB%betaMax,BB%betaMax,6)

    if(BB%NumOpenL.gt.0) then
       allocate(BB%R11(BB%NumOpenL,BB%NumOpenL),BB%R12(BB%NumOpenL,BB%NumOpenR))
       allocate(BB%R21(BB%NumOpenR,BB%NumOpenL),BB%R22(BB%NumOpenR,BB%NumOpenR))
       allocate(temp1(BB%NumOpenL,BB%NumOpenL),temp2(BB%NumOpenR,BB%NumOpenL))
       BB%R11=0d0
       BB%R12=0d0
       BB%R21=0d0
       BB%R22=0d0
       BB%RF=0d0

       do beta=1,BB%betaMax
          do i=1,BB%NumOpenL
             do j=1,BB%NumOpenL
                BB%R11(i,j) = BB%R11(i,j) - BB%Z(i,beta)*BB%Z(j,beta)/EIG%eval(ikeep(beta))
             enddo
          enddo
          do i=1,BB%NumOpenL
             do j=1,BB%NumOpenR
                BB%R12(i,j) = BB%R12(i,j) - BB%Z(i,beta)*BB%Z(j + BB%NumOpenL,beta)/EIG%eval(ikeep(beta)) ! 
             enddo
          enddo
          do i=1,BB%NumOpenR
             do j=1,BB%NumOpenL
                BB%R21(i,j) = BB%R21(i,j) - BB%Z(i + BB%NumOpenL,beta)*BB%Z(j,beta)/EIG%eval(ikeep(beta)) ! 
             enddo
          enddo
          do i=1,BB%NumOpenR
             do j=1,BB%NumOpenR
                BB%R22(i,j) = BB%R22(i,j) - BB%Z(i + BB%NumOpenL,beta)*BB%Z(j + BB%NumOpenL,beta)/EIG%eval(ikeep(beta)) ! 
             enddo
          enddo
       enddo
       
       !Flip the sign to account for the minus sign convention of Baluja
!!$       BB%R11 = -BB%R11
!!$       BB%R12 = -BB%R12
!!$       BB%R21 = -BB%R21
!!$       BB%R22 = -BB%R22
             
       temp1 = BB%R11 + x1*BA%RF
       
       allocate(tempR(BB%NumOpenR,BB%NumOpenR))
       call SqrMatInv(temp1,BB%NumOpenL)
       call dgemm('N','N',BB%NumOpenR,BB%NumOpenL,BB%NumOpenL,1.0d0,BB%R21,BB%NumOpenR,temp1,BB%NumOpenL,0.0d0,temp2,BB%NumOpenR) ! 
       call dgemm('N','N',BB%NumOpenR,BB%NumOpenR,BB%NumOpenL,1.0d0,temp2,BB%NumOpenR,BB%R12,BB%NumOpenL,0.0d0,tempR,BB%NumOpenR) !


       BB%RF = (BB%R22 - tempR)/x2
       deallocate(BB%R11,BB%R12,BB%R21,BB%R22)
       deallocate(temp1,temp2)
    else
       temp0 = BB%ZP
       call SqrMatInv(temp0,BB%betaMax)
       call dgemm('N','N',BB%betaMax,BB%betaMax,BB%betaMax,1.0d0,BB%Z,BB%betaMax,temp0,BB%betaMax,0.0d0,BB%RF,BB%betaMax) ! 
    end if

    deallocate(temp0)

    deallocate(ikeep)
    deallocate(tempnorm)
    deallocate(BB%NBeta,BB%Norm)
    deallocate(ZT)
  end subroutine RProp

  
!****************************************************************************************************
  subroutine BoxMatch(BA, BB, BPD, EIG)
    use DataStructures
    use GlobalVars
    implicit none
    type(BPData), intent(in) :: BPD
    type(GenEigVal), intent(in) :: EIG
    type(BoxData), intent(in) :: BA
    type(BoxData) BB
    double precision, allocatable :: tempnorm(:,:), temp1(:,:), temp2(:,:), temp0(:,:), ZT(:,:),tempR(:,:)
    double precision, allocatable :: BigZA(:,:), BigZB(:,:),Dvec(:,:),bfinal(:)
    integer, allocatable :: Dikeep(:)
    double precision x1, x2, NewNorm
    integer i,j,k,beta,betaprime,nch,mch!,betamax
    integer, allocatable :: ikeep(:)
    
    allocate(ikeep(BPD%MatrixDim))
!    call printmatrix(EIG%eval,BPD%MatrixDim,1,6)
    j=1
    ikeep = 0
    do i = 1, BPD%MatrixDim
       if(abs(EIG%eval(i)).ge.1e-12) then
          write(6,"(A,T20,I5,T30,E12.6)") 'eval',i, EIG%eval(i)
          ikeep(j)=i
          j = j+1
       endif
    enddo
    BB%betaMax=j-1
    allocate(BB%NBeta(BB%betaMax),BB%Norm(BB%betaMax,BB%betaMax))
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

    write(6,*) "Norm matrix:"
    call printmatrix(BB%Norm,BB%betaMax,BB%betaMax,6)

    x1=BPD%xPoints(1)
    x2=BPD%xPoints(BPD%xNumPoints)

    do beta = 1,BB%betaMax
       do i = 1,BB%NumOpenL
          BB%Z(i,beta) = EIG%evec((i-1)*BPD%xDim + 1, ikeep(beta))/BB%Nbeta(beta)! 
          BB%ZP(i,beta) = EIG%eval(ikeep(beta))*BB%Z(i,beta)
       enddo
       do i = 1, BB%NumOpenR
          BB%Z(i+BB%NumOpenL,beta) = EIG%evec((i-1)*BPD%xDim + BPD%xDim,ikeep(beta))/BB%Nbeta(beta) ! 
          BB%ZP(i+BB%NumOpenL,beta) = -EIG%eval(ikeep(beta))*BB%Z(i+BB%NumOpenL,beta)
       enddo
    enddo
    
    allocate(ZT(BB%betaMax,BB%betaMax))
    write(6,*) "Z:"
    call printmatrix(BB%Z,BB%betaMax,BB%betaMax,6)
!    call printmatrix(BB%ZP,BB%betaMax,BB%betaMax,6)

    allocate(temp0(BB%betaMax,BB%betaMax))
    ZT = BB%Z
    temp0 = 0d0
    call dgemm('T','N',BB%betaMax,BB%betaMax,BB%betaMax,1.0d0,ZT,BB%betaMax,BB%Z,BB%betaMax,0.0d0,temp0,BB%betaMax) !

    write(6,*) "Check norm:"
    call printmatrix(temp0,BB%betaMax,BB%betaMax,6)

    allocate(Dikeep(BB%betaMax + BB%NumOpenL))
    allocate(BigZA(BB%NumOpenL + BB%betaMax,BB%NumOpenL+BB%betaMax),BigZB(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax)) ! 
    allocate(Dvec(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax),bfinal(BB%NumOpenL+BB%betaMax)) ! 
    BigZA=0d0
    BigZB=0d0
    
    do i=1,BB%NumOpenL
       do beta=1,BB%NumOpenL
          BigZA(i,beta)=BA%Zf(i,beta)
       enddo
       do beta=1,BB%NumOpenL
          BigZA(i+BB%NumOpenL,beta)=BA%Zf(i,beta)*BA%bf(beta) ! 
       enddo
    enddo
    
    
    do i=1,BB%NumOpenL
       do beta=1,BB%betaMax
          BigZA(i, BB%NumOpenL+beta) = -BB%Z(i,beta)
          BigZA(i+BB%NumOpenL,beta+BB%NumOpenL) = BB%Z(i,beta)*EIG%eval(ikeep(beta))
       enddo
    enddo
    
    do i=1,BB%NumOpenR
       do beta=1,BB%betaMax
          BigZA(i+2*BB%NumOpenL, BB%NumOpenL+beta)=BB%Z(i+BB%NumOpenL,beta)*EIG%eval(ikeep(beta))
          BigZB(i+2*BB%NumOpenL, BB%NumOpenL+beta)=BB%Z(i+BB%NumOpenL,beta)
       enddo
    enddo
    
    print*, 'BigZA : '
    call printmatrix(BigZA, BB%NumOpenL+BB%betaMax, BB%NumOpenL+BB%betaMax,6)
    print*, 'BigZB : '
    call printmatrix(BigZB, BB%NumOpenL+BB%betaMax, BB%NumOpenL+BB%betaMax,6)

    call Mydggev(BB%NumOpenL+BB%betaMax,BigZA,BB%NumOpenL+BB%betaMax,BigZB,BB%NumOpenL+BB%betaMax,BB%b,Dvec) !
    print*, 'b : '
    call printmatrix(BB%b,BB%NumOpenL+BB%betaMax,1,6)
    
    print*, 'Dvec : '
    do j = 1,BB%betaMax
       write(6,*) (Dvec(j,k), k = 1,BB%betaMax)
    enddo
    print*, 'Dvec : '
    call printmatrix(Dvec,BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax,6)
    print*, 'Dvec 1: '
    call printmatrix(Dvec(:,1),BB%NumOpenL+BB%betaMax,1,6)
    print*, 'Dvec 2: '
    call printmatrix(Dvec(:,2),BB%NumOpenL+BB%betaMax,1,6)
    print*, 'Dvec 3: '
    call printmatrix(Dvec(:,3),BB%NumOpenL+BB%betaMax,1,6)
    write(6,*) "diagonal of D", Dvec(1,1), Dvec(2,2), Dvec(3,3)

    j=1
    do i = 1,BB%NumOpenL+BB%betaMax
       if(abs(BB%b(i)).ge.1e-12) then
          Dikeep(j)=i
          BB%bf(j)=BB%b(i)
          write(6,*) i, j, BB%bf(j), (Dvec(k,i), k=BB%NumOpenL+1,BB%NumOpenL+BB%betaMax)
          j=j+1
       endif
    enddo
    print*, 'Final Number of Channels : ',j-1
    do beta=1,BB%NumOpenR
       !for each beta, construct the final state with constant log-derivative at the right boundary
       write(6,*) "beta value being calculated now is: ", beta, BB%bf(beta)
       do i=1, BB%NumOpenR
          BB%Zf(i,beta) = 0.0d0
          do betaprime = 1,BB%betaMax

             BB%Zf(i,beta) =  BB%Zf(i,beta) + BB%Z(BB%NumOpenL+i,betaprime)*Dvec(BB%NumOpenL+betaprime,Dikeep(beta)) !
             write(6,*) "adding...",BB%NumOpenL+i, betaprime, BB%Z(BB%NumOpenL+i,betaprime),"*",&
                  Dvec(BB%NumOpenL+betaprime,Dikeep(beta))
          enddo
          write(6,*) i,beta,"Zf = ", BB%Zf(i,beta)
       enddo
       ! calculate the normalization of the final state
       NewNorm=0.0d0
       do i=1,BB%NumOpenR
          NewNorm=NewNorm+BB%Zf(i,beta)**2
       enddo
       NewNorm=dsqrt(NewNorm)
       
       ! normalize
       do i=1,BB%NumOpenR
          BB%Zf(i,beta) = BB%Zf(i,beta)/NewNorm

       enddo
    enddo
    print*, 'Dvec : '
    call printmatrix(Dvec,BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax,6)
    print*, 'Zf'
    call printmatrix(BB%Zf,BB%NumOpenR,BB%NumOpenR,6)
    deallocate(Dikeep,Dvec,BigZA,BigZB)         
    
    do i=1,BB%NumOpenR
         do j=1,BB%NumOpenR
            BB%RF(i,j)=0.0d0
            do beta = 1, BB%NumOpenR
               BB%RF(i,j) = BB%RF(i,j) - BB%Zf(i,beta)*BB%Zf(j,beta)/BB%bf(beta) ! 
            enddo
         enddo
      enddo

    deallocate(temp0)

    deallocate(ikeep)
    deallocate(tempnorm)
    deallocate(BB%NBeta,BB%Norm)
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
