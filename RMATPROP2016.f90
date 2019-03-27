!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module GlobalVars
  implicit none
  integer NumParticles, NumChannels, NumTot, Order, xDim, xNumPoints, Left, Right, xDimMin, n
  double precision reducedmass, alpha, xMin, xMax
  double precision SpatialDim, MultipoleCoupMat(3,3,2), threshold(3)
  double precision, allocatable :: mass(:)
  double precision RMatBaluja1(3,3)                                   ! this is the Baluja et al R-matrix at r=5 for their test case
  double precision RMatBaluja2(3,3)                                   ! this is the Baluja et al R-matrix at r=15 for their test case
  double precision, allocatable :: u(:,:,:), ux(:,:,:)                ! basis sets for the multi-sector domain
  double precision, allocatable :: su(:,:,:), sux(:,:,:)              ! basis sets for the elemental sector in the R-matrix propogator
  double precision, allocatable :: xPoints(:)                         ! the radial (x) grid in each sector
  double precision, allocatable :: xSector(:)                         ! the mesh for the x-sectors
  integer, allocatable :: xBounds(:)
  character*64 InputFile
  complex*16 II
  parameter(II=(0.0d0,1.0d0))
  
contains
  subroutine ReadInputFile
    ! Be sure to match the format of the input file to the format of the read statements below


    open(unit=7,file=InputFile(1:index(InputFile,' ')-1))
    read(7,*)
    read(7,*) NumParticles, NumChannels, NumTot, SpatialDim
    allocate(mass(NumParticles))
    read(7,*)
    read(7,*) (mass(n), n=1,NumParticles)
    
    close(unit=7)
  end subroutine ReadInputFile
  
end module GlobalVars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



program main
  use GlobalVars
  use Quadrature
  implicit none
  double precision, allocatable :: evec(:,:), eval(:)

  !--------------------------------------------------
  ! Set the total number of spatial dimensions
  ! (It's a double precision number for convenience)
  !--------------------------------------------------
  SpatialDim = 3.d0

  !--------------------------------------------------
  ! Read in the quadrature data stored in Quadrature module.
  !--------------------------------------------------
  LegendreFile = 'Legendre.dat'
  LegPoints = 10
  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)
  
  !--------------------------------------------------
  ! Initialize the Baluja et al data and print it to the screen
  !--------------------------------------------------
  call SetMultipoleCoup()
  write(6,*) 'The starting R-Matrix at r=5.d0:'
  call printmatrix(RMatBaluja1,3,3,6)

  !----------------------------------------------------------------------
  ! Read information from the input file
  !----------------------------------------------------------------------
  InputFile = 'RMATPROP.inp'
  call ReadInputFile()

  !--------------------------------------------------
  ! Find the eigenvalues and eigenvectors of the initial R matrix, and print to screen.
  !--------------------------------------------------
  allocate(evec(3,3),eval(3))
  evec = RMatBaluja1
  call Mydsyev(evec,3,eval,evec)
  write(6,*) 'The eigenvalues of the original R-Matrix:'
  call printmatrix(eval,3,1,6)
  write(6,*) 'The eigenvectors of the original R-Matrix:'
  call printmatrix(evec,3,3,6)
  deallocate(evec,eval)


  !------------------------------------------------------------
  ! Set the left/right boundary conditions
  ! Determine and set the dimensionality of the x-basis
  ! Set xMin, xMax
  !------------------------------------------------------------
  xMin = 5.0d0
  xMax = 15.0d0
  xNumPoints = 10
  Order = 5
  Left=2
  Right=2
  xDimMin=xNumPoints+order-3
  xDim=xDimMin
  if (Left .eq. 2) xDim = xDim + 1
  if (Right .eq. 2) xDim = xDim + 1
  print*, '#Left BC  Right BC = ', Left, Right
  print*, '# xDimMin = ',xDimMin,' xDim = ',xDim
  print*, '# xNumPoints = ', xNumPoints
  print*, '# Order = ', Order
  
  allocate(xPoints(xNumPoints))
  allocate(xBounds(xNumPoints+2*Order))   ! allocate xBounds before calling CalcBasisFuncs
  allocate(u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim)) ! allocate memory for basis functions
  
  call GridMaker(xNumPoints,xMin,xMax,xPoints)

  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,0,u)
  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,1,ux)
  
  call RMATPROP()



  deallocate(xBounds,u,ux, xPoints)
20 format(1P,100e14.8)
  
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RMATPROP
  use GlobalVars
  use Quadrature  
  implicit none

  

  
end subroutine RMATPROP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------------------------------------
! This subroutine sets the multipole coupling elements for the test case described in Baluja et al CPC (1982) paper
!----------------------------------------------------------------------------------------------------
subroutine SetMultipoleCoup()
  use GlobalVars
  implicit none

  
  threshold(1)=0.0d0
  threshold(2)=0.0d0
  threshold(3)=2.169d0
  
  MultipoleCoupMat(1,1,1)=0.0d0
  MultipoleCoupMat(1,2,1)=0.0d0
  MultipoleCoupMat(1,3,1)=-0.5682977d0
  MultipoleCoupMat(2,1,1)=0.0d0
  MultipoleCoupMat(2,2,1)=0.0d0
  MultipoleCoupMat(2,3,1)=-0.8036944d0
  MultipoleCoupMat(3,1,1)=-0.5682977d0
  MultipoleCoupMat(3,2,1)=-0.8036944d0
  MultipoleCoupMat(3,3,1)=0.0d0
  
  MultipoleCoupMat(1,1,2)=0.0d0
  MultipoleCoupMat(1,2,2)=-0.5554260d0
  MultipoleCoupMat(1,3,2)=0.0d0
  MultipoleCoupMat(2,1,2)=-0.5554260d0
  MultipoleCoupMat(2,2,2)=-0.3927455d0
  MultipoleCoupMat(2,3,2)=0.0d0
  MultipoleCoupMat(3,1,2)=0.0d0
  MultipoleCoupMat(3,2,2)=0.0d0
  MultipoleCoupMat(3,3,2)=0.0d0

  
  RMatBaluja1(1,1) = 0.145599d0
  RMatBaluja1(1,2) = 0.00559130d0
  RMatBaluja1(1,3) = 0.000514900d0
  RMatBaluja1(2,1) = 0.00559130d0
  RMatBaluja1(2,2) = -0.00811540d0
  RMatBaluja1(2,3) = -0.01735370d0
  RMatBaluja1(3,1) = 0.000514900d0
  RMatBaluja1(3,2) = -0.01735370d0
  RMatBaluja1(3,3) = 0.00829450d0
      
end subroutine SetMultipoleCoup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine printmatrix(M,nr,nc,file)
  implicit none
  integer nr,nc,file,j,k
  double precision M(nr,nc)
  
  do j = 1,nr
     write(file,20) (M(j,k), k = 1,nc)
  enddo

20 format(1P,100e14.5)
end subroutine printmatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$subroutine InitializePropBasis(Left,Right,xPoints,xNumPoints,xBounds)
!!$  use GlobalVars
!!$  use Quadrature
!!$  implicit none
!!$  integer Left, Right, xNumPoints, 
!!$  
!!$  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
!!$       xDim,xBounds,xNumPoints,0,u)
!!$ call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
!!$       xDim,xBounds,xNumPoints,0,u)
!!$  
!!$  
!!$end subroutine InitializePropBasis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GridMaker(xNumPoints,xMin,xMax,xPoints)
  implicit none
  integer xNumPoints
  double precision xMin,xMax,xPoints(xNumPoints)
  
  integer i,j,k
  double precision Pi
  double precision r0New
  double precision xRswitch
  double precision xDelt,x0,x1,x2
  
  
  Pi = 3.1415926535897932385d0
  
  x0 = xMin
  x1 = xMax
  k = 1
  xDelt = (x1-x0)/dble(xNumPoints-1)
  do i = 1,xNumPoints
     !xPoints(k) = x1*((i-1)*xDelt/(x1-x0))**2 + x0  ! use this for a quadratic grid that concentrates more points near x0.
     xPoints(k) = (i-1)*xDelt + x0 ! Simple linear grid
     k = k + 1
  enddo
  
15 format(6(1x,1pd12.5))

end subroutine GridMaker
