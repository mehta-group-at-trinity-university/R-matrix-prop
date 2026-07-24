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
  !--------------------------------------------------------------------
  ! Benchmark against Baluja, Burke & Morgan, Comput. Phys. Commun. 27
  ! (1982) 299-307, section 6 "Description of the test run": propagate
  ! the literature R-matrix RMatBaluja1 (given at r=5, from Berrington
  ! et al.) outward to r=15 and compare against RMatBaluja2 (the paper's
  ! own "R-MATRIX AT END OF OUTWARD PROPAGATION"). This is an outward-
  ! propagation-only check -- it does not exercise the regularity
  ! boundary condition at r=0, only CalcGamLam/BoxMatch/CalcK.
  !--------------------------------------------------------------------
  USE DataStructures
  USE GlobalVars
  USE BalujaParameters
  USE Quadrature
  USE scattering

  IMPLICIT NONE
  TYPE(BPData) BPD,BPD0
  TYPE(GenEigVal) EIG
  TYPE(BoxData) Bseed
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: EvecSeed(:,:), EvalSeed(:)
  DOUBLE PRECISION, ALLOCATABLE :: xprim(:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION xDelt
  INTEGER NumRealBoxes, beta, i, iBox,lx,kx, NumOpenChannels
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

  CALL SetMultipoleCoup()  ! fills lmom, BalujaEth, BalujaMCmat, Znet, RMatBaluja1/2

  !---------------------------------------------------------------------
  ! Fixed open-channel count for this run (mirrors CABA.EXT.f's
  ! NumOpenChannels: chosen once for the whole energy, not per-box).
  ! Channels are assumed ordered open-before-closed, same convention as
  ! CABA. k_1^2 = 2*Energy is the total energy relative to channel 1.
  !---------------------------------------------------------------------
  NumOpenChannels = 0
  DO i = 1,NumChannels
     IF (2d0*Energy.GE.BalujaEth(i)) NumOpenChannels = NumOpenChannels+1
  ENDDO
  WRITE(6,*) "NumOpenChannels = ",NumOpenChannels," (of ",NumChannels,") at k_1^2 = ",2d0*Energy

  !---------------------------------------------------------------------
  ! Seed a zero-width "box" at r=xStart directly from the literature
  ! R-matrix RMatBaluja1, instead of deriving it from a regularity
  ! condition at r=0. R = -Zf.diag(1/bf).Zf^T with Zf orthonormal is
  ! exactly the symmetric eigendecomposition of R (matches how CalcK
  ! builds SD%R from Zf/bf), so diagonalizing RMatBaluja1 gives Zf/bf/Zfp
  ! directly.
  !
  ! Convention mismatch: Baluja/Burke/Morgan (1982) Eq. (2) defines the
  ! R-matrix via u_i(a) = Sum_j R_ij(a)*[a*du_j/dr]_{r=a} -- note the
  ! explicit factor of "a" multiplying the derivative. CalcK's SD%R
  ! (following Burke's later thesis convention, Eq. 3.12 there) has no
  ! such factor: u_i(a) = Sum_j R_ij(a)*[du_j/dr]_{r=a}. So
  ! R_code(a) = a * R_paper(a); scale RMatBaluja1 by xStart on the way in
  ! and RMatBaluja2 by xEnd for the comparison below.
  !---------------------------------------------------------------------
  Bseed%NumOpenL = 0
  Bseed%NumOpenR = NumChannels
  CALL AllocateBox(Bseed)
  CALL InitZeroBox(Bseed)

  ALLOCATE(EvecSeed(NumChannels,NumChannels),EvalSeed(NumChannels))
  CALL MyDSYEV(xStart*RMatBaluja1,NumChannels,EvalSeed,EvecSeed)
  DO beta = 1,NumChannels
     Bseed%bf(beta) = -1.0d0/EvalSeed(beta)
     DO i = 1,NumChannels
        Bseed%Zf(i,beta) = EvecSeed(i,beta)
        Bseed%Zfp(i,beta) = -Bseed%Zf(i,beta)*Bseed%bf(beta)
     ENDDO
  ENDDO
  DEALLOCATE(EvecSeed,EvalSeed)

  !---------------------------------------------------------------------
  ! Real boxes spanning xStart (=5) to xEnd (=15)
  !---------------------------------------------------------------------
  NumRealBoxes = NumBoxes - 1
  ALLOCATE(Boxes(NumRealBoxes))
  xDelt = (xEnd-xStart)/DBLE(NumRealBoxes)
  DO i = 1,NumRealBoxes
     Boxes(i)%NumOpenL = NumChannels
     IF (i.EQ.NumRealBoxes) THEN
        Boxes(i)%NumOpenR = NumOpenChannels  ! final box: only physically-open channels retained
     ELSE
        Boxes(i)%NumOpenR = NumChannels      ! interior box-to-box boundary: all channels retained
     ENDIF
     IF (i.EQ.1) THEN
        Boxes(i)%xl = xStart
     ELSE
        Boxes(i)%xl = Boxes(i-1)%xr
     ENDIF
     Boxes(i)%xr = xStart + i*xDelt
     CALL AllocateBox(Boxes(i))
     CALL InitZeroBox(Boxes(i))
     WRITE(6,*) "Box",i,Boxes(i)%xl,Boxes(i)%xr
  ENDDO

  !---------------------------------------------------------------------
  ! "Primitive" basis set defined over [0,1], rescaled per box (avoids
  ! recomputing B-spline values from scratch for every equal-width box).
  !---------------------------------------------------------------------
  BPD0%NumChannels = NumChannels
  BPD0%Order = Order
  BPD0%Left = 2
  BPD0%Right = 2
  BPD0%xNumPoints = xTotNumPoints
  BPD0%kl = kStart
  BPD0%kr = kEnd
  CALL AllocateBPD(BPD0)
  ALLOCATE(xprim(xTotNumPoints))
  CALL GridMakerLinear(BPD0%xNumPoints,0d0,1d0,xprim)
  BPD0%xPoints = xprim
  CALL MakeBasis(BPD0)
  BPD0%lam(1:NumChannels) = lmom(1:NumChannels)

  BPD = BPD0  ! workhorse set actually used for propagation

  EIG%MatrixDim = BPD%MatrixDim
  CALL AllocateEIG(EIG)
  WRITE(6,*) "EIG%MatrixDim = ",EIG%MatrixDim

  CALL AllocateScat(SD,NumChannels)

  DO iBox = 1,NumRealBoxes
     BPD%xl = Boxes(iBox)%xl
     BPD%xr = Boxes(iBox)%xr
     BPD%xPoints = BPD%xl + (BPD%xr-BPD%xl)*xprim
     DO kx = 1,BPD%xNumPoints
        DO lx = 1,LegPoints
           BPD%ux(lx,kx,1:BPD%xDim) = BPD0%ux(lx,kx,1:BPD%xDim)/(BPD%xr-BPD%xl)
        ENDDO
     ENDDO
     CALL SetBalujaPotential(BPD)
     CALL CalcGamLam(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
     ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
     ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
     CALL PartitionAndEliminate(BPD,EIG,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR,evalRed,evecRed,LamooDiag)
     IF (iBox.EQ.1) THEN
        CALL BoxMatch(Bseed, Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
     ELSE
        CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPD, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
     ENDIF
     DEALLOCATE(evalRed,evecRed,LamooDiag)
  ENDDO

  CALL CalcK(Boxes(NumRealBoxes),BPD,SD,reducedmass,EffDim,AlphaFactor,Energy,BalujaEth)

  WRITE(6,*)
  WRITE(6,*) "k_1^2 = ",2d0*Energy," Ry   (Energy input = ",Energy,")"
  WRITE(6,*) "Thresholds (Ry):", BalujaEth(1:NumChannels)
  WRITE(6,*)
  WRITE(6,*) "Computed R-matrix at r=",Boxes(NumRealBoxes)%xr," (RMATPROP2016, full ",NumChannels,"x",NumChannels,"):"
  CALL printmatrix(SD%R,NumChannels,NumChannels,6)

  IF (2d0*Energy.LT.2.5d0) THEN
     WRITE(6,*) "(No literature reference R-matrix at this energy -- 2.5 Ry is the only Baluja/Burke/Morgan test point.)"
  ELSE
     WRITE(6,*) "Reference xEnd*RMatBaluja2 (Baluja/Burke/Morgan convention, u=a*R*u'):"
     CALL printmatrix(xEnd*RMatBaluja2,NumChannels,NumChannels,6)
     WRITE(6,*) "Difference:"
     CALL printmatrix(SD%R-xEnd*RMatBaluja2,NumChannels,NumChannels,6)
  ENDIF

  WRITE(6,*)
  WRITE(6,*) "Number of open channels = ",SIZE(SD%K,1)
  WRITE(6,*) "K-matrix (open channels only):"
  CALL printmatrix(SD%K,SIZE(SD%K,1),SIZE(SD%K,2),6)
  WRITE(6,*) "S-matrix, real part:"
  CALL printmatrix(DBLE(REAL(SD%S)),SIZE(SD%S,1),SIZE(SD%S,2),6)
  WRITE(6,*) "S-matrix, imag part:"
  CALL printmatrix(DBLE(AIMAG(SD%S)),SIZE(SD%S,1),SIZE(SD%S,2),6)

20 FORMAT(1P,100e14.8)
END PROGRAM main
!****************************************************************************************************

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

!C      SETS UP EXAMPLE POTENTIAL FOR USE IN RPROP

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
    double precision etest,mutest,ethtest,ktest
    DOUBLE PRECISION, ALLOCATABLE :: x(:)
    INTEGER d,file
    INTEGER n, npts

    mutest=0.5D0
    etest=2.5d0
    ethtest=0.0d0
    ktest = dsqrt(2d0*mutest*(etest-ethtest))
    ALLOCATE(x(npts))
    DO n = 1, npts
       x(n) = xmin + (n-1)*(xmax-xmin)/DBLE(npts-1)
    ENDDO
    DO n=1,npts
       !call hyperjy(d,lam,x(n),hypj,hypy,hypjp,hypyp)
       CALL hyperrjry(d,DBLE(0.5d0*(d-1d0)),lam,ktest*x(n),hypj,hypy,hypjp,hypyp)
       !call hyperik(d,lam,x(n),hypj,hypy,hypjp,hypyp)
       !call hyperrirk(d,dble(0.5d0*(d-1d0)),lam,x(n),hypj,hypy,hypjp,hypyp)
       WRITE(file,*) x(n), hypj, hypy, hypjp, hypyp
    ENDDO
20  FORMAT(1P,100E12.4)
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
