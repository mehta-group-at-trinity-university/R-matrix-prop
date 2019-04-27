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
