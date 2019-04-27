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
