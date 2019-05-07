MODULE DipoleDipole
  USE Quadrature
  USE DataStructures

  TYPE DPData
     INTEGER lmax,ml
     DOUBLE PRECISION, ALLOCATABLE :: cllp(:,:,:)
     DOUBLE PRECISION, ALLOCATABLE :: Eth(:)
     LOGICAL even
  END TYPE DPData

CONTAINS
  SUBROUTINE AllocateDP(DP,N)
    IMPLICIT NONE
    integer N ! Numchannels
    TYPE(DPData) DP
    ALLOCATE(DP%cllp(0:DP%lmax,0:DP%lmax,0:DP%lmax))
    ALLOCATE(DP%Eth(N))
  END SUBROUTINE AllocateDP

  SUBROUTINE MakeDipoleDipoleCouplingMatrix(DP)
    IMPLICIT NONE
    TYPE(DPData) :: DP
    INTEGER l,lp,m
    DOUBLE PRECISION, EXTERNAL :: THRJ
    DOUBLE PRECISION tj1,tj2,phase,prefact

    DP%Eth=0d0
    DP%cllp=0d0
    DO l=0,DP%lmax
       DO lp=MAX(l-2,0),MIN(DP%lmax,l+2),2

          prefact=DBLE((2*l+1)*(2*lp+1))
          prefact=dsqrt(prefact)
          DO m = 0,l
             phase=(-1)**m
             tj1 = THRJ(2*l,2*2,2*lp,-2*m,0,2*m)
             tj2 = THRJ(2*l,2*2,2*lp,0,0,0)
             DP%cllp(l,lp,m)=prefact*phase*tj1*tj2
             DP%cllp(lp,l,m)=DP%cllp(l,lp,m)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE MakeDipoleDipoleCouplingMatrix

  SUBROUTINE SetDipoleDipolePot(BPD,DP)

    IMPLICIT NONE
    TYPE(BPData) BPD
    TYPE(DPData) DP
    INTEGER m,l,lp
    INTEGER kx,lx,mch,nch,NChan
    DOUBLE PRECISION ax,bx,xScaledZero
    DOUBLE PRECISION xScale(BPD%xNumPoints)

    NChan = BPD%NumChannels
    ! IF(MOD(DP%lmax,2).EQ.0) THEN
    !    Nchan=DP%lmax/2
    ! ELSE
    !    Nchan=(DP%lmax-1)/2
    ! ENDIF

    BPD%Pot(:,:,:,:) = 0d0

    DO kx = 1,BPD%xNumPoints-1
       ax = BPD%xPoints(kx)
       bx = BPD%xPoints(kx+1)
       xScale(kx) = 0.5d0*(bx-ax)
       xScaledZero = 0.5d0*(bx+ax)
       DO lx = 1,LegPoints
          BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
          DO mch = 1,NChan
             DO nch = 1, mch
                IF(DP%even) THEN
                   l=2*(mch-1) ! l = 0, 2, 4, ...
                   lp=2*(nch-1) ! l' = 0, 2, 4, ...
                ELSE
                   l=2*(mch-1)+1 ! l = 1, 3, 5, ...
                   lp=2*(nch-1)+1 ! l' = 1, 3, 5, ...
                ENDIF
                BPD%Pot(mch,nch,lx,kx) = -2d0*DP%Cllp(l,lp,DP%ml)/BPD%x(lx,kx)**3
                IF(mch.EQ.nch) BPD%Pot(mch,nch,lx,kx) = BPD%Pot(mch,nch,lx,kx) + 0.5d0*l*(l+1)/BPD%x(lx,kx)**2
                BPD%Pot(nch,mch,lx,kx) = BPD%Pot(mch,nch,lx,kx)
                !write(6,*) BPD%x(lx,kx), BPD%Pot(mch,nch,lx,kx)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  

  END SUBROUTINE SetDipoleDipolePot
END MODULE DipoleDipole
