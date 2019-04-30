module DipoleDipole
  use Quadrature
  use DataStructures

  type DPData
    integer lmax,ml
    double precision, allocatable :: cllp(:,:,:)
    logical even
  end type DPData

CONTAINS
  subroutine AllocateDP(DP)
    implicit NONE
    type(DPData) DP
    allocate(DP%cllp(0:DP%lmax,0:DP%lmax,0:DP%lmax))

  end subroutine AllocateDP

  subroutine MakeDipoleDipoleCouplingMatrix(DP)
    implicit none
    type(DPData) :: DP
    integer l,lp,m
    double precision, external :: THRJ
    double precision tj1,tj2,phase,prefact


    DP%cllp=0d0
      do l=0,DP%lmax
        do lp=MAX(l-2,0),MIN(DP%lmax,l+2),2
        !do lp=0,lmax
          prefact=dble((2*l+1)*(2*lp+1))
          prefact=dsqrt(prefact)
          do m = 0,l
            phase=(-1)**m
            tj1 = THRJ(2*l,2*2,2*lp,-2*m,0,2*m)
            tj2 = THRJ(2*l,2*2,2*lp,0,0,0)
            DP%cllp(l,lp,m)=prefact*phase*tj1*tj2
            DP%cllp(lp,l,m)=DP%cllp(l,lp,m)
          ENDDO
        ENDDO
      ENDDO

  end subroutine MakeDipoleDipoleCouplingMatrix

  SUBROUTINE SetDipoleDipolePot(BPD,DP)

    IMPLICIT NONE
    TYPE(BPData) BPD
    type(DPData) DP
    integer m,l,lp
    INTEGER kx,lx,mch,nch,NChan
    DOUBLE PRECISION ax,bx,xScaledZero
    DOUBLE PRECISION xScale(BPD%xNumPoints)

    if(mod(DP%lmax,2).eq.0) then
      Nchan=DP%lmax/2
    ELSE
      Nchan=(DP%lmax-1)/2
    ENDIF

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
               if(DP%even) THEN
                 l=2*(mch-1)
                 lp=2*(nch-1)
               ELSE
                 l=2*(mch-1)+1
                 lp=2*(nch-1)+1
               ENDIF
                BPD%Pot(mch,nch,lx,kx) = -2d0*DP%Cllp(l,lp,DP%ml)/BPD%x(lx,kx)**3
                BPD%Pot(nch,mch,lx,kx) = BPD%Pot(mch,nch,lx,kx)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE SetDipoleDipolePot
end module DipoleDipole
