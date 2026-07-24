!****************************************************************************************************
! Imports the adiabatic hyperspherical P/Q-tilde/U data produced by
! Adiabatic-Scattering-BoundStates' adiabatic solvers (Uad.dat, Pmat.dat, Qmat.dat,
! masses.dat) plus the long-range fit outputs already required to run that repo's own
! CABA scattering pipeline (Fit.data for NumChannels/alpha/ddim, FitLeff.data for the
! per-channel threshold and effective angular momentum used for asymptotic Bessel
! matching in CalcK). Builds B-spline interpolants via the canonical interpolation
! library (~/Documents/GitHub/interpolation/Interpolation.f90, module InterpType) rather
! than CABS.f90's Akima tables, since R-matrix-prop's boxes need pointwise evaluation at
! arbitrary R rather than a fixed precomputed grid.
!****************************************************************************************************
MODULE AdiabaticPotential
  USE InterpType
  IMPLICIT NONE

CONTAINS
  !****************************************************************************************************
  SUBROUTINE ReadAdiabaticData(DataDir,NumChannels,NumDataPoints,muLocal,alpha,ddim, &
       Threshold,Leff,UInterp,PInterp,QInterp)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: DataDir
    INTEGER, INTENT(out) :: NumChannels, NumDataPoints
    DOUBLE PRECISION, INTENT(out) :: muLocal, alpha, ddim
    DOUBLE PRECISION, ALLOCATABLE, INTENT(out) :: Threshold(:), Leff(:)
    TYPE(InterpolatingFunction), ALLOCATABLE, INTENT(out) :: UInterp(:)
    TYPE(InterpolatingMatrix), INTENT(out) :: PInterp, QInterp

    INTEGER, PARAMETER :: SplineOrder = 4   ! cubic B-spline, matching CABA's Akima cubic
                                             ! interpolation -- testing whether the residual
                                             ! ~2-5x gap against CABA traces to interpolation
                                             ! order. (Was 2/piecewise-linear, matching the
                                             ! precedent set in Adiabatic-R-Mat-Prop's
                                             ! SetReadPotential, chosen there to stay well-behaved
                                             ! through Q-tilde's sharp features near
                                             ! adiabatic-curve avoided crossings -- revisit if
                                             ! cubic proves unstable there.)
    INTEGER NumChannelsUad, NumDataPointsUad, i, j, n, hashpos, NPts
    DOUBLE PRECISION dum
    DOUBLE PRECISION, ALLOCATABLE :: Rdata(:), Udata(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: PMatFull(:,:,:), QMatFull(:,:,:)
    CHARACTER(LEN=200) headerline

    !--- Fit.data header: NumChannels, NumDataPoints, alpha, ddim ---
    OPEN(unit=21,file=trim(DataDir)//'/Fit.data',status='old')
    READ(21,*) NumChannels, NumDataPoints, alpha, ddim
    CLOSE(21)

    !--- masses.dat: comment line, then mu (hyperradial reduced mass) ---
    OPEN(unit=22,file=trim(DataDir)//'/masses.dat',status='old')
    READ(22,*)
    READ(22,*) muLocal
    CLOSE(22)

    !--- Uad.dat: "# NumChannels NumDataPoints" header, then R + U(1:NumChannels) rows ---
    OPEN(unit=23,file=trim(DataDir)//'/Uad.dat',status='old')
    READ(23,'(A)') headerline
    hashpos = INDEX(headerline,'#')
    READ(headerline(hashpos+1:),*) NumChannelsUad, NumDataPointsUad
    IF (NumChannelsUad.NE.NumChannels .OR. NumDataPointsUad.NE.NumDataPoints) THEN
       WRITE(6,*) "ReadAdiabaticData: Uad.dat header (",NumChannelsUad,NumDataPointsUad, &
            ") disagrees with Fit.data header (",NumChannels,NumDataPoints,")"
       STOP
    ENDIF
    ALLOCATE(Rdata(NumDataPoints),Udata(NumDataPoints,NumChannels))
    DO i = 1,NumDataPoints
       READ(23,*) Rdata(i), Udata(i,1:NumChannels)
    ENDDO
    CLOSE(23)

    ! Drop row 1 (R~1e-6) before building any interpolant: it is numerically
    ! imprecise at that scale (off from the true near-origin -C/R asymptote by
    ! only ~1e-3 relative, but at R~1e-6 that is an O(1000) ABSOLUTE error in
    ! U/Q), and feeding it into a global interpolating spline causes Runge-type
    ! ringing that doesn't fully damp until R~0.005-0.01 and leaves a residual
    ! ~1% offset out to R~0.05 -- exactly the box-1 xMin range this matters at.
    ! Confirmed by direct comparison (with vs without this point) in Python.
    NPts = NumDataPoints - 1

    ALLOCATE(UInterp(NumChannels))
    DO n = 1,NumChannels
       CALL AllocateInterpolatingFunction(NPts,SplineOrder,UInterp(n))
       UInterp(n)%x = Rdata(2:NumDataPoints)
       UInterp(n)%y = Udata(2:NumDataPoints,n)
       CALL SetupInterpolatingFunction(UInterp(n))
    ENDDO

    !--- Pmat.dat / Qmat.dat: R-line, then NumChannels rows of NumChannels values ---
    ALLOCATE(PMatFull(NumDataPoints,NumChannels,NumChannels))
    ALLOCATE(QMatFull(NumDataPoints,NumChannels,NumChannels))
    OPEN(unit=24,file=trim(DataDir)//'/Pmat.dat',status='old')
    DO i = 1,NumDataPoints
       READ(24,*) dum
       DO j = 1,NumChannels
          READ(24,*) PMatFull(i,j,1:NumChannels)
       ENDDO
    ENDDO
    CLOSE(24)

    OPEN(unit=25,file=trim(DataDir)//'/Qmat.dat',status='old')
    DO i = 1,NumDataPoints
       READ(25,*) dum
       DO j = 1,NumChannels
          READ(25,*) QMatFull(i,j,1:NumChannels)
       ENDDO
    ENDDO
    CLOSE(25)

    CALL AllocateInterpolatingMatrix(NPts,SplineOrder,NumChannels,NumChannels,PInterp)
    CALL AllocateInterpolatingMatrix(NPts,SplineOrder,NumChannels,NumChannels,QInterp)
    PInterp%x = Rdata(2:NumDataPoints)
    QInterp%x = Rdata(2:NumDataPoints)
    PInterp%M = PMatFull(2:NumDataPoints,:,:)
    QInterp%M = QMatFull(2:NumDataPoints,:,:)
    DEALLOCATE(PMatFull,QMatFull)

    CALL SetupInterpolatingMatrix(PInterp)
    CALL SetupInterpolatingMatrix(QInterp)

    !--- FitLeff.data: 2 header lines + blank, then j, Threshold(j), <resid>, Leff(j), FitRangeMin(j) ---
    ALLOCATE(Threshold(NumChannels),Leff(NumChannels))
    OPEN(unit=26,file=trim(DataDir)//'/FitLeff.data',status='old')
    READ(26,*)
    READ(26,*)
    READ(26,*)
    DO n = 1,NumChannels
       READ(26,*) j, Threshold(n), dum, Leff(n)
    ENDDO
    CLOSE(26)

    DEALLOCATE(Rdata,Udata)
  END SUBROUTINE ReadAdiabaticData
  !****************************************************************************************************
  ! Fills BPD%Pot/BPD%Pcoup for one box by evaluating the interpolants at the box's own
  ! Gauss-Legendre quadrature points -- mirrors SetBalujaPotential's loop structure.
  ! BPD%Pot(mch,nch) = U_mch(R)*delta_mn + Qtilde_mn(R)/(2*mu): both pieces are symmetric,
  ! local, multiplicative couplings, so they share CalcGamLam's existing Pot mechanism.
  ! BPD%Pcoup(mch,nch) = P_mn(R): the antisymmetric derivative coupling CalcGamLam's new
  ! P-coupling loop consumes separately.
  SUBROUTINE SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp)
    USE DataStructures
    USE Quadrature
    IMPLICIT NONE
    TYPE(BPData) BPD
    DOUBLE PRECISION, INTENT(in) :: muLocal
    TYPE(InterpolatingFunction), INTENT(in) :: UInterp(BPD%NumChannels)
    TYPE(InterpolatingMatrix), INTENT(in) :: PInterp, QInterp

    INTEGER kx,lx,mch
    DOUBLE PRECISION ax,bx,xScaledZero,R
    DOUBLE PRECISION xScale(BPD%xNumPoints)

    DO kx = 1,BPD%xNumPoints-1
       ax = BPD%xPoints(kx)
       bx = BPD%xPoints(kx+1)
       xScale(kx) = 0.5d0*(bx-ax)
       xScaledZero = 0.5d0*(bx+ax)
       DO lx = 1,LegPoints
          BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
          R = BPD%x(lx,kx)
          BPD%Pcoup(:,:,lx,kx) = InterpolatedMatrix(R,PInterp)
          BPD%Pot(:,:,lx,kx) = InterpolatedMatrix(R,QInterp)/(2d0*muLocal)
          DO mch = 1,BPD%NumChannels
             BPD%Pot(mch,mch,lx,kx) = BPD%Pot(mch,mch,lx,kx) + Interpolated(R,UInterp(mch))
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE SetAdiabaticPotential
  !****************************************************************************************************
END MODULE AdiabaticPotential
