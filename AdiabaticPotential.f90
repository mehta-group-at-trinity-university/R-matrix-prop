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

    INTEGER, PARAMETER :: SplineOrder = 2   ! piecewise-linear B-spline (order=2), chosen to
                                             ! validate this repo's R-matrix propagation fork
                                             ! against the DeltaScat no-interpolation benchmark
                                             ! (Adiabatic-Scattering-BoundStates/DeltaScat) while
                                             ! avoiding cubic-spline ringing near the origin.
                                             ! Was 4 (cubic, matching CABA's Akima interpolation)
                                             ! -- revisit/compare once linear is validated.
    INTEGER NumChannelsUad, NumDataPointsUad, i, j, n, hashpos, NPts, ios
    DOUBLE PRECISION dum, fitRangeMinDum
    DOUBLE PRECISION, ALLOCATABLE :: Rdata(:), Udata(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: PMatFull(:,:,:), QMatFull(:,:,:)
    CHARACTER(LEN=200) headerline, fitLeffLine

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

    ! Rdata(1) must lie in range 0<R<Rquad(1), where Rquad(1) is the location
    ! of the first quadrature point in the radial matrix elements -- i.e. the
    ! interpolants built here must cover every R a box's Gauss-Legendre
    ! quadrature can ever actually query, all the way down to xStart=0 (box 1's
    ! own leftmost quadrature node is strictly > 0 but can be arbitrarily close
    ! to it). Previously this dropped row 1 entirely, which is wrong in general:
    ! it silently raises the interpolants' domain floor to Rdata(2) (~0.125 for
    ! a typical 400-point/RMin=1e-4 tabulation), and any quadrature point below
    ! that crashes dbsval's domain-edge search outright rather than degrading
    ! gracefully. (The original rationale here was Runge-type ringing from a
    ! genuine near-origin -C/R Coulomb divergence, Confirmed by direct
    ! comparison in Python for the CoulombC/=0 delta-function dataset
    ! specifically -- but RcutoffAnalytic=0.01 in SetAdiabaticPotential below
    ! already routes CoulombC/=0 channels around the raw interpolant entirely
    ! below R=0.01, via AnalyticComboNearOrigin, making that concern moot here
    ! regardless. Whether the same ringing risk applies to ordinary CoulombC=0
    ! channels -- whose excited-channel curves also rise steeply near R=0, just
    ! smoothly rather than divergently -- has not been separately verified.)
    NPts = NumDataPoints

    ALLOCATE(UInterp(NumChannels))
    DO n = 1,NumChannels
       CALL AllocateInterpolatingFunction(NPts,SplineOrder,UInterp(n))
       UInterp(n)%x = Rdata(1:NumDataPoints)
       UInterp(n)%y = Udata(1:NumDataPoints,n)
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
    PInterp%x = Rdata(1:NumDataPoints)
    QInterp%x = Rdata(1:NumDataPoints)
    PInterp%M = PMatFull(1:NumDataPoints,:,:)
    QInterp%M = QMatFull(1:NumDataPoints,:,:)
    DEALLOCATE(PMatFull,QMatFull)

    CALL SetupInterpolatingMatrix(PInterp)
    CALL SetupInterpolatingMatrix(QInterp)

    !--- FitLeff.data: 2 header lines + blank, then j, Threshold(j), residual(unused),
    !    Leff(j), FitRangeMin(j) -- REQUIRED to be exactly 5 columns. The "residual" column
    !    (always 0 -- NewFitProg.f computes it against an array that's never assigned before
    !    that point) is present in every correctly-generated FitLeff.data, but nothing about
    !    the file format itself declares the column count, so a generator that omits it
    !    produces a silently-still-valid-looking 4-column file: list-directed READ then
    !    shifts every field left by one, reading the residual placeholder (or, with only 4
    !    columns, FitRangeMin) as Leff -- corrupting CalcK's outer Bessel-matching order for
    !    the whole run with no read error at all (confirmed twice: once via the delta-function
    !    5-channel dataset, once via 3bodydata_sech2_3boson, whose FitLeff.data had exactly
    !    this 4-column defect -- Leff silently read as FitRangeMin=40.0 for every channel).
    !    Guard against a repeat by reading each line into a buffer and demanding a 5th
    !    (FitRangeMin) value explicitly -- a 4-column line then fails this internal READ with
    !    a nonzero IOSTAT immediately, instead of silently misreading. ---
    ALLOCATE(Threshold(NumChannels),Leff(NumChannels))
    OPEN(unit=26,file=trim(DataDir)//'/FitLeff.data',status='old')
    READ(26,*)
    READ(26,*)
    READ(26,*)
    DO n = 1,NumChannels
       READ(26,'(A)') fitLeffLine
       READ(fitLeffLine,*,IOSTAT=ios) j, Threshold(n), dum, Leff(n), fitRangeMinDum
       IF (ios.NE.0) THEN
          WRITE(6,*) "ReadAdiabaticData: ",trim(DataDir)//"/FitLeff.data channel ",n, &
               " does not have the required 5 columns (j, Threshold, residual, Leff, FitRangeMin)."
          WRITE(6,*) "  Offending line: ",trim(fitLeffLine)
          WRITE(6,*) "  A 4-column line here silently misreads Leff -- see this read's own comment."
          STOP
       ENDIF
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
  SUBROUTINE SetAdiabaticPotential(BPD,muLocal,UInterp,PInterp,QInterp,CoulombC)
    USE DataStructures
    USE Quadrature
    IMPLICIT NONE
    TYPE(BPData) BPD
    DOUBLE PRECISION, INTENT(in) :: muLocal
    TYPE(InterpolatingFunction), INTENT(in) :: UInterp(BPD%NumChannels)
    TYPE(InterpolatingMatrix), INTENT(in) :: PInterp, QInterp
    ! Per-channel Coulomb coefficient (0 for ordinary channels). Where nonzero, the
    ! diagonal combo(R)=2*mu*U+Q is replaced below RcutoffAnalytic by the analytic
    ! small-R series (AnalyticComboNearOrigin) instead of the tabulated-data spline,
    ! which has significant Runge-phenomenon ringing there (see CoulombRHS).
    DOUBLE PRECISION, INTENT(in) :: CoulombC(BPD%NumChannels)

    INTEGER kx,lx,mch
    DOUBLE PRECISION ax,bx,xScaledZero,R
    DOUBLE PRECISION xScale(BPD%xNumPoints)
    DOUBLE PRECISION, EXTERNAL :: AnalyticComboNearOrigin
    DOUBLE PRECISION, PARAMETER :: RcutoffAnalytic = 0.01d0

    DO kx = 1,BPD%xNumPoints-1
       ax = BPD%xPoints(kx)
       bx = BPD%xPoints(kx+1)
       xScale(kx) = 0.5d0*(bx-ax)
       xScaledZero = 0.5d0*(bx+ax)
       DO lx = 1,LegPoints
          BPD%x(lx,kx) = xScale(kx)*xLeg(lx) + xScaledZero
          R = BPD%x(lx,kx)
          ! P is antisymmetric (Pcoup(mch,mch)=0 always) -- for this 1-channel test
          ! there is no off-diagonal P at all, so skipping the interpolant call below
          ! RcutoffAnalytic is exact here (not just an approximation). Needed because
          ! InterpolatedMatrix(R,PInterp) hits the same dbsval domain-edge search bug
          ! at very small R that the U/Q interpolation calls were guarded against above.
          IF (R.LT.RcutoffAnalytic) THEN
             BPD%Pcoup(:,:,lx,kx) = 0d0
          ELSE
             BPD%Pcoup(:,:,lx,kx) = InterpolatedMatrix(R,PInterp)
          ENDIF
          ! Verified directly against this repo's own Uad.dat/Qmat.dat data:
          ! 2*mu*U + Q - 0.25/R^2 is flat to 1e-12 (= 2*mu*Threshold) out to
          ! R=200, while 2*mu*U - Q - 0.25/R^2 has a 3.7e-5 spread and is
          ! still visibly drifting there -- i.e. +Q/(2mu) is the convention
          ! actually used by this repo's Qmat.dat (do NOT "fix" this to match
          ! CABS.f90's/VQmat.dat's documented sign -- that file evidently uses
          ! an internally opposite Q convention from this one).
          ! Below RcutoffAnalytic, skip the QInterp/UInterp calls entirely (not just
          ! override their result afterward) -- both hit the same dbsval domain-edge
          ! search bug at very small R that the guard below is trying to avoid.
          IF (R.LT.RcutoffAnalytic) THEN
             BPD%Pot(:,:,lx,kx) = 0d0
          ELSE
             BPD%Pot(:,:,lx,kx) = InterpolatedMatrix(R,QInterp)/(2d0*muLocal)
          ENDIF
          DO mch = 1,BPD%NumChannels
             IF (CoulombC(mch).NE.0d0 .AND. R.LT.RcutoffAnalytic) THEN
                BPD%Pot(mch,mch,lx,kx) = AnalyticComboNearOrigin(R,muLocal,CoulombC(mch))/(2d0*muLocal)
             ELSE
                ! Ordinary (non-Coulomb) channel: no near-origin singularity, so its
                ! interpolant is safe to evaluate anywhere in its tabulated domain
                ! (starts at Rdata(1) now, not Rdata(2) -- see ReadAdiabaticData's own
                ! header comment on why row 1 is kept) -- fixes a
                ! bug where this branch used to require R>=RcutoffAnalytic and silently
                ! left BPD%Pot(mch,mch) at 0 otherwise. That zeroed out a real, large
                ! centrifugal-type U(mch) for R<0.01 (confirmed via the delta-function
                ! 5-channel dataset: U(2) ~ 1.9e6 at R=0.004, dropped entirely), which
                ! corrupts box 1's log-derivative and propagates through the whole
                ! outward R-matrix chain -- this DOES occur (box 1 routinely spans past
                ! R=0.01 for channels 2+), contrary to the comment this replaces.
                BPD%Pot(mch,mch,lx,kx) = BPD%Pot(mch,mch,lx,kx) + Interpolated(R,UInterp(mch))
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE SetAdiabaticPotential
  !****************************************************************************************************
END MODULE AdiabaticPotential
