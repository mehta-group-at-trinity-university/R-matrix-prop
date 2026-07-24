  !****************************************************************************************************
MODULE GlobalVars
  IMPLICIT NONE
  INTEGER NumParticles, NumChannels,  NumAllChan, Order
  INTEGER StartBC, EndBC, xTotNumPoints, NumBoxes
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION AlphaFactor ! This is the parameter that appears in the reduced wavefunction u(R) = R^(AlphaFactor) Psi(R)
  ! Typical choice is either AlphaFactor = 0 (reduced wavefunction = wavefunction), or AlphaFactor = (EffDim - 1)/2 (eliminates 1st derivative terms from KE)
  ! For this calculation, we ALWAYS set AlphaFactor = (EffDim-1)/2 and set StartBC = 0 (Left = 0).
  !----------------------------------------------------------------------------------------------------
  DOUBLE PRECISION reducedmass, xStart, xEnd, energy,kStart,kEnd
  DOUBLE PRECISION SpatialDim, EffDim, Pi
  DOUBLE PRECISION, ALLOCATABLE :: mass(:)
  CHARACTER*64 InputFile
  COMPLEX*16 II
  PARAMETER(II=(0.0d0,1.0d0))
!****************************************************************************************************
CONTAINS
  SUBROUTINE ReadGlobal
    USE Quadrature
    IMPLICIT NONE
    INTEGER n
    ! Be sure to match the format of the input file to the format of the read statements below
    OPEN(unit=7,file=InputFile(1:INDEX(InputFile,' ')-1),action='read')
    READ(7,*)
    READ(7,*) NumParticles, NumChannels, SpatialDim, NumAllChan, Order
    ALLOCATE(mass(NumParticles))
    READ(7,*)
    READ(7,*)
    READ(7,*) (mass(n), n=1,NumParticles)
    READ(7,*)
    READ(7,*)
    READ(7,*) xStart, xEnd, xTotNumPoints, LegPoints
    READ(7,*)
    READ(7,*)
    READ(7,*) NumBoxes, StartBC, EndBC, kStart, kEnd
    READ(7,*)
    READ(7,*)
    READ(7,*) Energy
    Pi = dacos(-1d0)
    CLOSE(unit=7)
    EffDim = NumParticles*SpatialDim - SpatialDim
    !AlphaFactor = 0d0
    AlphaFactor = (EffDim-1d0)/2d0

    IF (NumParticles.EQ.2) THEN
       reducedmass = mass(1)*mass(2)/(mass(1)+mass(2))
    ELSE
       WRITE(6,*) "Reduced mass not set. Must set reduced mass"
       STOP
    END IF


  END SUBROUTINE ReadGlobal
  !****************************************************************************************************
END MODULE GlobalVars
!****************************************************************************************************
!****************************************************************************************************
MODULE DataStructures
  IMPLICIT NONE
  !****************************************************************************************************
  TYPE BPData
    !This data type contains the array of basis functions and potential matrix
     INTEGER xNumPoints, xDim, Left, Right, Order, NumChannels, MatrixDim
     DOUBLE PRECISION, ALLOCATABLE :: u(:,:,:),ux(:,:,:),xPoints(:),x(:,:)
     DOUBLE PRECISION, ALLOCATABLE :: Pot(:,:,:,:),Pcoup(:,:,:,:),lam(:)
     DOUBLE PRECISION kl, kr, xl,xr
     INTEGER, ALLOCATABLE :: xBounds(:), kxMin(:,:), kxMax(:,:)
  END TYPE BPData
  !****************************************************************************************************
  TYPE GenEigVal
     INTEGER MatrixDim
     ! Gam0/Overlap: energy-independent pieces of Gam (Gam0 = kinetic+potential+
     ! centrifugal+P/Q-coupling, Overlap = bare channel-diagonal <u,u'> integral), built
     ! once per box and reused across an energy scan via CombineGam. Gam is the per-energy
     ! combination Gam0+Energy*Overlap actually fed to PartitionAndEliminate.
     DOUBLE PRECISION, ALLOCATABLE :: Gam(:,:), Gam0(:,:), Overlap(:,:), Lam(:,:), evec(:,:), eval(:)
  END TYPE GenEigVal
  !****************************************************************************************************
  TYPE BoxData
     INTEGER NumOpenL, NumOpenR, betaMax
     DOUBLE PRECISION xL, xR
     DOUBLE PRECISION, ALLOCATABLE :: Z(:,:), ZP(:,:), NBeta(:), Norm(:,:)
     DOUBLE PRECISION, ALLOCATABLE :: b(:), bf(:), Zf(:,:), Zfp(:,:)
     !double precision, allocatable :: RF(:,:),K(:,:)

  END TYPE BoxData
  TYPE ScatData
     DOUBLE PRECISION, ALLOCATABLE :: K(:,:), R(:,:), f(:,:), sigma(:,:)
     COMPLEX*8, ALLOCATABLE :: S(:,:), T(:,:)

  END TYPE ScatData

CONTAINS
  !****************************************************************************************************
  SUBROUTINE AllocateScat(SD,N)
    IMPLICIT NONE
    TYPE(ScatData) SD
    INTEGER N
    ALLOCATE(SD%K(N,N),SD%R(N,N),SD%f(N,N),SD%sigma(N,N))
    ALLOCATE(SD%S(N,N),SD%T(N,N))
    SD%K=0d0
    SD%R=0d0
    SD%f=0d0
    SD%sigma=0d0
    SD%S=(0d0,0d0)
    SD%T=(0d0,0d0)
  END SUBROUTINE AllocateScat
  !****************************************************************************************************
  SUBROUTINE DeAllocateScat(SD)
    IMPLICIT NONE
    TYPE(ScatData) SD

    DEALLOCATE(SD%K,SD%R,SD%f,SD%sigma,SD%S,SD%T)

  END SUBROUTINE DeAllocateScat
  !****************************************************************************************************
  SUBROUTINE AllocateEIG(EIG)
    IMPLICIT NONE
    TYPE(GenEigVal) :: EIG
    ALLOCATE(EIG%Gam(EIG%MatrixDim,EIG%MatrixDim),EIG%Lam(EIG%MatrixDim,EIG%MatrixDim))
    ALLOCATE(EIG%Gam0(EIG%MatrixDim,EIG%MatrixDim),EIG%Overlap(EIG%MatrixDim,EIG%MatrixDim))
    ALLOCATE(EIG%eval(EIG%MatrixDim),EIG%evec(EIG%MatrixDim,EIG%MatrixDim))
    EIG%Gam=0d0
    EIG%Gam0=0d0
    EIG%Overlap=0d0
    EIG%Lam=0d0
    EIG%eval=0d0
    EIG%Evec=0d0

  END SUBROUTINE AllocateEIG
  !****************************************************************************************************
  SUBROUTINE DeAllocateEIG(EIG)
    IMPLICIT NONE
    TYPE(GenEigVal) :: EIG
    DEALLOCATE(EIG%Gam,EIG%Gam0,EIG%Overlap,EIG%Lam,EIG%eval,EIG%evec)

  END SUBROUTINE DeAllocateEIG
  !****************************************************************************************************
  SUBROUTINE AllocateBox(B)
    IMPLICIT NONE
    TYPE(BoxData) :: B

    B%betaMax = B%NumOpenR + B%NumOpenL

    ALLOCATE(B%Z(B%betaMax,B%betaMax),B%ZP(B%betaMax,B%betaMax))
    ALLOCATE(B%b(B%betaMax))
    ALLOCATE(B%NBeta(B%betaMax))
    ALLOCATE(B%Norm(B%betaMax,B%betaMax))
    ALLOCATE(B%Zf(B%NumOpenR,B%NumOpenR),B%Zfp(B%NumOpenR,B%NumOpenR))
    ALLOCATE(B%bf(B%NumOpenR))

  END SUBROUTINE AllocateBox
  !****************************************************************************************************
  SUBROUTINE InitZeroBox(B)
    IMPLICIT NONE
    TYPE(BoxData) :: B
    B%Z=0d0
    B%ZP=0d0
    B%b=0d0
    B%NBeta=0d0
    B%Norm=0d0
    B%Zf=0d0
    B%Zfp=0d0
    B%bf=0d0

  END SUBROUTINE InitZeroBox
  !****************************************************************************************************
  SUBROUTINE DeAllocateBox(B)
    IMPLICIT NONE
    TYPE(BoxData) :: B

    DEALLOCATE(B%Z,B%ZP)
    DEALLOCATE(B%b)
    DEALLOCATE(B%NBeta,B%Norm,B%Zf,B%bf,B%Zfp)


  END SUBROUTINE DeAllocateBox
  !****************************************************************************************************
  SUBROUTINE AllocateBPD(BPD)
    USE Quadrature
    IMPLICIT NONE
    TYPE(BPData) BPD
    INTEGER xDimMin
    !Determine the basis dimension for these boundary conditions
    xDimMin=BPD%xNumPoints+BPD%Order-3
    BPD%xDim=xDimMin
    IF (BPD%Left .EQ. 2) BPD%xDim = BPD%xDim + 1
    IF (BPD%Right .EQ. 2) BPD%xDim = BPD%xDim + 1
    BPD%MatrixDim = BPD%NumChannels*BPD%xDim

    ALLOCATE(BPD%xPoints(BPD%xNumPoints))
    ALLOCATE(BPD%xBounds(BPD%xNumPoints + 2*BPD%Order))   ! allocate xBounds before calling CalcBasisFuncs
    ALLOCATE(BPD%x(LegPoints,BPD%xNumPoints-1))

    ALLOCATE(BPD%u(LegPoints,BPD%xNumPoints,BPD%xDim),BPD%ux(LegPoints,BPD%xNumPoints,BPD%xDim)) ! allocate memory for basis functions
    ALLOCATE(BPD%kxMin(BPD%xDim,BPD%xDim),BPD%kxMax(BPD%xDim,BPD%xDim)) !
    ALLOCATE(BPD%Pot(BPD%NumChannels,BPD%NumChannels,LegPoints,BPD%xNumPoints-1))
    ALLOCATE(BPD%Pcoup(BPD%NumChannels,BPD%NumChannels,LegPoints,BPD%xNumPoints-1))
    BPD%Pcoup = 0d0   ! zero unless a potential module (e.g. SetAdiabaticPotential) fills it in

    ALLOCATE(BPD%lam(BPD%NumChannels))

    BPD%lam=0.d0

  END SUBROUTINE AllocateBPD
  !****************************************************************************************************
  SUBROUTINE DeAllocateBPD(BPD)
    IMPLICIT NONE
    TYPE(BPData) BPD
    DEALLOCATE(BPD%xPoints,BPD%xBounds,BPD%x,BPD%u,BPD%ux,BPD%kxMin,BPD%kxMax,BPD%Pot,BPD%Pcoup,BPD%lam)
  END SUBROUTINE DeAllocateBPD
  !****************************************************************************************************
  SUBROUTINE MakeBasis(BPD)
    USE Quadrature
    IMPLICIT NONE

    TYPE(BPData) BPD
    INTEGER ix, ixp

    CALL CalcBasisFuncsBP(BPD%Left,BPD%Right,BPD%kl,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         BPD%xDim,BPD%xBounds,BPD%xNumPoints,0,BPD%u)
    CALL CalcBasisFuncsBP(BPD%Left,BPD%Right,BPD%kl,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         BPD%xDim,BPD%xBounds,BPD%xNumPoints,1,BPD%ux)

    ! Determine the bounds for integration of matrix elements
    DO ix = 1,BPD%xDim
       DO ixp = 1,BPD%xDim
          BPD%kxMin(ixp,ix) = MAX(BPD%xBounds(ix),BPD%xBounds(ixp))
          BPD%kxMax(ixp,ix) = MIN(BPD%xBounds(ix+BPD%Order+1),BPD%xBounds(ixp+BPD%Order+1))-1 !
       ENDDO
    ENDDO

  END SUBROUTINE MakeBasis
  !****************************************************************************************************
  SUBROUTINE checkpot(BPD,file)
    USE Quadrature
    IMPLICIT NONE
    INTEGER file,mch,nch,lx,kx
    TYPE(BPData) BPD

    DO mch = 1,BPD%NumChannels
       DO nch = 1,mch
          DO kx = 1, BPD%xNumPoints-1
             DO lx = 1,LegPoints
                WRITE(file,*) BPD%x(lx,kx), BPD%Pot(mch,nch,lx,kx)
             ENDDO
          ENDDO
          WRITE(file,*)
       ENDDO
       WRITE(file,*)
    ENDDO

  END SUBROUTINE checkpot

END MODULE DataStructures
!****************************************************************************************************
  SUBROUTINE CalcGamLam(BPD,EIG,NumOpenL,NumOpenR)
    USE DataStructures
    USE Quadrature
    USE GlobalVars
    IMPLICIT NONE
    TYPE(BPData), INTENT(in) :: BPD
    TYPE(GenEigVal) :: EIG
    ! NumOpenL/NumOpenR: number of channels (a prefix of the channel list --
    ! open channels must precede closed ones, same convention as CABA's
    ! NumberL/NumberR) that are boundary-retained ("R-matrix-theory open") at
    ! this box's left/right edge. Channels beyond these counts get NO Lambda
    ! entry, so the generalized eigenvalue problem Gam*c=b*Lam*c drives their
    ! eigenvalue to "infinite" and BoxMatch's finite-eigenvalue filter
    ! correctly folds them into the eliminated block instead of carrying them
    ! forward as retained boundary unknowns. This is the same elimination
    ! Burke's thesis Eq. 3.4 achieves via explicit Gcc/Gco/Goo partitioning
    ! (and what CABA.EXT.f's CalcClosed/CalcOpen do) -- just expressed as a
    ! singular Lambda instead of a partitioned linear solve.
    INTEGER, INTENT(in) :: NumOpenL, NumOpenR
    INTEGER ix,ixp,lx,kx,mch,nch
    DOUBLE PRECISION a, ax, bx, xIntScale(BPD%xNumPoints), TempG, TempS, xScaledZero
    EIG%Gam0=0d0
    EIG%Overlap=0d0
    EIG%Lam=0d0
    EIG%eval=0d0
    EIG%evec=0d0
    !Calculate the energy-independent Gamma0 matrix elements, plus the (also
    !energy-independent) channel-diagonal Overlap integral <u,u'> that Energy multiplies.
    !CombineGam assembles the actual Gam=Gam0+Energy*Overlap fed to PartitionAndEliminate.
    DO ix = 1,BPD%xDim
       DO ixp = MAX(1,ix-BPD%Order),MIN(BPD%xDim,ix+BPD%Order)
          DO mch = 1,BPD%NumChannels
             ! Do the channel-diagonal part first:
             DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                ax = BPD%xPoints(kx)
                bx = BPD%xPoints(kx+1)
                xIntScale(kx) = 0.5d0*(bx-ax)
                xScaledZero = 0.5d0*(bx+ax)
                TempG = 0.0d0
                TempS = 0.0d0
                DO lx = 1,LegPoints
                   a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                   ! The KE matrix elements
                   TempG = TempG - a*(BPD%ux(lx,kx,ix)*BPD%ux(lx,kx,ixp))
                   ! The diagonal "overlap" (bare integral -- Energy multiplies it later in
                   ! CombineGam) and additional term from reducing the wavefunction
                   a = a*BPD%u(lx,kx,ix)*BPD%u(lx,kx,ixp)*2d0*reducedmass
                   TempS = TempS + a
                   TempG = TempG - a*BPD%Pot(mch,mch,lx,kx)
                   TempG = TempG - a*AlphaFactor*(AlphaFactor-EffDim+2d0)/(2d0*reducedmass*BPD%x(lx,kx)**2)
                ENDDO
                EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) = EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) + TempG ! place values into Gamma0
                EIG%Overlap((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) = EIG%Overlap((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) + TempS
             ENDDO

             ! Now do the off-diagonal parts
             DO nch = 1, mch-1
                DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                   ax = BPD%xPoints(kx)
                   bx = BPD%xPoints(kx+1)
                   xIntScale(kx) = 0.5d0*(bx-ax)
                   xScaledZero = 0.5d0*(bx+ax)
                   TempG = 0.0d0
                   DO lx = 1,LegPoints
                      a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                      ! The potential matrix elements calculated here
                      TempG = TempG - a*BPD%u(lx,kx,ix)*2.0d0*reducedmass*BPD%Pot(mch,nch,lx,kx)*BPD%u(lx,kx,ixp)
                   ENDDO
                   EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) = EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) + TempG ! place values into EIG%Gam0
                   ! Accumulate (not overwrite) into the symmetric element mch<-->nch: Pot is
                   ! symmetric so both entries get the same per-kx contribution. Must be += here
                   ! (not a single post-loop assignment copying the final mch,nch total) because the
                   ! P-coupling loop below also writes into this same (nch,mch) entry independently,
                   ! during whichever earlier/later outer-mch pass has mch_outer=nch -- an overwrite
                   ! here would silently discard that contribution.
                   EIG%Gam0((ix-1)*BPD%NumChannels+nch,(ixp-1)*BPD%NumChannels+mch) = EIG%Gam0((ix-1)*BPD%NumChannels+nch,(ixp-1)*BPD%NumChannels+mch) + TempG
                ENDDO
             ENDDO

             ! Non-adiabatic P coupling (antisymmetric first-derivative coupling,
             ! Pcoup(mch,nch)=-Pcoup(nch,mch)): Gamma=2*mu*(E-H) turns CABS.f90's
             ! H-matrix term (-P_mn/2mu)*(Bi*dBj/dR-dBi/dR*Bj) into
             ! +P_mn*(Bi*dBj/dR-dBi/dR*Bj). Computed directly (no mirror-fill) for
             ! every (mch,nch) pair -- swapping (ix,mch)<->(ixp,nch) flips the sign of
             ! both Pcoup and the antisymmetric (u*ux'-ux*u') factor, so the two flips
             ! cancel and Gam(Row,Col)=Gam(Col,Row) comes out correct automatically.
             DO nch = 1,BPD%NumChannels
                IF (nch.NE.mch) THEN
                   DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                      ax = BPD%xPoints(kx)
                      bx = BPD%xPoints(kx+1)
                      xIntScale(kx) = 0.5d0*(bx-ax)
                      xScaledZero = 0.5d0*(bx+ax)
                      TempG = 0.0d0
                      DO lx = 1,LegPoints
                         a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                         TempG = TempG + a*BPD%Pcoup(mch,nch,lx,kx)* &
                              (BPD%u(lx,kx,ix)*BPD%ux(lx,kx,ixp) - BPD%ux(lx,kx,ix)*BPD%u(lx,kx,ixp))
                      ENDDO
                      EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) = &
                           EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) + TempG
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !Calculate the Lambda matrix elements
    EIG%Lam = 0.d0

    IF(BPD%Right.EQ.2) THEN
       DO mch = 1, NumOpenR
          EIG%Lam( (BPD%xDim-1)*BPD%NumChannels+mch, (BPD%xDim-1)*BPD%NumChannels+mch ) = BPD%xr**(EffDim-1d0-2d0*AlphaFactor)
       ENDDO
    ENDIF
    IF(BPD%Left.EQ.2) THEN
       DO mch = 1, NumOpenL
          EIG%Lam(mch,mch) = BPD%xl**(EffDim-1d0-2d0*AlphaFactor)  ! (1-1)*NumChannels+mch = mch
       ENDDO
    ENDIF

  END SUBROUTINE CalcGamLam
  !****************************************************************************************************
  ! Cheap per-energy recombination of the energy-independent Gam0/Overlap (built once per
  ! box by CalcGamLam) into the actual Gam PartitionAndEliminate needs -- avoids redoing
  ! the full quadrature/integration in CalcGamLam at every energy in a scan.
  SUBROUTINE CombineGam(EIG,Energy)
    USE DataStructures
    IMPLICIT NONE
    TYPE(GenEigVal) :: EIG
    DOUBLE PRECISION, INTENT(in) :: Energy
    EIG%Gam = EIG%Gam0 + Energy*EIG%Overlap
  END SUBROUTINE CombineGam
  !****************************************************************************************************
  ! Builds box 1's per-channel, energy-dependent left-boundary (Robin/log-derivative)
  ! basis, exactly matching the true regular solution's log-derivative at R=xMin via the
  ! same hyperrjry machinery CalcK already trusts for the far-field asymptotic match. BPD
  ! must already have been built with Left=2 (two independent boundary functions at
  ! spatial index ix=1,2 -- see CalcBasisFuncsBP's case(2)). For each channel, forms
  !   u1(:,mch) = BPD%u(:,1) + BPD%u(:,2)/constLeft(mch)
  ! -- exactly Bsplines.f90's own case(3) (Left=3) combination, but built pointwise from
  ! the already-evaluated Left=2 arrays instead of re-deriving from raw B-splines, which
  ! is what makes a per-channel (energy-dependent) combination coefficient practical.
  !
  ! kxMinC/kxMaxC: the combined function's true support is the UNION of BPD%u(:,1) and
  ! BPD%u(:,2)'s individual supports, not just BPD%u(:,1)'s -- using BPD%kxMin/kxMax
  ! (built for the uncombined Left=2 basis) would silently under-integrate it. Left=3's
  ! own xBounds construction already encodes exactly this union (verified directly
  ! against Bsplines.f90's source), and does not depend on kLeft, so it's fetched once
  ! via a throwaway Left=3 call and spliced into a corrected xBounds array.
  SUBROUTINE BuildBox1CombinedBasis(BPD,Leff,mu,d,alpha,EE,Threshold,u1,ux1,kLeftOut,CoulombC,UInterp,QInterp)
    ! BPD must already be built with Left=3 (ONE boundary function at spatial index
    ! ix=1, xDim reduced by one relative to Left=2 -- see CalcBasisFuncsBP's case(3)),
    ! NOT Left=2. This is essential, not cosmetic: box 1 always has NumOpenL=0, so its
    ! boundary DOF is *fully variationally eliminated* regardless of what basis spans
    ! it. If Left=2's two independent functions are both kept (one "combined" via
    ! relabeling, the other left as the plain complementary function), {u1,B2} spans
    ! the *identical* 2D subspace that {B1,B2} already did -- a change of basis within
    ! an already-eliminated subspace cannot change the variational answer. Confirmed
    ! empirically this session: varying kLeft by 30x left the final K-matrix completely
    ! unchanged under that (wrong) approach. Only actually *dropping* a trial-space
    ! dimension (Left=3's smaller xDim) imposes a real constraint.
    !
    ! Per channel, forms u1(:,mch) = B1(:) + B2(:)/constLeft(mch) -- exactly
    ! Bsplines.f90's own case(3) combination -- from a throwaway Left=2 evaluation
    ! (B1,B2 are NOT BPD%u(:,:,1:2); BPD itself only has ONE boundary function now).
    ! This is what replaces BPD%u(:,:,1)/BPD%ux(:,:,1) as position 1 in CalcGamLamBox1
    ! -- BPD's own %kxMin/%kxMax (from its own Left=3 MakeBasis call) are already
    ! correct for it, no separate union-of-support correction needed.
    USE DataStructures
    USE GlobalVars
    USE Quadrature
    USE InterpType
    IMPLICIT NONE
    TYPE(BPData), INTENT(in) :: BPD
    DOUBLE PRECISION, INTENT(in) :: Leff(BPD%NumChannels), mu, d, alpha, EE
    DOUBLE PRECISION, INTENT(in) :: Threshold(BPD%NumChannels)
    DOUBLE PRECISION, INTENT(out) :: u1(LegPoints,BPD%xNumPoints,BPD%NumChannels)
    DOUBLE PRECISION, INTENT(out) :: ux1(LegPoints,BPD%xNumPoints,BPD%NumChannels)
    DOUBLE PRECISION, INTENT(out) :: kLeftOut(BPD%NumChannels)
    ! Per-channel Coulomb-potential coefficient C0 (channel effective potential ->
    ! -CoulombC(mch)/R as R->0). CoulombC(mch)=0 (the default/ordinary case) keeps the
    ! original hyperrjry/Bessel-based log-derivative unchanged for that channel.
    DOUBLE PRECISION, INTENT(in) :: CoulombC(BPD%NumChannels)
    ! Interpolants for the real (not idealized) channel potential U+Q/2mu -- needed by
    ! the CoulombC/=0 branch's numerical near-origin integration (CoulombRegularNumeric).
    TYPE(InterpolatingFunction), INTENT(in) :: UInterp(BPD%NumChannels)
    TYPE(InterpolatingMatrix), INTENT(in) :: QInterp

    DOUBLE PRECISION, EXTERNAL :: MyBSpline
    DOUBLE PRECISION kLeft, constLeft, k, rj, ry, rjp, ryp, uval, uderiv, R0
    INTEGER mch, xDim2
    INTEGER xBounds2(BPD%xNumPoints+2*BPD%Order)
    DOUBLE PRECISION, ALLOCATABLE :: uAux(:,:,:), uxAux(:,:,:)

    xDim2 = BPD%xDim + 1   ! Left=2's own xDim (one more than Left=3's)
    ALLOCATE(uAux(LegPoints,BPD%xNumPoints,xDim2),uxAux(LegPoints,BPD%xNumPoints,xDim2))
    CALL CalcBasisFuncsBP(2,BPD%Right,1d0,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         xDim2,xBounds2,BPD%xNumPoints,0,uAux)
    CALL CalcBasisFuncsBP(2,BPD%Right,1d0,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         xDim2,xBounds2,BPD%xNumPoints,1,uxAux)

    DO mch = 1,BPD%NumChannels
       ! exact log-derivative of the true regular solution at R=xMin (same k(i)
       ! open/closed convention CalcK uses, so both boundaries stay self-consistent)
       IF (2d0*mu*EE.GE.Threshold(mch)) THEN
          k = DSQRT(2d0*mu*EE-Threshold(mch))
       ELSE
          k = DSQRT(Threshold(mch)-2d0*mu*EE)
       ENDIF
       IF (CoulombC(mch).EQ.0d0) THEN
          CALL hyperrjry(INT(d),alpha,Leff(mch),k*BPD%xl,rj,ry,rjp,ryp)
          ! kLeft = u'(x1)/u(x1) directly (verified by hand: substituting constLeft into
          ! u=B1+B2/constLeft and using B2(x1)=0 for a clamped knot vector gives
          ! u(x1)=B1(x1), u'(x1)=B1(x1)*kLeft, i.e. u'(x1)/u(x1)=kLeft exactly).
          ! rjp is d(rj)/d(kR) (derivative w.r.t. hyperrjry's dimensionless argument), NOT
          ! d(rj)/dR -- CalcK's own sp(i)=k(i)*...*rhypjp makes the same chain-rule
          ! correction; missing the factor of k gives kLeft ~30x too large.
          kLeft = k*rjp/rj
       ELSE
          ! This channel's effective potential has a genuine -CoulombC(mch)/R
          ! divergence at the origin (Coulomb-like, not the pure 1/R^2 centrifugal
          ! structure hyperrjry assumes). The naive idealized-potential ODE
          ! u''+(1/R)u'+[k^2+2*mu*C0/R]u=0 is NOT enough on its own: fitting the real
          ! combo(R)=U(R)+Q11(R)/(2mu) data to -C0/R+D shows a genuine finite additive
          ! constant D (~-0.33 for the delta-function test case) that is neither 0 nor
          ! the channel threshold -- so instead of guessing D, integrate the EXACT
          ! equation with the real tabulated potential (UInterp/QInterp) numerically
          ! (RK4) from a small R0, deep in the region where Pot(R)->-C0/R is essentially
          ! exact (so the leading-order regular series is a valid IC there), out to
          ! R1=xMin (see CoulombRegularNumeric).
          ! R0 previously used MAX(BPD%xl/100,UInterp(mch)%x(1)*2) to stay clear of the
          ! interpolant's own left domain edge -- but UInterp(mch)%x(1) is the SECOND
          ! tabulated point (row 1 is dropped in ReadAdiabaticData), 2*that = 1.021e-4,
          ! which exceeds BPD%xl for any xMin below ~0.0102, making R0>R1 and the RK4
          ! integration run backward over a near-empty, nearly-degenerate range (confirmed:
          ! R0=1.021e-4 > R1=xMin=1e-4 in that regime). CoulombRHS now uses the analytic
          ! series (not the interpolant) for all R<0.01 regardless, so R0 no longer needs
          ! to respect the interpolant's domain at all -- just stay safely below xMin.
          R0 = MIN(BPD%xl*0.01d0, 1d-6)
          CALL CoulombRegularNumeric(mu,CoulombC(mch),EE,mch,BPD%NumChannels,UInterp,QInterp,R0,BPD%xl,uval,uderiv)
          kLeft = uderiv/uval
       ENDIF
       kLeftOut(mch) = kLeft

       ! combination coefficient -- exactly Bsplines.f90's own case(3) formula
       constLeft = (MyBSpline(BPD%Order,1,BPD%xNumPoints,BPD%xPoints,2,BPD%xPoints(1)) &
            - MyBSpline(BPD%Order,0,BPD%xNumPoints,BPD%xPoints,2,BPD%xPoints(1))*kLeft) &
            / (MyBSpline(BPD%Order,0,BPD%xNumPoints,BPD%xPoints,1,BPD%xPoints(1))*kLeft &
            - MyBSpline(BPD%Order,1,BPD%xNumPoints,BPD%xPoints,1,BPD%xPoints(1)))

       u1(:,:,mch) = uAux(:,:,1) + uAux(:,:,2)/constLeft
       ux1(:,:,mch) = uxAux(:,:,1) + uxAux(:,:,2)/constLeft
    ENDDO

    DEALLOCATE(uAux,uxAux)

  END SUBROUTINE BuildBox1CombinedBasis
  !****************************************************************************************************
  ! Simple, energy-independent "Neumann sum" box-1 basis: u1 = B1 + B2 (equal weight, constLeft=1),
  ! mirroring the 2006/2007-era adrmatprop.f's own "Left=1" boundary treatment (Bsplines.f case(1):
  ! u = MyBSpline(...,1,x) + MyBSpline(...,2,x)) -- NOT an energy/channel-dependent Robin match to any
  ! Bessel/Coulomb log-derivative at all. Intended to be combined with a box starting very close to the
  ! true origin (near the smallest tabulated R in Uad.dat/Qmat.dat), relying on the weak form's own
  ! R^(d-1-2*alpha) weight to regularize any 1/R-type potential divergence there (confirmed: for
  ! alpha=0, d=2, this weight is R^1, which exactly cancels a -C/R potential tail in the integrand).
  SUBROUTINE BuildBox1NeumannBasis(BPD,u1,ux1,kLeftOut)
    USE DataStructures
    USE GlobalVars
    USE Quadrature
    IMPLICIT NONE
    TYPE(BPData), INTENT(in) :: BPD
    DOUBLE PRECISION, INTENT(out) :: u1(LegPoints,BPD%xNumPoints,BPD%NumChannels)
    DOUBLE PRECISION, INTENT(out) :: ux1(LegPoints,BPD%xNumPoints,BPD%NumChannels)
    DOUBLE PRECISION, INTENT(out) :: kLeftOut(BPD%NumChannels)

    DOUBLE PRECISION, EXTERNAL :: MyBSpline
    INTEGER mch, xDim2
    INTEGER xBounds2(BPD%xNumPoints+2*BPD%Order)
    DOUBLE PRECISION, ALLOCATABLE :: uAux(:,:,:), uxAux(:,:,:)
    DOUBLE PRECISION uval, uderiv

    xDim2 = BPD%xDim + 1   ! Left=2's own xDim (one more than Left=3's)
    ALLOCATE(uAux(LegPoints,BPD%xNumPoints,xDim2),uxAux(LegPoints,BPD%xNumPoints,xDim2))
    CALL CalcBasisFuncsBP(2,BPD%Right,1d0,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         xDim2,xBounds2,BPD%xNumPoints,0,uAux)
    CALL CalcBasisFuncsBP(2,BPD%Right,1d0,BPD%kr,BPD%Order,BPD%xPoints,LegPoints,xLeg, &
         xDim2,xBounds2,BPD%xNumPoints,1,uxAux)

    ! log-derivative of the FIXED combination u1=B1+B2 at R=xl (needed by CalcGamLamBox1's
    ! missing-surface-term correction; this combination is energy-independent, unlike the
    ! Robin-matched u1 in BuildBox1CombinedBasis).
    uval = MyBSpline(BPD%Order,0,BPD%xNumPoints,BPD%xPoints,1,BPD%xPoints(1)) &
         + MyBSpline(BPD%Order,0,BPD%xNumPoints,BPD%xPoints,2,BPD%xPoints(1))
    uderiv = MyBSpline(BPD%Order,1,BPD%xNumPoints,BPD%xPoints,1,BPD%xPoints(1)) &
         + MyBSpline(BPD%Order,1,BPD%xNumPoints,BPD%xPoints,2,BPD%xPoints(1))

    DO mch = 1,BPD%NumChannels
       u1(:,:,mch) = uAux(:,:,1) + uAux(:,:,2)
       ux1(:,:,mch) = uxAux(:,:,1) + uxAux(:,:,2)
       kLeftOut(mch) = uderiv/uval
    ENDDO

    DEALLOCATE(uAux,uxAux)
  END SUBROUTINE BuildBox1NeumannBasis
  !****************************************************************************************************
  ! Regular solution of the near-origin equation u''+(1/R)u'+2*mu*(EE-Pot(R))u=0 for a
  ! channel whose effective potential Pot(R)=U(R)+Q_diag(R)/(2*mu) has a genuine -C0/R
  ! (Coulomb-like) leading divergence, as opposed to the pure centrifugal 1/R^2 structure
  ! hyperrjry assumes.
  !
  ! An earlier version of this routine used the IDEALIZED potential -C0/R everywhere (a
  ! closed-form power series in a fictitious k^2). That is wrong at the finite xMin this
  ! is actually used at: fitting the real combo(R)=U(R)+Q11(R)/(2*mu) data (delta-function
  ! test case) to -C0/R+D over R in [1e-4,1e-2] gives D~-0.33 with tiny residuals (~1e-5
  ! relative) -- a genuine, non-negligible finite additive constant that is neither 0 nor
  ! the channel threshold, and not something to guess analytically. Rather than solve for
  ! D, integrate the EXACT equation with the REAL tabulated potential (UInterp/QInterp)
  ! numerically via RK4, starting from a small R0 deep in the region where Pot(R)->-C0/R
  ! is an excellent approximation (confirmed by the fit above), so the leading-order
  ! regular-series IC u(R0)=1-2*mu*C0*R0, u'(R0)=-2*mu*C0 is valid there, and integrating
  ! out to R1=xMin picks up whatever D and sub-leading structure the real potential has.
  SUBROUTINE CoulombRegularNumeric(mu,C0,EE,mch,NCh,UInterp,QInterp,R0,R1,uval,uderiv)
    USE InterpType
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: mu, C0, EE, R0, R1
    INTEGER, INTENT(in) :: mch, NCh
    TYPE(InterpolatingFunction), INTENT(in) :: UInterp(NCh)
    TYPE(InterpolatingMatrix), INTENT(in) :: QInterp
    DOUBLE PRECISION, INTENT(out) :: uval, uderiv

    INTEGER, PARAMETER :: nsteps = 4000
    DOUBLE PRECISION h, R, y(2), k1(2), k2(2), k3(2), k4(2)
    INTEGER istep

    y(1) = 1d0 - 2d0*mu*C0*R0
    y(2) = -2d0*mu*C0
    R = R0
    h = (R1-R0)/DBLE(nsteps)
    DO istep = 1,nsteps
       CALL CoulombRHS(R,        y,            mu,C0,EE,mch,NCh,UInterp,QInterp,k1)
       CALL CoulombRHS(R+0.5d0*h,y+0.5d0*h*k1, mu,C0,EE,mch,NCh,UInterp,QInterp,k2)
       CALL CoulombRHS(R+0.5d0*h,y+0.5d0*h*k2, mu,C0,EE,mch,NCh,UInterp,QInterp,k3)
       CALL CoulombRHS(R+h,      y+h*k3,       mu,C0,EE,mch,NCh,UInterp,QInterp,k4)
       y = y + (h/6d0)*(k1 + 2d0*k2 + 2d0*k3 + k4)
       R = R + h
    ENDDO
    uval = y(1)
    uderiv = y(2)
  END SUBROUTINE CoulombRegularNumeric
  !****************************************************************************************************
  ! Below RcutoffAnalytic, use the exact small-R series (AnalyticComboNearOrigin) instead of
  ! the tabulated-data spline interpolation, which has significant Runge-phenomenon ringing
  ! there (confirmed: 10-31% errors in combo(R)=2*mu*U+Q for R<0.005, vanishing by R~0.01).
  SUBROUTINE CoulombRHS(R,y,mu,C0,EE,mch,NCh,UInterp,QInterp,dydt)
    USE InterpType
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: R, y(2), mu, C0, EE
    INTEGER, INTENT(in) :: mch, NCh
    TYPE(InterpolatingFunction), INTENT(in) :: UInterp(NCh)
    TYPE(InterpolatingMatrix), INTENT(in) :: QInterp
    DOUBLE PRECISION, INTENT(out) :: dydt(2)
    DOUBLE PRECISION, EXTERNAL :: AnalyticComboNearOrigin
    DOUBLE PRECISION, PARAMETER :: RcutoffAnalytic = 0.01d0
    DOUBLE PRECISION Pot
    DOUBLE PRECISION QM(QInterp%nr,QInterp%nc)

    IF (R.LT.RcutoffAnalytic) THEN
       Pot = AnalyticComboNearOrigin(R,mu,C0)/(2d0*mu)
    ELSE
       QM = InterpolatedMatrix(R,QInterp)
       Pot = Interpolated(R,UInterp(mch)) + QM(mch,mch)/(2d0*mu)
    ENDIF
    dydt(1) = y(2)
    dydt(2) = -(1d0/R)*y(2) - 2d0*mu*(EE-Pot)*y(1)
  END SUBROUTINE CoulombRHS
  !****************************************************************************************************
  ! Box-1 analog of CalcGamLam: identical term-by-term derivation (see CalcGamLam), except
  ! spatial index ix=1 -- BPD's own single Left=3 boundary function, whatever reference
  ! basis it was built with -- is replaced, PER CHANNEL, by BuildBox1CombinedBasis's u1/ux1.
  ! BPD%kxMin/BPD%kxMax (from BPD's own Left=3 MakeBasis call) are already correct.
  SUBROUTINE CalcGamLamBox1(BPD,EIG,NumOpenL,NumOpenR,u1,ux1,kLeft)
    USE DataStructures
    USE Quadrature
    USE GlobalVars
    IMPLICIT NONE
    TYPE(BPData), INTENT(in) :: BPD
    TYPE(GenEigVal) :: EIG
    INTEGER, INTENT(in) :: NumOpenL, NumOpenR
    DOUBLE PRECISION, INTENT(in) :: u1(LegPoints,BPD%xNumPoints,BPD%NumChannels)
    DOUBLE PRECISION, INTENT(in) :: ux1(LegPoints,BPD%xNumPoints,BPD%NumChannels)
    DOUBLE PRECISION, INTENT(in) :: kLeft(BPD%NumChannels)
    INTEGER ix,ixp,lx,kx,mch,nch
    DOUBLE PRECISION a, ax, bx, xIntScale(BPD%xNumPoints), TempG, TempS, xScaledZero
    DOUBLE PRECISION uRow(LegPoints), uxRow(LegPoints), uCol(LegPoints), uxCol(LegPoints)
    EIG%Gam0=0d0
    EIG%Overlap=0d0
    EIG%Lam=0d0
    EIG%eval=0d0
    EIG%evec=0d0

    DO ix = 1,BPD%xDim
       DO ixp = MAX(1,ix-BPD%Order),MIN(BPD%xDim,ix+BPD%Order)
          DO mch = 1,BPD%NumChannels
             ! Do the channel-diagonal part first (row channel = col channel = mch):
             DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                ax = BPD%xPoints(kx)
                bx = BPD%xPoints(kx+1)
                xIntScale(kx) = 0.5d0*(bx-ax)
                xScaledZero = 0.5d0*(bx+ax)
                IF (ix.EQ.1) THEN
                   uRow = u1(1:LegPoints,kx,mch); uxRow = ux1(1:LegPoints,kx,mch)
                ELSE
                   uRow = BPD%u(1:LegPoints,kx,ix); uxRow = BPD%ux(1:LegPoints,kx,ix)
                ENDIF
                IF (ixp.EQ.1) THEN
                   uCol = u1(1:LegPoints,kx,mch); uxCol = ux1(1:LegPoints,kx,mch)
                ELSE
                   uCol = BPD%u(1:LegPoints,kx,ixp); uxCol = BPD%ux(1:LegPoints,kx,ixp)
                ENDIF
                TempG = 0.0d0
                TempS = 0.0d0
                DO lx = 1,LegPoints
                   a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                   TempG = TempG - a*(uxRow(lx)*uxCol(lx))
                   a = a*uRow(lx)*uCol(lx)*2d0*reducedmass
                   TempS = TempS + a
                   TempG = TempG - a*BPD%Pot(mch,mch,lx,kx)
                   TempG = TempG - a*AlphaFactor*(AlphaFactor-EffDim+2d0)/(2d0*reducedmass*BPD%x(lx,kx)**2)
                ENDDO
                EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) = EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) + TempG
                EIG%Overlap((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) = EIG%Overlap((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) + TempS
             ENDDO
             IF (ix.EQ.1 .AND. ixp.EQ.1) THEN
                ! Missing surface term from integrating the kinetic energy by parts:
                ! int w(R) u_i u_j'' dR = [w(R) u_i u_j']_boundary - int w(R) u_i'u_j' dR.
                ! The kx-loop above only ever assembles the bulk (second) piece. At the
                ! right edge the boundary piece is supplied by Lam (validated). At the
                ! left edge, standard Dirichlet/Neumann boundary functions have
                ! u(xl)*u'(xl)=0 identically (one factor vanishes), so this term has
                ! been silently absent with no consequence anywhere else in this code.
                ! For the Robin-combined u1, u1(xl)=1 and u1'(xl)=kLeft(mch) exactly by
                ! construction, so this term is w(xl)*kLeft(mch) -- genuinely nonzero.
                ! This is a BOUNDARY term (evaluated at the single point xl), so it must
                ! be added exactly ONCE per (ix,ixp,mch) -- NOT once per kx element of
                ! u1's support (u1 spans kxMin(1,1)=1 to kxMax(1,1)=2, i.e. 2 elements;
                ! placing this correction inside the kx-loop above added it twice,
                ! introducing a spurious extra -w(xl)*kLeft(mch) and giving an
                ! effective log-derivative roughly double the imposed value -- found via
                ! a direct dump-and-diff of EIG%Gam between a hand-rolled reference
                ! (fix applied once, outside any loop) and this subroutine.
                EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) = &
                     EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+mch) &
                     - BPD%xl**(EffDim-1d0-2d0*AlphaFactor)*kLeft(mch)
             ENDIF

             ! Off-diagonal Q-coupling (row channel mch, col channel nch):
             DO nch = 1, mch-1
                DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                   ax = BPD%xPoints(kx)
                   bx = BPD%xPoints(kx+1)
                   xIntScale(kx) = 0.5d0*(bx-ax)
                   xScaledZero = 0.5d0*(bx+ax)
                   IF (ix.EQ.1) THEN
                      uRow = u1(1:LegPoints,kx,mch)
                   ELSE
                      uRow = BPD%u(1:LegPoints,kx,ix)
                   ENDIF
                   IF (ixp.EQ.1) THEN
                      uCol = u1(1:LegPoints,kx,nch)
                   ELSE
                      uCol = BPD%u(1:LegPoints,kx,ixp)
                   ENDIF
                   TempG = 0.0d0
                   DO lx = 1,LegPoints
                      a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                      TempG = TempG - a*uRow(lx)*2.0d0*reducedmass*BPD%Pot(mch,nch,lx,kx)*uCol(lx)
                   ENDDO
                   EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) = EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) + TempG
                   EIG%Gam0((ix-1)*BPD%NumChannels+nch,(ixp-1)*BPD%NumChannels+mch) = EIG%Gam0((ix-1)*BPD%NumChannels+nch,(ixp-1)*BPD%NumChannels+mch) + TempG
                ENDDO
             ENDDO

             ! Non-adiabatic P coupling (row channel mch, col channel nch):
             DO nch = 1,BPD%NumChannels
                IF (nch.NE.mch) THEN
                   DO kx = BPD%kxMin(ixp,ix),BPD%kxMax(ixp,ix)
                      ax = BPD%xPoints(kx)
                      bx = BPD%xPoints(kx+1)
                      xIntScale(kx) = 0.5d0*(bx-ax)
                      xScaledZero = 0.5d0*(bx+ax)
                      IF (ix.EQ.1) THEN
                         uRow = u1(1:LegPoints,kx,mch); uxRow = ux1(1:LegPoints,kx,mch)
                      ELSE
                         uRow = BPD%u(1:LegPoints,kx,ix); uxRow = BPD%ux(1:LegPoints,kx,ix)
                      ENDIF
                      IF (ixp.EQ.1) THEN
                         uCol = u1(1:LegPoints,kx,nch); uxCol = ux1(1:LegPoints,kx,nch)
                      ELSE
                         uCol = BPD%u(1:LegPoints,kx,ixp); uxCol = BPD%ux(1:LegPoints,kx,ixp)
                      ENDIF
                      TempG = 0.0d0
                      DO lx = 1,LegPoints
                         a = wLeg(lx)*xIntScale(kx)*BPD%x(lx,kx)**(EffDim-1d0-2d0*AlphaFactor)
                         TempG = TempG + a*BPD%Pcoup(mch,nch,lx,kx)*(uRow(lx)*uxCol(lx) - uxRow(lx)*uCol(lx))
                      ENDDO
                      EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) = &
                           EIG%Gam0((ix-1)*BPD%NumChannels+mch,(ixp-1)*BPD%NumChannels+nch) + TempG
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !Calculate the Lambda matrix elements (identical to CalcGamLam -- box 1 always has
    !NumOpenL=0, so these loops are vacuous regardless of Left=3's ix=1 combination)
    EIG%Lam = 0.d0
    IF(BPD%Right.EQ.2) THEN
       DO mch = 1, NumOpenR
          EIG%Lam( (BPD%xDim-1)*BPD%NumChannels+mch, (BPD%xDim-1)*BPD%NumChannels+mch ) = BPD%xr**(EffDim-1d0-2d0*AlphaFactor)
       ENDDO
    ENDIF
    DO mch = 1, NumOpenL
       EIG%Lam(mch,mch) = BPD%xl**(EffDim-1d0-2d0*AlphaFactor)
    ENDDO

  END SUBROUTINE CalcGamLamBox1
  !****************************************************************************************************
  ! Burke thesis Eq. 3.4-3.6 partitioned closed/open solve, replacing a single dense
  ! Mydggev over the whole box with a banded linear solve on the (large) closed block
  ! plus a small generalized eigenproblem on the (tiny) open block. Reuses CalcGamLam's
  ! already-built, already-validated dense Gam/Lam unchanged -- only the linear-algebra
  ! strategy applied to them changes, not the matrix elements themselves.
  !
  ! Open DOF (global index (ix-1)*NumChannels+ic, spatial-major per CABS.f90's
  ! convention -- see CalcGamLam) are exactly: ic=1..NumOpenL at ix=1 (global index ic),
  ! and ic=1..NumOpenR at ix=xDim (global index (xDim-1)*NumChannels+ic). Removing this
  ! prefix+notch+suffix and renumbering the remaining ("closed") indices 1..Ncc can only
  ! *decrease* the separation between any two closed indices relative to their original
  ! global separation (removing indices between two points never increases the gap
  ! between them), so the original HalfBandWidth=(Order+1)*NumChannels-1 (CABS.f90's
  ! formula) remains a safe upper bound for Gcc's bandwidth after renumbering -- possibly
  ! slightly conservative near the last box's notch, but never wrong.
  SUBROUTINE PartitionAndEliminate(BPD,EIG,NumOpenL,NumOpenR,evalRed,evecRed,LamooDiag)
    USE DataStructures
    USE GlobalVars
    IMPLICIT NONE
    TYPE(BPData), INTENT(in) :: BPD
    TYPE(GenEigVal), INTENT(in) :: EIG
    INTEGER, INTENT(in) :: NumOpenL, NumOpenR
    DOUBLE PRECISION, INTENT(out) :: evalRed(NumOpenL+NumOpenR)
    DOUBLE PRECISION, INTENT(out) :: evecRed(NumOpenL+NumOpenR,NumOpenL+NumOpenR)
    DOUBLE PRECISION, INTENT(out) :: LamooDiag(NumOpenL+NumOpenR)

    INTEGER NumOpen, Ncc, hb, iband, i, j, ic, jc, info
    INTEGER, ALLOCATABLE :: OpenIdx(:), ClosedIdx(:), ipiv(:)
    LOGICAL, ALLOCATABLE :: isOpen(:)
    DOUBLE PRECISION, ALLOCATABLE :: Gcc(:,:), Gco(:,:), Goc(:,:), Goo(:,:), Lamoo(:,:)

    NumOpen = NumOpenL + NumOpenR
    Ncc = BPD%MatrixDim - NumOpen
    hb = (BPD%Order+1)*BPD%NumChannels - 1   ! HalfBandWidth, same formula as CABS.f90
    iband = 3*hb + 1                          ! dgbsv storage: 2*kl+ku+1 with kl=ku=hb

    !--- classify global DOF as open (boundary-retained) or closed (eliminated) ---
    ALLOCATE(isOpen(BPD%MatrixDim))
    isOpen = .FALSE.
    DO ic = 1,NumOpenL
       isOpen(ic) = .TRUE.                                    ! (1-1)*NumChannels+ic
    ENDDO
    DO ic = 1,NumOpenR
       isOpen((BPD%xDim-1)*BPD%NumChannels+ic) = .TRUE.
    ENDDO
    ALLOCATE(OpenIdx(NumOpen),ClosedIdx(Ncc))
    i = 0; j = 0
    DO ic = 1,BPD%MatrixDim
       IF (isOpen(ic)) THEN
          i = i+1; OpenIdx(i) = ic
       ELSE
          j = j+1; ClosedIdx(j) = ic
       ENDIF
    ENDDO
    DEALLOCATE(isOpen)
    ! OpenIdx is filled left-channels-first (small global index) then right-channels
    ! (large global index), matching what BoxMatch expects for BB%Z's index order.

    !--- gather Gcc (banded, dgbsv format), Gco, Goo, Lamoo from the dense Gam/Lam ---
    ALLOCATE(Gcc(iband,Ncc),Gco(Ncc,NumOpen),Goc(NumOpen,Ncc))
    ALLOCATE(Goo(NumOpen,NumOpen),Lamoo(NumOpen,NumOpen))
    Gcc = 0d0
    DO j = 1,Ncc
       DO i = MAX(1,j-hb),MIN(Ncc,j+hb)
          Gcc(2*hb+1+i-j,j) = EIG%Gam(ClosedIdx(i),ClosedIdx(j))
       ENDDO
    ENDDO
    DO jc = 1,NumOpen
       DO i = 1,Ncc
          Gco(i,jc) = EIG%Gam(ClosedIdx(i),OpenIdx(jc))
          Goc(jc,i) = Gco(i,jc)   ! save transpose before dgbsv overwrites Gco
       ENDDO
    ENDDO
    Goo = 0d0
    Lamoo = 0d0
    DO ic = 1,NumOpen
       DO jc = 1,NumOpen
          Goo(ic,jc) = EIG%Gam(OpenIdx(ic),OpenIdx(jc))
       ENDDO
       Lamoo(ic,ic) = EIG%Lam(OpenIdx(ic),OpenIdx(ic))
    ENDDO

    !--- eliminate the closed block: Gcc*X=Gco (banded LU), then Omega_oo=Goo-Goc*X ---
    ALLOCATE(ipiv(Ncc))
    CALL dgbsv(Ncc,hb,hb,NumOpen,Gcc,iband,ipiv,Gco,Ncc,info)
    IF (info.NE.0) THEN
       WRITE(6,*) "PartitionAndEliminate: dgbsv failed, info=",info
       STOP
    ENDIF
    CALL dgemm('N','N',NumOpen,NumOpen,Ncc,-1.0d0,Goc,NumOpen,Gco,Ncc,1.0d0,Goo,NumOpen)

    !--- small generalized eigenproblem on the reduced (open-DOF-only) system ---
    CALL Mydggev(NumOpen,Goo,NumOpen,Lamoo,NumOpen,evalRed,evecRed)
    DO ic = 1,NumOpen
       LamooDiag(ic) = Lamoo(ic,ic)
    ENDDO

    DEALLOCATE(OpenIdx,ClosedIdx,ipiv,Gcc,Gco,Goc,Goo,Lamoo)
  END SUBROUTINE PartitionAndEliminate
  !****************************************************************************************************
  ! Diagnostic companion to PartitionAndEliminate: identical computation, but also
  ! returns the closed-DOF back-substitution matrix X (Gcc*X=Gco, i.e. c_closed=-X*c_open
  ! for any open-DOF coefficient vector c_open) and the OpenIdx/ClosedIdx global-index
  ! maps, so a caller can reconstruct the FULL (open+closed) spline coefficient vector for
  ! a chosen eigen-solution branch -- not just the open/retained boundary amplitudes
  ! PartitionAndEliminate normally returns. Used to actually plot the propagated
  ! wavefunction, not just a single boundary basis function in isolation.
  SUBROUTINE PartitionAndEliminateFull(BPD,EIG,NumOpenL,NumOpenR,evalRed,evecRed,LamooDiag, &
       Xout,OpenIdxOut,ClosedIdxOut)
    USE DataStructures
    USE GlobalVars
    IMPLICIT NONE
    TYPE(BPData), INTENT(in) :: BPD
    TYPE(GenEigVal), INTENT(in) :: EIG
    INTEGER, INTENT(in) :: NumOpenL, NumOpenR
    DOUBLE PRECISION, INTENT(out) :: evalRed(NumOpenL+NumOpenR)
    DOUBLE PRECISION, INTENT(out) :: evecRed(NumOpenL+NumOpenR,NumOpenL+NumOpenR)
    DOUBLE PRECISION, INTENT(out) :: LamooDiag(NumOpenL+NumOpenR)
    DOUBLE PRECISION, INTENT(out) :: Xout(BPD%MatrixDim-NumOpenL-NumOpenR,NumOpenL+NumOpenR)
    INTEGER, INTENT(out) :: OpenIdxOut(NumOpenL+NumOpenR), ClosedIdxOut(BPD%MatrixDim-NumOpenL-NumOpenR)

    INTEGER NumOpen, Ncc, hb, iband, i, j, ic, jc, info
    INTEGER, ALLOCATABLE :: OpenIdx(:), ClosedIdx(:), ipiv(:)
    LOGICAL, ALLOCATABLE :: isOpen(:)
    DOUBLE PRECISION, ALLOCATABLE :: Gcc(:,:), Gco(:,:), Goc(:,:), Goo(:,:), Lamoo(:,:)

    NumOpen = NumOpenL + NumOpenR
    Ncc = BPD%MatrixDim - NumOpen
    hb = (BPD%Order+1)*BPD%NumChannels - 1
    iband = 3*hb + 1

    ALLOCATE(isOpen(BPD%MatrixDim))
    isOpen = .FALSE.
    DO ic = 1,NumOpenL
       isOpen(ic) = .TRUE.
    ENDDO
    DO ic = 1,NumOpenR
       isOpen((BPD%xDim-1)*BPD%NumChannels+ic) = .TRUE.
    ENDDO
    ALLOCATE(OpenIdx(NumOpen),ClosedIdx(Ncc))
    i = 0; j = 0
    DO ic = 1,BPD%MatrixDim
       IF (isOpen(ic)) THEN
          i = i+1; OpenIdx(i) = ic
       ELSE
          j = j+1; ClosedIdx(j) = ic
       ENDIF
    ENDDO
    DEALLOCATE(isOpen)

    ALLOCATE(Gcc(iband,Ncc),Gco(Ncc,NumOpen),Goc(NumOpen,Ncc))
    ALLOCATE(Goo(NumOpen,NumOpen),Lamoo(NumOpen,NumOpen))
    Gcc = 0d0
    DO j = 1,Ncc
       DO i = MAX(1,j-hb),MIN(Ncc,j+hb)
          Gcc(2*hb+1+i-j,j) = EIG%Gam(ClosedIdx(i),ClosedIdx(j))
       ENDDO
    ENDDO
    DO jc = 1,NumOpen
       DO i = 1,Ncc
          Gco(i,jc) = EIG%Gam(ClosedIdx(i),OpenIdx(jc))
          Goc(jc,i) = Gco(i,jc)
       ENDDO
    ENDDO
    Goo = 0d0
    Lamoo = 0d0
    DO ic = 1,NumOpen
       DO jc = 1,NumOpen
          Goo(ic,jc) = EIG%Gam(OpenIdx(ic),OpenIdx(jc))
       ENDDO
       Lamoo(ic,ic) = EIG%Lam(OpenIdx(ic),OpenIdx(ic))
    ENDDO

    ALLOCATE(ipiv(Ncc))
    CALL dgbsv(Ncc,hb,hb,NumOpen,Gcc,iband,ipiv,Gco,Ncc,info)
    IF (info.NE.0) THEN
       WRITE(6,*) "PartitionAndEliminateFull: dgbsv failed, info=",info
       STOP
    ENDIF
    Xout = Gco   ! Gcc*Gco(original)=Gco(post-solve)=X: c_closed = -X*c_open
    OpenIdxOut = OpenIdx
    ClosedIdxOut = ClosedIdx

    CALL dgemm('N','N',NumOpen,NumOpen,Ncc,-1.0d0,Goc,NumOpen,Gco,Ncc,1.0d0,Goo,NumOpen)

    CALL Mydggev(NumOpen,Goo,NumOpen,Lamoo,NumOpen,evalRed,evecRed)
    DO ic = 1,NumOpen
       LamooDiag(ic) = Lamoo(ic,ic)
    ENDDO

    DEALLOCATE(OpenIdx,ClosedIdx,ipiv,Gcc,Gco,Goc,Goo,Lamoo)
  END SUBROUTINE PartitionAndEliminateFull
  !****************************************************************************************************
  MODULE Scattering
    USE DataStructures
    USE GlobalVars
  CONTAINS
    SUBROUTINE CalcK(B,BPD,SD,mu,d,alpha,EE,Eth)
      IMPLICIT NONE
      ! TEMP diagnostic toggle: use order=0 (plain sin/cos, matching legacy
      ! adrmatprop.f's outer matching for the dimer-threshold channel) instead
      ! of BPD%lam(i)=Leff in the outer hyperrjry call only -- BuildBox1CombinedBasis's
      ! own Leff-based near-origin regularity matching is untouched.
      LOGICAL, PARAMETER :: UseOuterOrderZero = .FALSE.
      TYPE(BoxData), INTENT(in) :: B
      TYPE(BPData), INTENT(in) :: BPD
      TYPE(ScatData) :: SD
      DOUBLE PRECISION mu, EE,rm,d,alpha
      DOUBLE PRECISION, ALLOCATABLE :: s(:),c(:),sp(:),cp(:)
      DOUBLE PRECISION, ALLOCATABLE :: Rmat(:,:),tmp1(:,:),tmp2(:,:)
      DOUBLE PRECISION rhypj,rhypy,rhypjp,rhypyp
      DOUBLE PRECISION k(BPD%NumChannels),Eth(BPD%NumChannels)
      INTEGER i,j,no,nw,nc,beta

      rm=BPD%xr
      no=0
      nw=0
      nc=0

      DO i = 1,BPD%NumChannels
         IF (2d0*mu*EE.GE.Eth(i)) THEN   ! must match the k(i)^2=2*mu*EE-Eth(i) convention below
            k(i) = dsqrt(2d0*mu*(EE)-Eth(i)) ! k is real
            no=no+1
         ELSE
            k(i) = dsqrt(Eth(i)-2d0*mu*EE) ! k->kappa, kappa is real
            IF( (k(i)*rm).LT.10d0) nw = nw+1
            IF( (k(i)*rm).GE.10d0) nc = nc+1
         ENDIF
      ENDDO
!      write(6,*) "no = ", no
      IF((no+nw+nc).NE.BPD%NumChannels) THEN
         WRITE(6,*) "Channel miscount in calcK"
         STOP
      ENDIF

      ALLOCATE(s(no),c(no))
      ALLOCATE(sp(no),cp(no))
      DO i = 1,no
         IF (UseOuterOrderZero) THEN
            CALL hyperrjry(INT(d),alpha,0d0,k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
         ELSE
            CALL hyperrjry(INT(d),alpha,BPD%lam(i),k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
         ENDIF
         !s(i) = dsqrt(mu)*rhypj  ! the factor of sqrt(mu) is for energy normalization
         !c(i) = -dsqrt(mu)*rhypy ! the factor of sqrt(mu) is for energy normalization
         !sp(i) = k(i)*dsqrt(mu)*rhypjp
         !cp(i) = -k(i)*dsqrt(mu)*rhypyp
          s(i) = 1d0/dsqrt(Pi*k(i))*rhypj  ! energy-normalized: has the required 1/sqrt(k(i))
          c(i) = -1d0/dsqrt(Pi*k(i))*rhypy
          sp(i) = k(i)*1d0/dsqrt(Pi*k(i))*rhypjp
          cp(i) = -k(i)*1d0/dsqrt(Pi*k(i))*rhypyp

      ENDDO
      ! K-matrix construction now matches adrmatprop.f (2006/2007 legacy code) exactly:
      ! build the explicit surface R-matrix (PLUS sign: R(i,j)=+sum_beta Zf*Zf/bf,
      ! Baluja et al's own convention -- our earlier Imat/Jmat eigenchannel formula used
      ! an equivalent but structurally different combination), then combine with the
      ! outer functions via tmp1=f*delta+fp*R, tmp2=g*delta+gp*R, K=tmp1*tmp2^-1, where
      ! legacy's f (cosine-like, irregular) is our c, and legacy's g (sine-like, regular)
      ! is our s (verified against both codes' own outer-matching formulas: legacy's
      ! f~cos(kR), g~sin(kR); ours has s~sin-like via rhypj, c~cos-like via -rhypy).
      ! The relative PLUS sign on cp*R/sp*R (not MINUS) is deliberate -- adrmatprop.f's
      ! own comment flags this as the fix needed to reconcile Baluja et al's convention
      ! against the "Orange Review" (Aymar/Greene/Luc-Koenig)'s, and we're matching
      ! Baluja et al here.
      ALLOCATE(Rmat(no,no),tmp1(no,no),tmp2(no,no))
      Rmat=0d0
      DO i=1,no
         DO j=1,no
            DO beta=1,no
               Rmat(i,j) = Rmat(i,j) + B%Zf(i,beta)*B%Zf(j,beta)/B%bf(beta)
            ENDDO
         ENDDO
      ENDDO

      tmp1=0d0
      tmp2=0d0
      DO i=1,no
         DO j=1,no
            tmp1(i,j) = cp(i)*Rmat(i,j)
            tmp2(i,j) = sp(i)*Rmat(i,j)
         ENDDO
         tmp1(i,i) = tmp1(i,i) + c(i)
         tmp2(i,i) = tmp2(i,i) + s(i)
      ENDDO
      CALL SqrMatInv(tmp2, no)
      SD%K = MATMUL(tmp1,tmp2)

      DO i=1,B%NumOpenR
         DO j=1,B%NumOpenR
            SD%R(i,j) = Rmat(i,j)
         ENDDO
      ENDDO

      DEALLOCATE(s,c,sp,cp,Rmat,tmp1,tmp2)
    END SUBROUTINE CalcK

  END MODULE Scattering
!****************************************************************************************************
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)

  DO j = 1,nr
     WRITE(file,20) (M(j,k), k = 1,nc)
  ENDDO

20 FORMAT(1P,100D20.12)
30 FORMAT(100F12.6)
END SUBROUTINE printmatrix
!****************************************************************************************************
SUBROUTINE GridMakerLinear(xNumPoints,x1,x2,xPoints)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: xNumPoints
  DOUBLE PRECISION, INTENT(in) :: x1,x2
  DOUBLE PRECISION, INTENT(out) :: xPoints(xNumPoints)
  INTEGER i
  DOUBLE PRECISION xDelt
  xDelt = (x2-x1)/DBLE(xNumPoints-1)
  DO i = 1,xNumPoints
     xPoints(i) = (i-1)*xDelt + x1 ! Simple linear grid
  ENDDO
END SUBROUTINE GridMakerLinear
!****************************************************************************************************
! Same routine as ~/Documents/GitHub/AtomIon1D/AtomIon1D.f90 (and VeffAtomIon1D.f90,
! BackPropFT-AtomIon/milne_core.f90, TuneVeff.f90's GridMakerNew) -- not yet consolidated
! into ~/Documents/GitHub/lib/, copied here rather than linking cross-project.
SUBROUTINE GridMaker(grid,numpts,E1,E2,scale)
  IMPLICIT NONE
  DOUBLE PRECISION grid(numpts)
  DOUBLE PRECISION E1,E2,LE1,LE2,DE,LDE
  INTEGER numpts, iE
  CHARACTER(LEN=*), INTENT(IN) :: scale
  grid(1)=E1
  IF((scale.EQ."linear").and.(numpts.gt.1)) THEN
     DE=(E2-E1)/DBLE(numpts-1)
     DO iE=1,numpts
        grid(iE) = E1 + (iE-1)*DE
     ENDDO
  ENDIF
  IF((scale.EQ."log").and.(numpts.gt.1)) THEN
     LE1=dlog(E1)
     LE2=dlog(E2)
     LDE=(LE2-LE1)/DBLE(numpts-1d0)
     DO iE=1,numpts
        grid(iE) = dexp(LE1 + (iE-1)*LDE)
     ENDDO
  ENDIF
  IF((scale.EQ."quadratic").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**2*DE
     ENDDO
  ENDIF
  IF((scale.EQ."cubic").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**3*DE
     ENDDO
  ENDIF
  IF((scale.EQ."sqrroot").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**0.5d0*DE
     ENDDO
  ENDIF
END SUBROUTINE GridMaker
!****************************************************************************************************
SUBROUTINE BoxMatch(BA, BB, BPD, evalRed, evecRed, LamooDiag, dim, alphafact)
  ! Builds BB%Z/BB%ZP directly from PartitionAndEliminate's reduced (open-DOF-only)
  ! eval/evec, instead of extracting them from a full dense EIG%eval/EIG%evec via the
  ! old finite-eigenvalue filter (ikeep). evalRed/evecRed already contain exactly
  ! BB%betaMax=NumOpenL+NumOpenR entries (Lamoo is nonsingular by construction, so
  ! Mydggev never produces an "infinite" eigenvalue here), and evecRed's row index is
  ! already ordered left-channels-then-right-channels (see PartitionAndEliminate's
  ! OpenIdx), matching this loop's i=1..NumOpenL / i=1..NumOpenR structure directly.
  ! The normalization Nbeta is likewise just the Lamoo-weighted quadratic form, since
  ! the full Lam matrix is zero everywhere except at the open (boundary-retained) DOF.
  ! Everything below this extraction (the box-to-box matching via its own Mydggev call)
  ! is unchanged.
  USE DataStructures
  IMPLICIT NONE
  TYPE(BPData), INTENT(in) :: BPD
  TYPE(BoxData), INTENT(in) :: BA
  TYPE(BoxData) BB
  DOUBLE PRECISION, INTENT(in) :: evalRed(BB%betaMax), evecRed(BB%betaMax,BB%betaMax)
  DOUBLE PRECISION, INTENT(in) :: LamooDiag(BB%betaMax)

  DOUBLE PRECISION, ALLOCATABLE :: ZT(:,:), temp0(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: BigZA(:,:), BigZB(:,:),Dvec(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Dval(:)
  INTEGER, ALLOCATABLE :: Dikeep(:)
  DOUBLE PRECISION NewNorm, dim, alphafact, absd
  INTEGER i,j,k,beta,betaprime,nch,mch!,betamax

  BB%b = evalRed

  BB%Norm=0d0
  BB%Nbeta=0d0
  DO beta=1,BB%betaMax
    DO betaprime=1,BB%betaMax
      BB%Norm(beta,betaprime)=0d0
      DO i=1,BB%betaMax
        BB%Norm(beta,betaprime) = BB%Norm(beta,betaprime) + evecRed(i,beta)*LamooDiag(i)*evecRed(i,betaprime)
      ENDDO
    ENDDO
    BB%Nbeta(beta) = dsqrt(BB%Norm(beta,beta))
  ENDDO

  DO beta = 1,BB%betaMax
    DO i = 1,BB%NumOpenL
      BB%Z(i,beta) = BPD%xl**(0.5d0*(dim-1d0-2d0*alphafact))*&
      evecRed(i,beta)/BB%Nbeta(beta)!  left channels are OpenIdx(1..NumOpenL)
      BB%ZP(i,beta) = -BB%b(beta)*BB%Z(i,beta)
    ENDDO
    DO i = 1, BB%NumOpenR
      BB%Z(i+BB%NumOpenL,beta) = BPD%xr**(0.5d0*(dim-1d0-2d0*alphafact))*&
      evecRed(i+BB%NumOpenL,beta)/BB%Nbeta(beta) !  right channels are OpenIdx(NumOpenL+1..)
      BB%ZP(i+BB%NumOpenL,beta) = -BB%b(beta)*BB%Z(i+BB%NumOpenL,beta)
    ENDDO
  ENDDO

  ALLOCATE(ZT(BB%betaMax,BB%betaMax))
  ALLOCATE(temp0(BB%betaMax,BB%betaMax))
  ZT = BB%Z
  temp0 = 0d0
  temp0 = MATMUL(TRANSPOSE(ZT),BB%Z)

  ALLOCATE(Dikeep(BB%betaMax + BB%NumOpenL))
  ALLOCATE(BigZA(BB%NumOpenL + BB%betaMax,BB%NumOpenL+BB%betaMax),BigZB(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax)) !
  ALLOCATE(Dvec(BB%NumOpenL+BB%betaMax,BB%NumOpenL+BB%betaMax))
  ALLOCATE(Dval(BB%NumOpenL+BB%betaMax))

  Dvec=0d0
  BigZA=0d0
  BigZB=0d0

  DO i=1,BB%NumOpenL
    DO beta=1,BB%NumOpenL
      BigZA(i,beta)=BA%Zf(i,beta)
      BigZA(i+BB%NumOpenL,beta)=BA%Zfp(i,beta)
    ENDDO
  ENDDO

  DO i=1,BB%NumOpenL
    DO beta=1,BB%betaMax
      BigZA(i, BB%NumOpenL+beta) = -BB%Z(i,beta)
      BigZA(i+BB%NumOpenL,beta+BB%NumOpenL) = BB%ZP(i,beta)
    ENDDO
  ENDDO

  DO i=1,BB%NumOpenR
    DO beta=1,BB%betaMax
      BigZA(i+2*BB%NumOpenL, BB%NumOpenL+beta)=-BB%ZP(i+BB%NumOpenL,beta)
      BigZB(i+2*BB%NumOpenL, BB%NumOpenL+beta)=BB%Z(i+BB%NumOpenL,beta)
    ENDDO
  ENDDO

  CALL Mydggev(BB%NumOpenL+BB%betaMax,BigZA,BB%NumOpenL+BB%betaMax,BigZB,BB%NumOpenL+BB%betaMax,Dval,Dvec) !
  j=1
  DO i = 1,BB%NumOpenL+BB%betaMax
    absd=ABS(Dval(i))
    IF((absd.GE.1d-12).and.(absd.lt.1d12)) THEN
      Dikeep(j)=i
      BB%bf(j)=Dval(i)
      j=j+1
    ENDIF
  ENDDO

  DO beta=1,BB%NumOpenR
    !for each beta, construct the final state with constant log-derivative at the right boundary
    DO i=1, BB%NumOpenR
      BB%Zf(i,beta) = 0.0d0
      DO betaprime = 1,BB%betaMax
        BB%Zf(i,beta) =  BB%Zf(i,beta) + BB%Z(BB%NumOpenL+i,betaprime)*Dvec(BB%NumOpenL+betaprime,Dikeep(beta)) !
      ENDDO
    ENDDO
    ! calculate the normalization of the final state
    NewNorm=0.0d0
    DO i=1,BB%NumOpenR
      NewNorm=NewNorm+BB%Zf(i,beta)**2
    ENDDO
    NewNorm=dsqrt(NewNorm)
    ! normalize and set the derivative Zfp
    DO i=1,BB%NumOpenR
      BB%Zf(i,beta) = BB%Zf(i,beta)/NewNorm
      BB%Zfp(i,beta) = -BB%Zf(i,beta)*BB%bf(beta)
    ENDDO
  ENDDO

  DEALLOCATE(Dikeep,Dvec,Dval,BigZA,BigZB)
  DEALLOCATE(temp0)
  DEALLOCATE(ZT)
END SUBROUTINE BoxMatch
!****************************************************************************************************
! Analytic small-R series for combo(R) = 2*mu*U(R) + Q(R), for a genuine two-body
! bound-pair (Coulomb-like -C0/R) channel of the 3-body delta-function problem.
! Derived from the Kartavtsev-Malykh-Sofianos transcendental equation
! q*Tanh[phi*q] = b*R (a2=1 convention, phi=Pi/6), inverted as a series in R and
! substituted into the exact rational Qmat(1,1) formula (both done symbolically
! in Wolfram, not by hand). b relates to the already-computed CoulombC(mch)=C0
! via b = 2*mu*phi*C0 (verified: matches -C0/R exactly at leading order).
! Cross-checked against the tabulated Uad.dat/Qmat.dat data in the region
! (R>0.02) where that data is itself reliable -- matches to <0.001%. Exists
! because the cubic-spline interpolation of the tabulated data has significant
! Runge-phenomenon ringing for R below ~0.005-0.01 (confirmed directly against
! a fresh root-find of the transcendental equation), which corrupts
! CoulombRegularNumeric's RK4 integration since it starts as low as R0~2e-6.
DOUBLE PRECISION FUNCTION AnalyticComboNearOrigin(R,mu,C0)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(in) :: R, mu, C0
  DOUBLE PRECISION, PARAMETER :: phi = 0.523598775598298873d0  ! Pi/6
  DOUBLE PRECISION b

  b = 2d0*mu*phi*C0
  AnalyticComboNearOrigin = -b/(phi*R) + b**2*(phi**2-15d0)/45d0 &
       + b**3*phi*(phi**2-4d0)*R/45d0 + 4d0*b**4*phi**4*R**2/405d0
END FUNCTION AnalyticComboNearOrigin
