!****************************************************************************************************
! Generalized drop-in analog of RMATPROPSVDRobin.f90/RMATPROPSVDAnalytic.f90: an
! all-SVD R-matrix propagator, built as a REUSABLE SUBROUTINE (RunHybridSVDGeneric) rather than a
! PROGRAM, so any system with a PotentialProc/GridMakerProcI/RobinCoeffProc plugin
! (lib/AdiabaticInterfaces.f90) can drive a full box-chained energy scan without forking a new
! PROGRAM per system. A Fortran PROGRAM can't itself take procedure arguments -- only a subroutine
! can -- so genuine reuse across systems (delta-function benchmark, atom-ion, ...) needs the actual
! propagation logic factored out of the PROGRAM entry point; each system gets a thin driver PROGRAM
! that sets its own physics (NumChannels/Threshold/Leff/mu/AlphaFactor/EffDim) and plugins, then
! calls this subroutine. See RMATPROPHybridSVDGenericDelta.f90 for the delta-function regression
! test built this way (validates this file reproduces RMATPROPSVDRobin.x's own results
! exactly, before this machinery is trusted with new physics).
!
! Box structure and physics-independent parts of the pipeline are all unchanged from
! RMATPROPSVDRobin.f90 -- see that file's own header for the full rationale. The only real
! change is BuildSVDBoxRobin -> GenericSVDChannelBasis.f90's BuildSVDBoxGeneric, and
! physics/thresholds/plugins becoming subroutine arguments instead of hardcoded module-level
! constants.
!
! Box 1's grid convention (Box1GridType, 'Radau' or 'LobattoDirichlet') is now a caller-supplied
! argument rather than hardcoded to 'Radau' -- 'Radau' (a1 excluded entirely, open/unconstrained
! boundary) is correct only for AlphaFactor=0 (the delta-function/atom-ion regularity condition
! this file was originally validated against); this repo's ordinary AlphaFactor=(EffDim-1)/2
! convention needs true Dirichlet u(a1)=0 instead, which requires 'LobattoDirichlet' -- Radau
! structurally has no DOF at a1 to impose Dirichlet with (see SVDChannelBasis.f90's own BuildSVDBox
! header). BuildSVDBoxGeneric/GenericSVDChannelBasis.f90 already implemented both; only this
! subroutine's box-1 setup was hardcoded.
!
! GENUINE HYBRID (B-spline tail) capability: boxes 1..NumSVDBoxes are SVD/DVR (as always); when
! the caller-supplied NumBoxesIn > NumSVDBoxes, boxes NumSVDBoxes+1..NumBoxesIn are ordinary
! B-spline boxes reading tabulated Uad/Pmat/Qmat via AdiabaticPotential.f90's interpolants (the
! SAME Fit.data/FitLeff.data/Uad.dat/Pmat.dat/Qmat.dat fitting procedure RMATPROPAdiabatic.x uses)
! -- exactly RMATPROPHybridSVD.f90's original box structure (see that file's own header), just
! reimplemented on top of the pluggable SVD engine instead of the hardcoded OneDimChannels one.
! Rswitch = xStart + (xEnd-xStart)*NumSVDBoxes/NumBoxesIn falls out of the same uniform box-edge
! spacing formula used throughout this pipeline. NumBoxesIn == NumSVDBoxes (the default for the
! three existing thin drivers) reproduces the previous all-SVD behavior exactly -- confirmed by
! construction, since that branch never touches any of the new tail-only code paths below.
! BoxMatch/PartitionAndEliminate are basis-agnostic (already relied on by RMATPROPHybridSVD.f90
! to chain SVD boxes directly into B-spline boxes), so the only new work is building each region's
! own Gam0/Overlap/Lam the way its own driver already knows how.
!****************************************************************************************************
MODULE RMATPROPHybridSVDGenericMod
  IMPLICIT NONE
CONTAINS
  SUBROUTINE RunHybridSVDGeneric(NumChannelsIn,mu,AlphaFactorIn,EffDimIn,Threshold,Leff, &
       NumSVDBoxes,NumBoxesIn,xStartIn,xEndIn,BoxSpacingPower, &
       Emin,Emax,NumEnergies,EnergyGridType, &
       LSVD,Order1D,Left1D,Right1D,xNumPoints1D,xMin1D,xMax1D, &
       LegendreFile1D,LegPoints1D,ShiftInit,GridRebuildEveryR, &
       MyPotential,MyGridMaker,MyRobinCoeff,datadir, &
       PhaseShiftFile,KMatrixFile,RecombinationFile,Box1GridType, &
       OrderTail,xNumPointsTail,LegPointsTail,UInterp,PInterp,QInterp,CoulombC)
    USE DataStructures
    USE GlobalVars
    USE Quadrature
    USE scattering
    USE AdiabaticInterfaces
    USE AdiabaticPotential
    USE SVDChannelBasis, ONLY: RebuildSVDBoxLam, DeAllocateSVDBPD
    USE GenericSVDChannelBasis, ONLY: BuildSVDBoxGeneric
    IMPLICIT NONE
    ! NumChannelsIn/xStartIn/xEndIn (not NumChannels/xStart/xEnd) -- those bare names are
    ! GlobalVars' own public module variables (see RMatPropCore.f90's GlobalVars); reusing
    ! them as dummy-argument names here would shadow/collide with the module, corrupting
    ! unrelated later declarations (the exact "GlobalVars name-clash" gotcha documented in
    ! this repo's own history). Assigned into the real globals just below instead, since
    ! CombineGam/PartitionAndEliminate/BoxMatch/CalcK (RMatPropCore.f90) all read
    ! NumChannels/xStart/xEnd/NumBoxes/Energy as globals, not as passed arguments. NumBoxesIn
    ! is the same story -- NumBoxes is a GlobalVars name too.
    INTEGER, INTENT(in) :: NumChannelsIn
    INTEGER, INTENT(in) :: NumSVDBoxes, NumBoxesIn
    DOUBLE PRECISION, INTENT(in) :: mu,AlphaFactorIn,EffDimIn
    DOUBLE PRECISION, INTENT(in) :: Threshold(NumChannelsIn), Leff(NumChannelsIn)
    DOUBLE PRECISION, INTENT(in) :: xStartIn,xEndIn,BoxSpacingPower
    DOUBLE PRECISION, INTENT(in) :: Emin,Emax
    INTEGER, INTENT(in) :: NumEnergies
    CHARACTER(LEN=*), INTENT(in) :: EnergyGridType
    INTEGER, INTENT(in) :: LSVD,Order1D,Left1D,Right1D,xNumPoints1D,LegPoints1D
    DOUBLE PRECISION, INTENT(in) :: xMin1D,xMax1D,ShiftInit
    CHARACTER(LEN=*), INTENT(in) :: LegendreFile1D
    LOGICAL, INTENT(in) :: GridRebuildEveryR
    PROCEDURE(PotentialProc) :: MyPotential
    PROCEDURE(GridMakerProcI) :: MyGridMaker
    PROCEDURE(RobinCoeffProc) :: MyRobinCoeff
    CHARACTER(LEN=*), INTENT(in) :: datadir
    CHARACTER(LEN=*), INTENT(in) :: PhaseShiftFile,KMatrixFile,RecombinationFile
    ! Grid convention for box 1 only (box 2..NumSVDBoxes are always 'Lobatto', both
    ! edges open, NOpL=NumChannels). 'Radau' (a1 excluded entirely, open/unconstrained
    ! there -- correct for AlphaFactor=0, e.g. the delta-function/atom-ion regularity
    ! condition) or 'LobattoDirichlet' (a1 IS a formal node whose DOF is eliminated,
    ! true Dirichlet u(a1)=0 -- correct for AlphaFactor=(EffDim-1)/2, this repo's
    ! ordinary convention per the top-level CLAUDE.md). Passed straight through to
    ! BuildSVDBoxGeneric/SVDChannelBasis.f90, which already implements both -- see
    ! SVDChannelBasis.f90's own GridType header note for why these are NOT interchangeable
    ! (Radau structurally has no a1 DOF to impose Dirichlet with).
    CHARACTER(LEN=*), INTENT(in) :: Box1GridType
    ! B-spline TAIL region (boxes NumSVDBoxes+1..NumBoxesIn) -- all OPTIONAL, and all
    ! required together exactly when NumBoxesIn>NumSVDBoxes (checked below). Order/
    ! xNumPoints/LegPoints here are the ordinary B-spline radial-quadrature knobs (same
    ! meaning as RMATPROPAdiabatic.x's own Order/xNumPoints/LegPoints), NOT to be confused
    ! with Order1D/xNumPoints1D/LegPoints1D above (the SVD region's ANGULAR resolution).
    ! UInterp/PInterp/QInterp/CoulombC are exactly AdiabaticPotential.f90's own
    ! ReadAdiabaticData/SetAdiabaticPotential arguments -- the caller reads them once
    ! (typically from the same Fit.data/FitLeff.data/Uad.dat/Pmat.dat/Qmat.dat dataset that
    ! also supplied Threshold/Leff above) and passes them straight through.
    INTEGER, INTENT(in), OPTIONAL :: OrderTail,xNumPointsTail,LegPointsTail
    TYPE(InterpolatingFunction), INTENT(in), OPTIONAL :: UInterp(:)
    TYPE(InterpolatingMatrix), INTENT(in), OPTIONAL :: PInterp, QInterp
    DOUBLE PRECISION, INTENT(in), OPTIONAL :: CoulombC(:)

    TYPE(BPData), ALLOCATABLE :: BPDSVDBoxes(:)
    TYPE(GenEigVal) EIGSVD
    TYPE(BoxData) Bnull
    TYPE(BoxData), ALLOCATABLE :: Boxes(:)
    TYPE(ScatData) SD
    DOUBLE PRECISION, ALLOCATABLE :: Gam0SVDBoxes(:,:,:), OverlapSVDBoxes(:,:,:), LamSVDBoxes(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: ODiagLSVDBoxes(:,:), ODiagRSVDBoxes(:,:), wsq1SVDBoxes(:), wsqLSVDBoxes(:)
    DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
    DOUBLE PRECISION Shift, qOurs, Kelem, rawDelta, delta, prevDelta
    INTEGER NumOpenChannels
    INTEGER i, iBox, ie, nBranch, PhaseFile, KMatFile, RecombFile, kx, lx
    LOGICAL havePrevDelta, haveTail
    COMPLEX*16, ALLOCATABLE :: Amat(:,:), Smat(:,:)
    DOUBLE PRECISION recombProb
    INTEGER MatrixDimSVD

    ! --- B-spline tail region (only allocated/used when haveTail) ---
    TYPE(BPData) BPD0Tail, BPDTail
    TYPE(GenEigVal) EIGTail
    DOUBLE PRECISION, ALLOCATABLE :: xprimTail(:)
    DOUBLE PRECISION, ALLOCATABLE :: Gam0TailBoxes(:,:,:), OverlapTailBoxes(:,:,:), LamTailBoxes(:,:,:)

    Pi = dacos(-1d0)
    reducedmass = mu
    AlphaFactor = AlphaFactorIn
    EffDim = EffDimIn
    NumChannels = NumChannelsIn
    xStart = xStartIn
    xEnd = xEndIn
    NumBoxes = NumBoxesIn

    haveTail = (NumBoxesIn.GT.NumSVDBoxes)
    IF (haveTail) THEN
       IF (.NOT.(PRESENT(OrderTail).AND.PRESENT(xNumPointsTail).AND.PRESENT(LegPointsTail) &
            .AND.PRESENT(UInterp).AND.PRESENT(PInterp).AND.PRESENT(QInterp).AND.PRESENT(CoulombC))) THEN
          WRITE(6,*) "RunHybridSVDGeneric: NumBoxesIn>NumSVDBoxes (a B-spline tail) requires ", &
               "OrderTail/xNumPointsTail/LegPointsTail/UInterp/PInterp/QInterp/CoulombC."
          STOP
       ENDIF
    ENDIF

    ! SVD region's own angular Legendre factors -- BuildSVDBoxGeneric takes LegendreFile1D/
    ! LegPoints1D as explicit arguments and manages its own internal Gauss-Legendre setup;
    ! it never reads GlobalVars' xLeg/wLeg/LegPoints (confirmed: no reference to them anywhere
    ! in GenericSVDChannelBasis.f90). These globals are set here only because some of the
    ! generic infrastructure below expects them allocated; the B-spline tail region (if any)
    ! re-sets them to its own LegPointsTail further down, after all SVD box construction is
    ! complete and cached -- safe, since nothing SVD-related is re-integrated afterward.
    LegendreFile = LegendreFile1D
    LegPoints = LegPoints1D
    IF (ALLOCATED(xLeg)) DEALLOCATE(xLeg)
    IF (ALLOCATED(wLeg)) DEALLOCATE(wLeg)
    ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
    CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

    WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," EffDim = ",EffDim," mu = ",reducedmass
    WRITE(6,*) "Energy scan: Emin = ",Emin," Emax = ",Emax," NumEnergies = ",NumEnergies
    WRITE(6,*) "NumSVDBoxes = ",NumSVDBoxes," NumBoxes = ",NumBoxes," xStart = ",xStart," xEnd = ",xEnd, &
         " LSVD = ",LSVD
    IF (haveTail) THEN
       WRITE(6,*) "B-spline tail: boxes ",NumSVDBoxes+1," to ",NumBoxesIn," Rswitch = ", &
            xStart + (xEnd-xStart)*DBLE(NumSVDBoxes)/DBLE(NumBoxesIn)
    ENDIF

    PhaseFile = 50
    OPEN(unit=PhaseFile,file=PhaseShiftFile,status='replace')
    WRITE(PhaseFile,'(A)') '# q  Energy  K  delta(rad, unwrapped, no anchor/offset applied)'
    WRITE(PhaseFile,'(A)') '@    xaxis label "q"'
    WRITE(PhaseFile,'(A)') '@    yaxis label "delta (rad)"'

    KMatFile = 51
    OPEN(unit=KMatFile,file=KMatrixFile,status='replace')

    RecombFile = 52
    OPEN(unit=RecombFile,file=RecombinationFile,status='replace')
    WRITE(RecombFile,'(A)') '# Energy  NumOpenChannels  ||S^dagger S - I||_max  (general unitarity check)'
    WRITE(RecombFile,'(A)') '@    xaxis label "Energy"'
    WRITE(RecombFile,'(A)') '@    yaxis label "||S^dagger S - I||_max"'

    !---------------------------------------------------------------------
    ! ONE-TIME setup: SVD boxes 1..NumSVDBoxes -- box 1 uses whatever GridType/NumOpenL the
    ! caller's own AlphaFactor/regularity convention needs (Radau for AlphaFactor=0, matching
    ! DeltaScat.f90/RMATPROPAdiabatic.f90's own Left=2/NumOpenL=0 convention; LobattoDirichlet
    ! for AlphaFactor=(EffDim-1)/2 -- see SVDChannelBasis.f90's own GridType header note).
    ! This subroutine hardcodes 'Radau' for box 1 -- callers needing AlphaFactor=(EffDim-1)/2's
    ! LobattoDirichlet convention instead would need that made a further argument; not needed
    ! by either the delta-function or atom-ion use case (both use AlphaFactor=0).
    !
    ! Box edges use NumBoxesIn (the FULL chain length, including any tail) as the spacing
    ! denominator, not NumSVDBoxes -- this is exactly what makes Rswitch fall out of the
    ! ordinary uniform box-edge formula (RMATPROPHybridSVD.f90's own convention): the SVD
    ! region packs boxes 1..NumSVDBoxes into [xStart,Rswitch], not [xStart,xEnd].
    !---------------------------------------------------------------------
    ALLOCATE(BPDSVDBoxes(NumSVDBoxes))
    MatrixDimSVD = LSVD*NumChannels
    ALLOCATE(Gam0SVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
    ALLOCATE(OverlapSVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
    ALLOCATE(LamSVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
    ALLOCATE(ODiagLSVDBoxes(NumChannels,NumSVDBoxes),ODiagRSVDBoxes(NumChannels,NumSVDBoxes))
    ALLOCATE(wsq1SVDBoxes(NumSVDBoxes),wsqLSVDBoxes(NumSVDBoxes))

    Shift = ShiftInit
    DO iBox = 1,NumSVDBoxes
       BLOCK
         DOUBLE PRECISION :: a1,a2
         INTEGER :: NOpL,NOpR,LThis
         CHARACTER(LEN=16) :: GType
         a1 = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxesIn))**BoxSpacingPower
         a2 = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxesIn))**BoxSpacingPower
         IF (iBox.EQ.1) THEN
            GType = Box1GridType
            LThis = LSVD
            NOpL = 0
         ELSE
            GType = 'Lobatto'
            LThis = LSVD
            NOpL = NumChannels
         ENDIF
         NOpR = NumChannels
         ! Diagonalizes the fixed-R angular Hamiltonian at each of L DVR nodes R_i in [a1,a2]
         ! (ARPACK, via MyPotential), giving adiabatic eigenvalues eps_nu(R_i); Gam0=diag(eps_nu)
         ! (+ kinetic coupling) and Overlap=O(R_i,R_j) (cross-node channel-overlap matrix) replace
         ! CalcGamLam's Galerkin double integral entirely for this box.
         CALL BuildSVDBoxGeneric(a1,a2,LThis,NumChannels,Order1D,Left1D,Right1D, &
              xNumPoints1D,xMin1D,xMax1D,LegendreFile1D,LegPoints1D,Shift,datadir, &
              GridRebuildEveryR,MyPotential,MyGridMaker,MyRobinCoeff, &
              NOpL,NOpR,GType,BPDSVDBoxes(iBox),EIGSVD, &
              ODiagLSVDBoxes(:,iBox),ODiagRSVDBoxes(:,iBox),wsq1SVDBoxes(iBox),wsqLSVDBoxes(iBox))
         Gam0SVDBoxes(1:LThis*NumChannels,1:LThis*NumChannels,iBox) = EIGSVD%Gam0
         OverlapSVDBoxes(1:LThis*NumChannels,1:LThis*NumChannels,iBox) = EIGSVD%Overlap
         LamSVDBoxes(1:LThis*NumChannels,1:LThis*NumChannels,iBox) = EIGSVD%Lam
         CALL DeAllocateEIG(EIGSVD)
         WRITE(6,'(A,I4,A,1P,E14.6,A,E14.6,A,A,A,I4)') " SVD box ",iBox," xl=",a1," xr=",a2, &
              " grid=",TRIM(GType)," L=",LThis
       END BLOCK
    ENDDO

    ! Only when the SVD region reaches all the way to the outer edge (no B-spline tail) does
    ! the last SVD box also play the role CalcK's asymptotic Bessel matching needs Leff for.
    IF (.NOT.haveTail) THEN
       BPDSVDBoxes(NumSVDBoxes)%lam(1:NumChannels) = Leff(1:NumChannels)
    ENDIF

    EIGSVD%MatrixDim = MatrixDimSVD
    CALL AllocateEIG(EIGSVD)

    !---------------------------------------------------------------------
    ! ONE-TIME setup: B-spline tail boxes NumSVDBoxes+1..NumBoxesIn (only when haveTail).
    ! Exactly RMATPROPHybridSVD.f90's own B-spline-tail construction (same CalcGamLam/
    ! SetAdiabaticPotential calls, same primitive-basis-on-[0,1]-rescaled-per-box pattern),
    ! just parameterized by the caller-supplied OrderTail/xNumPointsTail/LegPointsTail/
    ! UInterp/PInterp/QInterp/CoulombC instead of reading them from a hardcoded DataDir.
    !---------------------------------------------------------------------
    IF (haveTail) THEN
       LegPoints = LegPointsTail
       IF (ALLOCATED(xLeg)) DEALLOCATE(xLeg)
       IF (ALLOCATED(wLeg)) DEALLOCATE(wLeg)
       ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
       CALL GetGaussFactors(LegendreFile1D,LegPoints,xLeg,wLeg)

       BPD0Tail%NumChannels = NumChannels
       BPD0Tail%Order = OrderTail
       BPD0Tail%Left = 2
       BPD0Tail%Right = 2
       BPD0Tail%xNumPoints = xNumPointsTail
       BPD0Tail%kl = 1d20
       BPD0Tail%kr = 1d20
       CALL AllocateBPD(BPD0Tail)
       ALLOCATE(xprimTail(xNumPointsTail))
       CALL GridMakerLinear(BPD0Tail%xNumPoints,0d0,1d0,xprimTail)
       BPD0Tail%xPoints = xprimTail
       CALL MakeBasis(BPD0Tail)
       BPD0Tail%lam(1:NumChannels) = Leff(1:NumChannels)

       BPDTail = BPD0Tail

       EIGTail%MatrixDim = BPDTail%MatrixDim
       CALL AllocateEIG(EIGTail)

       IF (NumBoxesIn.GE.NumSVDBoxes+2) THEN
          ALLOCATE(Gam0TailBoxes(BPDTail%MatrixDim,BPDTail%MatrixDim,NumSVDBoxes+1:NumBoxesIn-1))
          ALLOCATE(OverlapTailBoxes(BPDTail%MatrixDim,BPDTail%MatrixDim,NumSVDBoxes+1:NumBoxesIn-1))
          ALLOCATE(LamTailBoxes(BPDTail%MatrixDim,BPDTail%MatrixDim,NumSVDBoxes+1:NumBoxesIn-1))
          DO iBox = NumSVDBoxes+1,NumBoxesIn-1
             BPDTail%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxesIn))**BoxSpacingPower
             BPDTail%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxesIn))**BoxSpacingPower
             BPDTail%xPoints = BPDTail%xl + (BPDTail%xr-BPDTail%xl)*xprimTail
             DO kx = 1,BPDTail%xNumPoints
                DO lx = 1,LegPointsTail
                   BPDTail%ux(lx,kx,1:BPDTail%xDim) = BPD0Tail%ux(lx,kx,1:BPDTail%xDim)/(BPDTail%xr-BPDTail%xl)
                ENDDO
             ENDDO
             ! Tail region switches to the tabulated/interpolated potential: Pot(m,n)=
             ! Qtilde_mn(R)/(2*mu)+delta_mn*U_m(R), Pcoup=P(R) (same as RMATPROPAdiabatic.x).
             CALL SetAdiabaticPotential(BPDTail,reducedmass,UInterp,PInterp,QInterp,CoulombC)
             ! Ordinary B-spline Galerkin Gam0/Overlap (same formula as the SVD box's Gam0/Overlap
             ! above, just built by quadrature over B-spline basis functions instead of DVR nodes).
             CALL CalcGamLam(BPDTail,EIGTail,NumChannels,NumChannels)
             Gam0TailBoxes(:,:,iBox) = EIGTail%Gam0
             OverlapTailBoxes(:,:,iBox) = EIGTail%Overlap
             LamTailBoxes(:,:,iBox) = EIGTail%Lam
             WRITE(6,'(A,I4,A,1P,E14.6,A,E14.6)') " B-spline tail box ",iBox," xl=",BPDTail%xl," xr=",BPDTail%xr
          ENDDO
       ENDIF
    ENDIF

    havePrevDelta = .FALSE.

    DO ie = 1,NumEnergies
       IF (NumEnergies.GT.1) THEN
          IF (EnergyGridType.EQ."quadratic") THEN
             Energy = Emin + (Emax-Emin)*(DBLE(ie-1)/DBLE(NumEnergies-1))**2
          ELSE
             Energy = Emin + (Emax-Emin)*DBLE(ie-1)/DBLE(NumEnergies-1)
          ENDIF
       ELSE
          Energy = Emin
       ENDIF

       NumOpenChannels = 0
       DO i = 1,NumChannels
          IF (2d0*reducedmass*Energy.GE.Threshold(i)) NumOpenChannels = NumOpenChannels+1
       ENDDO

       ALLOCATE(Boxes(NumBoxesIn))
       DO i = 1,NumBoxesIn
          Boxes(i)%NumOpenL = MERGE(0,NumChannels,i.EQ.1)
          IF (i.EQ.NumBoxesIn) THEN
             Boxes(i)%NumOpenR = NumOpenChannels
          ELSE
             Boxes(i)%NumOpenR = NumChannels
          ENDIF
          IF (i.EQ.1) THEN
             Boxes(i)%xl = xStart
          ELSE
             Boxes(i)%xl = Boxes(i-1)%xr
          ENDIF
          Boxes(i)%xr = xStart + (xEnd-xStart)*(DBLE(i)/DBLE(NumBoxesIn))**BoxSpacingPower
          CALL AllocateBox(Boxes(i))
          CALL InitZeroBox(Boxes(i))
       ENDDO

       Bnull%NumOpenL = 0
       Bnull%NumOpenR = 0
       CALL AllocateBox(Bnull)
       CALL InitZeroBox(Bnull)

       CALL AllocateScat(SD,NumChannels)

       DO iBox = 1,NumSVDBoxes
          BLOCK
            INTEGER :: MDim
            MDim = BPDSVDBoxes(iBox)%MatrixDim
            EIGSVD%Gam0(1:MDim,1:MDim) = Gam0SVDBoxes(1:MDim,1:MDim,iBox)
            EIGSVD%Overlap(1:MDim,1:MDim) = OverlapSVDBoxes(1:MDim,1:MDim,iBox)
            IF (iBox.EQ.NumBoxesIn) THEN
               ! Only possible when NumBoxesIn==NumSVDBoxes (no tail): this SVD box IS the
               ! terminal box, so its NumOpenR tracks NumOpenChannels(Energy) and the cached
               ! Lam (built with NumOpenR=NumChannels at setup time) isn't valid here.
               ! Cheaply refills just Lam_nunu=wsq*a^(D-1-2alpha)*O(a,a) at the moving boundary
               ! (a=a1 or a2) for the new NumOpenR, without repeating the ARPACK diagonalization.
               CALL RebuildSVDBoxLam(EIGSVD,ODiagLSVDBoxes(:,iBox),ODiagRSVDBoxes(:,iBox), &
                    wsq1SVDBoxes(iBox),wsqLSVDBoxes(iBox),BPDSVDBoxes(iBox)%xl,BPDSVDBoxes(iBox)%xr, &
                    NumChannels,BPDSVDBoxes(iBox)%xDim,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
            ELSE
               EIGSVD%Lam(1:MDim,1:MDim) = LamSVDBoxes(1:MDim,1:MDim,iBox)
            ENDIF
          END BLOCK
          ! Gam(E) = Gam0 + E*Overlap.
          CALL CombineGam(EIGSVD,Energy)
          ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
          ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
          ! Schur-eliminates interior/closed DOF (Omega_oo=Gam_oo-Gam_oc*Gam_cc^-1*Gam_co), then
          ! solves the small generalized eigenproblem Omega_oo*c=beta*Lam_oo*c on the retained
          ! ("open") boundary DOF -- evalRed=beta, evecRed=c.
          CALL PartitionAndEliminate(BPDSVDBoxes(iBox),EIGSVD,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR, &
               evalRed,evecRed,LamooDiag)
          ! Builds Z,Z'=-beta*Z and matches this box's boundary amplitudes to the previous box's
          ! via a generalized eigenvalue match enforcing R-matrix/log-derivative continuity
          ! (Burke thesis Eqs. 3.13-3.14) -- basis-agnostic, so an SVD box chains into a B-spline
          ! tail box (below) exactly the same way it chains into another SVD box.
          IF (iBox.EQ.1) THEN
             CALL BoxMatch(Bnull, Boxes(1), BPDSVDBoxes(1), evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
          ELSE
             CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPDSVDBoxes(iBox), evalRed, evecRed, LamooDiag, &
                  EffDim, AlphaFactor)
          ENDIF
          DEALLOCATE(evalRed,evecRed,LamooDiag)
       ENDDO

       IF (haveTail) THEN
          !--- interior B-spline tail boxes NumSVDBoxes+1..NumBoxesIn-1 ---
          DO iBox = NumSVDBoxes+1,NumBoxesIn-1
             EIGTail%Gam0 = Gam0TailBoxes(:,:,iBox)
             EIGTail%Overlap = OverlapTailBoxes(:,:,iBox)
             EIGTail%Lam = LamTailBoxes(:,:,iBox)
             CALL CombineGam(EIGTail,Energy)
             ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
             ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
             CALL PartitionAndEliminate(BPDTail,EIGTail,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR, &
                  evalRed,evecRed,LamooDiag)
             CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPDTail, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
             DEALLOCATE(evalRed,evecRed,LamooDiag)
          ENDDO

          !--- last box (always the tail's own outermost box, rebuilt fresh every energy
          !    since NumOpenR tracks NumOpenChannels(Energy)) ---
          iBox = NumBoxesIn
          BPDTail%xl = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxesIn))**BoxSpacingPower
          BPDTail%xr = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxesIn))**BoxSpacingPower
          BPDTail%xPoints = BPDTail%xl + (BPDTail%xr-BPDTail%xl)*xprimTail
          DO kx = 1,BPDTail%xNumPoints
             DO lx = 1,LegPointsTail
                BPDTail%ux(lx,kx,1:BPDTail%xDim) = BPD0Tail%ux(lx,kx,1:BPDTail%xDim)/(BPDTail%xr-BPDTail%xl)
             ENDDO
          ENDDO
          CALL SetAdiabaticPotential(BPDTail,reducedmass,UInterp,PInterp,QInterp,CoulombC)
          CALL CalcGamLam(BPDTail,EIGTail,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
          CALL CombineGam(EIGTail,Energy)
          ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
          ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
          CALL PartitionAndEliminate(BPDTail,EIGTail,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR, &
               evalRed,evecRed,LamooDiag)
          CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPDTail, evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
          DEALLOCATE(evalRed,evecRed,LamooDiag)

          ! Rmat(i,j)=sum_beta Zf(i,beta)Zf(j,beta)/bf(beta) at xEnd, matched to asymptotic
          ! Riccati-Bessel-like reference functions s,c (order=Leff): K=(c+cp*Rmat)*(s+sp*Rmat)^-1.
          CALL CalcK(Boxes(NumBoxesIn),BPDTail,SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
       ELSE
          ! Same K=(c+cp*Rmat)*(s+sp*Rmat)^-1 matching, using the terminal SVD box's own Zf/bf.
          CALL CalcK(Boxes(NumBoxesIn),BPDSVDBoxes(NumSVDBoxes),SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)
       ENDIF

       IF (NumOpenChannels.EQ.NumChannels .AND. NumChannels.GT.1) THEN
          WRITE(KMatFile,'(A,1P,E20.12)') '# Energy = ',Energy
          DO i = 1,NumOpenChannels
             WRITE(KMatFile,'(1P,20E20.12)') SD%K(i,1:NumOpenChannels)
          ENDDO
       ENDIF

       BLOCK
         COMPLEX*16, ALLOCATABLE :: SdagS(:,:)
         DOUBLE PRECISION :: maxdev
         INTEGER :: ii,jj
         ALLOCATE(Amat(NumOpenChannels,NumOpenChannels),Smat(NumOpenChannels,NumOpenChannels))
         ALLOCATE(SdagS(NumOpenChannels,NumOpenChannels))
         Amat = (0d0,0d0)
         Smat = (0d0,0d0)
         DO i = 1,NumOpenChannels
            Amat(i,1:NumOpenChannels) = -(0d0,1d0)*SD%K(i,1:NumOpenChannels)
            Smat(i,1:NumOpenChannels) = (0d0,1d0)*SD%K(i,1:NumOpenChannels)
            Amat(i,i) = Amat(i,i) + (1d0,0d0)
            Smat(i,i) = Smat(i,i) + (1d0,0d0)
         ENDDO
         CALL CompSqrMatInv(Amat,NumOpenChannels)
         Smat = MATMUL(Amat,Smat)
         SdagS = MATMUL(CONJG(TRANSPOSE(Smat)),Smat)
         maxdev = 0d0
         DO ii = 1,NumOpenChannels
            DO jj = 1,NumOpenChannels
               maxdev = MAX(maxdev,ABS(SdagS(ii,jj)-MERGE((1d0,0d0),(0d0,0d0),ii.EQ.jj)))
            ENDDO
         ENDDO
         recombProb = 1d0 - ABS(Smat(1,1))**2
         WRITE(RecombFile,'(1P,E20.12,I5,E20.12,A,E20.12)') Energy, NumOpenChannels, maxdev, &
              '  (1-|S11|^2=',recombProb
         DEALLOCATE(Amat,Smat,SdagS)
       END BLOCK

       IF (NumOpenChannels.EQ.1) THEN
          Kelem = SD%K(1,1)
          rawDelta = DATAN(Kelem)
          IF (havePrevDelta) THEN
             nBranch = NINT((prevDelta-rawDelta)/Pi)
             delta = rawDelta + DBLE(nBranch)*Pi
          ELSE
             delta = rawDelta
          ENDIF
          prevDelta = delta
          havePrevDelta = .TRUE.
          qOurs = DSQRT(2d0*reducedmass*Energy-Threshold(1))
          WRITE(PhaseFile,'(1P,4E20.12)') qOurs, Energy, Kelem, delta
          WRITE(6,'(A,I5,A,1P,E14.6,A,E14.6,A,E14.6)') " ie=",ie," Energy=",Energy," q=",qOurs," delta=",delta
       ELSE
          WRITE(6,'(A,I5,A,1P,E14.6,A,I3,A)') " ie=",ie," Energy=",Energy," NumOpenChannels=",NumOpenChannels, &
               " (phase shift skipped -- see Recombination.dat-style output for the unitarity check)"
       ENDIF

       DO i = 1,NumBoxesIn
          CALL DeAllocateBox(Boxes(i))
       ENDDO
       DEALLOCATE(Boxes)
       CALL DeAllocateBox(Bnull)
       CALL DeAllocateScat(SD)
    ENDDO

    DO iBox = 1,NumSVDBoxes
       CALL DeAllocateSVDBPD(BPDSVDBoxes(iBox))
    ENDDO
    DEALLOCATE(BPDSVDBoxes,Gam0SVDBoxes,OverlapSVDBoxes,LamSVDBoxes)
    DEALLOCATE(ODiagLSVDBoxes,ODiagRSVDBoxes,wsq1SVDBoxes,wsqLSVDBoxes)
    CALL DeAllocateEIG(EIGSVD)

    IF (haveTail) THEN
       CALL DeAllocateBPD(BPD0Tail)
       CALL DeAllocateBPD(BPDTail)
       DEALLOCATE(xprimTail)
       CALL DeAllocateEIG(EIGTail)
       IF (ALLOCATED(Gam0TailBoxes)) DEALLOCATE(Gam0TailBoxes,OverlapTailBoxes,LamTailBoxes)
    ENDIF

    CLOSE(PhaseFile)
    CLOSE(KMatFile)
    CLOSE(RecombFile)
    WRITE(6,*)
    WRITE(6,*) "Phase shift written to ",TRIM(PhaseShiftFile)

  END SUBROUTINE RunHybridSVDGeneric
END MODULE RMATPROPHybridSVDGenericMod
