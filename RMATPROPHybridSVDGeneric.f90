!****************************************************************************************************
! Generalized drop-in analog of RMATPROPHybridSVDRobin.f90/RMATPROPHybridSVDAnalytic.f90: an
! all-SVD R-matrix propagator, built as a REUSABLE SUBROUTINE (RunHybridSVDGeneric) rather than a
! PROGRAM, so any system with a PotentialProc/GridMakerProcI/RobinCoeffProc plugin
! (lib/AdiabaticInterfaces.f90) can drive a full box-chained energy scan without forking a new
! PROGRAM per system. A Fortran PROGRAM can't itself take procedure arguments -- only a subroutine
! can -- so genuine reuse across systems (delta-function benchmark, atom-ion, ...) needs the actual
! propagation logic factored out of the PROGRAM entry point; each system gets a thin driver PROGRAM
! that sets its own physics (NumChannels/Threshold/Leff/mu/AlphaFactor/EffDim) and plugins, then
! calls this subroutine. See RMATPROPHybridSVDGenericDelta.f90 for the delta-function regression
! test built this way (validates this file reproduces RMATPROPHybridSVDRobin.x's own results
! exactly, before this machinery is trusted with new physics).
!
! Box structure and physics-independent parts of the pipeline are all unchanged from
! RMATPROPHybridSVDRobin.f90 -- see that file's own header for the full rationale. The only real
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
!****************************************************************************************************
MODULE RMATPROPHybridSVDGenericMod
  IMPLICIT NONE
CONTAINS
  SUBROUTINE RunHybridSVDGeneric(NumChannelsIn,mu,AlphaFactorIn,EffDimIn,Threshold,Leff, &
       NumSVDBoxes,xStartIn,xEndIn,BoxSpacingPower, &
       Emin,Emax,NumEnergies,EnergyGridType, &
       LSVD,Order1D,Left1D,Right1D,xNumPoints1D,xMin1D,xMax1D, &
       LegendreFile1D,LegPoints1D,ShiftInit,GridRebuildEveryR, &
       MyPotential,MyGridMaker,MyRobinCoeff,datadir, &
       PhaseShiftFile,KMatrixFile,RecombinationFile,Box1GridType)
    USE DataStructures
    USE GlobalVars
    USE Quadrature
    USE scattering
    USE AdiabaticInterfaces
    USE SVDChannelBasis, ONLY: RebuildSVDBoxLam, DeAllocateSVDBPD
    USE GenericSVDChannelBasis, ONLY: BuildSVDBoxGeneric
    IMPLICIT NONE
    ! NumChannelsIn/xStartIn/xEndIn (not NumChannels/xStart/xEnd) -- those bare names are
    ! GlobalVars' own public module variables (see RMatPropCore.f90's GlobalVars); reusing
    ! them as dummy-argument names here would shadow/collide with the module, corrupting
    ! unrelated later declarations (the exact "GlobalVars name-clash" gotcha documented in
    ! this repo's own history). Assigned into the real globals just below instead, since
    ! CombineGam/PartitionAndEliminate/BoxMatch/CalcK (RMatPropCore.f90) all read
    ! NumChannels/xStart/xEnd/NumBoxes/Energy as globals, not as passed arguments.
    INTEGER, INTENT(in) :: NumChannelsIn
    INTEGER, INTENT(in) :: NumSVDBoxes
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
    ! SVDChannelBasis.f90's own BuildSVDBox header for why these are NOT interchangeable
    ! (Radau structurally has no a1 DOF to impose Dirichlet with).
    CHARACTER(LEN=*), INTENT(in) :: Box1GridType

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
    INTEGER i, iBox, ie, nBranch, PhaseFile, KMatFile, RecombFile
    LOGICAL havePrevDelta
    COMPLEX*16, ALLOCATABLE :: Amat(:,:), Smat(:,:)
    DOUBLE PRECISION recombProb
    INTEGER MatrixDimSVD

    Pi = dacos(-1d0)
    LegendreFile = LegendreFile1D
    LegPoints = LegPoints1D
    IF (ALLOCATED(xLeg)) DEALLOCATE(xLeg)
    IF (ALLOCATED(wLeg)) DEALLOCATE(wLeg)
    ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
    CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)
    reducedmass = mu
    AlphaFactor = AlphaFactorIn
    EffDim = EffDimIn
    NumChannels = NumChannelsIn
    xStart = xStartIn
    xEnd = xEndIn

    WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," ddim = ",EffDim," mu = ",reducedmass
    WRITE(6,*) "Energy scan: Emin = ",Emin," Emax = ",Emax," NumEnergies = ",NumEnergies
    WRITE(6,*) "NumSVDBoxes = ",NumSVDBoxes," xStart = ",xStart," xEnd = ",xEnd," LSVD = ",LSVD

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
         a1 = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumSVDBoxes))**BoxSpacingPower
         a2 = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumSVDBoxes))**BoxSpacingPower
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

    ! Every box is SVD and the SVD region reaches the outer edge -- the last box's own BPD%lam
    ! carries Leff for CalcK's asymptotic Bessel matching.
    BPDSVDBoxes(NumSVDBoxes)%lam(1:NumChannels) = Leff(1:NumChannels)

    EIGSVD%MatrixDim = MatrixDimSVD
    CALL AllocateEIG(EIGSVD)

    havePrevDelta = .FALSE.
    NumBoxes = NumSVDBoxes

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

       ALLOCATE(Boxes(NumSVDBoxes))
       DO i = 1,NumSVDBoxes
          Boxes(i)%NumOpenL = MERGE(0,NumChannels,i.EQ.1)
          IF (i.EQ.NumSVDBoxes) THEN
             Boxes(i)%NumOpenR = NumOpenChannels
          ELSE
             Boxes(i)%NumOpenR = NumChannels
          ENDIF
          IF (i.EQ.1) THEN
             Boxes(i)%xl = xStart
          ELSE
             Boxes(i)%xl = Boxes(i-1)%xr
          ENDIF
          Boxes(i)%xr = xStart + (xEnd-xStart)*(DBLE(i)/DBLE(NumSVDBoxes))**BoxSpacingPower
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
            IF (iBox.EQ.NumSVDBoxes) THEN
               CALL RebuildSVDBoxLam(EIGSVD,ODiagLSVDBoxes(:,iBox),ODiagRSVDBoxes(:,iBox), &
                    wsq1SVDBoxes(iBox),wsqLSVDBoxes(iBox),BPDSVDBoxes(iBox)%xl,BPDSVDBoxes(iBox)%xr, &
                    NumChannels,BPDSVDBoxes(iBox)%xDim,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR)
            ELSE
               EIGSVD%Lam(1:MDim,1:MDim) = LamSVDBoxes(1:MDim,1:MDim,iBox)
            ENDIF
          END BLOCK
          CALL CombineGam(EIGSVD,Energy)
          ALLOCATE(evalRed(Boxes(iBox)%betaMax),evecRed(Boxes(iBox)%betaMax,Boxes(iBox)%betaMax))
          ALLOCATE(LamooDiag(Boxes(iBox)%betaMax))
          CALL PartitionAndEliminate(BPDSVDBoxes(iBox),EIGSVD,Boxes(iBox)%NumOpenL,Boxes(iBox)%NumOpenR, &
               evalRed,evecRed,LamooDiag)
          IF (iBox.EQ.1) THEN
             CALL BoxMatch(Bnull, Boxes(1), BPDSVDBoxes(1), evalRed, evecRed, LamooDiag, EffDim, AlphaFactor)
          ELSE
             CALL BoxMatch(Boxes(iBox-1), Boxes(iBox), BPDSVDBoxes(iBox), evalRed, evecRed, LamooDiag, &
                  EffDim, AlphaFactor)
          ENDIF
          DEALLOCATE(evalRed,evecRed,LamooDiag)
       ENDDO

       CALL CalcK(Boxes(NumSVDBoxes),BPDSVDBoxes(NumSVDBoxes),SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)

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

       DO i = 1,NumSVDBoxes
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

    CLOSE(PhaseFile)
    CLOSE(KMatFile)
    CLOSE(RecombFile)
    WRITE(6,*)
    WRITE(6,*) "Phase shift written to ",TRIM(PhaseShiftFile)

  END SUBROUTINE RunHybridSVDGeneric
END MODULE RMATPROPHybridSVDGenericMod
