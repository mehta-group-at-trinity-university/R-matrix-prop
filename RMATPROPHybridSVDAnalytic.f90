!****************************************************************************************************
! All-analytic-SVD R-matrix propagator for the equal-mass three-boson delta-function problem of
! Mehta, Esry & Greene, PRA 76, 022711 (2007) -- the SAME system DeltaScat.f90 solves exactly.
! EVERY box is an SVD (slow-variable-discretization) box built by AnalyticSVDChannelBasis.f90's
! BuildSVDBoxAnalytic, which gets its per-R channel data (Uad, and the channel-function overlap
! needed to relate different R's) directly from KMSFormulas.f's closed-form SolveQ/EvalU/
! EvalOverlapAt -- no B-spline basis, no OneDimChannels/ARPACK diagonalization, no interpolation,
! no basis-set error anywhere. Physics (mu/alpha/ddim/Threshold/Leff) is hardcoded to the exact
! delta-function values, matching DeltaScat.f90's own convention exactly -- NOT read from
! Fit.data/FitLeff.data/masses.dat (there is no tabulated dataset involved at all, so the
! FitLeff.data column-misalignment bug fixed in AdiabaticPotential.f90 is structurally impossible
! here: this driver never touches that file).
!
! Structurally a simplified NumSVDBoxes==NumBoxes special case of RMATPROPHybridSVD.f90 (box 1:
! Gauss-Lobatto DVR with the R=a1 node's DOF eliminated -- true Dirichlet regularity, matching
! AlphaFactor=(EffDim-1)/2=0.5 for this problem; boxes 2..NumBoxes: ordinary two-sided Lobatto DVR)
! with the entire B-spline tail removed, since every box here is SVD.
!
! Purpose: test the hybrid SVD/adiabatic propagator against the known-exact DeltaScat/Mehta-
! Shepard benchmark, the same way RMATPROPAdiabatic_seeded.f90 validated the B-spline pipeline.
!****************************************************************************************************
PROGRAM main
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE scattering
  USE SVDChannelBasis, ONLY: RebuildSVDBoxLam, DeAllocateSVDBPD
  USE AnalyticSVDChannelBasis, ONLY: BuildSVDBoxAnalytic

  IMPLICIT NONE
  CHARACTER(LEN=16) EnergyGridType
  INTEGER ios
  TYPE(BPData), ALLOCATABLE :: BPDSVDBoxes(:)
  TYPE(GenEigVal) EIGSVD
  TYPE(BoxData) Bnull
  TYPE(BoxData), ALLOCATABLE :: Boxes(:)
  TYPE(ScatData) SD
  DOUBLE PRECISION, ALLOCATABLE :: Gam0SVDBoxes(:,:,:), OverlapSVDBoxes(:,:,:), LamSVDBoxes(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: ODiagLSVDBoxes(:,:), ODiagRSVDBoxes(:,:), wsq1SVDBoxes(:), wsqLSVDBoxes(:)
  DOUBLE PRECISION, ALLOCATABLE :: evalRed(:), evecRed(:,:), LamooDiag(:)
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  DOUBLE PRECISION Emin, Emax, qOurs, Kelem, rawDelta, delta, prevDelta
  INTEGER NumOpenChannels
  INTEGER i, iBox, NumEnergies, ie, nBranch, PhaseFile, KMatFile, RecombFile
  LOGICAL havePrevDelta
  COMPLEX*16, ALLOCATABLE :: Amat(:,:), Smat(:,:)
  DOUBLE PRECISION recombProb
  DOUBLE PRECISION, PARAMETER :: BoxSpacingPower = 1.0d0
  DOUBLE PRECISION, PARAMETER :: DeltaMu = 0.57735026918962584d0
  DOUBLE PRECISION, PARAMETER :: DeltaAlpha = 0d0
  DOUBLE PRECISION, PARAMETER :: DeltaDdim = 2d0

  INTEGER NumSVDBoxes, LSVD, MatrixDimSVD

  !----------------------------------------------------------------------
  ! Read the run's numerical parameters. Physics (mu/alpha/ddim/Threshold/Leff) is hardcoded
  ! below to the exact delta-function values, not read from any file.
  !----------------------------------------------------------------------
  OPEN(unit=7,file='RMATPROPHybridSVDAnalytic.inp',status='old')
  READ(7,*)
  READ(7,*) NumChannels
  READ(7,*)
  READ(7,*)
  READ(7,*) NumSVDBoxes, xStart, xEnd
  READ(7,*)
  READ(7,*)
  READ(7,*) Emin, Emax, NumEnergies
  READ(7,*,IOSTAT=ios) EnergyGridType
  IF (ios.NE.0) EnergyGridType = "linear"
  READ(7,*)
  READ(7,*)
  READ(7,*) LSVD
  CLOSE(7)
  NumBoxes = NumSVDBoxes

  Pi = dacos(-1d0)
  reducedmass = DeltaMu
  AlphaFactor = DeltaAlpha
  EffDim = DeltaDdim

  ALLOCATE(Threshold(NumChannels),Leff(NumChannels))
  Threshold(1) = -1d0
  Leff(1) = 0.5d0
  DO i = 2,NumChannels
     Threshold(i) = 0d0
     Leff(i) = 6d0*DBLE(i) - 9d0
  ENDDO
  Threshold = 2d0*reducedmass*Threshold

  WRITE(6,*) "NumChannels = ",NumChannels," alpha = ",AlphaFactor," ddim = ",EffDim," mu = ",reducedmass
  WRITE(6,*) "Energy scan: Emin = ",Emin," Emax = ",Emax," NumEnergies = ",NumEnergies
  WRITE(6,*) "NumSVDBoxes = ",NumSVDBoxes," xStart = ",xStart," xEnd = ",xEnd," LSVD = ",LSVD

  PhaseFile = 50
  OPEN(unit=PhaseFile,file='PhaseShift.dat',status='replace')
  WRITE(PhaseFile,'(A)') '# q  Energy  K  delta(rad, unwrapped, no anchor/offset applied)'
  WRITE(PhaseFile,'(A)') '@    xaxis label "q"'
  WRITE(PhaseFile,'(A)') '@    yaxis label "delta (rad)"'

  KMatFile = 51
  OPEN(unit=KMatFile,file='KMatrixFull.dat',status='replace')

  RecombFile = 52
  OPEN(unit=RecombFile,file='Recombination.dat',status='replace')
  WRITE(RecombFile,'(A)') '# Energy  NumOpenChannels  ||S^dagger S - I||_max  (general unitarity check)'
  WRITE(RecombFile,'(A)') '@    xaxis label "Energy"'
  WRITE(RecombFile,'(A)') '@    yaxis label "||S^dagger S - I||_max"'

  !---------------------------------------------------------------------
  ! ONE-TIME setup: SVD boxes 1..NumSVDBoxes -- box 1: Gauss-Lobatto DVR with the a1=xStart
  ! node's DOF eliminated (true Dirichlet u(xStart)=0, matching AlphaFactor=(EffDim-1)/2=0.5
  ! for this problem -- see SVDChannelBasis.f90's GridType header note); boxes 2..NumSVDBoxes:
  ! ordinary two-sided Lobatto DVR, all channels retained at both edges.
  !---------------------------------------------------------------------
  ALLOCATE(BPDSVDBoxes(NumSVDBoxes))
  MatrixDimSVD = LSVD*NumChannels
  ALLOCATE(Gam0SVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
  ALLOCATE(OverlapSVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
  ALLOCATE(LamSVDBoxes(MatrixDimSVD,MatrixDimSVD,NumSVDBoxes))
  ALLOCATE(ODiagLSVDBoxes(NumChannels,NumSVDBoxes),ODiagRSVDBoxes(NumChannels,NumSVDBoxes))
  ALLOCATE(wsq1SVDBoxes(NumSVDBoxes),wsqLSVDBoxes(NumSVDBoxes))

  DO iBox = 1,NumSVDBoxes
     BLOCK
       DOUBLE PRECISION :: a1,a2
       INTEGER :: NOpL,NOpR,LThis
       CHARACTER(LEN=16) :: GType
       a1 = xStart + (xEnd-xStart)*(DBLE(iBox-1)/DBLE(NumBoxes))**BoxSpacingPower
       a2 = xStart + (xEnd-xStart)*(DBLE(iBox)/DBLE(NumBoxes))**BoxSpacingPower
       IF (iBox.EQ.1) THEN
          ! 'Radau' (a1 excluded as a DVR node, left OPEN/unconstrained there -- NOT
          ! Dirichlet) is the correct regularity condition for this problem's
          ! AlphaFactor=0 -- see SVDChannelBasis.f90's own GridType header note:
          ! 'LobattoDirichlet' is only right at AlphaFactor=(EffDim-1)/2=0.5.
          ! Matches the Left=2/NumOpenL=0 convention DeltaScat.f90/RMATPROPAdiabatic.f90
          ! already use for box 1's regularity at R=0.
          GType = 'Radau'
          LThis = LSVD
          NOpL = 0
       ELSE
          GType = 'Lobatto'
          LThis = LSVD
          NOpL = NumChannels
       ENDIF
       NOpR = NumChannels
       CALL BuildSVDBoxAnalytic(a1,a2,LThis,NumChannels,NOpL,NOpR,GType,BPDSVDBoxes(iBox),EIGSVD, &
            ODiagLSVDBoxes(:,iBox),ODiagRSVDBoxes(:,iBox),wsq1SVDBoxes(iBox),wsqLSVDBoxes(iBox))
       ! BuildSVDBoxAnalytic sizes EIGSVD%Gam0/Overlap/Lam to LThis*NumChannels -- box 1's
       ! LThis=LSVD-1 is smaller than MatrixDimSVD=LSVD*NumChannels, so only fill the
       ! matching leading block; PartitionAndEliminate/BoxMatch below always index by
       ! BPDSVDBoxes(iBox)%MatrixDim, never by the padded MatrixDimSVD cache array size.
       Gam0SVDBoxes(1:LThis*NumChannels,1:LThis*NumChannels,iBox) = EIGSVD%Gam0
       OverlapSVDBoxes(1:LThis*NumChannels,1:LThis*NumChannels,iBox) = EIGSVD%Overlap
       LamSVDBoxes(1:LThis*NumChannels,1:LThis*NumChannels,iBox) = EIGSVD%Lam
       CALL DeAllocateEIG(EIGSVD)
       WRITE(6,'(A,I4,A,1P,E14.6,A,E14.6,A,A,A,I4)') " SVD box ",iBox," xl=",a1," xr=",a2, &
            " grid=",TRIM(GType)," L=",LThis
     END BLOCK
  ENDDO

  ! Every box is SVD and the SVD region reaches the outer edge -- the last box's own BPD%lam
  ! carries Leff for CalcK's asymptotic Bessel matching (BuildSVDBoxAnalytic otherwise leaves
  ! it at its placeholder zero).
  BPDSVDBoxes(NumSVDBoxes)%lam(1:NumChannels) = Leff(1:NumChannels)

  EIGSVD%MatrixDim = MatrixDimSVD
  CALL AllocateEIG(EIGSVD)

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

     ALLOCATE(Boxes(NumBoxes))
     DO i = 1,NumBoxes
        Boxes(i)%NumOpenL = MERGE(0,NumChannels,i.EQ.1)
        IF (i.EQ.NumBoxes) THEN
           Boxes(i)%NumOpenR = NumOpenChannels
        ELSE
           Boxes(i)%NumOpenR = NumChannels
        ENDIF
        IF (i.EQ.1) THEN
           Boxes(i)%xl = xStart
        ELSE
           Boxes(i)%xl = Boxes(i-1)%xr
        ENDIF
        Boxes(i)%xr = xStart + (xEnd-xStart)*(DBLE(i)/DBLE(NumBoxes))**BoxSpacingPower
        CALL AllocateBox(Boxes(i))
        CALL InitZeroBox(Boxes(i))
     ENDDO

     Bnull%NumOpenL = 0
     Bnull%NumOpenR = 0
     CALL AllocateBox(Bnull)
     CALL InitZeroBox(Bnull)

     CALL AllocateScat(SD,NumChannels)

     !--- SVD boxes 1..NumSVDBoxes ---
     DO iBox = 1,NumSVDBoxes
        BLOCK
          INTEGER :: MDim
          MDim = BPDSVDBoxes(iBox)%MatrixDim
          EIGSVD%Gam0(1:MDim,1:MDim) = Gam0SVDBoxes(1:MDim,1:MDim,iBox)
          EIGSVD%Overlap(1:MDim,1:MDim) = OverlapSVDBoxes(1:MDim,1:MDim,iBox)
          IF (iBox.EQ.NumBoxes) THEN
             ! Terminal box: NumOpenR tracks NumOpenChannels(Energy), so the cached Lam
             ! (built with NumOpenR=NumChannels at setup) isn't valid here -- rebuild
             ! cheaply instead of redoing the analytic per-R solve.
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

     CALL CalcK(Boxes(NumBoxes),BPDSVDBoxes(NumBoxes),SD,reducedmass,EffDim,AlphaFactor,Energy,Threshold)

     IF (NumOpenChannels.EQ.NumChannels .AND. NumChannels.GT.1) THEN
        WRITE(KMatFile,'(A,1P,E20.12)') '# Energy = ',Energy
        DO i = 1,NumOpenChannels
           WRITE(KMatFile,'(1P,20E20.12)') SD%K(i,1:NumOpenChannels)
        ENDDO
     ENDIF

     !---------------------------------------------------------------------
     ! General multichannel unitarity check: ||S^dagger S - I||_max over ALL open channels.
     !---------------------------------------------------------------------
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
             " (phase shift skipped -- see Recombination.dat for the multichannel unitarity check)"
     ENDIF

     DO i = 1,NumBoxes
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
  WRITE(6,*) "Phase shift written to PhaseShift.dat (xmgrace: xmgrace PhaseShift.dat)"

END PROGRAM main
