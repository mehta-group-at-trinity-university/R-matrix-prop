!****************************************************************************************************
! Bound-state SVD calculation, following 4BodySVD.f90's structure (single Gauss-Lobatto DVR
! grid over the WHOLE radial domain, Dirichlet at both ends via interior-point restriction,
! one direct Mydsyev diagonalization -- no R-matrix box chaining, no eigenchannel/energy-
! dependent re-diagonalization, no BoxMatch/CalcK reference-function matching at all), but
! generalized to call Adiabatic-Scattering-BoundStates/adiabaticSolver1D.f90's OneDimChannels
! for the adiabatic angular eigenfunctions (Uad/Psi/S_out) instead of 4BodySVD's own local
! 2D angular solver -- the one piece of that template genuinely specific to the 4-body problem.
!
! Radial grid: R=[0,20], LobattoPoints total Gauss-Lobatto nodes including both endpoints;
! both endpoints are DROPPED (Dirichlet elimination, matching AlphaFactor=(EffDim-1)/2's
! regularity requirement at R=0 -- see delta_function_validation.md -- and ordinary bound-
! state decay/normalizability at R=Rmax) -- only the LobattoPoints-2 interior nodes are
! independent DOF, exactly 4BodySVD.f90's own sRc=sR-2 convention. Eigenfunctions are
! reported ONLY at those interior DVR nodes (u(R_l)=c_l/sqrt(w_l), the normalized DVR
! basis's own Kronecker-delta property -- no interpolation needed, matching
! PlotWavefunctionSVD.f90's plotting convention, not 4BodySVD.f90's own Plotwfn, whose
! "*sqrt(weight)" (rather than "/sqrt(weight)") looks like an unvalidated diagnostic quirk).
!
! H is assembled and diagonalized exactly as 4BodySVD.f90 does: TmatR(i,j)*O(inu,jmu) off-
! diagonal (dense in the DVR spatial index, Suno Eq. 18's kinetic term), +Uad(i,nu,1) added
! on the DVR-and-channel diagonal, then a single SIMPLE (non-generalized) Mydsyev call. This
! is exact for AlphaFactor=(EffDim-1)/2 specifically: RPowExp=EffDim-1-2*AlphaFactor=0
! identically for that alpha choice (regardless of EffDim), so the flat (unweighted)
! quadrature 4BodySVD.f90 uses is already the correct measure -- no R^RPowExp weighting is
! needed in TmatR/O the way BuildSVDBox's R-matrix construction needed it for general alpha.
! Deliberately does NOT include CalcGamLam's AlphaFactor*(AlphaFactor-EffDim+2)/R^2
! correction term -- and validation confirms it does NOT need to.
!
! VALIDATION (free-particle dataset, V2Depth=0, Leff=1,3,5,7,9, R=[0,20], Dirichlet-
! Dirichlet, LobattoPoints=120 or 200 -- fully converged, identical to 8 digits at both
! resolutions): the resulting spectrum matches the EXACT eigenvalues of
! -u''/(2mu) + L^2/(2mu R^2) u = E u on [0,20] with u(0)=u(20)=0 to 10-13 significant
! figures. The exact solution is u(R)=sqrt(k R)*J_nu(k R) with nu=sqrt(L^2+1/4), NOT
! nu=L -- from chi=sqrt(x)*J_nu(x) satisfying chi''+[1-(nu^2-1/4)/x^2]chi=0, so
! representing a BARE L^2/R^2 potential (exactly what Uad tabulates, no extra term)
! requires nu^2-1/4=L^2, i.e. nu=sqrt(L^2+1/4). This is the opposite mapping from what
! hyperrjry (besselnew.f, used by CalcK for the R-matrix/scattering reference functions)
! currently does: it's called with lam=Leff=L directly, giving order=L via
! order=halfd+lam-1 (halfd=EffDim/2), which represents (L^2-1/4)/R^2 instead --
! consistent with the CalcGamLam/CalcGamLamBox1 correction term's role (compensating
! for reducing the wavefunction, u=R^alpha*Psi), but NOT consistent with a bare-Uad SVD
! box that has no such correction. This SVDBound.f90 result is the first clean,
! fully-converged confirmation of which convention (nu=L vs nu=sqrt(L^2+1/4)) is
! actually correct for a GIVEN box's own potential -- resolving (retroactively) the
! session-long free-particle R-matrix K-matrix investigation: CalcK must pass
! sqrt(Leff^2+1/4), not Leff, to hyperrjry when matching against a bare-Uad box.
!****************************************************************************************************
PROGRAM SVDBound
  USE DataStructures
  USE GlobalVars
  USE Quadrature
  USE SVDChannelBasis, ONLY: SVDCalcDX
  IMPLICIT NONE

  CHARACTER(LEN=64) DataDir
  CHARACTER*64 LobattoFile64
  INTEGER :: LobattoPoints, sR, sRc, NumStates
  INTEGER :: Order1D, Left1D, Right1D, xNumPoints1D, LegPoints1D
  DOUBLE PRECISION :: xMin1D, xMax1D, Shift1D, Rmin, Rmax, muLocal
  INTEGER :: NumEigsToPlot

  DOUBLE PRECISION, ALLOCATABLE :: nodesR(:), weightsR(:)
  DOUBLE PRECISION, ALLOCATABLE :: dXr(:,:), TmatR(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Uad(:,:,:), Psi(:,:,:), S_out(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: O(:,:), TempPsi(:)
  DOUBLE PRECISION, ALLOCATABLE :: H(:,:), EVals(:), EVecs(:,:)
  INTEGER, ALLOCATABLE :: indexOf(:,:)
  DOUBLE PRECISION, EXTERNAL :: ddot

  INTEGER :: PsiDim, HalfBandWidth1D, Ufile
  INTEGER :: i, j, iL, jL, nu, mu2, inu, jmu, probSize, ivec
  CHARACTER(LEN=64) :: wfnfile
  CHARACTER(LEN=3) :: wfnsuffix

  OPEN(unit=7,file='SVDBound.inp',status='old')
  READ(7,*)
  READ(7,*) DataDir
  READ(7,*)
  READ(7,*)
  READ(7,*) LobattoPoints, NumStates, Order1D, Left1D, Right1D
  READ(7,*)
  READ(7,*)
  READ(7,*) xNumPoints1D, xMin1D, xMax1D, LegPoints1D
  READ(7,*)
  READ(7,*)
  READ(7,*) Shift1D
  READ(7,*)
  READ(7,*)
  READ(7,*) NumEigsToPlot
  CLOSE(7)

  ! Radial grid fixed at R=[0,20] per this tool's spec -- not a tunable .inp parameter.
  Rmin = 0.d0
  Rmax = 20.d0

  Pi = dacos(-1d0)
  LegendreFile = 'Legendre.dat'
  LegPoints = LegPoints1D
  ALLOCATE(xLeg(LegPoints),wLeg(LegPoints))
  CALL GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  ! muLocal read straight from masses.dat via the SAME small helper ReadAdiabaticData uses
  ! internally would require the full Fit.data/Uad.dat table this tool doesn't need --
  ! instead just read masses.dat directly (mu is the only quantity needed from DataDir here,
  ! since Uad/Psi/S_out all come fresh from OneDimChannels, not tabulated data).
  OPEN(unit=8,file=TRIM(DataDir)//'/masses.dat',status='old')
  READ(8,*)
  READ(8,*) muLocal
  CLOSE(8)
  reducedmass = muLocal
  WRITE(6,*) "DataDir=",TRIM(DataDir)," mu=",muLocal," NumStates=",NumStates

  sR = LobattoPoints
  sRc = sR - 2
  ALLOCATE(nodesR(sR),weightsR(sR))
  LobattoFile64 = 'GaussLobatto.dat'
  CALL GetGaussLobattoFactors(LobattoFile64,sR,nodesR,weightsR)
  nodesR = 0.5d0*(Rmin+Rmax) + 0.5d0*(Rmax-Rmin)*nodesR
  weightsR = 0.5d0*(Rmax-Rmin)*weightsR
  WRITE(6,*) "Radial grid: Rmin=",Rmin," Rmax=",Rmax," sR=",sR," sRc(interior)=",sRc

  !--- DVR derivative matrix: interior basis functions (full index j+1, j=1..sRc),
  !    evaluated at ALL sR nodes (Lagrange polynomial through the whole grid) -- row
  !    index is the FULL quadrature point (1..sR), column is the interior basis index
  !    (1..sRc), matching BuildSVDBox's own dXr(LQuad,L) convention (SVDChannelBasis.f90) ---
  ALLOCATE(dXr(sR,sRc))
  DO iL = 1,sR
     DO jL = 1,sRc
        CALL SVDCalcDX(nodesR,weightsR,jL+1,iL,sR,dXr(iL,jL))
     ENDDO
  ENDDO

  !--- bulk kinetic matrix: quadrature sum over ALL sR nodes (both eliminated endpoints'
  !    weights still contribute), (1/2mu) prefactor baked in directly since this is a
  !    plain (non-generalized) energy-unit Hamiltonian, matching 4BodySVD.f90 ---
  ALLOCATE(TmatR(sRc,sRc))
  DO iL = 1,sRc
     DO jL = 1,sRc
        TmatR(iL,jL) = (0.5d0/reducedmass) * SUM(weightsR(:)*dXr(:,iL)*dXr(:,jL))
     ENDDO
  ENDDO

  !--- per-R adiabatic channel data at the sRc interior nodes only (CouplingFlag=0) ---
  PsiDim = xNumPoints1D + Order1D - 3
  IF (Left1D.EQ.2)  PsiDim = PsiDim + 1
  IF (Right1D.EQ.2) PsiDim = PsiDim + 1
  HalfBandWidth1D = Order1D
  ALLOCATE(Uad(sRc,NumStates,2),Psi(sRc,PsiDim,NumStates),S_out(HalfBandWidth1D+1,PsiDim))
  Ufile = 71
  OPEN(unit=Ufile,file=TRIM(DataDir)//'/SVDBoundUad.dat',status='replace')
  CALL OneDimChannels(NumStates,PsiDim,0,0,LegendreFile,LegPoints1D,Shift1D, &
       Order1D,Left1D,Right1D,reducedmass,xNumPoints1D,xMin1D,xMax1D, &
       nodesR(2:sR-1),sRc,Ufile,72,73,74,75,Uad,Psi,S_out,DataDir)
  CLOSE(Ufile)
  WRITE(6,*) "OneDimChannels done."

  !--- channel-overlap matrix O(inu,jmu) for ALL interior-node pairs (dense, Suno Eq. 20) ---
  probSize = sRc*NumStates
  ALLOCATE(indexOf(sRc,NumStates))
  i = 1
  DO iL = 1,sRc
     DO nu = 1,NumStates
        indexOf(iL,nu) = i
        i = i+1
     ENDDO
  ENDDO

  ALLOCATE(O(probSize,probSize),TempPsi(PsiDim))
  DO iL = 1,sRc
     DO nu = 1,NumStates
        CALL dsbmv('U',PsiDim,HalfBandWidth1D,1.0d0,S_out,HalfBandWidth1D+1, &
             Psi(iL,:,nu),1,0.0d0,TempPsi,1)
        inu = indexOf(iL,nu)
        DO jL = 1,sRc
           DO mu2 = 1,NumStates
              jmu = indexOf(jL,mu2)
              O(inu,jmu) = ddot(PsiDim,TempPsi,1,Psi(jL,:,mu2),1)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE(TempPsi)
  WRITE(6,*) "Channel-overlap matrix O done."

  !--- assemble H (4BodySVD.f90's own convention: plain, non-generalized) ---
  ALLOCATE(H(probSize,probSize))
  H = 0d0
  DO iL = 1,sRc
     DO nu = 1,NumStates
        inu = indexOf(iL,nu)
        DO jL = 1,sRc
           DO mu2 = 1,NumStates
              jmu = indexOf(jL,mu2)
              H(inu,jmu) = TmatR(iL,jL)*O(inu,jmu)
           ENDDO
        ENDDO
        H(inu,inu) = H(inu,inu) + Uad(iL,nu,1)
     ENDDO
  ENDDO
  WRITE(6,*) "Hamiltonian assembled, probSize=",probSize

  ALLOCATE(EVals(probSize),EVecs(probSize,probSize))
  WRITE(6,*) "Diagonalizing..."
  CALL Mydsyev(H,probSize,EVals,EVecs)
  WRITE(6,*) "Done. Lowest ",MIN(10,probSize)," eigenvalues:"
  DO i = 1,MIN(10,probSize)
     WRITE(6,*) i,EVals(i)
  ENDDO

  OPEN(unit=10,file='SVDBound_Evals.dat',status='replace')
  DO i = 1,probSize
     WRITE(10,*) i, EVals(i)
  ENDDO
  CLOSE(10)

  !--- eigenfunctions at the interior DVR nodes only, u(R_l)=EVecs(indexOf(l,nu),ivec)/
  !    sqrt(weightsR(l+1)) -- the normalized DVR basis's own Kronecker-delta property, no
  !    interpolation needed (see header note) ---
  DO ivec = 1,MIN(NumEigsToPlot,probSize)
     WRITE(wfnsuffix,'(I3.3)') ivec
     wfnfile = 'SVDBound_wfn'//wfnsuffix//'.dat'
     OPEN(unit=20,file=TRIM(wfnfile),status='replace')
     DO iL = 1,sRc
        WRITE(20,*) nodesR(iL+1), (EVecs(indexOf(iL,nu),ivec)/DSQRT(weightsR(iL+1)), nu=1,NumStates)
     ENDDO
     CLOSE(20)
  ENDDO
  WRITE(6,*) "Wrote SVDBound_Evals.dat and ",MIN(NumEigsToPlot,probSize)," SVDBound_wfn###.dat files"

END PROGRAM SVDBound
