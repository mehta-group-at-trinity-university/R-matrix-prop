!****************************************************************************************************
! No-interpolation analog of AdiabaticPotential.f90's SetAdiabaticPotential: fills BPD%Pot/BPD%Pcoup
! by running the real per-R adiabatic diagonalization (OneDimChannelsGeneric, ARPACK shift-invert,
! via SechFunctionPlugins.f90's sech^2-3-identical-boson plugins) DIRECTLY at each of the box's own
! Gauss-Legendre quadrature points -- exactly mirroring how DeltaScat.f90 evaluates its (closed-form)
! delta-function potential directly at quadrature points, with no interpolant/tabulated-Rdata-array
! domain-floor question ever arising, since there is no tabulated array at all.
!
! Physics/domain hardcoded to the equal-mass 3-identical-boson sech^2 problem (domain [0,pi/6],
! Left=Right=1, per this repo's 3boson dataset generation -- DriveAdiabaticSolver1D_3boson.par).
! Cost: one ARPACK diagonalization per quadrature point (~0.05s each at xNumPointsPhi=200, confirmed
! this session via SechCurvesSVD.x) -- call this ONCE per box and cache BPD%Pot/BPD%Pcoup (mirroring
! how RMATPROPAdiabatic.f90 already caches Gam0/Overlap/Lam for the energy-independent boxes), never
! per-energy -- re-evaluating per energy would multiply the cost by NumEnergies for no reason, since
! the physics here doesn't depend on Energy at all.
!****************************************************************************************************
MODULE SechAdiabaticPotentialDirect
  IMPLICIT NONE
  INTEGER, PARAMETER :: DirectOrder1D=5, DirectLeft1D=1, DirectRight1D=1
  INTEGER, PARAMETER :: DirectxNumPoints1D=200, DirectLegPoints1D=10
  DOUBLE PRECISION, PARAMETER :: DirectShiftInit=-5.0d0
CONTAINS
  SUBROUTINE SetAdiabaticPotentialDirect(BPD,muLocal)
    USE DataStructures
    USE Quadrature, ONLY: LegPoints, xLeg
    USE SechFunctionPlugins
    USE AdiabaticSolverGeneric, ONLY: OneDimChannelsGeneric
    IMPLICIT NONE
    TYPE(BPData) BPD
    DOUBLE PRECISION, INTENT(in) :: muLocal

    INTEGER :: kx,lx,idx,NumQuad,i,j,PsiDim,NumStates
    DOUBLE PRECISION :: ax,bx,xIntScale,xScaledZero,Shift,dum
    DOUBLE PRECISION, ALLOCATABLE :: Rquad(:)
    DOUBLE PRECISION, ALLOCATABLE :: UadOut(:,:,:), PsiOut(:,:,:), S_out(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: PMatFull(:,:,:), QMatFull(:,:,:)
    CHARACTER*64 :: LegendreFile64, datadirLoc

    NumStates = BPD%NumChannels

    !--- box's own quadrature R array, in the EXACT same (kx,lx) nested order
    !    SetAdiabaticPotential's own BPD%x(lx,kx) fill uses, so BPD%Pot/BPD%Pcoup end up
    !    indexed identically regardless of which routine filled them. ---
    NumQuad = (BPD%xNumPoints-1)*LegPoints
    ALLOCATE(Rquad(NumQuad))
    idx = 0
    DO kx = 1,BPD%xNumPoints-1
       ax = BPD%xPoints(kx)
       bx = BPD%xPoints(kx+1)
       xIntScale = 0.5d0*(bx-ax)
       xScaledZero = 0.5d0*(bx+ax)
       DO lx = 1,LegPoints
          idx = idx+1
          Rquad(idx) = xIntScale*xLeg(lx) + xScaledZero
          BPD%x(lx,kx) = Rquad(idx)
       ENDDO
    ENDDO

    LegendreFile64 = 'Legendre.dat'
    datadirLoc = '.'
    PsiDim = DirectxNumPoints1D + DirectOrder1D - 3
    ALLOCATE(UadOut(NumQuad,NumStates,2),PsiOut(NumQuad,PsiDim,NumStates),S_out(DirectOrder1D+1,PsiDim))

    !--- run the real diagonalization at every one of this box's quadrature points, in one
    !    call (CouplingFlag=1 so P/Q come out via Hellmann-Feynman, using ADJACENT Rquad
    !    points for the dP/dR finite difference -- correct here since Rquad is genuinely
    !    monotonic: increasing within each sub-interval, and sub-intervals are
    !    non-overlapping and increasing, exactly matching how SetAdiabaticPotential's own
    !    (kx,lx) loop builds BPD%x). Scratch P/Q/dP/VQ files: written then immediately read
    !    back using AdiabaticPotential.f90's own Pmat.dat/Qmat.dat reading convention, so the
    !    sign/normalization convention matches exactly (OneDimChannelsGeneric's own Q output
    !    is P^2's negative, i.e. Qtilde, the SAME quantity AdiabaticPotential.f90 already
    !    assumes when it does +Q/(2mu) -- confirmed via that file's own header comment). ---
    Shift = DirectShiftInit
    OPEN(unit=81,file='SechDirectUad.scratch',status='replace')
    OPEN(unit=82,file='SechDirectPmat.scratch',status='replace')
    OPEN(unit=83,file='SechDirectQmat.scratch',status='replace')
    OPEN(unit=84,file='SechDirectdPmat.scratch',status='replace')
    OPEN(unit=85,file='SechDirectVQmat.scratch',status='replace')
    CALL OneDimChannelsGeneric(NumStates,PsiDim,0,1,LegendreFile64,DirectLegPoints1D,Shift, &
         DirectOrder1D,DirectLeft1D,DirectRight1D,muLocal,DirectxNumPoints1D,0d0,dacos(-1d0)/6.0d0, &
         Rquad,NumQuad,.TRUE.,0d0, &
         SechPotential,SechGridMaker,SechRobinCoeff, &
         81,82,83,84,85,UadOut,PsiOut,S_out,datadirLoc)
    CLOSE(81)
    CLOSE(82)
    CLOSE(83)
    CLOSE(84)
    CLOSE(85)

    !--- read Pmat.dat/Qmat.dat scratch files back exactly as AdiabaticPotential.f90's
    !    ReadAdiabaticData does (same R-line-then-NumStates-rows-of-NumStates-values format). ---
    ALLOCATE(PMatFull(NumQuad,NumStates,NumStates),QMatFull(NumQuad,NumStates,NumStates))
    OPEN(unit=82,file='SechDirectPmat.scratch',status='old')
    DO i = 1,NumQuad
       READ(82,*) dum
       DO j = 1,NumStates
          READ(82,*) PMatFull(i,j,1:NumStates)
       ENDDO
    ENDDO
    CLOSE(82)
    OPEN(unit=83,file='SechDirectQmat.scratch',status='old')
    DO i = 1,NumQuad
       READ(83,*) dum
       DO j = 1,NumStates
          READ(83,*) QMatFull(i,j,1:NumStates)
       ENDDO
    ENDDO
    CLOSE(83)

    !--- fill BPD%Pot/BPD%Pcoup in the SAME convention SetAdiabaticPotential uses:
    !    Pot(mch,nch) = Qtilde(mch,nch)/(2*mu) + delta_{mch,nch}*U(mch); Pcoup = P. No
    !    RcutoffAnalytic/CoulombC special-casing here at all -- there is no near-origin
    !    interpolant-domain-floor concern to guard against in the first place. ---
    idx = 0
    DO kx = 1,BPD%xNumPoints-1
       DO lx = 1,LegPoints
          idx = idx+1
          BPD%Pcoup(:,:,lx,kx) = PMatFull(idx,:,:)
          BPD%Pot(:,:,lx,kx) = QMatFull(idx,:,:)/(2d0*muLocal)
          DO i = 1,NumStates
             BPD%Pot(i,i,lx,kx) = BPD%Pot(i,i,lx,kx) + UadOut(idx,i,1)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(Rquad,UadOut,PsiOut,S_out,PMatFull,QMatFull)
  END SUBROUTINE SetAdiabaticPotentialDirect
END MODULE SechAdiabaticPotentialDirect
