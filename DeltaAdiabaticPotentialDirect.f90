!****************************************************************************************************
! No-interpolation, ordinary-adiabatic-box analog of AdiabaticPotential.f90's SetAdiabaticPotential
! for the equal-mass 3-boson delta-function problem, mirroring SechAdiabaticPotentialDirect.f90
! exactly but wired to DeltaFunctionPlugins.f90 (domain [0,pi/6], Left=1 Neumann at phi=0, Right=3
! Robin at phi=pi/6 with kRight=sqrt(2*mu)*R -- the contact interaction).
!
! Purpose: isolate whether OneDimChannelsGeneric's ARPACK-shift-invert diagonalization + FixPhase
! is itself trustworthy for this Robin-BC problem, independent of GenericSVDChannelBasis.f90's own
! cross-R channel-overlap construction (which reuses a single box-rightmost-R overlap matrix S_out
! for every DVR-node pair in the box -- invalid here since the Robin-BC-satisfying combined basis
! functions are genuinely R-dependent, see this session's own diagnosis). This module needs NO
! cross-R overlap matrix at all -- CouplingFlag=1 gets P/Q directly via Hellmann-Feynman using
! ADJACENT quadrature R points within one OneDimChannelsGeneric call, exactly as ordinary adiabatic
! box propagation (MakeBasis/CalcGamLam) already expects -- so if this driver reproduces the exact
! result while RMATPROPHybridSVDGenericDelta.x does not, that pins the bug specifically on
! GenericSVDChannelBasis.f90's shared-S_out trick, not on ARPACK/FixPhase/OneDimChannelsGeneric.
!****************************************************************************************************
MODULE DeltaAdiabaticPotentialDirect
  IMPLICIT NONE
  INTEGER, PARAMETER :: DirectOrder1D=6, DirectLeft1D=1, DirectRight1D=3
  INTEGER, PARAMETER :: DirectxNumPoints1D=20, DirectLegPoints1D=10
  DOUBLE PRECISION, PARAMETER :: DirectShiftInit=-1.0d0
  ! Number of angular states carried in the INTERMEDIATE sum that builds Q.
  ! OneDimChannelsGeneric forms Q as -P.P (dgemm), which is the closure relation
  ! Qtilde_mn = <dPhi_m/dR|dPhi_n/dR> = sum_k P_km P_kn -- EXACT only when k runs over a
  ! complete set.  Truncating that sum at the NumChannels physical channels is a real error:
  ! feeding RMATPROPAdiabatic.x a dataset copy with Qmat.dat replaced by -P.P over 5 states
  ! degrades it from 0.0076 rad to 0.046 rad against the exact Mehta-Shepard result.  So solve
  ! for DirectNumSolve >= NumChannels states, build P/Q over all of them, and only THEN truncate
  ! to the leading NumChannels block -- keeping the 5-channel PHYSICAL truncation (which sets
  ! the ~0.008 rad benchmark floor) while converging the closure sum.  Capped against PsiDim
  ! below since ARPACK needs ncv=MIN(MatrixDim,10*NumStates) comfortably above NumStates.
  INTEGER, PARAMETER :: DirectNumSolve=20
CONTAINS
  SUBROUTINE SetAdiabaticPotentialDirect(BPD,muLocal)
    USE DataStructures
    USE Quadrature, ONLY: LegPoints, xLeg
    USE DeltaFunctionPlugins
    USE AdiabaticSolverGeneric, ONLY: OneDimChannelsGeneric
    IMPLICIT NONE
    TYPE(BPData) BPD
    DOUBLE PRECISION, INTENT(in) :: muLocal

    INTEGER :: kx,lx,idx,NumQuad,i,j,PsiDim,NumStates,NumSolve
    DOUBLE PRECISION :: ax,bx,xIntScale,xScaledZero,Shift,dum
    DOUBLE PRECISION, ALLOCATABLE :: Rquad(:)
    DOUBLE PRECISION, ALLOCATABLE :: UadOut(:,:,:), PsiOut(:,:,:), S_out(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: PMatFull(:,:,:), QMatFull(:,:,:)
    CHARACTER*64 :: LegendreFile64, datadirLoc

    NumStates = BPD%NumChannels
    NumSolve = MIN(MAX(DirectNumSolve,NumStates), DirectxNumPoints1D+DirectOrder1D-6)

    !--- box's own quadrature R array, same (kx,lx) nested order SetAdiabaticPotential uses ---
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
    ALLOCATE(UadOut(NumQuad,NumSolve,2),PsiOut(NumQuad,PsiDim,NumSolve),S_out(DirectOrder1D+1,PsiDim))

    !--- real diagonalization at every quadrature point in this box, one call, CouplingFlag=1
    !    so P/Q come out via Hellmann-Feynman using adjacent Rquad points -- no cross-R overlap
    !    matrix construction anywhere in this file. ---
    Shift = DirectShiftInit
    OPEN(unit=81,file='DeltaDirectUad.scratch',status='replace')
    OPEN(unit=82,file='DeltaDirectPmat.scratch',status='replace')
    OPEN(unit=83,file='DeltaDirectQmat.scratch',status='replace')
    OPEN(unit=84,file='DeltaDirectdPmat.scratch',status='replace')
    OPEN(unit=85,file='DeltaDirectVQmat.scratch',status='replace')
    CALL OneDimChannelsGeneric(NumSolve,PsiDim,0,1,LegendreFile64,DirectLegPoints1D,Shift, &
         DirectOrder1D,DirectLeft1D,DirectRight1D,muLocal,DirectxNumPoints1D,0d0,dacos(-1d0)/6.0d0, &
         Rquad,NumQuad,.TRUE.,0d0, &
         DeltaZeroPotential,DeltaWedgeGridMaker,DeltaRobinCoeff, &
         81,82,83,84,85,UadOut,PsiOut,S_out,datadirLoc)
    CLOSE(81)
    CLOSE(82)
    CLOSE(83)
    CLOSE(84)
    CLOSE(85)

    !--- read Pmat.dat/Qmat.dat scratch files back exactly as AdiabaticPotential.f90's
    !    ReadAdiabaticData does. ---
    ALLOCATE(PMatFull(NumQuad,NumSolve,NumSolve),QMatFull(NumQuad,NumSolve,NumSolve))
    OPEN(unit=82,file='DeltaDirectPmat.scratch',status='old')
    DO i = 1,NumQuad
       READ(82,*) dum
       DO j = 1,NumSolve
          READ(82,*) PMatFull(i,j,1:NumSolve)
       ENDDO
    ENDDO
    CLOSE(82)
    OPEN(unit=83,file='DeltaDirectQmat.scratch',status='old')
    DO i = 1,NumQuad
       READ(83,*) dum
       DO j = 1,NumSolve
          READ(83,*) QMatFull(i,j,1:NumSolve)
       ENDDO
    ENDDO
    CLOSE(83)

    !--- fill BPD%Pot/BPD%Pcoup in SetAdiabaticPotential's own convention:
    !    Pot(mch,nch) = Qtilde(mch,nch)/(2*mu) + delta_{mch,nch}*U(mch); Pcoup = P. ---
    idx = 0
    DO kx = 1,BPD%xNumPoints-1
       DO lx = 1,LegPoints
          idx = idx+1
          ! Truncate to the leading NumChannels block only HERE -- P/Q above were built with
          ! NumSolve intermediate states so the -P.P closure sum is converged (see DirectNumSolve).
          BPD%Pcoup(:,:,lx,kx) = PMatFull(idx,1:NumStates,1:NumStates)
          BPD%Pot(:,:,lx,kx) = QMatFull(idx,1:NumStates,1:NumStates)/(2d0*muLocal)
          DO i = 1,NumStates
             BPD%Pot(i,i,lx,kx) = BPD%Pot(i,i,lx,kx) + UadOut(idx,i,1)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(Rquad,UadOut,PsiOut,S_out,PMatFull,QMatFull)
  END SUBROUTINE SetAdiabaticPotentialDirect
END MODULE DeltaAdiabaticPotentialDirect
