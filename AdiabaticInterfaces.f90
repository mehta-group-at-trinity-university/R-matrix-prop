!****************************************************************************************************
! Abstract interfaces making the per-system pieces of the adiabaticSolver1D.f90/adiabaticSolver2D.f90
! family (currently forked across ~5 repos -- 4BodySVD, Rmatrix1D, VeffAtomIon1D, R-matrix-prop, plus
! the canonical Adiabatic-Scattering-BoundStates copies) pluggable instead of copy-pasted. Three axes:
!
!   (1) Hamiltonian/potential matrix elements -- PotentialProc. Deliberately abstracted to a plain
!       (R,phi)->(V,dV/dR) evaluation at the quadrature points CalcHamiltonianGeneric already built,
!       NOT the pairwise-distance-conversion machinery canonical's own CalcHamiltonian bakes in
!       (r12/r13/r23 for 3-body, up to 6 pairwise separations for 4-body) -- that conversion is
!       system-specific and belongs inside each concrete PotentialProc implementation, not the
!       shared engine. System-specific scalar knobs (V2Depth, PotRange, alpha, masses, ...) are NOT
!       threaded through this interface -- own them as module-level state in the concrete
!       implementation's own module (the same pattern Adiabatic-Scattering-BoundStates/lib/
!       Potential.f90 already uses for V2Depth/PotRange/alpha), set by the driver before calling
!       OneDimChannelsGeneric, rather than widening this interface per caller.
!
!   (2) Angular grid generation -- GridMakerProc. Mirrors GridMakerHHL's own signature (R,xNumPoints,
!       xMin,xMax,xPoints,CalcNewBasisFunc) minus its module-owned r0/PotRange parameter (same
!       module-state-ownership reasoning as (1)). Grid REBUILD CADENCE (every R vs. once, anchored
!       at a caller-chosen R) is a separate, explicit control on OneDimChannelsGeneric itself, NOT
!       something GridMakerProc encodes -- see that subroutine's own header for why (grid maker has
!       no notion of "box", it just returns a grid for whatever R it's given).
!
!   (3) Hyperangular boundary conditions -- RobinCoeffProc. CalcBasisFuncsBP's Left/Right BC codes
!       (0/1/2/3) are already an ordinary input, but the Robin/mixed code (3) needs a per-R kLeft/
!       kRight pair -- e.g. this repo's own delta-function benchmark needs kRight=sqrt(2*mu)*R,
!       genuinely R-dependent even though the angular GRID itself is completely fixed for that
!       problem. Grid rebuild cadence and Robin-coefficient update are therefore independent axes:
!       RobinCoeffProc is called every R regardless of GridMakerProc's cadence. For non-Robin
!       Left/Right in {0,1,2}, a concrete RobinCoeffProc can just return placeholder values (e.g.
!       1d20, matching this codebase's existing "unused" sentinel for BPD%kl/BPD%kr) -- they're
!       never read by CalcBasisFuncsBP's non-Robin branches.
!
! Radial grid (already a plain input array, R(RSteps)) is NOT generalized further here -- already
! parameterized in every existing fork.
!****************************************************************************************************
MODULE AdiabaticInterfaces
  IMPLICIT NONE

  ABSTRACT INTERFACE
     !--- (1) potential value + R-derivative at every (LegPoints x (xNumPoints-1)) quadrature point.
     !    x0(lx,kx) are the already-mapped phi-quadrature points (CalcHamiltonianGeneric's own x0
     !    array); Pot/dPot follow that file's existing (LegPoints,xNumPoints) allocation convention
     !    (only kx=1..xNumPoints-1 columns are ever read -- the last column is unused padding,
     !    preserved for drop-in compatibility with the existing xT/xV/xdV assembly loop).
     SUBROUTINE PotentialProc(R,LegPoints,xNumPoints,x0,Pot,dPot)
       INTEGER, INTENT(in) :: LegPoints,xNumPoints
       DOUBLE PRECISION, INTENT(in) :: R
       DOUBLE PRECISION, INTENT(in) :: x0(LegPoints,xNumPoints)
       DOUBLE PRECISION, INTENT(out) :: Pot(LegPoints,xNumPoints), dPot(LegPoints,xNumPoints)
     END SUBROUTINE PotentialProc

     !--- (2) angular grid at radius R -- see this module's header for why r0/mu/PotRange are
     !    deliberately NOT parameters here (own them as the concrete implementation's own
     !    module-level state instead).
     SUBROUTINE GridMakerProcI(R,xNumPoints,xMin,xMax,xPoints,CalcNewBasisFunc)
       INTEGER, INTENT(in) :: xNumPoints
       INTEGER, INTENT(out) :: CalcNewBasisFunc
       DOUBLE PRECISION, INTENT(in) :: R,xMin,xMax
       DOUBLE PRECISION, INTENT(out) :: xPoints(xNumPoints)
     END SUBROUTINE GridMakerProcI

     !--- (3) Robin/mixed boundary-condition coefficients at radius R (CalcBasisFuncsBP's kLeft/
     !    kRight). Called every R regardless of grid-rebuild cadence -- see this module's header.
     SUBROUTINE RobinCoeffProc(R,kLeft,kRight)
       DOUBLE PRECISION, INTENT(in) :: R
       DOUBLE PRECISION, INTENT(out) :: kLeft,kRight
     END SUBROUTINE RobinCoeffProc
  END INTERFACE

END MODULE AdiabaticInterfaces
