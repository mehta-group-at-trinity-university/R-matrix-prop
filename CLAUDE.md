# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Overhauled multichannel R-matrix propagation code for few-body quantum scattering (per
README.md). It propagates a log-derivative/R-matrix outward in hyperradius across a set of
"boxes" spanning `[xStart, xEnd]`, using a B-spline (or DVR/SVD) basis per box, and extracts
scattering K/S-matrices and phase shifts by matching to asymptotic Bessel reference functions
at the outer edge. See the top-of-repo `CLAUDE.md` at `~/Documents/GitHub/CLAUDE.md` for the
conventions shared across the whole `~/Documents/GitHub/` workspace (shared library migration
policy, general build pattern) — this file covers only what's specific to this repo.

There is no `.gitignore`; build artifacts (`*.o`/`*.mod`/`*.x`), run logs, and generated
`.dat` output routinely show up as untracked/modified in `git status` and are not meaningful
diffs on their own.

## Build

```bash
make            # builds RMATPROP2016.x and RMATPROPAdiabatic.x (the `all` target)
make clean      # rm *.mod *.o *.x
```

Every other executable is a standalone target not built by `all` — build them individually:

```bash
make RMATPROPHybridSVD.x
make SVDBound.x
make SVDRmat.x
make BSplineFree.x
make PlotWavefunction.x
make PlotWavefunctionFree.x
make PlotWavefunctionSVD.x
```

`makefile` and `Makefile` are the same file (case-insensitive filesystem, same inode) — don't
edit one expecting the other to diverge.

All targets compile with `${DEBUG} = -fcheck=all` unconditionally (bounds/shape checking) —
this is intentional for a code base with no test suite; don't strip it to "optimize" a build
without checking whether the user wants a perf run.

Executables read parameters from stdin:

```bash
./RMATPROPAdiabatic.x < RMATPROPAdiabatic.inp > run.log
```

## Architecture

### RMatPropCore.f90 — the shared propagation engine

All drivers link this module and reuse its pipeline unchanged; only how a box's *basis* is
built differs between drivers. Key pieces, roughly in the order a driver calls them:

- `GlobalVars` — run-wide scalars (`NumChannels`, `xStart`/`xEnd`, `reducedmass`, `AlphaFactor`,
  ...). **`AlphaFactor` is always set to `(EffDim-1)/2`** (never 0) with `StartBC=0` — this is a
  fixed convention for this codebase, not a free parameter, per the comment at the top of the
  module.
- `DataStructures` — `BPData` (one box's basis + potential matrices), `GenEigVal` (`Gam0`/
  `Overlap` built once per box and reused across an energy scan via `CombineGam`; `Gam =
  Gam0 + Energy*Overlap` is what actually gets fed to `PartitionAndEliminate`), `BoxData`
  (per-box propagated log-derivative amplitudes `Z`/`ZP`), `ScatData` (K/R/S/T matrices).
- `MakeBasis` / `CalcGamLam` (or a box-1-specific variant, see below) — build the per-box
  kinetic+potential+coupling matrix.
- `BuildBox1CombinedBasis` / `CalcGamLamBox1` / `CoulombRegularNumeric` — box 1 (innermost,
  touching the origin) is special-cased: it needs a regularity condition instead of an ordinary
  two-sided basis, and can absorb a near-origin `-C/R`-type Coulomb divergence (`CoulombC`) via
  numeric RK integration of the regular solution rather than the ordinary spline basis.
- `PartitionAndEliminate` / `PartitionAndEliminateFull` — eliminate closed-channel/interior DOF
  per box, reducing to the retained boundary amplitudes. The `Full` variant also returns the
  full eigenvector data (`Xout`/`OpenIdxOut`/`ClosedIdxOut`) needed to reconstruct the actual
  spline coefficients (used only by the `PlotWavefunction*` diagnostics, not the scattering
  drivers).
- `BoxMatch` (via `BoxMatchEigensystem`) — propagates the R-matrix/log-derivative from one box
  to the next.
- `CalcK` — final asymptotic Bessel-function matching at `xEnd` to produce the K-matrix (and
  phase shift when exactly one channel is open).

### Driver programs

| Program | Basis / data source | Notes |
|---|---|---|
| `RMATPROP2016.x` | `BalujaParameters` module — synthetic multipole-coupling test potential | Reproduces the Baluja et al. CPC (1982) benchmark; uses `RMATPROP.inp`-style input (`NumParticles`/`NumChannels`/masses/... — format documented in the `.inp` file's own trailer comment), read via `RMatPropCore`'s `ReadGlobal`, not the `DataDir`-based format below. |
| `RMATPROPAdiabatic.x` | `AdiabaticPotential.f90`, reading tabulated `Uad.dat`/`Pmat.dat`/`Qmat.dat` from a `DataDir` (see below) | Real adiabatic hyperspherical data in place of the synthetic Baluja potential; energy-scan driver (Emin/Emax/NumEnergies, linear or quadratic spacing) validated against Amaya-Tapia, Larsen & Popiel, Few-Body Systems 23, 87-109 (1997). |
| `RMATPROPHybridSVD.x` | Box 1..`NumSVDBoxes`: SVD/DVR boxes via `SVDChannelBasis.f90` (interfacing `adiabaticSolver1D.f90`'s `OneDimChannels`); remaining boxes: same B-spline basis as `RMATPROPAdiabatic.x` | `Rswitch` between the two basis types falls out of the uniform box-edge spacing formula. |
| `RMATPROPSVDAnalytic.x` | Every box SVD, built by `AnalyticSVDChannelBasis.f90`'s `BuildSVDBoxAnalytic` — per-R channel data (`Uad`, cross-R overlap) from `KMSFormulas.f`'s closed-form `SolveQ`/`EvalU`/`EvalOverlapAt` | Zero-basis-error benchmark driver for the equal-mass three-boson delta-function problem (Mehta, Esry & Greene, PRA 76, 022711 (2007)); physics hardcoded, no `Fit.data`/`FitLeff.data`. Not part of `all`. |
| `RMATPROPSVDRobin.x` | Every box SVD, built by `RobinSVDChannelBasis.f90`'s `BuildSVDBoxRobin` — per-R channel data from a genuine numerical diagonalization (`CalcBasisFuncsBP`'s Robin/`Left=1,Right=3` BC, fed the exact delta-function log-derivative `sqrt(2*mu)*R` at `phi=pi/6`) | Same benchmark as above but exercising the SVD machinery's generic (non-closed-form) numerical path. See "Delta-function SVD-box benchmark" below — two real bugs were found and fixed building this. Not part of `all`. |
| `RMATPROPHybridSVDGenericDelta.x` | Every box SVD via `GenericSVDChannelBasis.f90`'s `BuildSVDBoxGeneric`, driven by `DeltaFunctionPlugins.f90` through `lib/AdiabaticSolverGeneric.f90`'s `OneDimChannelsGeneric` | **BROKEN for this problem — 8.822 rad off the exact result, do not use.** The generic pipeline sandwiches ONE overlap matrix (`S_out`, captured at the box's last radial node) for every cross-R pair, which is only valid when the angular basis is R-independent. The `Right=3` Robin BC (`kRight=sqrt(2*mu)*R`) rebuilds the basis at every R, so `O`'s same-R block is off the identity by up to 3.3e-2. Not fixable here — see the driver's own header banner. Unrelated to the gauge bug below, which is fixed. |
| `RMATPROPHybridSVDGenericSech.x` / `RMATPROPHybridSVDGenericSech3Boson.x` | Same generic SVD pipeline via `SechFunctionPlugins.f90`. Both support a B-spline tail beyond `NumSVDBoxes`: `Sech3Boson` reads `NumSVDBoxes, NumBoxes` on one line, `Sech` takes an **optional** trailer (`NumBoxesTotal` + `OrderTail/xNumPointsTail/LegPointsTail`) so an `.inp` without it still runs pure-SVD unchanged | The angular BC is R-independent (`Right=0` and `Right=1` respectively, fixed grid), so the basis never varies with R and the shared-`S_out` trick above is EXACT — these are unaffected by that defect. |
| `RMATPROPAdiabaticDirectDelta.x` / `RMATPROPAdiabaticDirect3Boson.x` | Ordinary B-spline boxes (`RMatPropCore`), but `Uad`/`P`/`Q` computed directly at every radial quadrature point via `OneDimChannelsGeneric` — no interpolation, no `DataDir` for the delta variant | Needs no cross-R overlap matrix, so it sidesteps the defect above entirely. Delta variant reaches ~7.9e-3 rad, the 5-channel floor. Boxes MUST be built in ascending-R order (phase chaining) — see the note at each build site. |
| `SVDBound.x` | Single Gauss-Lobatto DVR grid, Dirichlet-Dirichlet | Bound-state solver — one direct diagonalization, no R-matrix box-chaining at all. |
| `SVDRmat.x` | Same T+V construction as `SVDBound.x` | Wigner-Eisenbud-form single-box R-matrix: one energy-independent diagonalization, then `R(E)` from the pole sum — deliberately *not* the log-derivative/`PartitionAndEliminate`/`BoxMatch` eigenchannel approach used elsewhere in this repo. |
| `BSplineFree.x` | Single free-particle (`V=0`) box, 2019-era box-1 recipe | Sanity check against the exact free-particle solution; explicitly avoids the newer `BuildBox1CombinedBasis`/`CalcGamLamBox1` machinery. |
| `PlotWavefunction.x` / `PlotWavefunctionFree.x` / `PlotWavefunctionSVD.x` | n/a (diagnostics) | Standalone, single-energy, single-box wavefunction reconstructions with hard-coded `DataDir`/energy — not general-purpose tools, used to validate against independent RK45 integration / exact analytic forms. |

**Naming convention — `Hybrid` means the driver can switch box type mid-chain.** A driver is named
`RMATPROPHybridSVD*` only if it can propagate SVD boxes *and* ordinary adiabatic (B-spline) boxes in
one run, i.e. it accepts a total `NumBoxes` exceeding `NumSVDBoxes` and builds the remainder from
tabulated `Uad`/`Pmat`/`Qmat`. Drivers that only ever build SVD boxes are named `RMATPROPSVD*`
(`RMATPROPSVDAnalytic.x`, `RMATPROPSVDRobin.x`) — these two have no tail machinery at all: no
`ReadAdiabaticData`, no `BPData`, no `NumSVDBoxes+1..NumBoxes` loop anywhere, and `NumBoxes` is
*assigned* from `NumSVDBoxes` rather than read, so no input can reach a tail. Their `NumBoxes` exists
only as chain bookkeeping (box-edge spacing denominator, `Boxes(:)` allocation, `CalcK`'s outermost
box), not as a switch.

### Data directories (`3bodydata_*_5ch/`)

Each is a self-contained dataset read by `AdiabaticPotential.f90`'s `ReadAdiabaticData`, produced
externally by `Adiabatic-Scattering-BoundStates`'s adiabatic solvers:

- `Fit.data` — `NumChannels, NumDataPoints, alpha, ddim`
- `masses.dat` — comment line, then `mu` (hyperradial reduced mass)
- `Uad.dat` — `# NumChannels NumDataPoints` header, then `R, U(1:NumChannels)` rows
- `Pmat.dat`, `Qmat.dat` — coupling matrices
- `FitLeff.data` — per-channel threshold + effective angular momentum for `CalcK`'s asymptotic
  Bessel matching

## Known gotchas

- **The generic SVD pipeline's shared overlap matrix is only valid for an R-INDEPENDENT angular
  basis.** `GenericSVDChannelBasis.f90` builds a box's cross-R channel-overlap `O` by sandwiching a
  single `S_out`, and `lib/AdiabaticSolverGeneric.f90` assigns `S_out = S` only *after* its R loop —
  so it is the overlap at the box's **last** radial node. With `Left/Right` in `{0,1,2}` the basis
  never changes with R and this is exact. With a Robin BC (`Left/Right = 3`) whose coefficient
  depends on R — the delta-function problem's `kRight = sqrt(2*mu)*R` — `CalcBasisFuncsBP` is
  rebuilt at every R, and the trick silently breaks: `O`'s same-R block, which must be the identity
  at every node, is off by up to 3.3e-2 on the diagonal at the far end of a box, decaying linearly
  toward the node where `S` was captured. Returning `S` per R would not fix it either — the correct
  cross-R overlap needs a cross-basis Gram `int u_i^(Ra) u_j^(Rb)`, obtainable only by pointwise
  real-space evaluation (what `RobinSVDChannelBasis.f90` does). This is what makes
  `RMATPROPHybridSVDGenericDelta.x` unusable for the delta benchmark.
- **`OneDimChannelsGeneric`'s `FixPhase` only chains phases WITHIN one call.** Callers that drive an
  R-matrix box chain invoke it once per box, so without help each box inherits whatever arbitrary
  sign ARPACK returned at its own first R — an eigenvector-gauge jump at every box boundary, exactly
  where `BoxMatch` matches amplitudes assuming a shared channel basis. A gauge held constant over the
  whole domain is unobservable; one that jumps per box is not. Fixed by threading `PhaseSeed` between
  calls, which **requires callers to build boxes in ascending-R order** — the ordinary-box drivers
  previously built `2..N-1, then 1, then N` (a caching grouping, harmless while
  `SetAdiabaticPotential` was stateless) and had to be reordered. Symptom to recognize: per-box data
  correct in isolation, but the answer degrading non-monotonically under refinement.
- **`Threshold(i) < 0` is treated as "genuine Coulomb-diverging bound-pair channel"** in
  `RMATPROPAdiabatic.f90`'s `CoulombC` setup (and mirrored in `RMATPROPHybridSVD.f90`), using an
  equal-mass delta-function-specific `gamma0` formula (Mehta, Esry & Greene, PRA 76, 022711
  (2007)). This is **not general** — a different potential, unequal masses, or numerical noise
  producing a spuriously negative threshold (see `3bodydata_free_5ch`, where `Threshold` is
  ~1e-13 negative from ARPACK noise, not physical) will silently get the wrong near-origin
  physics unless this block is adapted or zeroed out for that dataset.
- `PlotWavefunctionFree.f90` targets `xMaxLocal=9.99`, not `10.0`, because `3bodydata_free_5ch`'s
  fit range tops out exactly at `R=10.0` and a quadrature point landing on `10.000000000000002`
  crashes the interpolator.
- `BSplineFree.x`'s box starts at exactly `xStart=0.0` (safe only because Gauss-Legendre
  quadrature points never land exactly on the box edge) — don't copy this pattern into a driver
  whose quadrature rule could include an endpoint.
- `adiabaticSolver1D.f90` here is a **local copy** of `Adiabatic-Scattering-BoundStates`'s file
  with one dead subroutine removed to resolve a link-time name collision with
  `RMatPropCore.f90`'s own `GridMaker` — not a physics fork; keep it in sync with upstream
  otherwise.
- The `Potential.o`/`adiabaticSolver1D.o` build rules intentionally link against *this* repo's
  own `Bsplines.o`/`matrix_stuff.o`/`Quadrature.o`, not `Adiabatic-Scattering-BoundStates`'s
  local copies, to avoid linking two competing copies of the same routines into one executable —
  preserve this if editing the makefile.
- **`DSYGV`'s per-R eigenvector sign is arbitrary** and must be pinned to a channel-appropriate,
  never-vanishing reference point before using it in any cross-R overlap (`RobinSVDChannelBasis.f90`).
  Left uncorrected, adjacent DVR nodes can pick up an independent sign flip, corrupting the
  cross-R channel overlap with a magnitude-correct-but-wrong-sign entry. Pin via `Phi(0)` for
  `cos(kappa*phi)`-branch channels (never zero there) but via `Phi(pi/6)` for the `cosh(q*phi)`
  channel specifically — `Phi(0)` for that branch is exponentially suppressed relative to
  `Phi(pi/6)` at large `kr` and decays into roundoff noise (confirmed: `evecPhi(1,1)~3e-5` vs
  `~2` for other channels by `R~65`), reintroducing the same sign-ambiguity failure if used as
  the reference there. See `RobinSVDChannelBasis.f90`'s own comment at the sign-fix for the full
  derivation of why a single reference point doesn't work for every channel.
- **A fixed angular B-spline resolution (`OrderPhi`/`xNumPointsPhi`) does not stay adequate as R
  grows** in `RobinSVDChannelBasis.f90` — the Robin condition's `kr=sqrt(2*mu)*R` makes the
  channel-1 (`cosh`) eigenfunction increasingly steep near `phi=pi/6` as R increases, and a
  resolution adequate at small R can leave `RMATPROPSVDRobin.x`'s phase shift off by up to
  ~0.2 rad at low q even though the per-box `Gam0`/`Overlap`/`Lam`/cross-R-overlap all look
  correct in isolation (low-q/near-threshold states have long healing lengths and are
  disproportionately sensitive to this large-R truncation). `xNumPointsPhi=12` was not enough
  for the `xEnd=100`, `NumBoxes=20` benchmark config; `xNumPointsPhi=20` brought the discrepancy
  down to the same ~0.008 rad floor set by ordinary 5-channel truncation (see validation below).

## Validation

No automated test suite. Validation is against known benchmarks: the Baluja et al. CPC (1982)
multichannel test case (`RMATPROP2016.x`), the Amaya-Tapia/Larsen/Popiel delta-function
three-body phase shifts (`RMATPROPAdiabatic.x` + `3bodydata_delta_5ch`), and exact free-particle
solutions (`BSplineFree.x`, `PlotWavefunctionFree.x`). `verify_AT_MS_tandelta.nb` is a
Mathematica notebook doing this comparison against literature phase shifts.

### Delta-function SVD-box benchmark (`RMATPROPSVDAnalytic.x` / `RMATPROPSVDRobin.x`)

Both drivers solve the equal-mass three-boson delta-function problem end-to-end using an
all-SVD box chain (no B-spline tail, no `Fit.data`/`FitLeff.data`), compared against
`Adiabatic-Scattering-BoundStates/DeltaScat`'s independent no-interpolation benchmark and the
exact closed-form Mehta-Shepard atom-dimer K-matrix. Compare branch-continued phase shift `delta`
(`atan(K)`, unwrapped), not raw `K` — `K=tan(delta)` has a mathematical pole wherever `delta`
crosses an odd multiple of `pi/2`, which can dominate a raw-`K` diff with an artifact even where
`delta` itself agrees closely.

- **Analytic** (`NumBoxes=20, xStart=0, xEnd=100, LSVD=24`, `Emin/Emax/NumEnergies =
  -0.9999999635/-0.02/150` quadratic): matches `DeltaScat` to `max|delta-delta_DeltaScat|=9e-4`
  rad, both matching exact to `~0.008` rad (ordinary 5-channel truncation, not a numerical
  artifact — see `Adiabatic-Scattering-BoundStates/CLAUDE.md`'s own note on this).
- **Robin** (same box grid, `OrderPhi=6, xNumPointsPhi=20`): matches exact to the same `~0.008`
  rad floor once `xNumPointsPhi` is large enough (see the angular-resolution gotcha above) —
  confirming the SVD box-chaining/`BoxMatch`/`CalcK` machinery is correct via a fully
  independent, non-closed-form code path (genuine per-R diagonalization instead of
  `KMSFormulas.f`'s analytic shortcut).
