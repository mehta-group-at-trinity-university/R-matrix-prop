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
| `SVDBound.x` | Single Gauss-Lobatto DVR grid, Dirichlet-Dirichlet | Bound-state solver — one direct diagonalization, no R-matrix box-chaining at all. |
| `SVDRmat.x` | Same T+V construction as `SVDBound.x` | Wigner-Eisenbud-form single-box R-matrix: one energy-independent diagonalization, then `R(E)` from the pole sum — deliberately *not* the log-derivative/`PartitionAndEliminate`/`BoxMatch` eigenchannel approach used elsewhere in this repo. |
| `BSplineFree.x` | Single free-particle (`V=0`) box, 2019-era box-1 recipe | Sanity check against the exact free-particle solution; explicitly avoids the newer `BuildBox1CombinedBasis`/`CalcGamLamBox1` machinery. |
| `PlotWavefunction.x` / `PlotWavefunctionFree.x` / `PlotWavefunctionSVD.x` | n/a (diagnostics) | Standalone, single-energy, single-box wavefunction reconstructions with hard-coded `DataDir`/energy — not general-purpose tools, used to validate against independent RK45 integration / exact analytic forms. |

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

## Validation

No automated test suite. Validation is against known benchmarks: the Baluja et al. CPC (1982)
multichannel test case (`RMATPROP2016.x`), the Amaya-Tapia/Larsen/Popiel delta-function
three-body phase shifts (`RMATPROPAdiabatic.x` + `3bodydata_delta_5ch`), and exact free-particle
solutions (`BSplineFree.x`, `PlotWavefunctionFree.x`). `verify_AT_MS_tandelta.nb` is a
Mathematica notebook doing this comparison against literature phase shifts.
