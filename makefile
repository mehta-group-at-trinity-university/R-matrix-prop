CMP     = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
ARPACK =  -L/opt/homebrew/lib/ -larpack
INCLUDE =  -I/opt/homebrew/include
LAPACK =  -framework Accelerate
INTERP_DIR = $(HOME)/Documents/GitHub/interpolation
LIB_DIR = $(HOME)/Documents/GitHub/lib
ASBS_DIR = $(HOME)/Documents/GitHub/Adiabatic-Scattering-BoundStates
OBJS  = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o RMATPROP2016.o Quadrature.o Interpolation.o
ADIABATIC_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o RMATPROPAdiabatic.o
ADIABATIC_SEEDED_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o RMATPROPAdiabatic_seeded.o
SVD_OBJS = besselnew.o Bsplines.o matrix_stuff.o Quadrature.o Interpolation.o Potential.o adiabaticSolver1D.o RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o

all: RMATPROP2016.x RMATPROPAdiabatic.x

RMATPROP2016.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK}  ${CMPFLAGS} ${FORCEDP} -o RMATPROP2016.x

RMATPROPAdiabatic.x: ${ADIABATIC_OBJS}
	${CMP} ${DEBUG} ${ADIABATIC_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPAdiabatic.x

# Diagnostic: box 1 seeded from DeltaScat's exact R-matrix (Box1Rmatrix.dat) instead of the
# interpolated pipeline -- see RMATPROPAdiabatic_seeded.f90's own header. Not part of `all`;
# build/run separately with `make RMATPROPAdiabatic_seeded.x && ./RMATPROPAdiabatic_seeded.x`.
RMATPROPAdiabatic_seeded.x: ${ADIABATIC_SEEDED_OBJS}
	${CMP} ${DEBUG} ${ADIABATIC_SEEDED_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPAdiabatic_seeded.x

RMATPROPAdiabatic_seeded.o: RMATPROPAdiabatic_seeded.f90 RMatPropCore.o AdiabaticPotential.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPAdiabatic_seeded.f90

RMatPropCore.o: RMatPropCore.f90 Quadrature.o Interpolation.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMatPropCore.f90

RMATPROP2016.o: RMATPROP2016.f90 RMatPropCore.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROP2016.f90

Interpolation.o: $(INTERP_DIR)/Interpolation.f90
	${CMP} ${FORCEDP} -c $(INTERP_DIR)/Interpolation.f90

AdiabaticPotential.o: AdiabaticPotential.f90 Interpolation.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c AdiabaticPotential.f90

RMATPROPAdiabatic.o: RMATPROPAdiabatic.f90 RMatPropCore.o AdiabaticPotential.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPAdiabatic.f90

# Standalone diagnostic (see PlotWavefunction.f90 header) -- not part of `all`,
# build/run separately with `make PlotWavefunction.x && ./PlotWavefunction.x`.
PLOTWAVEFUNCTION_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o PlotWavefunction.o

PlotWavefunction.x: ${PLOTWAVEFUNCTION_OBJS}
	${CMP} ${DEBUG} ${PLOTWAVEFUNCTION_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o PlotWavefunction.x

PlotWavefunction.o: PlotWavefunction.f90 RMatPropCore.o AdiabaticPotential.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c PlotWavefunction.f90

# Standalone diagnostic (see PlotWavefunctionCh1Robin.f90 header) -- Left=3/Robin box 1 with
# kLeft fixed to the true near-origin log-derivative, single-channel 3boson dataset.
PLOTWAVEFUNCTIONCH1ROBIN_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o PlotWavefunctionCh1Robin.o

PlotWavefunctionCh1Robin.x: ${PLOTWAVEFUNCTIONCH1ROBIN_OBJS}
	${CMP} ${DEBUG} ${PLOTWAVEFUNCTIONCH1ROBIN_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o PlotWavefunctionCh1Robin.x

PlotWavefunctionCh1Robin.o: PlotWavefunctionCh1Robin.f90 RMatPropCore.o AdiabaticPotential.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c PlotWavefunctionCh1Robin.f90

TESTFREEPARTICLE_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o TestFreeParticle.o

TestFreeParticle.x: ${TESTFREEPARTICLE_OBJS}
	${CMP} ${DEBUG} ${TESTFREEPARTICLE_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o TestFreeParticle.x

TestFreeParticle.o: TestFreeParticle.f90 RMatPropCore.o AdiabaticPotential.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c TestFreeParticle.f90

TESTALPHAZERO_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o TestAlphaZero.o

TestAlphaZero.x: ${TESTALPHAZERO_OBJS}
	${CMP} ${DEBUG} ${TESTALPHAZERO_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o TestAlphaZero.x

TestAlphaZero.o: TestAlphaZero.f90 RMatPropCore.o AdiabaticPotential.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c TestAlphaZero.f90

# Free-particle variant (see PlotWavefunctionFree.f90 header) -- same pattern as
# PlotWavefunction.x, targets 3bodydata_free_5ch instead of 3bodydata_delta_5ch.
PLOTWAVEFUNCTIONFREE_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o PlotWavefunctionFree.o

PlotWavefunctionFree.x: ${PLOTWAVEFUNCTIONFREE_OBJS}
	${CMP} ${DEBUG} ${PLOTWAVEFUNCTIONFREE_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o PlotWavefunctionFree.x

PlotWavefunctionFree.o: PlotWavefunctionFree.f90 RMatPropCore.o AdiabaticPotential.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c PlotWavefunctionFree.f90

# SVD analog (see PlotWavefunctionSVD.f90 header) -- single-box SVD wavefunction
# reconstruction, overlaid against PlotWavefunctionFree.x's B-spline result.
PLOTWAVEFUNCTIONSVD_OBJS = besselnew.o Bsplines.o matrix_stuff.o Quadrature.o Interpolation.o Potential.o adiabaticSolver1D.o RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o PlotWavefunctionSVD.o

PlotWavefunctionSVD.x: ${PLOTWAVEFUNCTIONSVD_OBJS}
	${CMP} ${DEBUG} ${PLOTWAVEFUNCTIONSVD_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o PlotWavefunctionSVD.x

PlotWavefunctionSVD.o: PlotWavefunctionSVD.f90 RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c PlotWavefunctionSVD.f90

# Bound-state SVD calculation (see SVDBound.f90 header) -- single Gauss-Lobatto DVR grid,
# Dirichlet-Dirichlet, one direct Mydsyev diagonalization; no R-matrix machinery involved.
SVDBOUND_OBJS = besselnew.o Bsplines.o matrix_stuff.o Quadrature.o Interpolation.o Potential.o adiabaticSolver1D.o RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o SVDBound.o

SVDBound.x: ${SVDBOUND_OBJS}
	${CMP} ${DEBUG} ${SVDBOUND_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o SVDBound.x

SVDBound.o: SVDBound.f90 RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c SVDBound.f90

# Wigner-Eisenbud-form single-box R-matrix (see SVDRmat.f90 header) -- one energy-
# independent diagonalization of SVDBound.f90's own T+V, R(E) from the pole sum.
SVDRMAT_OBJS = besselnew.o Bsplines.o matrix_stuff.o Quadrature.o Interpolation.o Potential.o adiabaticSolver1D.o RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o SVDRmat.o

SVDRmat.x: ${SVDRMAT_OBJS}
	${CMP} ${DEBUG} ${SVDRMAT_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o SVDRmat.x

SVDRmat.o: SVDRmat.f90 RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c SVDRmat.f90

# B-spline single-box free-particle test, 2019/DeltaScat-style simple boundary condition
# (see BSplineFree.f90 header) -- no BuildBox1CombinedBasis, no tabulated interpolation.
BSPLINEFREE_OBJS = besselnew.o Bsplines.o matrix_stuff.o Quadrature.o Interpolation.o Potential.o adiabaticSolver1D.o RMatPropCore.o BSplineFree.o

BSplineFree.x: ${BSPLINEFREE_OBJS}
	${CMP} ${DEBUG} ${BSPLINEFREE_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o BSplineFree.x

BSplineFree.o: BSplineFree.f90 RMatPropCore.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c BSplineFree.f90

matrix_stuff.o: $(LIB_DIR)/matrix_stuff.f90
	${CMP} ${FORCEDP} ${CMPFLAGS} -c $(LIB_DIR)/matrix_stuff.f90

Bsplines.o: $(LIB_DIR)/Bsplines.f90
	${CMP} ${FORCEDP} -c $(LIB_DIR)/Bsplines.f90

nrtype.mod: modules_qd.o
	${CMP} ${FORCEDP} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${FORCEDP} -c modules_qd.f90

besselnew.o:	besselnew.f
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c besselnew.f

Quadrature.mod: Quadrature.o
		${CMP} ${FORCEDP} Quadrature.o

Quadrature.o: $(LIB_DIR)/Quadrature.f90
	${CMP} ${FORCEDP} -c $(LIB_DIR)/Quadrature.f90

# --- SVD R-matrix propagator: interfaces with Adiabatic-Scattering-BoundStates'
# OneDimChannels (adiabaticSolver1D.f90) for per-R adiabatic channel data. Links
# against THIS repo's own Bsplines.o/matrix_stuff.o/Quadrature.o (confirmed to
# already provide everything adiabaticSolver1D.f90 needs -- CalcBasisFuncs,
# MyDsband, GetGaussFactors/GetGaussLobattoFactors -- since those are plain
# external subroutines/a shared MODULE Quadrature, not per-repo module forks),
# NOT Adiabatic-Scattering-BoundStates's own local lib/ copies, to avoid linking
# two competing copies of the same routines into one executable.
Potential.o: $(ASBS_DIR)/lib/Potential.f90
	${CMP} ${FORCEDP} -c $(ASBS_DIR)/lib/Potential.f90

# adiabaticSolver1D.f90 is a LOCAL copy of $(ASBS_DIR)/adiabaticSolver1D.f90 with one
# dead subroutine removed (see the file's own header note) to resolve a link-time name
# collision with RMatPropCore.f90's own GridMaker -- not a physics fork.
adiabaticSolver1D.o: adiabaticSolver1D.f90 Quadrature.o Potential.o
	${CMP} ${FORCEDP} ${CMPFLAGS} -c adiabaticSolver1D.f90

SVDChannelBasis.o: SVDChannelBasis.f90 RMatPropCore.o Quadrature.o adiabaticSolver1D.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c SVDChannelBasis.f90

# Analytic delta-function overlap formulas (SolveQ/EvalU/EvalOverlapAt), copied from
# Adiabatic-Scattering-BoundStates/DeltaScat/KMSFormulas.f -- keep in sync with that file.
KMSFormulas.o: KMSFormulas.f
	${CMP} ${CMPFLAGS} -c KMSFormulas.f

# Analytic drop-in replacement for BuildSVDBox (see AnalyticSVDChannelBasis.f90 header) --
# no B-spline/OneDimChannels dependency, only needs SVDChannelBasis.o (for
# GetGaussRadauFactors/SVDCalcDX reuse), Quadrature.o, and KMSFormulas.o.
AnalyticSVDChannelBasis.o: AnalyticSVDChannelBasis.f90 SVDChannelBasis.o Quadrature.o KMSFormulas.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c AnalyticSVDChannelBasis.f90

# Third variant: genuine numerical diagonalization (CalcBasisFuncsBP's Robin/Left=3-Right=3 BC,
# fed the exact delta-function log-derivative) instead of either OneDimChannels' Poschl-Teller
# well or KMSFormulas.f's closed-form shortcut -- see RobinSVDChannelBasis.f90's own header.
RobinSVDChannelBasis.o: RobinSVDChannelBasis.f90 SVDChannelBasis.o RMatPropCore.o Quadrature.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RobinSVDChannelBasis.f90

HYBRIDSVD_OBJS = $(SVD_OBJS) RMATPROPHybridSVD.o

RMATPROPHybridSVD.x: ${HYBRIDSVD_OBJS}
	${CMP} ${DEBUG} ${HYBRIDSVD_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPHybridSVD.x

RMATPROPHybridSVD.o: RMATPROPHybridSVD.f90 RMatPropCore.o AdiabaticPotential.o SVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPHybridSVD.f90

# All-analytic-SVD delta-function benchmark driver (see RMATPROPHybridSVDAnalytic.f90 header) --
# every box built by AnalyticSVDChannelBasis.f90, no B-spline/OneDimChannels/FitLeff.data
# dependency at all. Not part of `all`; build/run separately.
HYBRIDSVDANALYTIC_OBJS = $(SVD_OBJS) KMSFormulas.o AnalyticSVDChannelBasis.o RMATPROPHybridSVDAnalytic.o

RMATPROPHybridSVDAnalytic.x: ${HYBRIDSVDANALYTIC_OBJS}
	${CMP} ${DEBUG} ${HYBRIDSVDANALYTIC_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPHybridSVDAnalytic.x

RMATPROPHybridSVDAnalytic.o: RMATPROPHybridSVDAnalytic.f90 RMatPropCore.o SVDChannelBasis.o AnalyticSVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPHybridSVDAnalytic.f90

# All-Robin-SVD delta-function benchmark driver (see RMATPROPHybridSVDRobin.f90 header) -- every
# box built by RobinSVDChannelBasis.f90's genuine numerical diagonalization (CalcBasisFuncsBP's
# Robin BC, fed the exact log-derivative). Not part of `all`; build/run separately.
HYBRIDSVDROBIN_OBJS = $(SVD_OBJS) KMSFormulas.o RobinSVDChannelBasis.o RMATPROPHybridSVDRobin.o

RMATPROPHybridSVDRobin.x: ${HYBRIDSVDROBIN_OBJS}
	${CMP} ${DEBUG} ${HYBRIDSVDROBIN_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPHybridSVDRobin.x

RMATPROPHybridSVDRobin.o: RMATPROPHybridSVDRobin.f90 RMatPropCore.o SVDChannelBasis.o RobinSVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPHybridSVDRobin.f90

# Canonical shared home is $(LIB_DIR)/AdiabaticInterfaces.f90 -- promoted there from this repo
# once VeffAtomIon1D became a second real consumer; see AdiabaticInterfaces.f90's own header.
AdiabaticInterfaces.o: $(LIB_DIR)/AdiabaticInterfaces.f90
	${CMP} ${FORCEDP} ${CMPFLAGS} -c $(LIB_DIR)/AdiabaticInterfaces.f90

AdiabaticSolverGeneric.o: $(LIB_DIR)/AdiabaticSolverGeneric.f90 AdiabaticInterfaces.o
	${CMP} ${FORCEDP} ${CMPFLAGS} -c $(LIB_DIR)/AdiabaticSolverGeneric.f90

# Generic drop-in replacement for SVDChannelBasis.f90's BuildSVDBox (see file header) --
# builds an SVD box from ANY AdiabaticInterfaces-conforming PotentialProc/GridMakerProcI/
# RobinCoeffProc plugin instead of hardcoding a call to OneDimChannels. This is the engine
# VeffAtomIon1D's atom-ion pipeline plugs into via AtomIonPotential.f90's own procs.
GENERICSVD_OBJS = $(SVD_OBJS) AdiabaticInterfaces.o AdiabaticSolverGeneric.o GenericSVDChannelBasis.o

GenericSVDChannelBasis.o: GenericSVDChannelBasis.f90 RMatPropCore.o Quadrature.o AdiabaticInterfaces.o AdiabaticSolverGeneric.o SVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c GenericSVDChannelBasis.f90

# Delta-function plugin (PotentialProc/GridMakerProcI/RobinCoeffProc) for the generic engine --
# the AdiabaticInterfaces-conforming counterpart to adiabaticSolver1D.f90's hardcoded physics.
DeltaFunctionPlugins.o: DeltaFunctionPlugins.f90 RMatPropCore.o AdiabaticInterfaces.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c DeltaFunctionPlugins.f90

# Generic hybrid-SVD driver engine (RunHybridSVDGeneric), physics-agnostic -- see
# RMATPROPHybridSVDGenericDelta.f90 for the delta-function plugin instantiation.
RMATPROPHybridSVDGeneric.o: RMATPROPHybridSVDGeneric.f90 RMatPropCore.o AdiabaticInterfaces.o SVDChannelBasis.o GenericSVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPHybridSVDGeneric.f90

HYBRIDSVDGENERICDELTA_OBJS = ${GENERICSVD_OBJS} DeltaFunctionPlugins.o RMATPROPHybridSVDGeneric.o RMATPROPHybridSVDGenericDelta.o

RMATPROPHybridSVDGenericDelta.x: ${HYBRIDSVDGENERICDELTA_OBJS}
	${CMP} ${DEBUG} ${HYBRIDSVDGENERICDELTA_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPHybridSVDGenericDelta.x

RMATPROPHybridSVDGenericDelta.o: RMATPROPHybridSVDGenericDelta.f90 RMatPropCore.o RMATPROPHybridSVDGeneric.o DeltaFunctionPlugins.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPHybridSVDGenericDelta.f90

# sech^2-well plugin (PotentialProc/GridMakerProcI/RobinCoeffProc) for the generic engine --
# independent numerical cross-check of RMATPROPAdiabatic.x's tabulated-data baseline on the
# 3bodydata_sech2_5ch/_fine datasets. See SechFunctionPlugins.f90's own header.
SechFunctionPlugins.o: SechFunctionPlugins.f90 RMatPropCore.o AdiabaticInterfaces.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c SechFunctionPlugins.f90

HYBRIDSVDGENERICSECH_OBJS = ${GENERICSVD_OBJS} SechFunctionPlugins.o RMATPROPHybridSVDGeneric.o RMATPROPHybridSVDGenericSech.o

RMATPROPHybridSVDGenericSech.x: ${HYBRIDSVDGENERICSECH_OBJS}
	${CMP} ${DEBUG} ${HYBRIDSVDGENERICSECH_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPHybridSVDGenericSech.x

RMATPROPHybridSVDGenericSech.o: RMATPROPHybridSVDGenericSech.f90 RMatPropCore.o RMATPROPHybridSVDGeneric.o SechFunctionPlugins.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPHybridSVDGenericSech.f90

# Same SechFunctionPlugins.f90 physics as above, but for the THREE-IDENTICAL-BOSON case (domain
# [0,pi/6], Left/Right=1/1) instead of the 2-fermion+1-distinguishable case (domain [0,pi/2],
# Left/Right=1/0) -- cross-checks RMATPROPAdiabatic.x's 3bodydata_sech2_3boson* baseline.
HYBRIDSVDGENERICSECH3BOSON_OBJS = ${GENERICSVD_OBJS} SechFunctionPlugins.o RMATPROPHybridSVDGeneric.o RMATPROPHybridSVDGenericSech3Boson.o

RMATPROPHybridSVDGenericSech3Boson.x: ${HYBRIDSVDGENERICSECH3BOSON_OBJS}
	${CMP} ${DEBUG} ${HYBRIDSVDGENERICSECH3BOSON_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPHybridSVDGenericSech3Boson.x

RMATPROPHybridSVDGenericSech3Boson.o: RMATPROPHybridSVDGenericSech3Boson.f90 RMatPropCore.o RMATPROPHybridSVDGeneric.o SechFunctionPlugins.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPHybridSVDGenericSech3Boson.f90

# Standalone A/B comparator (Validation Playbook technique): single-channel, single-box R-matrix,
# tabulated/adiabatic path vs SVD/direct path, for the SAME 3-identical-boson sech^2 physics --
# isolates box construction from CalcK's asymptotic matching and from the full driver's cost.
COMPARERMATRIX_OBJS = ${GENERICSVD_OBJS} SechFunctionPlugins.o CompareRmatrixSVDvsAdiabatic.o

CompareRmatrixSVDvsAdiabatic.x: ${COMPARERMATRIX_OBJS}
	${CMP} ${DEBUG} ${COMPARERMATRIX_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o CompareRmatrixSVDvsAdiabatic.x

CompareRmatrixSVDvsAdiabatic.o: CompareRmatrixSVDvsAdiabatic.f90 RMatPropCore.o AdiabaticPotential.o SechFunctionPlugins.o GenericSVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c CompareRmatrixSVDvsAdiabatic.f90

# Box-by-box wavefunction reconstruction, single channel, adiabatic B-spline chain vs SVD/direct
# chain -- extends CompareRmatrixSVDvsAdiabatic.x's single-box test to a genuine multi-box chain.
PLOTWAVEFUNCTIONBOXBYBOX_OBJS = ${GENERICSVD_OBJS} SechFunctionPlugins.o PlotWavefunctionBoxByBox.o

PlotWavefunctionBoxByBox.x: ${PLOTWAVEFUNCTIONBOXBYBOX_OBJS}
	${CMP} ${DEBUG} ${PLOTWAVEFUNCTIONBOXBYBOX_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o PlotWavefunctionBoxByBox.x

PlotWavefunctionBoxByBox.o: PlotWavefunctionBoxByBox.f90 RMatPropCore.o AdiabaticPotential.o SechFunctionPlugins.o GenericSVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c PlotWavefunctionBoxByBox.f90

# Sanity check: does the final propagated Rmat depend on box count (1/2/3 boxes spanning the
# same total range)? Should not -- isolates BoxMatch chaining correctness from the separate
# wavefunction-reconstruction-only kink issue.
TESTBOXCOUNTRMAT_OBJS = ${GENERICSVD_OBJS} SechFunctionPlugins.o TestBoxCountRmat.o

TestBoxCountRmat.x: ${TESTBOXCOUNTRMAT_OBJS}
	${CMP} ${DEBUG} ${TESTBOXCOUNTRMAT_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o TestBoxCountRmat.x

TestBoxCountRmat.o: TestBoxCountRmat.f90 RMatPropCore.o AdiabaticPotential.o SechFunctionPlugins.o GenericSVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c TestBoxCountRmat.f90

# Direct test: does embedding channel 1 in a NumStates=5 SVD box construction (matching what
# the production driver actually calls) change its own R-matrix relative to building it alone
# (NumStates=1, already validated), at E=-0.999 where channels 2-5 are deep closed and should
# contribute negligibly if the physics genuinely factors as F_1(R)*Phi_1(R,phi)?
COMPARERMATRIXMULTICHANNEL_OBJS = ${GENERICSVD_OBJS} SechFunctionPlugins.o CompareRmatrixMultichannel.o

CompareRmatrixMultichannel.x: ${COMPARERMATRIXMULTICHANNEL_OBJS}
	${CMP} ${DEBUG} ${COMPARERMATRIXMULTICHANNEL_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o CompareRmatrixMultichannel.x

CompareRmatrixMultichannel.o: CompareRmatrixMultichannel.f90 RMatPropCore.o AdiabaticPotential.o SechFunctionPlugins.o GenericSVDChannelBasis.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c CompareRmatrixMultichannel.f90

# No-interpolation B-spline-box adiabatic potential (DeltaScat-style: real ARPACK diagonalization
# directly at each box's own Gauss-Legendre quadrature points via SechFunctionPlugins, not a
# tabulated-then-interpolated Uad.dat/Pmat.dat/Qmat.dat) -- see this file's own header.
SechAdiabaticPotentialDirect.o: SechAdiabaticPotentialDirect.f90 RMatPropCore.o AdiabaticInterfaces.o AdiabaticSolverGeneric.o SechFunctionPlugins.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c SechAdiabaticPotentialDirect.f90

RMATPROPADIABATICDIRECT3BOSON_OBJS = besselnew.o Bsplines.o matrix_stuff.o Quadrature.o Interpolation.o Potential.o adiabaticSolver1D.o RMatPropCore.o AdiabaticInterfaces.o AdiabaticSolverGeneric.o SechFunctionPlugins.o SechAdiabaticPotentialDirect.o RMATPROPAdiabaticDirect3Boson.o

RMATPROPAdiabaticDirect3Boson.x: ${RMATPROPADIABATICDIRECT3BOSON_OBJS}
	${CMP} ${DEBUG} ${RMATPROPADIABATICDIRECT3BOSON_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPAdiabaticDirect3Boson.x

RMATPROPAdiabaticDirect3Boson.o: RMATPROPAdiabaticDirect3Boson.f90 RMatPropCore.o SechAdiabaticPotentialDirect.o
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROPAdiabaticDirect3Boson.f90

clean:
	rm -f *.mod *.o *.x
