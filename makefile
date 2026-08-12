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

Quadrature.o: Quadrature.f90
	${CMP} ${FORCEDP} -c Quadrature.f90

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

clean:
	rm -f *.mod *.o *.x
