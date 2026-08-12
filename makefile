CMP     = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
ARPACK =  -L/opt/homebrew/lib/ -larpack
INCLUDE =  -I/opt/homebrew/include
LAPACK =  -framework Accelerate
INTERP_DIR = $(HOME)/Documents/GitHub/interpolation
LIB_DIR = $(HOME)/Documents/GitHub/lib
OBJS  = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o RMATPROP2016.o Quadrature.o Interpolation.o
ADIABATIC_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o RMATPROPAdiabatic.o
ADIABATIC_SEEDED_OBJS = besselnew.o Bsplines.o matrix_stuff.o RMatPropCore.o Quadrature.o Interpolation.o AdiabaticPotential.o RMATPROPAdiabatic_seeded.o

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

clean:
	rm -f *.mod *.o *.x
