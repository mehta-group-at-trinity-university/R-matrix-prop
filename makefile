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

all: RMATPROP2016.x RMATPROPAdiabatic.x

RMATPROP2016.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK}  ${CMPFLAGS} ${FORCEDP} -o RMATPROP2016.x

RMATPROPAdiabatic.x: ${ADIABATIC_OBJS}
	${CMP} ${DEBUG} ${ADIABATIC_OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROPAdiabatic.x

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
