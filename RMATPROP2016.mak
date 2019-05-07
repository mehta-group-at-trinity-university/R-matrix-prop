CMP = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3 -fdec
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE = -I/usr/local/opt/lapack/include
LAPACK =  -framework accelerate
ARPACK = -L/Users/mehtan/Code/ARPACK/ARPACK -larpack_OSX
OBJS = Quadrature.o DataStructures.o besselnew.o Bsplines.o matrix_stuff.o zgensub.o MorsePotential.o BalujaParameters.o DipoleDipole.o RMATPROP2016.o

RMATPROP2016.x:	   ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${CMPFLAGS} ${FORCEDP} -o RMATPROP2016.x

RMATPROP2016.o: RMATPROP2016.f90
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c RMATPROP2016.f90

matrix_stuff.o: matrix_stuff.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c Bsplines.f

nrtype.mod: modules_qd.o
	${CMP} ${FORCEDP} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${FORCEDP} -c modules_qd.f90

MorsePotential.mod: MorsePotential.o
	${CMP} ${FORCEDP} MorsePotential.o

MorsePotential.o: MorsePotential.f90
	${CMP} ${FORCEDP} -c MorsePotential.f90

besselnew.o:	besselnew.f
	${CMP} ${DEBUG} ${FORCEDP} ${CMPFLAGS} -c besselnew.f

Quadrature.mod: Quadrature.o
	${CMP} ${FORCEDP} Quadrature.o

Quadrature.o: Quadrature.f90
	${CMP} ${FORCEDP} -c Quadrature.f90

zgensub.o: zgensub.f
	${CMP} ${FORCEDP} -c zgensub.f

BalujaParameters.mod: BalujaParameters.o
	${CMP} ${FORCEDP} BalujaParameters.o

BalujaParameters.o: BalujaParameters.f90
	${CMP} ${FORCEDP} -c BalujaParameters.f90

DipoleDipole.mod: DipoleDipole.o
	${CMP} ${FORCEDP} DipoleDipole.o

DipoleDipole.o: DipoleDipole.f90
	${CMP} ${FORCEDP} -c DipoleDipole.f90

DataStructures.mod: DataStructures.o
	${CMP} ${FORCEDP} DataStructures.o

DataStructures.o: DataStructures.f90
	${CMP} ${FORCEDP} -c DataStructures.f90
clean:
	rm -f *.mod *.o *.x
