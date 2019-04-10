CMP     = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3
DEBUG   =
LAPACK = -framework accelerate
ARPACK = -L/Users/mehtan/Code/ARPACK/ARPACK -larpack_OSX
OBJS  = besselnew.o Bsplines.o matrix_stuff.o RMATPROP2016.o

RMATPROP2016.x:	   ${OBJS}
	${CMP} ${CMPFLAGS} ${OBJS} ${LAPACK} ${ARPACK}  -o RMATPROP2016.x 

RMATPROP2016.o: RMATPROP2016.f90
	${CMP} ${CMPFLAGS} -c RMATPROP2016.f90

matrix_stuff.o: matrix_stuff.f	
	${CMP} ${CMPFLAGS} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${CMPFLAGS} -c Bsplines.f	

nrtype.mod: modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} -c modules_qd.f90	

besselnew.o:	besselnew.f
	${CMP} ${CMPFLAGS} -c besselnew.f

clean:
	rm -f *.mod *.o *.x
