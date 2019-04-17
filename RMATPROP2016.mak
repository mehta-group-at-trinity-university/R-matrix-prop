CMP     = gfortran
CMPFLAGS = -ffixed-line-length-132 -O3
DEBUG   =
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE = -I/usr/local/opt/lapack/include
LAPACK =  -framework accelerate
ARPACK = -L/Users/mehtan/Code/ARPACK/ARPACK -larpack_OSX
OBJS  = besselnew.o Bsplines.o matrix_stuff.o RMATPROP2016.o

RMATPROP2016.x:	   ${OBJS}
	${CMP} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK}  ${CMPFLAGS} ${FORCEDP} -o RMATPROP2016.x 

RMATPROP2016.o: RMATPROP2016.f90
	${CMP} ${FORCEDP} ${CMPFLAGS} -c RMATPROP2016.f90

matrix_stuff.o: matrix_stuff.f	
	${CMP} ${FORCEDP} ${CMPFLAGS} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c Bsplines.f	

nrtype.mod: modules_qd.o
	${CMP} ${FORCEDP} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${FORCEDP} -c modules_qd.f90	

besselnew.o:	besselnew.f
	${CMP} ${FORCEDP} ${CMPFLAGS} -c besselnew.f

clean:
	rm -f *.mod *.o *.x
