MODULE Quadrature
  IMPLICIT NONE
  INTEGER LegPoints
  DOUBLE PRECISION, ALLOCATABLE :: xLeg(:),wLeg(:)
  CHARACTER*64 LegendreFile

CONTAINS
  SUBROUTINE GetGaussFactors(File,Points,x,w)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c     This subroutine retrieves the points and weights for
    !c      Gaussian Quadrature from a given file
    !c
    !c     Variables:
    !c      File		name of file containing points and
    !c			 weights
    !c      Points		number of points to be used for
    !c			quadrature
    !c      x		array of nodes
    !c      w		array of weights
    !c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    IMPLICIT NONE
    INTEGER Points
    DOUBLE PRECISION x(Points),w(Points)
    CHARACTER*64 File
    INTEGER i,tempPoints

    OPEN(unit=7,file=File(1:INDEX(File,' ')-1))
    DO i = 1,18
       READ(7,*)
    ENDDO
    READ(7,*) tempPoints
    DO WHILE (tempPoints .NE. Points)
       DO i = 1,tempPoints
          READ(7,*)
       ENDDO
       READ(7,*) tempPoints
    ENDDO
    DO i = 1,Points
       READ(7,*) x(i),w(i)
    ENDDO
    CLOSE(unit=7)
    RETURN
  END SUBROUTINE GetGaussFactors

END MODULE Quadrature
