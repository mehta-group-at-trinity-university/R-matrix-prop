! Standalone diagnostic: dump the actual interpolated combo(R) = 2*mu*U(R) + Q11(R)
! at a fine R grid, using the exact same ReadAdiabaticData/Interpolated/InterpolatedMatrix
! pipeline RMATPROPAdiabatic.x uses -- so "interpolated" here means precisely what the
! propagator's SetAdiabaticPotential sees, not a hand-rolled re-implementation.
PROGRAM DumpCombo
  USE InterpType
  USE AdiabaticPotential
  IMPLICIT NONE
  CHARACTER(LEN=200) DataDir
  INTEGER NumChannels, NumDataPoints
  DOUBLE PRECISION muLocal, alpha, EffDimLocal
  DOUBLE PRECISION, ALLOCATABLE :: Threshold(:), Leff(:)
  TYPE(InterpolatingFunction), ALLOCATABLE :: UInterp(:)
  TYPE(InterpolatingMatrix) :: PInterp, QInterp
  DOUBLE PRECISION QM(1,1), R, combo
  INTEGER i, NumPts
  DOUBLE PRECISION Rmin, Rmax

  CALL GET_COMMAND_ARGUMENT(1,DataDir)
  CALL ReadAdiabaticData(TRIM(DataDir),NumChannels,NumDataPoints,muLocal,alpha,EffDimLocal, &
       Threshold,Leff,UInterp,PInterp,QInterp)

  Rmin = 6d-5
  Rmax = 0.5d0
  NumPts = 400
  OPEN(unit=30,file='combo_interpolated.dat',status='replace')
  DO i = 1,NumPts
     R = Rmin*(Rmax/Rmin)**(DBLE(i-1)/DBLE(NumPts-1))   ! log-spaced
     QM = InterpolatedMatrix(R,QInterp)
     combo = 2d0*muLocal*Interpolated(R,UInterp(1)) + QM(1,1)
     WRITE(30,'(1P,3E20.12)') R, combo, QM(1,1)
  ENDDO
  CLOSE(30)
  WRITE(6,*) 'wrote combo_interpolated.dat, mu=',muLocal
END PROGRAM DumpCombo
