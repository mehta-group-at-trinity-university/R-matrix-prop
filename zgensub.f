********************************************************************************
C      gensub.f  -- standard subroutines, to calculate wigner coeffs,
c			interpolation, root-finding, numerical integration
c			Coulomb function programs also.
C**********************************************************C
C
C WIGNER 3J, 6J, AND 9J FUNCTION SUBPROGRAMS.
C NOTE:  ALL ANGULAR MOMENTUM QUANTUM NUMBERS SUPPLIED TO THESE
C        FUNCTIONS ARE INTEGERS WHICH ARE TWICE THE VALUE OF THE
C        ACTUAL ANGULAR MOMENTUM DESIRED.  (THIS ALLOWS FOR HALF-
C        INTEGER VALUES CONVENIENTLY.)  ALSO YOU MUST CALL SETUPAM
C        ONE TIME BEFORE CALLING ANY OF THESE SO THAT THE RELE-
C        VANT FACTORIALS CAN BE CALCULATED ONCE AND FOR ALL AND
C        STORED, AS IN THE ABOVE EXAMPLE.
C
C**********************************************************C
C
      function coeff(l1,l2,l1p,l2p,l,k)
	    implicit real*8(a-h,o-z)
			front=(2*l1+1)*(2*l2+1)*(2*l1p+1)*(2*l2p+1)
			front=dsqrt(front)*(-1)**(l1+l1p+l)
			l1d=2*l1
			l2d=2*l2
			l1pd=2*l1p
			l2pd=2*l2p
			ld=2*l
			kd=2*k
			iz=0
			t1=thrj(l1d,kd,l1pd,iz,iz,iz)
			t2=thrj(l2d,kd,l2pd,iz,iz,iz)
			s1=sixj(l1d,l2d,ld,l2pd,l1pd,kd)
			coeff=front*t1*t2*s1
			return
   		end
c
      FUNCTION XNINEJ(J11,J12,J13,J21,J22,J23,J31,J32,J33)
      IMPLICIT REAL*8(A-H,O-Z)
      KMIN1 = IABS(J11-J33)
      KMIN2 = IABS(J32-J21)
      KMIN3 = IABS(J23-J12)
      KMAX1 = J11+J33
      KMAX2 = J32+J21
      KMAX3 = J23+J12
C
      IF(KMIN2.GT.KMIN1) KMIN1=KMIN2
       IF(KMIN3.GT.KMIN1) KMIN1=KMIN3
      IF(KMAX2.LT.KMAX1) KMAX1=KMAX2
      IF(KMAX3.LT.KMAX1) KMAX1=KMAX3
C
      KMIN1 = KMIN1 + 1
      KMAX1 = KMAX1 + 1
       XNINEJ = 0.D0
      IF(KMIN1.GT.KMAX1) GO TO 1000
      DO 100 K1 = KMIN1,KMAX1,2
      K = K1 - 1
      S1 = SIXJ(J11,J21,J31,J32,J33,K)
      S2 = SIXJ(J12,J22,J32,J21,K,J23)
      S3 = SIXJ(J13,J23,J33,K,J11,J12)
      P = (K+1)*((-1)**K)
      XNINEJ = XNINEJ + P*S1*S2*S3
  100 CONTINUE
 1000 CONTINUE
      RETURN
      END
C
      FUNCTION THRJ(J1D,J2D,J3D,M1D,M2D,M3D)
      IMPLICIT REAL*8(A-H,O-Z)
      X1 = J1D/2.D0
      X2 = J2D/2.D0
      X3 = J3D/2.D0
      Y1 = M1D/2.D0
      Y2 = M2D/2.D0
      Y3 = M3D/2.D0
C
C -- NEXT COME THE TRIANGULARITY CHECKS:
C
      IF(J1D+J2D-J3D.LT.0) GO TO 9998
      IF(J2D+J3D-J1D.LT.0) GO TO 9998
      IF(J3D+J1D-J2D.LT.0) GO TO 9998
      IF(J3D.LT.IABS(J1D-J2D)) GO TO 9998
      IF(M1D+M2D+M3D.NE.0) GO TO 9998
      LLL = J1D+J2D+J3D
      IF(2*(LLL/2) .NE. LLL) GO TO 9998
C
      KMIN = (J3D-J1D-M2D)/2
      KMIN1 = KMIN
      KMIN2 = (J3D-J2D+M1D)/2
      IF(KMIN2.LT.KMIN) KMIN=KMIN2
      KMIN = (-1)*KMIN
      KMAX = X1+X2-X3 +0.1D0
      KMAX1 = KMAX
      KMAX2 = X1-Y1
      KMAX3 = X2+Y2
      IF(KMAX2.LT.KMAX) KMAX=KMAX2
      IF(KMAX3.LT.KMAX) KMAX=KMAX3
      IF(KMIN.LT.0) KMIN = 0
      IF(KMIN.GT.KMAX) GO TO 9998
C
      JMIN = KMIN+1
      JMAX = KMAX+1
      TERM1 = FRONTL(X1,X2,X3,Y1,Y2,Y3)
	iphase=iabs((j1d-j2d-m3d)/2)
	msign=(-1)**iphase
cg     MSIGN = (-1)**((J1D-J2D-M3D)/2)
      SUM = 0.D0
      DO 10 I1 = JMIN,JMAX
      I = I1 - 1
      TERM2 = FL(I1)+FL(KMIN1+I1)+FL(KMIN2+I1)
      TERM2 = TERM2+FL(KMAX1-I+1)+FL(KMAX2-I+1)+FL(KMAX3-I+1)
      TERM= DEXP(TERM1-TERM2)
      TERM = TERM*MSIGN*((-1)**I)
      SUM = SUM + TERM
  10  CONTINUE
      THRJ = SUM
      GO TO 9999
 9998 THRJ = 0.D0
 9999 CONTINUE
      RETURN
      END
C
      FUNCTION FL(I)
       IMPLICIT REAL*8(A-H,O-Z)
CCC   DIMENSION FACL(60)
      COMMON/FACTOR/FACL(200)
      FL = FACL(I)
      RETURN
      END
C
C** ** **
C-- THIS SUBROUTINE INITIALIZES BY FINDING THE LOGARITHM
C---OF THE FIRST 199 FACTORIALS AND STORING THEM.
C
      SUBROUTINE SETUPAM
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/FACTOR/FACL(200)
      N = 199
      FACL(1) = 0.D0
      DO 100 I = 2,N
      I1 = I - 1
      FACL(I) = FACL(I1) + DLOG(I1*1.D0)
  100 CONTINUE
      RETURN
      END
C** ** **
C
      FUNCTION FRONTL(X1,X2,X3,Y1,Y2,Y3)
      IMPLICIT REAL*8(A-H,O-Z)
      L1 = X1+X2-X3 +1.1D0
      L2 = X2+X3-X1 +1.1D0
      L3 = X3+X1-X2 +1.1D0
      L4 = X1+X2+X3+1+1.1D0
      L5 = X1+Y1+1.1D0
      L6 = X1-Y1+1.1D0
      L7 = X2+Y2+1.1D0
      L8 = X2-Y2+1.1D0
      L9 = X3+Y3+1.1D0
      L10 = X3-Y3+1.1D0
      FRONTL = FL(L1)+FL(L2)+FL(L3)-FL(L4)+FL(L5)+FL(L6)
      FRONTL = FRONTL +FL(L7)+FL(L8)+FL(L9)+FL(L10)
      FRONTL = FRONTL/2.D0
      RETURN
      END
C
      FUNCTION SIXJ(J1D,J2D,J3D,J4D,J5D,J6D)
      IMPLICIT REAL*8(A-H,O-Z)
C
C -- CHECK THAT TRIANGULARITY CONDITIONS ARE SATISFIED.
C
      IF(J3D.LT.IABS(J1D-J2D) .OR. J3D.GT.J1D+J2D) GO TO 9998
      IF(J6D.LT.IABS(J4D-J2D) .OR. J6D.GT.J4D+J2D) GO TO 9998
      IF(J6D.LT.IABS(J1D-J5D) .OR. J6D.GT.J1D+J5D) GO TO 9998
      IF(J3D.LT.IABS(J4D-J5D) .OR. J3D.GT.J4D+J5D) GO TO 9998
      K1=J1D+J2D+J3D
      K2=J4D+J2D+J6D
      K3=J6D+J1D+J5D
      K4=J3D+J4D+J5D
      IF(2*(K1/2).NE.K1 .OR. 2*(K2/2).NE.K2) GO TO 9998
      IF(2*(K3/2).NE.K3 .OR. 2*(K4/2).NE.K4) GO TO 9998
C
C -- NOW GO AHEAD AND CALCULATE THE SIXJ.
C
      JM1 = (J1D+J2D+J3D)/2
      JM2 = (J1D+J5D+J6D)/2
      JM3 = (J4D+J2D+J6D)/2
      JM4 = (J4D+J5D+J3D)/2
      JX1 = (J1D+J2D+J4D+J5D)/2
      JX2 = (J2D+J3D+J5D+J6D)/2
      JX3 = (J3D+J1D+J6D+J4D)/2
C
      JM = JM1
      IF(JM2.GT.JM) JM = JM2
      IF(JM3.GT.JM) JM = JM3
      IF(JM4.GT.JM) JM = JM4
      JX = JX1
      IF(JX2.LT.JX) JX = JX2
      IF(JX3.LT.JX) JX = JX3
      KM = JM+1
      KX = JX+1
      IF(KM.GT.KX) GO TO 9998
      TERM1 = FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)
      SIXJ = 0.D0
      DO 10 I1 = KM,KX
      I = I1 - 1
      TERM2 = FL(I+2)-FL(I+1-JM1)-FL(I+1-JM2)-FL(I+1-JM3)
      TERM2 = TERM2-FL(I+1-JM4)-FL(JX1-I+1)-FL(JX2-I+1)
      TERM2 = TERM2 -FL(JX3-I+1)
      TERM = DEXP(TERM1+TERM2) * ((-1)**I)
      SIXJ = SIXJ + TERM
   10 CONTINUE
      GO TO 9999
 9998 CONTINUE
      SIXJ = 0.D0
 9999 CONTINUE
      RETURN
      END
C
      FUNCTION FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)
      IMPLICIT REAL*8(A-H,O-Z)
      FRTSJL = DL(J1D,J2D,J3D) + DL(J1D,J5D,J6D)
      FRTSJL = FRTSJL + DL(J4D,J2D,J6D) + DL(J4D,J5D,J3D)
      RETURN
      END
C
      FUNCTION DL(J1D,J2D,J3D)
      IMPLICIT REAL*8(A-H,O-Z)
      L1 = (J1D+J2D-J3D)/2
      L2 = (J2D+J3D-J1D)/2
      L3 = (J3D+J1D-J2D)/2
      L4 = (J1D+J2D+J3D)/2 + 1
      DL = FL(L1+1)+FL(L2+1)+FL(L3+1)-FL(L4+1)
      DL = DL/2.D0
      RETURN
      END
C
      FUNCTION CLEBSH(J1D,J2D,JD,M1D,M2D,MD)
      IMPLICIT REAL*8(A-H,O-Z)
      CLEBSH = 0.D0
      IF(M1D+M2D.NE.MD) GO TO 100
      MMD = -MD
      Q = THRJ(J1D,J2D,JD,M1D,M2D,MMD)
      PHASE = ((-1)**(20 + (J2D-J1D-MD)/2))
      CLEBSH = Q*PHASE*DSQRT(JD+1.D0)
  100 CONTINUE
      RETURN
      END
************************************************************************
        function yfix(x,vx,v,ndata,igo)
        implicit real*8(a-h,o-z)
c
c --- interpolation routine yfix.
c --- vx and v are input data arrays for the x and y data points,
c -----  ndata is the number of points in these arrays, and x is the
c -----  final x-value at which y=yfix is desired.  Note that the
c -----  vx-values must be strictly increasing.  The parameter igo
c -----  should be set to 0 on the first call at a given x, and
c -----  different from 0 on subsequent calls at the same x.
c
        dimension vx(ndata),v(ndata)
        nx = ndata - 1
        if(x.ge.vx(nx)) go to 75
        if(x.le.vx(2)) go to 65
        if(igo.ne.0) go to 20
        n1 = 2
        n2 = ndata - 1
        i = 0
   10   continue
        i =   i + 1
        j = (n1+n2)/2
        if(vx(j).le.x) go to 15
        if(n2-n1 .le. 1) go to 20
        n2 = j
        go  to 10
   15   continue
        n1 = j
        if(n2-n1 .le. 1) go to 20
        go to 10
   20   continue
        n = n2-1
        np = n2
        np1 = np + 1
        n1 = n - 1
        x0 = vx(n1)
        y0 = v(n1)
        x1 = vx(n)
        y1 = v(n)
        x2 = vx(np)
        y2 = v(np)
        x3 = vx(np1)
        y3 = v(np1)
        q0 = (x-x1)*(x-x2)*(x-x3)/( (x0-x1)*(x0-x2)*(x0-x3) )
        q1 = (x-x0)*(x-x2)*(x-x3)/( (x1-x0)*(x1-x2)*(x1-x3) )
        q2 = (x-x0)*(x-x1)*(x-x3)/( (x2-x0)*(x2-x1)*(x2-x3) )
        q3 = (x-x0)*(x-x1)*(x-x2)/( (x3-x0)*(x3-x1)*(x3-x2) )
        y = q0*y0 + q1*y1 + q2*y2 + q3*y3
        go to 100
   65   n1 = 2
        n2 = 3
        go to 20
   75   nx = ndata - 1
        slope = (v(ndata)-v(nx))/(vx(ndata)-vx(nx))
        y = v(ndata) - slope*(vx(ndata)-x)
  100   yfix = y
        return
        end
***********************************************************************
**********************************************************************
	function rint(f,na,nb,nq,h)
c
c -- function rint performs a numerical integration over the function f,
c -----  assuming an equally spaced x-axis grid with step size h and a
c -----  total of nb mesh points.  Use na=1, nq=10, and nb>20 always.
c
	implicit real*8(a-h,o-z)
	dimension c(55),d(10),f(nb)
	data c/1.d0,
     1  2.d0,1.d0,
     2  23.d0, 28.d0, 9.d0,
     3  25.d0, 20.d0, 31.d0, 8.d0,
     4  1413.d0, 1586.d0, 1104.d0, 1902.d0, 475.d0,
     5  1456.d0, 1333.d0, 1746.d0, 944.d0,1982.d0,459.d0,
     6  119585.d0, 130936.d0, 89437.d0, 177984.d0, 54851.d0,
     7  176648.d0, 36799.d0,
     8  122175.d0, 111080.d0, 156451.d0,46912.d0,220509.d0,
     9  29336.d0, 185153.d0, 35584.d0,
     a  7200319.d0, 7783754.d0, 5095890.d0,12489922.d0,-1020160.d0,
     b  16263486.d0, 261166.d0, 11532470.d0, 2082753.d0,
     c  7305728.d0, 6767167.d0, 9516362.d0,1053138.d0, 18554050.d0,
     d   -7084288.d0, 20306238.d0, -1471442.d0, 11965622.d0,
     e  2034625.d0/
	data d/  2.d0,    2.d0, 24.d0, 24.d0,1440.d0,1440.d0,120960.d0,
     &     120960.d0, 7257600.d0, 7257600.d0/
	a=0.d0
	l=na
	m=nb
	i=nq*(nq+1)/2
	do 10 j=1,nq
	a=a+c(i)*(f(l)+f(m))
	l=l+1
	m=m-1
   10   i=i-1
	a=a/d(nq)
	do 20 n=l,m
   20   a=a+f(n)
	rint=a*h
	return
	end
c***********************************************************************
c***********************************************************************
	function zrint(zf,na,nb,nq,h)
c
c -- function rint performs a numerical integration over the function f,
c -----  assuming an equally spaced x-axis grid with step size h and a
c -----  total of nb mesh points.  Use na=1, nq=10, and nb>20 always.
c
	implicit real*8(a-h,o-y),complex*16(z)
	dimension c(55),d(10),zf(nb)
	data c/1.d0,
     1  2.d0,1.d0,
     2  23.d0, 28.d0, 9.d0,
     3  25.d0, 20.d0, 31.d0, 8.d0,
     4  1413.d0, 1586.d0, 1104.d0, 1902.d0, 475.d0,
     5  1456.d0, 1333.d0, 1746.d0, 944.d0,1982.d0,459.d0,
     6  119585.d0, 130936.d0, 89437.d0, 177984.d0, 54851.d0,
     7  176648.d0, 36799.d0,
     8  122175.d0, 111080.d0, 156451.d0,46912.d0,220509.d0,
     9  29336.d0, 185153.d0, 35584.d0,
     a  7200319.d0, 7783754.d0, 5095890.d0,12489922.d0,-1020160.d0,
     b  16263486.d0, 261166.d0, 11532470.d0, 2082753.d0,
     c  7305728.d0, 6767167.d0, 9516362.d0,1053138.d0, 18554050.d0,
     d   -7084288.d0, 20306238.d0, -1471442.d0, 11965622.d0,
     e  2034625.d0/
	data d/  2.d0,    2.d0, 24.d0, 24.d0,1440.d0,1440.d0,120960.d0,
     &     120960.d0, 7257600.d0, 7257600.d0/
	za=dcmplx(0.d0,0.d0)
	l=na
	m=nb
	i=nq*(nq+1)/2
	do 10 j=1,nq
	za=za+c(i)*(zf(l)+zf(m))
	l=l+1
	m=m-1
   10   i=i-1
	za=za/d(nq)
	do 20 n=l,m
   20   za=za+zf(n)
	zrint=za*h
	return
	end
c***********************************************************************
c*****************************************************************************
        double precision function delt(i,j)
        implicit double precision (a-h,o-z)
        delt=0.d0
        if(i.eq.j) delt=1.d0
        return
        end
c***********************************************************************
        subroutine root(xtry,hh,iroute,ep,xest,yest,yfcn)
        implicit real*8(a-h,o-z)
		external yfcn


c***  subroutine root locates a root of yfcn(x,y=0),
c***  returning finally (xest,yest)
c***  xtry is initial guess,hh is step size for the iteration.
c***  iroute=0 means hh is kept constant in the initial (coarse)
c***  root location loop.
c***  iroute=1 means hh can change sign but not magnitude.
c***  iroute=2 means both magnitude & sign of hh can vary.
c***  after the initial loop, a 3-point method is used to
c***  speed convergence.
c***
        iq = 3
        u1 = xtry
	f1=yfcn(u1)
        mtest = 0
        do 420 ki = 1,30
        if(ki.eq.2 .and. iroute.eq.2) go to 421
        if(mtest.ne.0) go to 920
        u2 = u1 + hh
  421   f2=yfcn(u2)
        if(f1*f2 .le. 0.d0) mtest = 1000
        if(mtest.ne.0)  go to 3055
        q1 = dabs(f1)
        q2 = dabs(f2)
        if(ki.eq.1 .and. iroute.eq.2) go to 3047
        if(iroute-1) 3033,3044,3044
 3044   if(q1.ge.q2 .and. mtest.eq.0) u1=u2
        if(q1.ge.q2 .and. mtest.eq.0) f1=f2
        if(q2.gt.q1 .and. mtest.eq.0) hh = (-1)*hh
        go to 3055
 3033   if(mtest.eq.0) u1 = u2
        if(mtest.eq.0) f1 = f2
        go to 3055
 3047   slp = (f2-f1)/(u2-u1)
        uest = u1 - f1/slp
        if(q2.lt.q1) f1 = f2
        if(q2.lt.q1) u1 = u2
        u2 = uest
 3055   continue
  420   continue
        if(mtest.ne.1000)  go to 8500
c
  920   continue
ccccccccccccccccccccccccccccccccc      ep = 1.d-06
        mchk = 0
        do 440 i = 1,50
        if(mchk.ne.0) go to 8600
        if(i.gt.iq) go to 429
        slp = (f2-f1)/(u2-u1)
        uest = u1 - f1/slp
        if(i.lt.iq) uest = (u1+u2)/2.d0
        u3 = uest
        q3 = u3
        q2 = u2
        q1 = u1
  429   fest=yfcn(uest)
        f3 = fest
        u3 = uest
        if(i.ge.iq) go to 439
        if(fest*f1 .lt. 0.d0) f2=fest
        if(fest*f1 .lt. 0.d0) u2=uest
        if(fest*f1 .ge. 0.d0) f1=fest
        if(fest*f1 .ge. 0.d0) u1=uest
        go to 440
  439   continue
c
        aa = (u2-u1)*f3+(u1-u3)*f2+(u3-u2)*f1
        denom = (u3-u2)*(u2-u1)*(u2-u1)
        aa = aa/denom
        bb = (u2-u1)*(2.d0*u3-u2-u1)*f3
        bb = bb- (u3-u1)*(u3-u1)*f2+(u3-u2)*(u3-u2)*f1
        bb = bb/denom
        cc = (u3-u1)*f3/(u2-u1)
        q0q = dabs(bb*bb - 4.d0*aa*cc)
        q0q = dsqrt(q0q)
        den1 = bb+q0q
        den2 = bb-q0q
        qden1 = dabs(den1)
        qden2 = dabs(den2)
        if(qden2.ge.qden1) den1 = den2
        u1 = u2
        u2 = u3
        u3 = u2 - 2.d0*cc/den1
        f1 = f2
        f2 = f3
        tst = dabs(u3-u2)
        if(tst.lt.ep) mchk = 1000
        uest = u3
  440   continue
        if(mchk.ne.1000) go to 8500
        go to 8600
 8500   write(6,8501)
 8501   format(3x,'no convergence was achieved')
 8600   continue
        xest = uest
        yest = fest
        return
        end
***********************************************************************
c ********
c	C. Greene, 6-8-87   -- modified 11-9-87
c  The following is a main program which shows how to use subroutine seaton.
c
c       The form of the call is:
c
c	     call seaton(l,eryd,r,zion,f,fp,g,gp)
c
c	Here l=orbital angular momentum
c	     eryd = electron energy in rydbergs
c            r=radial distance
c 	     zion=nuclear charge
c	(CAUTION:  I have only tested my transformation program for zion=1,
c 		   and I suspect it must be modified for other cases.      )
c	     (f,g)=(regular,irregular) energy normalized Coulomb functions.
c	     (fp,gp)= (d/dr)(f,g).
c	     Check:  W(f,g)=f*gp-fp*g  should be 2/pi if the program works.
c
c *********
c
	subroutine dummy
	implicit real*8(a-h,o-z)
	data pi/3.14159265358979d0/
	l=0
	eryd=-1.d0/(1.5d0**2) + 1.d-10
	zion=1.d0
	write(6,*) 'e=',eryd
	do 100 ir=1,10
		r=ir*0.5d0
		call seaton(l,eryd,r,zion,f,fp,g,gp)
		wfgtst=f*gp-fp*g - 2.d0/pi
		write(6,*) r,g,wfgtst
  100	continue
	epos=-eryd
	write(6,*) 'e=',epos
	do 200 ir=1,10
		r=ir*0.5d0
		call seaton(l,epos,r,zion,f,fp,g,gp)
		wfgtst=f*gp-fp*g - 2.d0/pi
		write(6,*) r,g,wfgtst
  200 	continue
	stop
	end
c
      subroutine seaton(l,eryd,r,zion,f,fp,g,gp)
      implicit real*8(a-h,o-z)
      acc=1.d-10
      rl=l
      rr=r*zion
      eps=eryd/(zion**2)
      call coulfg(l,eps,rr,acc,f0,f0p,g0,g0p,k,ierr,actacc)
      if(.not.(eryd.lt.0))goto 23000
         ea=dabs(eryd)
         call ganda(a,gg,l,ea,zion,999)
         goto 23001
c     else
23000    continue
         gam=1.d0/dsqrt(eps)
         call gcmplx(a,gg,rl,gam)
23001 continue
      a5=dsqrt(dabs(a))
      f=a5*f0
      fp=a5*f0p
      g=(g0+gg*f0)/a5
      gp=(g0p+gg*f0p)/a5
c
c ** the next five lines changed on 1-22-88 by c.greene thanks to h. gao
c
	factor = dsqrt(zion)
	f=f/factor
	g=g/factor
	fp=fp*factor
	gp=gp*factor
c
      return
      end
c
c
      SUBROUTINE GANDA(A,G,L,E,ZION,NOPT)                               00000020
      IMPLICIT REAL*8(A-H,O-Z)                                          00000030
      DATA PI/3.1415926535897932D0/                                     00000040
      DPI = PI*2.D0                                                     00000050
C*** THIS PROGRAM RETURNS THE QDT TRANSFORMATION PARAMETERS A&G.        00000060
      IF(ZION.EQ.0) WRITE(6,90)                                         00000070
   90 FORMAT(1X,'***** ZION = 0')                                       00000080
      IF(ZION.EQ.0) RETURN                                              00000090
      E = DABS(E)                                                       00000100
      XNUE = ZION/DSQRT(E)                                              00000110
C*** EVALUATE A(K,L) FIRST.                                             00000120
      A = 1.D0                                                          00000130
      IF(L.EQ.0) GO TO 109                                              00000140
      DO 100 I = 1,L                                                    00000150
      A = A* (1.D0 -I*I*E/(ZION**2) )                                   00000160
  100 CONTINUE                                                          00000170
  109 continue
C*** GIVE WARNINGS IN CASE A < OR = 0 .                                 00000180
      IF(A.LE.0) WRITE(6,5555)                                          00000190
 5555 FORMAT(1X,'****** A < OR = 0')                                    00000200
      IF(NOPT.EQ.1) RETURN                                              00000210
C*** CHECK WHETHER XNUE = INTEGER.                                      00000220
      N = XNUE + 1.D-02                                                 00000230
      Z = XNUE - N                                                      00000240
      IF(Z.EQ.0) G = 0.D0                                               00000250
      IF(Z.EQ.0) RETURN                                                 00000260
C*** G(K,L) IS NOW EVALUATED USING THE DIGAMMA PROGRAM.                 00000270
      G = A*(DIGAM(L+1+XNUE) + DIGAM(-L+XNUE)                           00000280
     1  - 2.D0 * DLOG(XNUE) )/DPI                                       00000290
      RETURN                                                            00000300
      END                                                               00000310
      SUBROUTINE GCMPLX(B,G,RL,GAM)                                     00000040
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)                            00000050
CCCC  COMPLEX*16 ZQ1,ZQ2,ZQ3,ZSUM,ZG                                    00000060
      COMPLEX*16 CDLOG
      DIMENSION GG(2)                                                   00000070
ccc   EQUIVALENCE (GG(1),ZG)                                            00000080
      DATA PI/3.1415926535897932D0/                                     00000090
C                                                                       00000100
      ZQ1 = DCMPLX(RL+1.D0,GAM)                                         00000110
      ZQ2 = DCMPLX(-RL,GAM)                                             00000120
      ZQ3 = DCMPLX(0.D0,GAM)                                            00000130
      ZSUM = ZDIGAM(ZQ1) + ZDIGAM(ZQ2) + CDLOG(ZQ3)*(-2.D0)             00000140
      PROD = 1.D0                                                       00000150
      L = RL                                                            00000160
      DO 100 I = 1,L                                                    00000170
      IF(L.EQ.0) GO TO 100                                              00000180
      PROD = PROD*( 1.D0 + I*I/(GAM*GAM) )                              00000190
  100 CONTINUE                                                          00000200
      A1 = PROD                                                         00000210
      ZG = ZSUM*A1/(2.D0*PI)                                            00000220
ccc   G = GG(1)                                                         00000230
	g = dreal(zg)
      QQ = 1.D0 - DEXP(-2.D0*PI*GAM)                                    00000240
      B = A1 / QQ                                                       00000250
      RETURN                                                            00000260
      END                                                               00000270
c
c      function cdabs(z)
c      implicit real*8(a-h,o-y),complex*16(z)
c      cdabs = cabs(z)
c      return
c      end
ccc
c	function cdlog(z)
c	implicit complex*16(a-h,o-z)
c	cdlog=clog(z)
c	return
c	end
c
c      function cdexp(z)
c      implicit complex*16(a-h,o-z)
c      cdexp=cexp(z)
c      return
c      end
c
      FUNCTION DIGAM(ARGG)                                              00000030
      IMPLICIT REAL*8(A-H,O-Z)                                          00000040
      DIGAM = 0.D0                                                      00000050
      ARG5 =-0.57721566490153286D0                                      00000060
      EN = 1.D0                                                         00000070
      ARG  = ARGG                                                       00000080
      ARG2 = ARGG                                                       00000090
    1 IF(ARG2-40.D0) 2,3,3                                              00000100
    2 DIGAM = DIGAM - 1.D0/ARG                                          00000110
      ARG = 1.D0+ ARG                                                   00000120
      ARG2 = 1.D0 + ARG2                                                00000130
      GO TO 1                                                           00000140
    3 PSI = DLOG(ARG) - 1.D0/(2.D0*ARG) - 1.D0/(12.D0*ARG**2)           00000150
      PSI = PSI + 1.D0/(120.D0*ARG**4) - 1.D0/(252.D0*ARG**6)           00000160
      DIGAM = DIGAM + PSI                                               00000170
      RETURN                                                            00000180
      END                                                               00000190
c     COMPLEX FUNCTION ZDIGAM*16(ARG)                                   00000040
	function zdigam(arg)
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)                            00000050
      COMPLEX*16 ARG,ZDIGAM,cdlog                                       00000060
      COMPLEX*8 ARGG                                                    00000070
      REAL*8 INC                                                        00000080
      ZDIGAM = (0.D0,0.D0)                                              00000090
      ARG5 =-0.57721566490153286D0                                      00000100
      PI = 3.1415926535897932D0                                         00000110
      EN = 1.D0                                                         00000120
      ARGG = ARG                                                        00000130
      ARG2 = REAL(ARGG)                                                 00000140
      ARG3 = AIMAG(ARGG)                                                00000150
      IF(ARG3) 4,1,4                                                    00000160
    1 IF(ARG2-40.D0) 2,3,3                                              00000170
    2 ZDIGAM = ZDIGAM - 1.D0/ARG                                        00000180
      ARG = 1.D0+ ARG                                                   00000190
      ARG2 = 1.D0 + ARG2                                                00000200
      GO TO 1                                                           00000210
    3 PSI = CDLOG(ARG)-1.D0/(2.D0*ARG)-1.D0/(12.D0*ARG**2)              00000220
      PSI=PSI +1.D0/(120.D0*ARG**4)-1.D0/(252.D0*ARG**6)                00000230
      ZDIGAM = ZDIGAM + PSI                                             00000240
      GO TO 12                                                          00000250
    4 IF(ARG2) 5,7,6                                                    00000260
    5 ZDIGAM = ZDIGAM - 1.D0/ARG                                        00000270
      ARG = ARG + 1.D0                                                  00000280
      ARG2 = ARG2 + 1.D0                                                00000290
      GO TO 4                                                           00000300
    6 ARG = ARG - 1.D0                                                  00000310
      ARG2 = ARG2 - 1.D0                                                00000320
      ZDIGAM = ZDIGAM + 1.D0/ARG                                        00000330
      GO TO 4                                                           00000340
    7 Y = CDABS(ARG)                                                    00000350
      ARG7 = PI*Y                                                       00000360
      ARG4 = 0.5D0/Y + (PI/2.D0)/DTANH(ARG7)                            00000370
      IF(Y-20.D0) 8,10,10                                               00000380
    8 INC = Y*Y/(EN*(EN*EN+Y*Y))                                        00000390
      ARG5 = ARG5 + INC                                                 00000400
      IF(INC - 1.D-12) 11,11,9                                          00000410
    9 EN = EN + 1.D0                                                    00000420
      GO TO 8                                                           00000430
   10 ARG5 = 1.D0/(12.D0*Y**2) + 1.D0/(120.D0*Y**4)                     00000440
      ARG5 = ARG5 + 1.D0/(252.D0*Y**6)+ DLOG(Y)                         00000450
   11 ZDIGAM = DCMPLX(ARG5,ARG4) + ZDIGAM                               00000460
C     XQ1 = REAL(ZDIGAM)                                                00000470
C     XQ2 = AIMAG(ZDIGAM)                                               00000480
   12 continue                                                          00000490
	return
      END                                                               00000500
c***********************************************************************
c
c
      subroutine coulfg(ll,eps,rho,acc,f,fp,g,gp,k,ierr,actacc)
c
c  calculates coulomb functions f and g and their derivatives
c
c  input -
c        ll=angular momentum quantum number
c        eps=z-scaled energy in rydbergs
c        rho=z-scaled radial variable in atomic units
c        acc=accuracy required
c
c  output -
c        f=regular function
c        fp=derivative of f
c        g=irregular function
c        gp=derivative of g
c        k=number of terms needed in expansion
c        ierr=error code
c        actacc=accuracy actually achieved
c
c  convergence criterion -
c        value of wronskian converged to accuracy of 0.5*acc
c
c  error codes -
c        ierr=0, converged with actacc.lt.acc
c        ierr=1, converged with actacc.gt.acc
c        ierr=2, not converged with 101 terms in main summation
c
c  initialization
c
c
      implicit real*8 (a-h,o-z)
cd    delete previous card for single precision
c
      data r2pi,ps0/.159154943,-.154431330/
      ierr=0
      lp1=ll+1
      l2=2*ll
      l2p1=l2+1
      fl=ll
      flp1=lp1
      fl2p1=l2p1
      e2=.5*eps
      r2=2.*rho
      acc2=2.*acc
c
c     initialize fa=factorial(2*ll+1)
c     and ps=psi(2*ll+2)+psi(1)
c
      fa=1.
      ps=ps0
c
c
c  calculate alpha(n) and beta(n) and initialize s and sp
c  continue calculation of fa and ps
c
c     s and sp for n=0
      x3=-l2
      x2=l2p1
      x1=-2.*r2**(-lp1)
      sp=x3*x1
      x1=r2*x1
      s=x1
c
c     initialize for coefficients in recursion formulae
      p1=fl*e2
      p2=p1
      q1=-e2
c
c     initialize alpha and beta
      alp1=1.
      alp2=1.+p2
      bet1=0.
      bet2=q1
c
      if(ll.eq.0)goto 20
c
c     s and sp for n=1
      x3=x3+2.
      x2=x2-1.
      x1=x1/x2
      sp=sp+x3*x1
      x1=r2*x1
      s=s+x1
c
c     loop for n=2 to 2*ll
      do 10 n=2,l2
c
c     continue calculation of fa and psi
      fn=n
      fa=fn*fa
      ps=ps+1./fn
c
c     continue calculation of s and sp
      x3=x3+2.
      x2=x2-1.
      x1=x1/(x2*fn)
      sp=sp+x3*x1*alp2
      x1=r2*x1
      s=s+x1*alp2
c
c     compute coefficients in recursion formulae
      p1=p1-e2
      p2=p2+p1
      q1=q1-e2
c     now have p2=-n*(n-2*ll-1)*eps/4
c     and q1=-n*eps/2
c
c     new alpha and beta
      alp0=alp1
      alp1=alp2
      alp2=alp1+p2*alp0
      bet0=bet1
      bet1=bet2
   10 bet2=bet1+p2*bet0+q1*alp0
c
c     normalize s and sp, complete calculation of fa and ps
      s=s*fa
      sp=sp*fa
      fa=fl2p1*fa
      ps=ps+1./fl2p1
c
c     complete calculation of alpha and beta
      p1=p1-e2
      p2=p2+p1
      q1=q1-e2
      alp0=alp1
      alp1=alp2
      bet0=bet1
      bet1=bet2
      bet2=bet1+p2*bet0+q1*alp0
c
   20 continue
c     now have alp1=alpha(2*ll+1)
c     and bet1=beta(2*ll+1), bet2=beta(2*ll+2)
c
c     value of a=a(eps,ll)
      a=alp1
      a4=4.*a
      cl=2.*a*dlog(dabs(r2))
cd    for single precision replace dlog by alog and dabs by abs
      clp=2.*a/rho
c
c  calculate a(n) and d(n), f and fp and
c  complete calculation of s and sp
c
c     calculate a0,a1,d0,d1
      a0=(2.**lp1)/fa
      a1=-a0/flp1
      ps=2.*ps*a
      d0=(bet1-ps)*a0
      d1=(bet2-ps-(2.+1./flp1)*a)*a1
c
c     initialize f,fp, continue calculation of s,sp
c     - values for n=0
      fnplp1=flp1
      c1=rho**ll
      c1p=fnplp1*c1
      fp=c1p*a0
      sp=sp+c1p*d0
      c1=c1*rho
      f=c1*a0
      s=s+c1*d0
      w1=f*(clp*f+sp)-fp*s
c
c     - values for n=1
      fnplp1=fnplp1+1.
      c1p=fnplp1*c1
      fp=fp+c1p*a1
      sp=sp+c1p*d1
      c1=c1*rho
      f=f+c1*a1
      s=s+c1*d1
      w2=f*(clp*f+sp)-fp*s
      dw2=dabs(w2-w1)
cd    for single precision replace dabs by abs
c
c     initialize for coefficients in recursion formulae
      p1=-2.*flp1
      p2=p1
      q1=a4+2.*a*fl2p1
c
c     loop for n=2 to 100
      do 40 n=2,100
c
c     compute coefficients in recursion formulae
      p1=p1-2.
      p2=p2+p1
      q1=q1+a4
c     now have p2=-n*(n+2*ll+1)
c     and q1=2*a*(2*n+2*ll+1)
c
c     compute a2=a(n) and d2=d(n)
      a2=(2.*a1+eps*a0)/p2
      d2=(2.*d1+eps*d0+q1*a2)/p2
c
c     increment fp and sp
      fnplp1=fnplp1+1.
      fp=fp+a2*fnplp1*rho**(ll+2)
      sp=sp+d2*fnplp1*rho**(ll+2)
c
c     increment f and s
      f=f+a2*rho**(ll+3)
      s=s+d2*rho**(ll+3)
c
c     calculate wronskian
      w1=w2
      dw1=dw2
      w2=f*(clp*f+sp)-fp*s
      dw2=dabs(w2-w1)
cd    for single precision replace dabs by abs
c
c     convergence test
      k=n+1
      if(dw1.gt.acc2)goto 30
      if(dw2.gt.acc2)goto 30
      goto 50
c
c     new a0,a1,do,d1
   30 a0=a1*rho
      a1=a2*rho
      d0=d1*rho
      d1=d2*rho
c
   40 continue
c
c  not converged
c
      ierr=2
      actacc=dabs(0.25*w2-1.)
cd    for single precision replace dabs by abs
      goto 60
c
c  converged
c
   50 actacc=dabs(0.25*w2-1.)
cd    for single precision replace dabs by abs
      if(actacc.gt.acc)ierr=1
c
c  complete calculation of g and gp
c
   60 g=(s+cl*f)*r2pi
      gp=(sp+cl*fp+clp*f)*r2pi
c
      return
      end
