c****************************************************************************************************
c KMS analytic delta-function formulas, ported from
c Adiabatic-Scattering-BoundStates/DeltaBoundValidation/WriteDeltaCurves.f (itself copied from
c DeltaBound.f). Rewritten to evaluate at a SINGLE arbitrary R (scalar q1/qm/qn/R arguments)
c instead of a tabulated array indexed by irho, so they can be called directly at whatever R
c the R-matrix box quadrature needs -- no interpolation anywhere. Physical constants (a2=1,
c mass=1, mu=sqrt(1/3) for three equal unit masses) are the same fixed values WriteDeltaCurves.f
c itself hardcodes; this module is specific to that exact system, not a generic potential.
c****************************************************************************************************
      double precision function transeq(resold,n,R,mu)
      implicit none
      double precision resold,R,mu
      integer n
      double precision PIOVER6
      PIOVER6=0.52359877559829887307710723054658381403286156656252d0
      if(n.ne.0) transeq = resold*tan(resold*PIOVER6) + sqrt(mu*2.0d0)*R
      if(n.eq.0) transeq = resold*tanh(resold*PIOVER6) - sqrt(mu*2.0d0)*R
      return
      end
c****************************************************************************************************
c Newton-solve for channel l's q(R) (l=0 is the genuine two-body bound-pair channel; l=m-1
c for channel index m otherwise), warm-started from qguess -- exactly WriteDeltaCurves.f's own
c Calcq inner loop, just for one (l,R) pair at a time instead of a whole tabulated sweep.
c nEscalate counts how many times the 1d-14 convergence tolerance had to be relaxed
c (x10 each time) because Newton hadn't converged within 1000 iterations -- 0 means the
c root was accepted at the original 1d-14 tolerance. niter is the total Newton iteration
c count. Diagnostic-only outputs, added to test whether 1d-14 (near the machine-precision
c floor) is actually being hit routinely via the escalation fallback, per user request.
      subroutine SolveQ(l,R,mu,qguess,qout,nEscalate,niter)
      implicit none
      double precision, external :: transeq
      integer l,n,ncond,nEscalate,niter
      double precision R,mu,qguess,qout
      double precision f1,f2,slope,resold,resnew,tol,dq,condition
      dq=1d-12
      resold=qguess
      condition=1d-14
      ncond=1000
      n=0
      nEscalate=0
      f1=transeq(resold,l,R,mu)
      f2=transeq(resold+dq,l,R,mu)
      slope=(f2-f1)/dq
      resnew=resold-f1/slope
      tol=abs(resnew-resold)
      do while(tol.gt.condition)
         f1=transeq(resold,l,R,mu)
         f2=transeq(resold+dq,l,R,mu)
         slope=(f2-f1)/dq
         resnew=resold-f1/slope
         tol=abs(resnew-resold)
         n=n+1
         resold=resnew
         if(n.ge.ncond) then
            condition=condition*10
            ncond=ncond+1000
            nEscalate=nEscalate+1
         endif
      enddo
      qout=resnew
      niter=n
      end
c****************************************************************************************************
c Diagonal adiabatic potential U(m,R) = +q(m,R)^2/(2 mu R^2) for m>1, -q(1,R)^2/(2 mu R^2) for
c the bound-pair channel m=1 -- exact port of CalcCurves' own U(m,irho) assignment.
      double precision function EvalU(m,R,mu,qm)
      implicit none
      integer m
      double precision R,mu,qm
      if (m.eq.1) then
         EvalU = -(qm**2)/(2.0d0*mu*R**2)
      else
         EvalU = qm**2/(2.0d0*mu*R**2)
      endif
      end
c****************************************************************************************************
c Exact port of WriteDeltaCurves.f's CalcPmat (antisymmetric P-matrix element), rewritten for a
c single R: q(1,irho)->q1, q(m,irho)->qm, q(n,irho)->qn, rho(irho)->R.
      subroutine EvalPmatAt(m,n,q1,qm,qn,R,result)
      implicit none
      integer m,n
      double precision, external :: sec, sech
      double precision tempnormM,tempnormN,numer,denom,b,phi
      double precision q1,qm,qn,R,Pi,a2,result
      Pi = 3.14159265358979324d0
      a2 = 1.0d0
      b=sqrt(2.0d0)/3.d0**0.25d0
      phi=Pi/6.0d0
      if(m.eq.n) then 
         result =0.0d0
         goto 99
      endif 
      if((m.eq.1).and.(n.gt.1)) then
         result =              (12*Sqrt((a2**2*b**2*
     -        q1**2*qn**2)/
     -      ((a2**2*Pi*q1**2 + 
     -          b*R*
     -           (6*a2 - b*Pi*R))
     -         *(a2**2*Pi*qn**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          ))))/
     -  (q1**2 + qn**2)

      endif
      if((m.gt.1).and.(n.eq.1)) then
         result =          (-12*Sqrt((a2**2*b**2*
     -        q1**2*qm**2)/
     -      ((a2**2*Pi*q1**2 + 
     -          b*R*
     -           (6*a2 - b*Pi*R))
     -         *(a2**2*Pi*qm**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          ))))/
     -  (q1**2 + qm**2)

      endif
      
      if((m.gt.1).and.(n.gt.1)) then
         result =         (-12*Sqrt((a2**2*b**2*
     -        qm**2*qn**2)/
     -      ((a2**2*Pi*qm**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          )*
     -        (a2**2*Pi*qn**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          ))))/
     -  (qm**2 - qn**2)
         
      endif
      
 99   continue      
      end
c****************************************************************************************************
c Exact port of WriteDeltaCurves.f's CalcQmat (symmetric Q-matrix element), rewritten for a
c single R: q(1,irho)->q1, q(m,irho)->qm, q(n,irho)->qn, rho(irho)->R.
      subroutine EvalQmatAt(m,n,q1,qm,qn,R,result)
      implicit none
      integer n,m
      double precision, external :: sec, sech
      double precision temp1, temp2, numer, denom,E,b,phi
      double precision q1,qm,qn,R,Pi,a2,result
      Pi = 3.14159265358979324d0
      a2 = 1.0d0
      E=exp(1.0d0)
      b=sqrt(2.0d0)/3.d0**(0.25d0)
      phi=Pi/6.0d0
c     NEW EXPRESSIONS CALCULATED ON 7/15/13 NPM via KARTAVTSEV and MALYKH expressions.  
C     THE NOTATION, HOWEVER, HOLDS TO OUR USUAL CONVENTION, AND NOT THE CONVENTION OF KARTAVTSEV AND MALYKH (WHO REVERSE THE ROLES OF Q AND P).
c     NOTE THAT THIS QMAT IS THE SYMMETRIC Q, EQUAL TO KARTAVTSEV AND MALYKH'S P-MATRIX

c     FIRST CALCULATE THE DIAGONAL ELEMENTS
C     Q(1,1)
      if((m.eq.1).and.(n.eq.1)) then
         result =         (a2**2*b**2*
     -    (-4*a2**4*phi**4*
     -       q1**6 + 
     -      3*b**2*R**2*
     -       (a2 - b*phi*R)**2 - 
     -      2*b*phi*q1**2*
     -       R*
     -       (-a2 + b*phi*R)*
     -       (3*a2**2 + 
     -         6*a2*b*phi*R + 
     -         2*b**2*phi**2*R**2
     -         ) + 
     -      a2**2*phi**2*q1**4*
     -       (-9*a2**2 + 
     -         8*a2*b*phi*R + 
     -         8*b**2*phi**2*R**2
     -         )))/
     -  (12.d0*(a2**2*phi*q1**2 + 
     -       b*R*
     -        (a2 - b*phi*R))**4)

      endif

C     Q(n,n) with n > 1
      if((m.gt.1).and.(n.gt.1).and.(m.eq.n)) then
         result =         (a2**2*b**2*
     -    (4*a2**4*phi**4*qn**6 + 
     -      3*b**2*R**2*
     -       (a2 - b*phi*R)**2 + 
     -      2*b*phi*qn**2*
     -       R*
     -       (-a2 + b*phi*R)*
     -       (3*a2**2 + 
     -         6*a2*b*phi*R + 
     -         2*b**2*phi**2*R**2
     -         ) + 
     -      a2**2*phi**2*qn**4*
     -       (-9*a2**2 + 
     -         8*a2*b*phi*R + 
     -         8*b**2*phi**2*R**2
     -         )))/
     -  (12.*(a2**2*phi*qn**2 + 
     -       b*R*
     -        (-a2 + b*phi*R))**4
     -    )

      endif

c     Q(1,n) with n > 1
      if((m.eq.1).and.(n.gt.1)) then
         result =       (6*b*Sqrt((a2**2*b**2*
     -        q1**2*qn**2)/
     -      ((a2**2*Pi*q1**2 + 
     -          b*R*
     -           (6*a2 - b*Pi*R))
     -         *(a2**2*Pi*qn**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          )))*
     -    ((2*phi*
     -         (-a2 + b*phi*R)*
     -         (a2**2*q1**2 - 
     -           b**2*R**2))/
     -       (a2**2*phi*q1**2 + 
     -          b*R*
     -           (a2 - b*phi*R))
     -         **2 + 
     -      a2/
     -       (a2**2*phi*q1**2 + 
     -         b*R*
     -          (a2 - b*phi*R))
     -       - (4*a2*q1**2)/
     -       ((q1**2 + 
     -           qn**2)*
     -         (a2**2*phi*q1**2 + 
     -           b*R*
     -            (a2 - b*phi*R))
     -         ) - 
     -      (2*phi*
     -         (a2 - b*phi*R)*
     -         (a2**2*qn**2 + 
     -           b**2*R**2))/
     -       (a2**2*phi*qn**2 + 
     -          b*R*
     -           (-a2 + b*phi*R))
     -         **2 + 
     -      a2/
     -       (a2**2*phi*qn**2 + 
     -         b*R*
     -          (-a2 + b*phi*R))
     -       - (4*a2*qn**2)/
     -       ((q1**2 + 
     -           qn**2)*
     -         (a2**2*phi*qn**2 + 
     -           b*R*
     -            (-a2 + b*phi*R)
     -           ))))/
     -  (q1**2 + qn**2)

      endif

c     Q(m,1) with m > 1
      if((m.gt.1).and.(n.eq.1)) then
         result =         (-6*b*Sqrt((a2**2*b**2*
     -        q1**2*qm**2)/
     -      ((a2**2*Pi*q1**2 + 
     -          b*R*
     -           (6*a2 - b*Pi*R))
     -         *(a2**2*Pi*qm**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          )))*
     -    ((-2*phi*
     -         (-a2 + b*phi*R)*
     -         (a2**2*q1**2 - 
     -           b**2*R**2))/
     -       (a2**2*phi*q1**2 + 
     -          b*R*
     -           (a2 - b*phi*R))
     -         **2 - 
     -      a2/
     -       (a2**2*phi*q1**2 + 
     -         b*R*
     -          (a2 - b*phi*R))
     -       + (4*a2*q1**2)/
     -       ((q1**2 + 
     -           qm**2)*
     -         (a2**2*phi*q1**2 + 
     -           b*R*
     -            (a2 - b*phi*R))
     -         ) + 
     -      (2*phi*
     -         (a2 - b*phi*R)*
     -         (a2**2*qm**2 + 
     -           b**2*R**2))/
     -       (a2**2*phi*qm**2 + 
     -          b*R*
     -           (-a2 + b*phi*R))
     -         **2 - 
     -      a2/
     -       (a2**2*phi*qm**2 + 
     -         b*R*
     -          (-a2 + b*phi*R))
     -       + (4*a2*qm**2)/
     -       ((q1**2 + 
     -           qm**2)*
     -         (a2**2*phi*qm**2 + 
     -           b*R*
     -            (-a2 + b*phi*R)
     -           ))))/
     -  (q1**2 + qm**2)

      endif

c     Q(m,n) with m != n
      if((m.gt.1).and.(n.gt.1).and.(m.ne.n)) then

         result =       (-6*b*Sqrt((a2**2*b**2*
     -        qm**2*qn**2)/
     -      ((a2**2*Pi*qm**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          )*
     -        (a2**2*Pi*qn**2 + 
     -          b*R*
     -           (-6*a2 + b*Pi*R)
     -          )))*
     -    ((2*phi*(a2 - b*phi*R)*
     -         (a2**2*qm**2 + 
     -           b**2*R**2))/
     -       (a2**2*phi*qm**2 + 
     -          b*R*
     -           (-a2 + b*phi*R))
     -         **2 - 
     -      a2/
     -       (a2**2*phi*qm**2 + 
     -         b*R*
     -          (-a2 + b*phi*R))
     -       + (4*a2*qm**2)/
     -       ((qm**2 - 
     -           qn**2)*
     -         (a2**2*phi*qm**2 + 
     -           b*R*
     -            (-a2 + b*phi*R)
     -           )) - 
     -      (2*phi*
     -         (a2 - b*phi*R)*
     -         (a2**2*qn**2 + 
     -           b**2*R**2))/
     -       (a2**2*phi*qn**2 + 
     -          b*R*
     -           (-a2 + b*phi*R))
     -         **2 + 
     -      a2/
     -       (a2**2*phi*qn**2 + 
     -         b*R*
     -          (-a2 + b*phi*R))
     -       + (4*a2*qn**2)/
     -       ((qm**2 - 
     -           qn**2)*
     -         (a2**2*phi*qn**2 + 
     -           b*R*
     -            (-a2 + b*phi*R)
     -           ))))/
     -  (qm**2 - qn**2)

      endif


      end
c****************************************************************************************************
c Analytic hyperangular channel-function overlap <Phi_m(R_i)|Phi_n(R_j)>, for the equal-mass
c three-boson delta-function problem of Mehta, Esry & Greene, PRA 76, 022711 (2007) ("MEG").
c Uses the SAME q_m(R)/qn(R) already solved by SolveQ -- no new transcendental solve needed,
c since only the q/kappa values (not R itself) enter the formula.
c
c Derivation: in the reduced wedge 0<=phi<=pi/6 (MEG Appendix A: the full 2*pi domain is built
c from this wedge by the boson symmetry group S=1+P12+P23+P31+P12P23+P12P31, a 12-fold
c reduction), the even-parity boundary condition Phi'(0)=0 kills the odd term in both MEG's
c general solutions (Eqs. 16, 20), leaving
c   channel 1        (bound-pair, MEG Eq. 16): Phi(R;phi) = B1(R)*cosh(q1(R)*phi)
c   channel m>1  (continuum, MEG Eq. 20):  Phi(R;phi) = Bm(R)*cos(qm(R)*phi)
c with q1/qm exactly what SolveQ(0,...)/SolveQ(m-1,...) return (verified: SolveQ's transeq
c is exactly MEG Eqs. 19/22 in this code's own sign convention). Unit-normalizing each channel
c at its own R (integral over the FULL domain = 12x the wedge integral) gives, via the
c elementary integrals int cosh^2 = phi/2+sinh(2*q*phi)/(4*q) and int cos^2 = phi/2+sin(2*q*phi)/(4*q):
c   B1(R) = 1/sqrt(Pi + 3*sinh(q1*Pi/3)/q1)
c   Bm(R) = 1/sqrt(Pi + 3*sin(qm*Pi/3)/qm)
c The cross-R overlap is then 12*Bm(Ri)*Bn(Rj) times one of three elementary wedge integrals
c (cosh*cosh, cos*cos, or the mixed cosh*cos when exactly one of m,n is channel 1), each a
c standard closed-form product integral -- see the three branches below. Self-consistency
c checks (verified symbolically by hand before coding): m=n at the same R collapses to exactly
c 1 in all three cases (unit normalization), and the a=b limit of the same-type branches is
c continuous (no removable-singularity issue) with the a!=b formula.
      subroutine EvalOverlapAt(m,qm,n,qn,result)
      implicit none
      integer m,n
      double precision qm,qn,result
      double precision Pi,phi6,a,b,intgrl,Bm,Bn
      Pi = 3.14159265358979324d0
      phi6 = Pi/6.0d0

      if (m.eq.1) then
         Bm = 1.0d0/sqrt(Pi + 3.0d0*sinh(qm*Pi/3.0d0)/qm)
      else
         Bm = 1.0d0/sqrt(Pi + 3.0d0*sin(qm*Pi/3.0d0)/qm)
      endif
      if (n.eq.1) then
         Bn = 1.0d0/sqrt(Pi + 3.0d0*sinh(qn*Pi/3.0d0)/qn)
      else
         Bn = 1.0d0/sqrt(Pi + 3.0d0*sin(qn*Pi/3.0d0)/qn)
      endif

      if ((m.eq.1).and.(n.eq.1)) then
c        cosh(a*phi)*cosh(b*phi), a=qm, b=qn
         a = qm
         b = qn
         if (abs(a-b).lt.1d-10) then
            intgrl = phi6/2.0d0 + sinh(a*Pi/3.0d0)/(4.0d0*a)
         else
            intgrl = 0.5d0*(sinh((a+b)*phi6)/(a+b)
     -           + sinh((a-b)*phi6)/(a-b))
         endif
      elseif ((m.gt.1).and.(n.gt.1)) then
c        cos(a*phi)*cos(b*phi), a=qm, b=qn
         a = qm
         b = qn
         if (abs(a-b).lt.1d-10) then
            intgrl = phi6/2.0d0 + sin(a*Pi/3.0d0)/(4.0d0*a)
         else
            intgrl = 0.5d0*(sin((a-b)*phi6)/(a-b)
     -           + sin((a+b)*phi6)/(a+b))
         endif
      elseif (m.eq.1) then
c        cosh(a*phi)*cos(b*phi), a=qm (channel 1), b=qn (channel n>1)
         a = qm
         b = qn
         intgrl = (a*sinh(a*phi6)*cos(b*phi6)
     -        + b*cosh(a*phi6)*sin(b*phi6))/(a**2+b**2)
      else
c        cos(a*phi)*cosh(b*phi), a=qm (channel m>1), b=qn (channel 1)
         a = qn
         b = qm
         intgrl = (a*sinh(a*phi6)*cos(b*phi6)
     -        + b*cosh(a*phi6)*sin(b*phi6))/(a**2+b**2)
      endif

      result = 12.0d0*Bm*Bn*intgrl
      end
