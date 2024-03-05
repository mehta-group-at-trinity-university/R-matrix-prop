c234567890
      program CoupledAdiabaticBoxAvg
      
      integer LegPoints,xNumPoints
      integer Order,ishift,NumRmatch,NumNewSectors
      integer NumInpChannels,NumChannels,NumOpenChannels
      integer NumOpenL, NumOpenR, NumTotStates,NumSectors
      integer NumDataPoints,NumBoxes,Loops_E
      integer TotalAngMomentum
      integer, allocatable :: MEMap(:,:),NumberL(:),NumberR(:),NumSectorsBox(:)
      double precision, allocatable :: Thresholds(:),leff(:),xPoints(:)

      double precision rStart,rEnd,xStart,xEnd,xStep,rmsError,rmsErrorP
      double precision rmsErrorMax,rmsErrorPMax,WhereStart,BoxAvg,CrossTot
      double precision, allocatable :: xSector(:),logderiv(:),xSectorTot(:)
      double precision, allocatable :: xData(:),xSave(:)
      double precision, allocatable :: VQmat(:,:,:),Pmat(:,:,:)
      double precision, allocatable :: ub(:,:),up(:,:)
      double precision, allocatable :: PFit(:,:,:),VQFit(:,:,:)
      double precision, allocatable :: PExponents(:,:,:),VQExponents(:,:,:)

      integer i,j,k,ix,i1,i2,HalfBandWidth,lwork,info
      integer NumClosedChannels
      integer iEnergy,NumEnergySteps
      integer xDim,xOpenDim,xClosedDim,LeadDim
      integer ClosedMatrixDim,OpenMatrixDim
      integer, allocatable :: xBounds(:)
      integer, allocatable :: ipiv(:)
      double precision Tol,Einc
      double precision TotalMemory,ddot,Mass,kelvinPERau
      double precision m1,m2,m3
      double precision Energy,Einput,RMatch,RMatch1,f,fp,g,gp,Pi,Wronskian
      double precision, allocatable :: work(:)
      double precision, allocatable :: xLeg(:),wLeg(:)
      double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:)
      double precision, allocatable :: uTemp(:),uxTemp(:)
      double precision, allocatable :: RData(:),SData(:,:,:,:),HData(:,:,:,:,:)
      double precision, allocatable :: Scc(:,:),Sco(:,:),Soc(:,:),Soo(:,:)
      double precision, allocatable :: Hcc(:,:),Hco(:,:),Hoc(:,:),Hoo(:,:)
      double precision, allocatable :: Gcc(:,:),Gco(:,:),Goc(:,:),Goo(:,:)
      double precision, allocatable :: b(:)

      double precision, allocatable :: psi_in(:,:),solution(:,:),deriv(:)

      double precision, allocatable :: Tmatrix(:,:)
      double precision, allocatable :: TmatrixAvg(:,:)


      double precision, allocatable :: CrossSections(:,:)
      double precision, allocatable :: BoxAvgPartial(:,:)
      double complex, allocatable :: Smatrix(:,:)
      character*64 LegendreFile,SFile,HFile,PFile,QFile,FitFile

      double precision muK
      double precision sec,time,secp,timep
      double precision step

      double precision EinputAux
      double precision Ri,Rf,Xi,Xf,StepX,R(1:100000000)
      double precision Emax,wlenght,rstep,XP
      integer NPoints,npwave,NPointsW,NPS
      double precision URmin,wlenghti,rstepi,dumb
      integer ncount,ndumb
      character*64 GridName

C     TIME INICIALIZATION
      call cpu_time(time)
      SEC = time
      secp = time


      Pi = 3.1415926535897932385d0
      kelvinPERau=315774.65d0
      muK = 3.166815355d-12     !jpdincao ! muK=1.d-6/kelvinPERau

c     read in Gauss-Legendre info
      read(5,*)
      read(5,1002) LegendreFile
      read(5,*)

      read(5,*)
      read(5,*) LegPoints, WhereStart
      read(5,*)

      read(5,*)
      read(5,*) TotalAngMomentum,Einput,Einc,Loops_E,ithresh
      read(5,*)

      read(5,*)
      read(5,*) m1,m2,m3
      read(5,*)

      Mass = dsqrt(m1*m2*m3/(m1+m2+m3))

      read(5,*)
      read(5,*) NumChannels,NumOpenChannels,NumBoxes
      read(5,*)

      read(5,*)
      allocate(NumberL(NumBoxes),NumberR(NumBoxes))
      do i = 1,NumBoxes
         read(5,*) NumberL(i),NumberR(i)
      enddo 
      read(5,*)
      read(5,*)
      read(5,*)NumRmatch,NumNewSectors
c     allocate(leff(NumberL(NumBoxes)),Thresholds(NumberR(NumBoxes)))
      allocate(leff(NumberL(NumBoxes)),Thresholds(NumberL(NumBoxes)))
      read(5,*)
      read(5,*)  
      do i = 1,NumberR(NumBoxes)
         read(5,*)leff(i),Thresholds(i)
      enddo

      read(5,*) 
      read(5,*) 
      read(5,1002) SFile
      read(5,*)
      read(5,*)
      read(5,1002) HFile
      read(5,*)
      read(5,*)
      read(5,1002) PFile
      read(5,*)
      read(5,*)
      read(5,1002) QFile  
      read(5,*)
      read(5,*)
      read(5,1002) FitFile

c     Thresholds(i) = 2B Energies from the Extrapolation Code
      open(18,file=FitFile,status='old')
      read(18,*)NumTotStates,NumDataPoints
      close(18)
      open(19,file='FitLeff.data',status='old')
      read(19,*)
      read(19,*)
      read(19,*)
      do i=1,NumTotStates
         read(19,*)ndumb,Thresholds(i),dumb,dumb,dumb
         if (Thresholds(i).gt.0.d0) then
            Thresholds(i) = 0.d0
         endif
      enddo
      close(19)

      Energy = Einput + Thresholds(ithresh)

c     Checking WhereStart
      open(18,file=FitFile,status='old')
      read(18,*)NumTotStates,NumDataPoints
      close(18)
      allocate(VQmat(NumDataPoints,Numchannels,Numchannels))
      allocate(xdata(NumDataPoints))
      open(300,file=QFile,status='old')
      do i=1,NumDataPoints
         read(300,*)xdata(i)
         do j=1,Numchannels
            read(300,301)(VQmat(i,j,k),k=1,Numchannels)
         enddo
         do j=Numchannels+1,NumTotStates
            read(300,301)
         enddo
      enddo
      close(300)
 301  format(100(e18.12,1x))
c     301  format(100(e20.12,1x))

      if (xdata(NumDataPoints).ne.WhereStart) then
         WhereStart = xdata(NumDataPoints)
      endif

      deallocate(VQmat,xdata)

      read(5,*)
      read(5,*)
      read(5,*)Ri,Rf
      read(5,*)
      read(5,*)
      read(5,*)NPoints
      read(5,*)
      read(5,*)
      read(5,*)npwave

c     RADIAL GRID

      Xi = Ri**(1.0d0/3.0d0)
      Xf = WhereStart**(1.0d0/3.0d0)
      StepX = (Xf-Xi)/dfloat(NPoints-1)
      do iR = 1,NPoints
         R(iR) = (Xi+(iR-1)*StepX)**3
      enddo

      Pi = dacos(-1.d0)
      Emax = 3.166815355d-08    !(10^4 muK - Max income Energy)
      wlenght = 2.d0*Pi/dsqrt(2.d0*Mass*(Emax-Thresholds(1)))
      rstep = wlenght/dfloat(npwave) 

      NPointsW = int((Rf-WhereStart)/rstep)+1
      do iR = 1,NPointsW
         R(iR+NPoints) = WhereStart+(iR)*rstep
      enddo

      if (R(NPoints)-R(Npoints-1).gt.R(NPoints+1)-R(Npoints)) then
         write(*,*)'Step size bigger than the extrapolation step'
         write(*,*)'... Increase NPoints ...'
         write(*,*)'Program STOPED'
         stop
      endif   

      NumSectorsTot = NPoints+NPointsW

      allocate(xSectorTot(NumSectorsTot))
      do iR = 1,NumSectorsTot
         xSectorTot(iR) = R(iR)
      enddo

      allocate(NumSectorsBox(NumBoxes))
      XP = 1.d0/(NumBoxes-1)
      NPS = int(XP*NumSectorsTot)
      do k = 1,NumBoxes
         NumSectorsBox(k) = NPS
         if (k.eq.(NumBoxes-1)) NumSectorsBox(k) = 0.5*NPS
         if (k.eq.(NumBoxes))   NumSectorsBox(k) = 2
      enddo


c     GridName = 'grid.used.data'
c     open(unit=17,file=GridName)
c     write(17,*)
c     write(17,*)'----------------------'      
c     write(17,1705)URmin
c     write(17,1700)NumSectorsTot
c     write(17,1704)NPoints,NPointsW
c     write(17,1701)wlenghti,wlenght
c     write(17,1703)rstepi,rstep
c     write(17,1702)npwave
c     write(17,*)'----------------------'      
c     do k = 1,NumBoxes
c     write(17,171)k,NumSectorsBox(k)   
c     enddo
c     write(17,*)'----------------------'      
c     do iR = 1,NumSectorsTot
c     write(17,*)iR,xSectorTot(iR)
c     enddo
c     write(17,*)'----------------------'      
c     close(17)
c     
c     1700 format(1x,'NumSectorsTot =',I8)
c     1701 format(1x,'Wave Lenght =',f14.8',',f14.8)
c     1702 format(1x,'Num. Points per Wave Lenght =',I3)
c     1703 format(1x,'Step =',f14.8',',f14.8)
c     1704 format(1x,'NPoints =',I6,',',I6)
c     1705 format(1x,'Potential Minimum =',f14.8)
c     171  format(1x,'Box(',I3,')=',I8)

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      allocate(xLeg(LegPoints),wLeg(LegPoints))
      call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)
      
      goto 54 
      
 54   allocate(VQExponents(Numchannels,Numchannels,3))
      allocate(VQFit(Numchannels,Numchannels,3))
      allocate(PExponents(Numchannels,Numchannels,3))
      allocate(PFit(Numchannels,Numchannels,3))
      
      open(18,file=FitFile,status='old')
      read(18,*)NumTotStates,NumDataPoints
      rmsErrorMax=0.d0
      rmsErrorPMax=0.d0
      do j=1,NumChannels
         do k=j,NumChannels
            read(18,*) jj,kk,(VQExponents(jj,kk,ii),ii=1,3)
            read(18,*) (VQFit(jj,kk,ii),ii=1,3),rmsError
            IF(j.NE.k) THEN
               read(18,*) jj,kk,(PExponents(jj,kk,ii),ii=1,3)
               read(18,*) (PFit(jj,kk,ii),ii=1,3),rmsErrorP
            ENDIF
            IF(jj.NE.j .OR. kk.NE.k) THEN
               write(6,*) 'index mismatch!',j,jj,k,kk
            ENDIF
            DO i=1,3
               ! enforce symmetry in VQ and antisymmetry in P
               VQFit(k,j,i)=VQFit(j,k,i)
               VQExponents(k,j,i)=VQExponents(j,k,i)
               PFit(k,j,i)=-PFit(j,k,i)
               PExponents(k,j,i)=PExponents(j,k,i)
               IF(j.EQ.k) THEN
                  PFit(k,k,i)=0.d0
                  PExponents(k,k,i)=1
               ENDIF
            ENDDO
            IF(rmsError .GT. rmsErrorMax) THEN
               rmsErrorMax=dabs(rmsError)
            ENDIF
            IF(rmsErrorP.GT. rmsErrorPMax) THEN
               rmsErrorPMax=dabs(rmsErrorP)
            ENDIF
         ENDDO
         do k=Numchannels+1,NumTotStates
            read(18,*)
            read(18,*)
            IF(j.NE.k) THEN
               read(18,*)
               read(18,*)
            ENDIF
         enddo
      ENDDO
      close(18)

c     XXX
c     do j = 1,NumChannels
c     do k = j+1,NumChannels
c     j=1; k=3;
c     do i = 1,3   
c     VQFit(j,k,i) = 0.d0
c     VQFit(k,j,i) = 0.d0
c     PFit(j,k,i)  = 0.d0
c     PFit(k,j,i)  = 0.d0
c     c       VQFit(j,k,i) = VQFit(k,j,i) 
c     c        PFit(j,k,i) = -PFit(k,j,i) 
c     enddo
c     enddo
c     enddo

      allocate(xData(NumDataPoints))
      allocate(VQmat(NumDataPoints,Numchannels,Numchannels))
      allocate(Pmat(NumDataPoints,Numchannels,Numchannels))
      
      call potential(NumChannels,NumDataPoints,NumTotStates,
     >     PFile,QFile,xdata,VQmat,Pmat) ! read in the P, VQ matrix data to be interpolated later 



      do ie=1,Loops_E

c     if (energy.gt.Thresholds(ithresh+1)) goto 1234 
c     This condition is valid for the relaxation problem. 
c     It means that for energy.gt.Thresholds(ithresh+1) a 
c     new channel will be open. So, re-edit the input file 
c     to include those channels.

         ishift = 0
         do ibox = 1,NumBoxes-1

            NumSectors = NumSectorsBox(ibox)
            NumOpenL = NumberL(ibox)
            NumOpenR = NumberR(ibox)
            NumOpen = NumOpenL + NumOpenR
            MatrixDim = Numchannels*NumSectors*4
            kl = 6*Numchannels - 1
            ku = kl
            iband = 3*kl + 1
            allocate(xSector(NumSectors+1))
            
            if (ibox .gt. 1) ishift =  ishift + NumSectorsBox(ibox-1)
            do is=1,NumSectors+1
               xSector(is) = xSectorTot(ishift+is)
            enddo
            Rstart = xSector(1)
            Rmatch = xSector(NumSectors+1)

!     write(6,*)'Adiabatic Box:         ',ibox
!     write(6,*)'Box Boundaries:        ',Rstart,Rmatch
!     write(6,*)'NumOpenChannels Left:  ',NumOpenL
!     write(6,*)'NumOpenChannels Right: ',NumOpenR
!     write(6,*)'OpenMatrixDim:         ',NumOpen
!     write(6,*)'ClosedMatrixDim:       ',MatrixDim
!     write(6,*)'HalfBandwidth:         ',kl
!     write(6,*)

            allocate(MEMap(Numchannels*NumSectors,6))
            call GenMEMap(Numchannels,Numsectors,MEMap)
            
            allocate(Gcc(iband,MatrixDim))
            
!     write(6,*)'starting CalcClosed'
            call CalcClosed(Numchannels,NumSectors,MatrixDim,kl,iband,MEMap,Mass,energy,
     >           xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,PFit,VQFit,
     >           WhereStart,PExponents,VQExponents,Gcc)
!     write(6,*)'done'
!     write(6,*)
            
            allocate(Gco(MatrixDim,NumOpen),Goo(NumOpen,NumOpen))
            
!     write(6,*)'starting CalcOpen'
            call CalcOpen(NumOpenL,NumOpenR,Numchannels,NumSectors,MatrixDim,MEMap,Mass,
     >           energy,xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,
     >           PFit,VQFit,WhereStart,PExponents,VQExponents,Gco,Goo)
!     write(6,*)'done'
!     write(6,*)
            
            allocate(Goc(NumOpen,MatrixDim))
            
            do i=1,NumOpen
               do j=1,MatrixDim
                  Goc(i,j) = Gco(j,i)
               enddo
            enddo
            
c     Lapack Banded Linear Equation Solver
            
            allocate(ipiv(MatrixDim))
            
            call dgbsv(MatrixDim,kl,ku,NumOpen,Gcc,iband,ipiv,Gco,
     >           MatrixDim,info)
            
            call dgemm('N','N',NumOpen,NumOpen,MatrixDim,-1.0d0,Goc,NumOpen,
     >           Gco,MatrixDim,1.0d0,Goo,NumOpen)
            
            lwork = 3*NumOpen
            allocate(work(lwork)) 
            allocate(logderiv(NumOpen))
            
            call dsyev('V','U',NumOpen,Goo,NumOpen,logderiv,work,lwork,info)
            
c     write(406,*)'log derivs in'
c     write(406,407)b(1:NumOpenL)
c     write(406,*)
c     write(406,*)'solution'
c     do i=1,NumOpenL
c     write(406,407)psi_in(i,1:NumOpenL)
c     enddo 
c     write(406,*)
c     write(406,*)'log derivs out'
c     write(406,407)logderiv(1:NumOpen)
c     write(406,*)
c     write(406,*)'solution'
c     do i=1,NumOpen
c     write(406,407)Goo(i,1:NumOpen)
c     enddo

            deallocate(xSector,work,ipiv)
            deallocate(MEMap,Gcc,Gco,Goc)
            
c     starting matching section
            if (ibox .eq. 1)then
               allocate(psi_in(NumOpenR,NumOpenR),b(NumOpenR))
               psi_in = Goo
               b = logderiv
               deallocate(Goo,logderiv)
            else
!     write(6,*)'starting matching section'

               allocate(solution(NumOpenR,NumOpenR))
               allocate(deriv(NumOpenR))
               call BoxMatching(NumOpenL,NumOpenR,psi_in,b,Goo,logderiv,
     >              solution,deriv)

               deallocate(psi_in,b,Goo,logderiv)
               allocate(psi_in(NumOpenR,NumOpenR),b(NumOpenR))
               psi_in = solution
               b = deriv
               deallocate(solution,deriv)
            endif

         enddo

         allocate(BoxAvgPartial(NumberR(NumBoxes),NumberR(NumBoxes)))
         allocate(TmatrixAvg(NumOpenR,NumOpenR))
         BoxAvg = 0.0d0
         BoxAvgPartial = 0.0d0
         TmatrixAvg = 0.0d0

         do ir=1,NumRmatch

            NumSectors = NumSectorsBox(NumBoxes) + (ir-1)*NumNewSectors
            allocate(xSector(NumSectors+1))
            if(ir .eq. 1 .and. NumBoxes .eq. 1) ishift = 0 
            if(ir .eq. 1 .and. NumBoxes .gt. 1) ishift = ishift + NumSectorsBox(NumBoxes-1)
            do is=1,NumSectors+1
               xSector(is) = xSectorTot(ishift+is)
            enddo 
            Rstart = xSector(1)
c     xSector(NumSectors+1) = xSector(NumSectors+1)+100.d0
            Rmatch = xSector(NumSectors+1)

c     c     jpdincao 08-09-04
c     NumSectors = 100
c     allocate(xSector(NumSectors+1))
c     do is=1,NumSectors+1
c     xSector(is) = xSectorTot(NumSectorsTot) +0.1d0*(is-1) +1.d0*(ir-1)
c     enddo 
c     Rstart = xSector(1)
c     Rmatch = xSector(NumSectors+1)


c     write(*,*)ir,NumSectors,NumBoxes,Rmatch

            NumOpenL = NumberL(NumBoxes)
            NumOpenR = NumberR(NumBoxes)
            NumOpen = NumOpenL + NumOpenR
            MatrixDim = Numchannels*NumSectors*4
            kl = 6*Numchannels - 1
            ku = kl
            iband = 3*kl + 1

!     write(6,*)'Adiabatic Box:         ',NumBoxes
!     write(6,*)'Box Boundaries:        ',Rstart,Rmatch
!     write(6,*)'NumOpenChannels Left:  ',NumOpenL
!     write(6,*)'NumOpenChannels Right: ',NumOpenR
!     write(6,*)'OpenMatrixDim:         ',NumOpen
!     write(6,*)'ClosedMatrixDim:       ',MatrixDim
!     write(6,*)'HalfBandwidth:         ',kl
!     write(6,*) 

            
            allocate(MEMap(Numchannels*NumSectors,6))
            call GenMEMap(Numchannels,Numsectors,MEMap)
            
            allocate(Gcc(iband,MatrixDim))
            
!     write(6,*)'starting CalcClosed'
            call CalcClosed(Numchannels,NumSectors,MatrixDim,kl,iband,MEMap,Mass,energy,
     >           xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,PFit,VQFit,
     >           WhereStart,PExponents,VQExponents,Gcc)
!     write(6,*)'done'
!     write(6,*)
            
            allocate(Gco(MatrixDim,NumOpen),Goo(NumOpen,NumOpen))
            
!     write(6,*)'starting CalcOpen'
            call CalcOpen(NumOpenL,NumOpenR,Numchannels,NumSectors,MatrixDim,MEMap,Mass,
     >           energy,xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,
     >           PFit,VQFit,WhereStart,PExponents,VQExponents,Gco,Goo)
!     write(6,*)'done'
!     write(6,*)
            
            allocate(Goc(NumOpen,MatrixDim))
            
            do i=1,NumOpen
               do j=1,MatrixDim
                  Goc(i,j) = Gco(j,i)
               enddo
            enddo
            
c     Lapack Banded Linear Equation Solver
            
            allocate(ipiv(MatrixDim))
            
            call dgbsv(MatrixDim,kl,ku,NumOpen,Gcc,iband,ipiv,Gco,
     >           MatrixDim,info)
            
            call dgemm('N','N',NumOpen,NumOpen,MatrixDim,-1.0d0,Goc,NumOpen,
     >           Gco,MatrixDim,1.0d0,Goo,NumOpen)
            
            lwork = 3*NumOpen
            allocate(work(lwork))
            allocate(logderiv(NumOpen))
            
            call dsyev('V','U',NumOpen,Goo,NumOpen,logderiv,work,lwork,info)
            
            deallocate(xSector,work,ipiv)
            deallocate(MEMap,Gcc,Gco,Goc)
            
c     starting matching section
!     write(6,*)'starting matching section'
            
            allocate(solution(NumOpenR,NumOpenR))
            allocate(deriv(NumOpenR))
            call BoxMatching(NumOpenL,NumOpenR,psi_in,b,Goo,logderiv,
     >           solution,deriv)

            deallocate(Goo,logderiv)

c     calculate S-matrix and scattering cross sections
c     CrossSections is actually the K_3 (or VR) rate constant

            allocate(Smatrix(NumOpenR,NumOpenR),Tmatrix(NumOpenR,NumOpenR))
            allocate(CrossSections(NumOpenR,NumOpenR))

            call CalcSmatrix(NumOpenR,leff,Thresholds,TotalAngMomentum,Mass,Energy,RMatch,
     >           solution,deriv,Smatrix,Tmatrix,CrossSections)

            do i = 1,NumOpenR
               do j = 1,NumOpenR
                  BoxAvgPartial(i,j) = BoxAvgPartial(i,j) + CrossSections(i,j)
                  TmatrixAvg(i,j) = TmatrixAvg(i,j) + Tmatrix(i,j)
               enddo
            enddo

c     Recombination Rate
            Crosstot = 0.0d0
            do i=1,ithresh-1
               do j=ithresh,NumOpenR
                  Crosstot = Crosstot + CrossSections(i,j)
               enddo
            enddo

c     c     Relaxation Rate 
c     Crosstot = 0.0d0
c     do i=1,ithresh-1  ! ithresh must indicate the entrance channel
c     do j=ithresh,ithresh 
c     Crosstot = Crosstot + CrossSections(i,j)
c     enddo
c     enddo


c     do ic=1,ithresh-1
c     Crosstot = Crosstot + CrossSections(ic,ithresh)
c     enddo
            Boxavg = Boxavg + Crosstot

!     write(6,*)'Rate Constant K_3(cm^6/sec)'
!     write(6,*) Crosstot

c     write(29+Loops_E,28)Rmatch,Crosstot,CrossSections(1:ithresh-1,ithresh)
c     write(*,20)(energy-Thresholds(ithresh))/muK,1.d0*ir,((CrossSections(i,j),j=1,NumOpenR),i=1,NumOpenR)

            deallocate (Smatrix,Tmatrix,CrossSections)
            deallocate(solution,deriv)

         enddo

!     write(6,*)'Done'
         do i = 1,NumberR(NumBoxes)
            do j = 1,NumberR(NumBoxes)
               BoxAvgPartial(i,j) = BoxAvgPartial(i,j)/dfloat(NumRmatch)
               TmatrixAvg(i,j) = TmatrixAvg(i,j)/dfloat(NumRmatch)
            enddo
         enddo

         do ic = 1,ithresh-1
c     Partial Vibrational Relaxation Rate
            write(9+ic,20)(energy-Thresholds(ithresh))/muK,(BoxAvgPartial(ic,j),j=ithresh,ithresh) 
c     Partial Recombination Rate
c     write(9+ic,20)(energy-Thresholds(ithresh))/muK,BoxAvgPartial(ic,ithresh:NumberR(NumBoxes))
         enddo

         BoxAvg = BoxAvg/float(NumRmatch)
!     write(6,*)
!     write(6,*)'Box average results'
!     write(6,*)BoxAvg
c     write(25,*)energy/3.165226467522434d-12,BoxAvg
         write(25,*)(energy-Thresholds(ithresh))/muK,BoxAvg

c     Tmatrix outPut      
         write(1000,20)(energy-Thresholds(ithresh))/muK,((TmatrixAvg(i,j),j=1,NumOpenR),i=1,NumOpenR)

         call cpu_time(timep)
c     Total and Partial Vibrational Relaxation Rate
c     write(*,20)(energy-Thresholds(ithresh))/muK,
c     .           ((BoxAvgPartial(i,j),j=ithresh,ithresh),i=1,ithresh-1),BoxAvg,timep-secp
         write(*,20)(energy-Thresholds(ithresh))/muK,BoxAvg,
     .        ((BoxAvgPartial(i,j),j=ithresh,NumOpenR),i=1,ithresh-1),timep-secp


         call cpu_time(timep)
         secp = timep 

         deallocate(BoxAvgPartial,TmatrixAvg)
         deallocate(psi_in,b)

c     energy = Einput + float(ie)*Einc + Thresholds(ithresh)
c     energy = Einput*(10.d0**(float(ie)*0.1d0)) 

         EinputAux = Einput*(10.d0**(float(ie)*0.1d0)) 
         energy = (EinputAux+Thresholds(ithresh))


      enddo

 1234 continue

      deallocate(xData,VQmat,Pmat,Pfit,VQfit)
      deallocate(VQexponents,Pexponents)
      deallocate(xLeg,wLeg)

 20   format(1P,100e16.8)
 28   format(f10.2,1P,100e14.5)
 407  format(10(e16.8,1x))
 1002 format(a64)

      call cpu_time(time)
      WRITE(*,*)'>>> EXECUTION TIME =',time-SEC

      stop
      end

c******************************************************************
      subroutine GetGaussFactors(File,Points,x,w)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine retrieves the points and weights for
c      Gaussian Quadrature from a given file
c
c     Variables:
c      File		name of file containing points and 
c			 weights
c      Points		number of points to be used for 
c			quadrature
c      x		array of nodes
c      w		array of weights
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      integer Points
      double precision x(Points),w(Points)
      character*64 File
      
      integer i,tempPoints
      
      open(unit=7,file=File(1:index(File,' ')-1))

       do i = 1,18
        read(7,*)
       enddo
       read(7,*) tempPoints
       do while (tempPoints .ne. Points)
        do i = 1,tempPoints
         read(7,*)
        enddo
        read(7,*) tempPoints
       enddo
 
       do i = 1,Points
        read(7,*) x(i),w(i)
       enddo
      
      close(unit=7)

      return
      end
      
      double precision function BSpline(Order,Deriv,xNumPoints,xPoints,n,t,b,x)
      
      integer Order,Deriv,xNumPoints,n
      double precision xPoints(*),t(*),b(xNumPoints+Order),x
      
      integer i
      double precision bvalue

      b = 0.0d0
      b(n) = 1.0d0

      BSpline = bvalue(t,b,n,Order+1,x,Deriv)

      return
      end

      subroutine CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,MatrixDim,xBounds,xNumPoints,Deriv,u)

      integer Left,Right,Order,LegPoints,MatrixDim,xBounds(*),xNumPoints,Deriv
      double precision xPoints(*),xLeg(*)
      double precision u(LegPoints,xNumPoints,MatrixDim)

      integer i,k,l,Count
      integer, allocatable :: t(:)
      double precision BSpline
      double precision ax,bx,xIntScale,xScaledZero
      double precision, allocatable :: x(:,:),tx(:),b(:)

      allocate(x(LegPoints,xNumPoints))
      allocate(t(xNumPoints+2*Order),tx(xNumPoints+2*Order),b(xNumPoints+Order))

      do i = 1,Order
       t(i) = 1
       tx(i) = xPoints(1)
      enddo
      do i = 1,xNumPoints
       t(i+Order) = i
       tx(i+Order) = xPoints(i)
      enddo
      do i = 1,Order
       t(i+Order+xNumPoints) = xNumPoints
       tx(i+Order+xNumPoints) = xPoints(xNumPoints)
      enddo

      select case (Left)
       case (0:1)
        select case (Right)
         case (0:1)
          do i = 2,xNumPoints+2*Order-1
           xBounds(i-1) = t(i)
          enddo
         case (2)
          do i = 2,xNumPoints+2*Order
           xBounds(i-1) = t(i)
          enddo
        end select
       case (2)
        select case (Right)
         case (0:1)
          do i = 1,xNumPoints+2*Order-1
           xBounds(i) = t(i)
          enddo
         case (2)
          do i = 1,xNumPoints+2*Order
           xBounds(i) = t(i)
          enddo
        end select
      end select

      deallocate(t)

      u = 0.0d0

      do k = 1,xNumPoints-1
       ax = xPoints(k)
       bx = xPoints(k+1)
       xIntScale = 0.5d0*(bx-ax)
       xScaledZero = 0.5d0*(bx+ax)
       do l = 1,LegPoints
        x(l,k) = xIntScale*xLeg(l)+xScaledZero
       enddo
      enddo

      Count = 1
      select case (Left)
       case (0)
        do k = xBounds(Count),xBounds(Count+Order+1)-1
         do l = 1,LegPoints
          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(l,k))
         enddo
        enddo
        Count = Count + 1
       case (1)
        do k = xBounds(Count),xBounds(Count+Order+1)-1
         do l = 1,LegPoints
          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,1,tx,b,x(l,k))+BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(l,k))
         enddo
        enddo
        Count = Count + 1
       case(2)
        do i = 1,2
         do k = xBounds(Count),xBounds(Count+Order+1)-1
          do l = 1,LegPoints
           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(l,k))
          enddo
         enddo
         Count = Count + 1
        enddo
      end select

      do i = 3,xNumPoints+Order-3
       do k = xBounds(Count),xBounds(Count+Order+1)-1
        do l = 1,LegPoints
         u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(l,k))
        enddo
       enddo
       Count = Count + 1
      enddo

      select case (Right)
       case (0)
        do k = xBounds(Count),xBounds(Count+Order+1)-1
         do l = 1,LegPoints
          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(l,k))
         enddo
        enddo
       case (1)
        do k = xBounds(Count),xBounds(Count+Order+1)-1
         do l = 1,LegPoints
          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(l,k))+
     >                   BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,tx,b,x(l,k))
         enddo
        enddo
       case(2)
        do i = xNumPoints+Order-2,xNumPoints+Order-1
         do k = xBounds(Count),xBounds(Count+Order+1)-1
          do l = 1,LegPoints
           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(l,k))
          enddo
         enddo
         Count = Count + 1
        enddo
      end select

      deallocate(x,tx,b)

      return
      end

      subroutine CalcOverlapMatrix(NumChannels,NumOpenChannels,Order,xPoints,LegPoints,xLeg,wLeg,xClosedDim,xOpenDim,xDim,
     >                          xNumPoints,u,xBounds,SData,HalfBandWidth,ClosedMatrixDim,OpenMatrixDim,Scc,Sco,Soc,Soo)

      integer NumChannels,NumOpenChannels,Order,LegPoints,xClosedDim,xOpenDim,xDim,xNumPoints,xBounds(*)
      integer HalfBandWidth,OpenMatrixDim,ClosedMatrixDim,MatrixDim
      double precision xPoints(*),xLeg(*),wLeg(*)
      double precision u(LegPoints,xNumPoints,xDim)
      double precision SData(Order+1,NumChannels,NumChannels,xDim)
c      double precision Scc(ClosedMatrixDim,ClosedMatrixDim),Sco(NumChannels*Order,OpenMatrixDim)
      double precision Scc(2*HalfBandWidth+1,ClosedMatrixDim),Sco(NumChannels*Order,OpenMatrixDim)
      double precision Soc(OpenMatrixDim,NumChannels*Order),Soo(OpenMatrixDim,OpenMatrixDim)

      integer ic,icp,ix,ixp,kx,lx
      integer NonZero,Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      double precision ax,bx,xTempS
      double precision, allocatable :: xIntScale(:),xS(:,:)

      allocate(xIntScale(xNumPoints),xS(xDim,xDim))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))

      do kx = 1,xNumPoints-1
       ax = xPoints(kx)
       bx = xPoints(kx+1)
       xIntScale(kx) = 0.5d0*(bx-ax)
      enddo

      do ix = 1,xDim
       do ixp = 1,xDim
        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
       enddo
      enddo

      do ix = 1,xDim
       do ixp = max(1,ix-Order),min(xDim,ix+Order)
        xS(ixp,ix) = 0.0d0
        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
         xTempS = 0.0d0
         do lx = 1,LegPoints
          xTempS = xTempS + wLeg(lx)*u(lx,kx,ix)*u(lx,kx,ixp)
         enddo
         xS(ixp,ix) = xS(ixp,ix) + xIntScale(kx)*xTempS
        enddo
       enddo
      enddo

c cc
      Scc = 0.0d0
      do ix = 1,xClosedDim
       do ic = 1,NumChannels
        Row = (ix-1)*NumChannels+ic
        do ixp = max(1,ix-Order),min(xClosedDim,ix+Order)
         do icp = 1,NumChannels
          Col = (ixp-1)*NumChannels+icp
          NewRow = HalfBandWidth+1+Row-Col
          if (ix .le. ixp) then
c           Scc(Row,Col) = xS(ixp,ix)*SData(Order+1-(ixp-ix),icp,ic,ixp)
           Scc(NewRow,Col) = xS(ixp,ix)*SData(Order+1-(ixp-ix),icp,ic,ixp)
          else
c           Scc(Row,Col) = xS(ixp,ix)*SData(Order+1-(ix-ixp),ic,icp,ix)
           Scc(NewRow,Col) = xS(ixp,ix)*SData(Order+1-(ix-ixp),ic,icp,ix)
          endif
         enddo
        enddo
       enddo
      enddo

c oo

      Soo = 0.0d0
      do ic = 1,NumOpenChannels
       Row = (xDim-1)*NumChannels+ic-ClosedMatrixDim
       Soo(Row,Row) = xS(xDim,xDim)
      enddo

c co, oc
      Sco = 0.0d0
      do ix = xClosedDim-Order+1,xClosedDim
       do ic = 1,NumChannels
        Row = (ix-1)*NumChannels+ic-(xClosedDim-Order)*NumChannels
        do icp = 1,NumOpenChannels
         Col = (xDim-1)*NumChannels+icp-ClosedMatrixDim
         Sco(Row,Col) = xS(xDim,ix)*SData(Order+1-(xDim-ix),icp,ic,xDim)
        enddo
       enddo
      enddo

      do Row = 1,OpenMatrixDim
       do Col = 1,NumChannels*Order
        Soc(Row,Col) = Sco(Col,Row)
       enddo
      enddo

      deallocate(xIntScale,xS)
      deallocate(kxMin,kxMax)

      return
      end

      subroutine CalcHamiltonian(Mass,NumChannels,NumOpenChannels,Order,xPoints,LegPoints,xLeg,wLeg,xClosedDim,xOpenDim,xDim,
     >             xNumPoints,u,ux,uxx,xBounds,HalfBandWidth,OpenMatrixDim,ClosedMatrixDim,RData,SData,HData,Hcc,Hco,Hoc,Hoo)

      integer NumChannels,NumOpenChannels,Order,LegPoints,xClosedDim,xOpenDim,xDim,xNumPoints,xBounds(*),HalfBandWidth
      integer OpenMatrixDim,ClosedMatrixDim
      double precision Mass,xPoints(*),xLeg(*),wLeg(*)
      double precision u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)
      double precision SData(Order+1,NumChannels,NumChannels,xDim)
      double precision RData(xDim),HData(Order+6,Order+1,NumChannels,NumChannels,xDim)
c      double precision Hcc(ClosedMatrixDim,ClosedMatrixDim),Hco(NumChannels*Order,OpenMatrixDim)
      double precision Hcc(2*HalfBandWidth+1,ClosedMatrixDim),Hco(NumChannels*Order,OpenMatrixDim)
      double precision Hoc(OpenMatrixDim,NumChannels*Order),Hoo(OpenMatrixDim,OpenMatrixDim)

      integer i,j,k,ic,icp,ix,ixp,kx,lx
      integer HalfOrder,NumDataPoints,NumInterpPoints
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:)
      double precision m,a,b
      double precision ax,bx,xScaledZero
      double precision xTempT,xTempH,xH
      double precision, allocatable :: x(:,:),x2(:,:),xT(:,:),RTempData(:),HTempData(:),Hpp(:),HInterp(:),xIntScale(:)

      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
      allocate(x(LegPoints,xNumPoints),x2(LegPoints,xNumPoints),xIntScale(xNumPoints))
      allocate(xT(xDim,xDim),RTempData(Order+6),HTempData(Order+6),Hpp(Order+6),HInterp((Order+1)*LegPoints))

      m = -0.5d0/Mass

      do kx = 1,xNumPoints-1
       ax = xPoints(kx)
       bx = xPoints(kx+1)
       xIntScale(kx) = 0.5d0*(bx-ax)
       xScaledZero = 0.5d0*(bx+ax)
       do lx = 1,LegPoints
        x(lx,kx) = xIntScale(kx)*xLeg(lx)+xScaledZero
        x2(lx,kx) = 1.0d0/x(lx,kx)**2
       enddo
      enddo

      do ix = 1,xDim
       do ixp = 1,xDim
        kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
        kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
       enddo
      enddo

      do ix = 1,xDim
       do ixp = max(1,ix-Order),min(xDim,ix+Order)
        xT(ixp,ix) = 0.0d0
        do kx = kxMin(ixp,ix),kxMax(ixp,ix)
         xTempT = 0.0d0
         do lx = 1,LegPoints
          xTempT = xTempT + wLeg(lx)*u(lx,kx,ix)*(uxx(lx,kx,ixp)-3.75d0*u(lx,kx,ixp)*x2(lx,kx)) ! The 3.75=15/4 is for 3-body case
         enddo
         xT(ixp,ix) = xT(ixp,ix) + m*xIntScale(kx)*xTempT
        enddo
       enddo
      enddo

      HalfOrder = (Order+1)/2.0+3.5

c cc
      Hcc = 0.0d0

      do ix = 1,xClosedDim
       do ic = 1,NumChannels
        Row = (ix-1)*NumChannels+ic
        do ixp = max(1,ix-Order),min(xClosedDim,ix+Order)
         do icp = 1,NumChannels
          Col = (ixp-1)*NumChannels+icp
          NewRow = HalfBandWidth+1+Row-Col
          if (ix .le. ixp) then
c           Hcc(Row,Col) = xT(ixp,ix)*SData(Order+1-(ixp-ix),icp,ic,ixp)
           Hcc(NewRow,Col) = xT(ixp,ix)*SData(Order+1-(ixp-ix),icp,ic,ixp)
          else
c           Hcc(Row,Col) = xT(ixp,ix)*SData(Order+1-(ix-ixp),ic,icp,ix)
           Hcc(NewRow,Col) = xT(ixp,ix)*SData(Order+1-(ix-ixp),ic,icp,ix)
          endif
         enddo
        enddo
       enddo
      enddo

      do ix = 1,xClosedDim
       do ic = 1,NumChannels
        Row = (ix-1)*NumChannels+ic
        do ixp = max(1,ix-Order),min(xClosedDim,ix+Order)
         do icp = 1,NumChannels
          Col = (ixp-1)*NumChannels+icp
          if (Col .ge. Row) then
           if (ix .le. ixp) then
            j = 0
            do i = 1,Order+6-(ixp-ix)
             k = ixp+i-HalfOrder
             if ((k .ge. 1) .AND. (k .le. xClosedDim)) then
              j = j + 1
              RTempData(j) = RData(k)
              HTempData(j) = HData(i,Order+1-(ixp-ix),icp,ic,ixp)
             endif
            enddo
            NumDataPoints = j
           else
            j = 0
            do i = 1,Order+6-(ix-ixp)
             k = ix+i-HalfOrder
             if ((k .ge. 1) .AND. (k .le. xClosedDim)) then
              j = j + 1
              RTempData(j) = RData(k)
              HTempData(j) = HData(i,Order+1-(ix-ixp),ic,icp,ix)
             endif
            enddo
            NumDataPoints = j
           endif
           call CubicSpline(RTempData,HTempData,NumDataPoints,Hpp)
           NumInterpPoints = (kxMax(ixp,ix)-kxMin(ixp,ix)+1)*LegPoints
           call CubicSplineInterp(NumDataPoints,RTempData,HTempData,Hpp,x(1,kxMin(ixp,ix)),HInterp,NumInterpPoints)
           xH = 0.0d0
           i = 1
           do kx = kxMin(ixp,ix),kxMax(ixp,ix)
            xTempH = 0.0d0
            do lx = 1,LegPoints
             xTempH = xTempH + wLeg(lx)*u(lx,kx,ix)*HInterp(i)*u(lx,kx,ixp)
             i = i + 1
            enddo
            xH = xH + xIntScale(kx)*xTempH
           enddo
           NewRow = HalfBandWidth+1+Row-Col
c           Hcc(Row,Col) = Hcc(Row,Col) + xH
           Hcc(NewRow,Col) = Hcc(NewRow,Col) + xH
           Hcc(HalfBandWidth+1-Row+Col,Row) = Hcc(NewRow,Col)
          endif
         enddo
        enddo
       enddo
      enddo

c oo

      Hoo = 0.0d0

      do ic = 1,NumOpenChannels
       Row = (xDim-1)*NumChannels+ic-ClosedMatrixDim
       Hoo(Row,Row) = xT(xDim,xDim) - m*u(1,xNumPoints,xDim)*ux(1,xNumPoints,xDim)
      enddo

      do ic = 1,NumOpenChannels
       Row = (xDim-1)*NumChannels+ic-ClosedMatrixDim
       do icp = 1,NumOpenChannels
        Col = (xDim-1)*NumChannels+icp-ClosedMatrixDim
        j = 0
        do i = 1,Order+6
         k = xDim+i-HalfOrder
         if ((k .ge. 1) .AND. (k .le. xClosedDim)) then
          j = j + 1
          RTempData(j) = RData(k)
          HTempData(j) = HData(i,Order+1,icp,ic,xDim)
         endif
        enddo
        NumDataPoints = j
        call CubicSpline(RTempData,HTempData,NumDataPoints,Hpp)
        NumInterpPoints = (kxMax(xDim,xDim)-kxMin(xDim,xDim)+1)*LegPoints
        call CubicSplineInterp(NumDataPoints,RTempData,HTempData,Hpp,x(1,kxMin(xDim,xDim)),HInterp,NumInterpPoints)
        xH = 0.0d0
        i = 1
        do kx = kxMin(xDim,xDim),kxMax(xDim,xDim)
         xTempH = 0.0d0
         do lx = 1,LegPoints
          xTempH = xTempH + wLeg(lx)*u(lx,kx,xDim)*HInterp(i)*u(lx,kx,xDim)
          i = i + 1
         enddo
         xH = xH + xIntScale(kx)*xTempH
        enddo
        Hoo(Row,Col) = Hoo(Row,Col) + xH
       enddo
      enddo

c co

      Hco = 0.0d0

      do ix = xClosedDim-Order+1,xClosedDim
       do ic = 1,NumChannels
        Row = (ix-1)*NumChannels+ic-(xClosedDim-Order)*NumChannels
        do icp = 1,NumOpenChannels
         Col = (xDim-1)*NumChannels+icp-ClosedMatrixDim
         Hco(Row,Col) = (xT(xDim,ix)-m*u(1,xNumPoints,ix)*ux(1,xNumPoints,xDim))*SData(Order+1-(xDim-ix),icp,ic,xDim)
        enddo
       enddo
      enddo

      do ix = xClosedDim-Order+1,xClosedDim
       do ic = 1,NumChannels
        Row = (ix-1)*NumChannels+ic-(xClosedDim-Order)*NumChannels
        do icp = 1,NumOpenChannels
         Col = (xDim-1)*NumChannels+icp-ClosedMatrixDim
         j = 0
         do i = 1,Order+6-(xDim-ix)
          k = xDim+i-HalfOrder
          if ((k .ge. 1) .AND. (k .le. xClosedDim)) then
           j = j + 1
           RTempData(j) = RData(k)
           HTempData(j) = HData(i,Order+1-(xDim-ix),icp,ic,xDim)
          endif
         enddo
         NumDataPoints = j
         call CubicSpline(RTempData,HTempData,NumDataPoints,Hpp)
         NumInterpPoints = (kxMax(xDim,ix)-kxMin(xDim,ix)+1)*LegPoints
         call CubicSplineInterp(NumDataPoints,RTempData,HTempData,Hpp,x(1,kxMin(xDim,ix)),HInterp,NumInterpPoints)
         xH = 0.0d0
         i = 1
         do kx = kxMin(xDim,ix),kxMax(xDim,ix)
          xTempH = 0.0d0
          do lx = 1,LegPoints
           xTempH = xTempH + wLeg(lx)*u(lx,kx,ix)*HInterp(i)*u(lx,kx,xDim)
           i = i + 1
          enddo
          xH = xH + xIntScale(kx)*xTempH
         enddo
         Hco(Row,Col) = Hco(Row,Col) + xH
        enddo
       enddo
      enddo

      do Row = 1,OpenMatrixDim
       do Col = 1,NumChannels*Order
        Hoc(Row,Col) = Hco(Col,Row)
       enddo
      enddo

      deallocate(kxMin,kxMax)
      deallocate(x,x2,xIntScale)
      deallocate(xT,RTempData,HTempData,Hpp,HInterp)

      return
      end

      double precision function bvalue ( t, bcoef, n, k, x, jderiv )
c  from  * a practical guide to splines *  by c. de boor    
calls  interv
c
calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef(.),                     ,  j .eq. 0
c                   /
c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t(.+k-j) - t(.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c    a(.,x)  =    ---------------------------------------
c                 (x - t(.))      + (t(.+m-1) - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c
      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1
     *                     ,mflag,nmi,jdrvp1
      parameter (kmax = 20)
C     double precision bcoef(n),t(1),x,   aj(20),dl(20),dr(20),fkmj
      double precision bcoef(*),t(*),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
c      dimension t(n+k)
c former fortran standard made it impossible to specify the length of  t
c  precisely without the introduction of otherwise superfluous addition-
c  al arguments.
      bvalue = 0.0d0
      if (jderiv .ge. k)                go to 99
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c      at  t(n+k) where it is leftcontinuous.)
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do j=1,i
         dl(j) = x - t(i+1-j)
      enddo
      do j=i,km1
         aj(k-j) = 0.0d0
         dl(j) = dl(i)
      enddo
                                        go to 10
    8 do j=1,km1
         dl(j) = x - t(i+1-j)
      enddo
c
   10 jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do j=1,jcmax
         dr(j) = t(i+j) - x
      enddo
      do j=jcmax,km1
         aj(j+1) = 0.0d0
         dr(j) = dr(jcmax)
      enddo
                                        go to 20
   18 do j=1,km1
         dr(j) = t(i+j) - x
      enddo
c
   20 do jc=jcmin,jcmax
         aj(jc) = bcoef(imk + jc)
      enddo
c
c               *** difference the coefficients  jderiv  times.
      if (jderiv .eq. 0)                go to 30
      do j=1,jderiv
         kmj = k-j
         fkmj = float(kmj)
         ilo = kmj
         do jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
            ilo = ilo - 1
         enddo
      enddo
c
c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
   30 if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1     
      do j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
            ilo = ilo - 1
         enddo
      enddo
   39 bvalue = aj(1)
c
   99                                   return
      end
      subroutine interv ( xt, lxt, x, left, mflag )
c  from  * a practical guide to splines *  by C. de Boor    
computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a double precision sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the left at  xt(lxt) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
      data ilo /1/
      save ilo  
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
c
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
										go to 111
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ratint(xa,ya,n,x,y,dy)

      INTEGER n
      double precision dy,x,y,xa(n),ya(n),TINY
      PARAMETER (TINY=1.0d-25)

      INTEGER i,m,ns
      double precision dd,h,hh,t,w
      double precision, allocatable :: c(:),d(:)

      allocate(c(n),d(n))

      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.)pause 'failure in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue

      deallocate(c,d)

      return
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CubicSpline(x,y,N,ypp)

      integer N
      double precision x(*),y(*),ypp(*)

      integer NMax,i,k
      double precision Sig,qN,uN,p
      double precision, allocatable :: u(:)

      allocate(u(N))

      ypp(1) = 0.0d0
      u(1) = 0.0d0
      do i = 2,N-1
       Sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p = Sig * ypp(i-1)+2.0d0
       ypp(i) = (Sig-1.0d0)/p
       u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-Sig*u(i-1))/p
      enddo
      qN=0.0d0
      uN=0.0d0
      ypp(N)=(uN-qN*u(N-1))/(qN*ypp(N-1)+1.0d0)
      do k = N-1,1,-1
       ypp(k)= ypp(k)*ypp(k+1)+u(k)
      enddo

      deallocate(u)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CubicSplineDerivInterp(N,x,y,ypp,xInterp,ypInterp,NumPoints)

      integer N,NumPoints
      double precision x(*),y(*),ypp(*),xInterp(*),ypInterp(*)

      integer i,j
      double precision xDelt,yDelt,a,b,c,d

      j = 1
      do i = 1,NumPoints
       do while ((x(j+1) .lt. xInterp(i)) .AND. (j .lt. N))
        j = j + 1
       enddo
       xDelt = x(j+1)-x(j)
       yDelt = y(j+1)-y(j)
       a = (x(j+1)-xInterp(i))/xDelt
       b = 1.0d0-a
       ypInterp(i) = yDelt/xDelt+((3.0d0*b*b-1.0d0)*ypp(j+1)-(3.0d0*a*a-1.0d0)*ypp(j))*xDelt/6.0d0
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CubicSplineInterp(N,x,y,ypp,xInterp,yInterp,NumPoints)

      integer N,NumPoints
      double precision x(*),y(*),ypp(*),xInterp(*),yInterp(*)

      integer i,j
      double precision Delt,c,d

      j = 1
      do i = 1,NumPoints
       do while ((x(j+1) .lt. xInterp(i)) .AND. (j .lt. N))
        j = j + 1
       enddo
       Delt = x(j+1)-x(j)
       c = (x(j+1)-xInterp(i))/Delt
       d = (xInterp(i)-x(j))/Delt
       yInterp(i) = c*y(j)+d*y(j+1)+((c**3-c)*ypp(j)+(d**3-d)*ypp(j+1))*Delt**2/6.0d0
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GetData(SFile,HFile,Order,xDim,NumChannels,RData,SData,HData)

      integer Order,xDim,NumChannels
      double precision RData(xDim)
      double precision SData(Order+1,NumChannels,NumChannels,xDim)
      double precision HData(Order+6,Order+1,NumChannels,NumChannels,xDim)
      character*64 SFile,HFile

      integer i,j,k,l
      integer NumInpChannels,i1,i2,i3
      integer xNumPoints,HalfOrder
      integer, allocatable :: Count(:)
      double precision xIntScale,xScaledZero,Sum

      open(unit=7,file=SFile)
       read(7,*)
       read(7,*) NumInpChannels
       read(7,*)
       read(7,*)
       read(7,*)
       read(7,*)
       read(7,*)
       read(7,*) xNumPoints
       do i = 1,xNumPoints
        read(7,*)
       enddo
       read(7,*)
       read(7,*)
       read(7,*)
       do i = 1,xDim
        read(7,*) RData(i)
       enddo
       SData = 0.0d0
       read(7,*)
       do i = 1,xDim
        do j = 1,min(i,Order+1)
         read(7,*) i1,i2
         do k = 1,NumChannels
          read(7,*) SData(Order+1-(i-i2),1:NumChannels,k,i)
         enddo
         do k = NumChannels+1,NumInpChannels
          read(7,*)
         enddo
        enddo
        do j = 1,NumChannels
         do k = 1,NumChannels
          SData(Order+1,j,k,i) = 0.0d0
         enddo
         SData(Order+1,j,j,i) = 1.0d0
        enddo
       enddo
      close(unit=7)

      HalfOrder = (Order+1)/2.0+3.5

c HData - first index: R; second: spline

      HData = 0.0d0
      open(unit=7,file=HFile,form='unformatted')
       do i = 1,xDim
        do j = 1,min(i,Order+1)
         do k = max(1,j-HalfOrder+1),min(i,HalfOrder)
          read(7) i1,i2,i3
          do l = 1,NumChannels
           read(7) HData(HalfOrder-(i2-i1),Order+1-(i2-i3),1:NumChannels,l,i2)
          enddo
          do l = NumChannels+1,NumInpChannels
           read(7)
          enddo
         enddo
        enddo
        do j = 2,min(i,Order+1)
         do k = j,min(i,HalfOrder)
          read(7) i1,i2,i3
          do l = 1,NumChannels
           read(7) HData(HalfOrder-(i2-i1),Order+1-(i2-i3),1:NumChannels,l,i2)
          enddo
          do l = NumChannels+1,NumInpChannels
           read(7)
          enddo
         enddo
        enddo
       enddo
      close(unit=7)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcPlotBasisFuncs(Left,Right,Order,xPoints,x,xNumPlotPoints,MatrixDim,xBounds,xNumPoints,Deriv,u)

      integer Left,Right,Order,LegPoints,xNumPlotPoints,MatrixDim,xBounds(*),xNumPoints,Deriv
      double precision xPoints(*),x(*)
      double precision u(MatrixDim,xNumPlotPoints)

      integer i,k,l,Count
      integer, allocatable :: t(:)
      double precision BSpline
      double precision, allocatable :: tx(:),b(:)

      allocate(t(xNumPoints+2*Order),tx(xNumPoints+2*Order),b(xNumPoints+Order))

      do i = 1,Order
       t(i) = 1
       tx(i) = xPoints(1)
      enddo
      do i = 1,xNumPoints
       t(i+Order) = i
       tx(i+Order) = xPoints(i)
      enddo
      do i = 1,Order
       t(i+Order+xNumPoints) = xNumPoints
       tx(i+Order+xNumPoints) = xPoints(xNumPoints)
      enddo

      do i = 1,xNumPoints+2*Order
       xBounds(i) = 0.0d0
      enddo

      select case (Left)
       case (0:1)
        select case (Right)
         case (0:1)
          do i = 2,xNumPoints+2*Order-1
           xBounds(i-1) = t(i)
          enddo
         case (2)
          do i = 2,xNumPoints+2*Order
           xBounds(i-1) = t(i)
          enddo
        end select
       case (2)
        select case (Right)
         case (0:1)
          do i = 1,xNumPoints+2*Order-1
           xBounds(i) = t(i)
          enddo
         case (2)
          do i = 1,xNumPoints+2*Order
           xBounds(i) = t(i)
          enddo
        end select
      end select

      deallocate(t)

      Count = 1
      select case (Left)
       case (0)
        do k = 1,xNumPlotPoints
         u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(k))
        enddo
        Count = Count + 1
       case (1)
        do k = 1,xNumPlotPoints
         u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,1,tx,b,x(k))+BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(k))
        enddo
        Count = Count + 1
       case(2)
        do i = 1,2
         do k = 1,xNumPlotPoints
          u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(k))
         enddo
         Count = Count + 1
        enddo
      end select

      do i = 3,xNumPoints+Order-3
       do k = 1,xNumPlotPoints
        u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(k))
       enddo
       Count = Count + 1
      enddo

      select case (Right)
       case (0)
        do k = 1,xNumPlotPoints
         u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(k))
        enddo
       case (1)
        do k = 1,xNumPlotPoints
         u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(k))+
     >                BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,tx,b,x(k))
        enddo
       case(2)
        do i = xNumPoints+Order-2,xNumPoints+Order-1
         do k = 1,xNumPlotPoints
          u(Count,k) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(k))
         enddo
         Count = Count + 1
        enddo
      end select

      deallocate(tx,b)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BesselBasePair(k,r,l,f,fp,g,gp)

      double precision l,k,r,f,fp,g,gp

      integer lInt,info
      double precision x,a,b,alpha,alphaComp,Pi
      double precision pj,py,mj,my
      double precision, allocatable :: j(:),y(:)

      Pi = 3.1415926535897932385d0

      x = k*r
      a = dsqrt(r)
      b = 1.0d0/a
      lInt = int(l+0.5d0)
      alpha = l+0.5d0-dfloat(lInt)
      alphaComp = dabs(alpha-1.0d0)

      if (lInt .eq. 0) then

       allocate(j(2),y(2))

       call RJBESL(x,alpha,2,j,info)
       call RYBESL(x,alpha,2,y,info)
       pj = j(1)
       py = y(1)

       call RJBESL(x,alphaComp,2,j,info)
       call RYBESL(x,alphaComp,2,y,info)
       mj = dcos(Pi*alphaComp)*j(1)-dsin(Pi*alphaComp)*y(1)
       my = dsin(Pi*alphaComp)*j(1)+dcos(Pi*alphaComp)*y(1)

       f = a*pj
       g = a*py

       fp = b*(x*mj-l*pj)
       gp = b*(x*my-l*py)

       deallocate(j,y)

      else

       allocate(j(lInt+1),y(lInt+1))

       call RJBESL(x,alpha,lInt+1,j,info)
       call RYBESL(x,alpha,lInt+1,y,info)

       f = a*j(lInt+1)
       g = a*y(lInt+1)

       fp = b*(x*j(lInt)-l*j(lInt+1))
       gp = b*(x*y(lInt)-l*y(lInt+1))

       deallocate(j,y)

      endif

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE RJBESL(X, ALPHA, NB, B, NCALC)
C---------------------------------------------------------------------
C This routine calculates Bessel functions J sub(N+ALPHA) (X)
C   for non-negative argument X, and non-negative order N+ALPHA.
C
C
C  Explanation of variables in the calling sequence.
C
C   X     - working precision non-negative real argument for which
C           J's are to be calculated.
C   ALPHA - working precision fractional part of order for which
C           J's or exponentially scaled J'r (J*exp(X)) are
C           to be calculated.  0 <= ALPHA < 1.0.
C   NB  - integer number of functions to be calculated, NB > 0.
C           The first function calculated is of order ALPHA, and the
C           last is of order (NB - 1 + ALPHA).
C   B  - working precision output vector of length NB.  If RJBESL
C           terminates normally (NCALC=NB), the vector B contains the
C           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
C           corresponding exponentially scaled functions.
C   NCALC - integer output variable indicating possible errors.
C           Before using the vector B, the user should check that
C           NCALC=NB, i.e., all orders have been calculated to
C           the desired accuracy.  See Error Returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C  Explanation of machine-dependent constants
C
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
C            then at least N iterations of the backward recursion
C            will be executed.  The value of 10.0 ** 4 is used on
C            every machine.
C
C
C     Approximate values for some important machines are:
C
C
C                            it    NSIG    ENTEN       ENSIG
C
C   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15
C   Cyber 180/855
C     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16
C   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5
C   VAX           (S.P.)     24      8    1.0E+38     1.0E+8
C   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17
C   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16
C
C
C                           RTNSIG      ENMTEN      XLARGE
C
C   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4
C   Cyber 180/855
C     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4
C   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4
C   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4
C   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4
C   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4
C
C*******************************************************************
C*******************************************************************
C
C  Error returns
C
C    In case of an error,  NCALC .NE. NB, and not all J's are
C    calculated to the desired accuracy.
C
C    NCALC .LT. 0:  An argument is out of range. For example,
C       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
C       In this case, B(1) is set to zero, the remainder of the
C       B-vector is not calculated, and NCALC is set to
C       MIN(NB,0)-1 so that NCALC .NE. NB.
C
C    NB .GT. NCALC .GT. 0: Not all requested function values could
C       be calculated accurately.  This usually occurs because NB is
C       much larger than ABS(X).  In this case, B(N) is calculated
C       to the desired accuracy for N .LE. NCALC, but precision
C       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C       for N .GT. NCALC (because it is too small to be represented),
C       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C       significant figures of B(N) can be trusted.
C
C
C  Intrinsic and other functions required are:
C
C     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,
C
C     REAL, SIN, SQRT
C
C
C  Acknowledgement
C
C   This program is based on a program written by David J. Sookne
C   (2) that computes values of the Bessel functions J or I of real
C   argument and integer order.  Modifications include the restriction
C   of the computation to the J Bessel function of non-negative real
C   argument, the extension of the computation to arbitrary positive
C   order, and the elimination of most underflow.
C
C  References: "A Note on Backward Recurrence Algorithms," Olver,
C               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C               pp 941-947.
C
C              "Bessel Functions of Real Argument and Integer Order,"
C               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C               125-132.
C
C  Latest modification: March 19, 1990
C
C  Author: W. J. Cody
C          Applied Mathematics Division
C          Argonne National Laboratory
C          Argonne, IL  60439
C
C---------------------------------------------------------------------
      INTEGER I,J,K,L,M,MAGX,N,NB,NBMX,NCALC,NEND,NSTART
CS    REAL               GAMMA,
      DOUBLE PRECISION  DGAMMA,
     1 ALPHA,ALPEM,ALP2EM,B,CAPP,CAPQ,CONV,EIGHTH,EM,EN,ENMTEN,ENSIG,
     2 ENTEN,FACT,FOUR,FUNC,GNU,HALF,HALFX,ONE,ONE30,P,PI2,PLAST,
     3 POLD,PSAVE,PSAVEL,RTNSIG,S,SUM,T,T1,TEMPA,TEMPB,TEMPC,TEST,
     4 THREE,THREE5,TOVER,TWO,TWOFIV,TWOPI1,TWOPI2,X,XC,XIN,XK,XLARGE,
     5 XM,VCOS,VSIN,Z,ZERO
      DIMENSION B(NB), FACT(25)
C---------------------------------------------------------------------
C  Mathematical constants
C
C   PI2    - 2 / PI
C   TWOPI1 - first few significant digits of 2 * PI
C   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
C            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
C---------------------------------------------------------------------
CS    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0,
CS   1 1.935307179586476925286767E-3/
CS    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/
CS    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/
CS    DATA ONE30, THREE5 /130.0E0,35.0E0/
      DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535D0,6.28125D0,
     1 1.935307179586476925286767D-3/
      DATA ZERO, EIGHTH, HALF, ONE /0.0D0,0.125D0,0.5D0,1.0D0/
      DATA TWO, THREE, FOUR, TWOFIV /2.0D0,3.0D0,4.0D0,25.0D0/
      DATA ONE30, THREE5 /130.0D0,35.0D0/
C---------------------------------------------------------------------
C  Machine-dependent parameters
C---------------------------------------------------------------------
CS    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN, XLARGE /1.2E-37,1.0E4/
      DATA ENTEN, ENSIG, RTNSIG /1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN, XLARGE /8.9D-308,1.0D4/
C---------------------------------------------------------------------
C     Factorial(N)
C---------------------------------------------------------------------
CS    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3,
CS   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9,
CS   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14,
CS   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18,
CS   4 5.109094217170944E19,1.12400072777760768E21,
CS   5 2.585201673888497664E22,6.2044840173323943936E23/
      DATA FACT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,1.2D2,7.2D2,5.04D3,
     1 4.032D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,
     2 8.71782912D10,1.307674368D12,2.0922789888D13,3.55687428096D14,
     3 6.402373705728D15,1.21645100408832D17,2.43290200817664D18,
     4 5.109094217170944D19,1.12400072777760768D21,
     5 2.585201673888497664D22,6.2044840173323943936D23/
C---------------------------------------------------------------------
C Statement functions for conversion and the gamma function.
C---------------------------------------------------------------------
CS    CONV(I) = REAL(I)
CS    FUNC(X) = GAMMA(X)
      CONV(I) = DBLE(I)
      FUNC(X) = DGAMMA(X)
C---------------------------------------------------------------------
C Check for out of range arguments.
C---------------------------------------------------------------------
      MAGX = INT(X)
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (X.LE.XLARGE) 
     1       .AND. (ALPHA.GE.ZERO) .AND. (ALPHA.LT.ONE))  
     2   THEN
C---------------------------------------------------------------------
C Initialize result array to zero.
C---------------------------------------------------------------------
            NCALC = NB
            DO 20 I=1,NB
              B(I) = ZERO
   20       CONTINUE
C---------------------------------------------------------------------
C Branch to use 2-term ascending series for small X and asymptotic
C form for large X when NB is not too large.
C---------------------------------------------------------------------
            IF (X.LT.RTNSIG) THEN
C---------------------------------------------------------------------
C Two-term ascending series for small X.
C---------------------------------------------------------------------
               TEMPA = ONE
               ALPEM = ONE + ALPHA
               HALFX = ZERO
               IF (X.GT.ENMTEN) HALFX = HALF*X
               IF (ALPHA.NE.ZERO)
     1            TEMPA = HALFX**ALPHA/(ALPHA*FUNC(ALPHA))
               TEMPB = ZERO
               IF ((X+ONE).GT.ONE) TEMPB = -HALFX*HALFX
               B(1) = TEMPA + TEMPA*TEMPB/ALPEM
               IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0
               IF (NB .NE. 1) THEN
                  IF (X .LE. ZERO) THEN
                        DO 30 N=2,NB
                          B(N) = ZERO
   30                   CONTINUE
                     ELSE
C---------------------------------------------------------------------
C Calculate higher order functions.
C---------------------------------------------------------------------
                        TEMPC = HALFX
                        TOVER = (ENMTEN+ENMTEN)/X
                        IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB
                        DO 50 N=2,NB
                          TEMPA = TEMPA/ALPEM
                          ALPEM = ALPEM + ONE
                          TEMPA = TEMPA*TEMPC
                          IF (TEMPA.LE.TOVER*ALPEM) TEMPA = ZERO
                          B(N) = TEMPA + TEMPA*TEMPB/ALPEM
                          IF ((B(N).EQ.ZERO) .AND. (NCALC.GT.N))
     1                       NCALC = N-1
   50                   CONTINUE
                  END IF
               END IF
            ELSE IF ((X.GT.TWOFIV) .AND. (NB.LE.MAGX+1)) THEN
C---------------------------------------------------------------------
C Asymptotic series for X .GT. 21.0.
C---------------------------------------------------------------------
               XC = SQRT(PI2/X)
               XIN = (EIGHTH/X)**2
               M = 11
               IF (X.GE.THREE5) M = 8
               IF (X.GE.ONE30) M = 4
               XM = FOUR*CONV(M)
C---------------------------------------------------------------------
C Argument reduction for SIN and COS routines.
C---------------------------------------------------------------------
               T = AINT(X/(TWOPI1+TWOPI2)+HALF)
               Z = ((X-T*TWOPI1)-T*TWOPI2) - (ALPHA+HALF)/PI2
               VSIN = SIN(Z)
               VCOS = COS(Z)
               GNU = ALPHA + ALPHA
               DO 80 I=1,2
                 S = ((XM-ONE)-GNU)*((XM-ONE)+GNU)*XIN*HALF
                 T = (GNU-(XM-THREE))*(GNU+(XM-THREE))
                 CAPP = S*T/FACT(2*M+1)
                 T1 = (GNU-(XM+ONE))*(GNU+(XM+ONE))
                 CAPQ = S*T1/FACT(2*M+2)
                 XK = XM
                 K = M + M
                 T1 = T
                 DO 70 J=2,M
                   XK = XK - FOUR
                   S = ((XK-ONE)-GNU)*((XK-ONE)+GNU)
                   T = (GNU-(XK-THREE))*(GNU+(XK-THREE))
                   CAPP = (CAPP+ONE/FACT(K-1))*S*T*XIN
                   CAPQ = (CAPQ+ONE/FACT(K))*S*T1*XIN
                   K = K - 2
                   T1 = T
   70            CONTINUE
                 CAPP = CAPP + ONE
                 CAPQ = (CAPQ+ONE)*(GNU*GNU-ONE)*(EIGHTH/X)
                 B(I) = XC*(CAPP*VCOS-CAPQ*VSIN)
                 IF (NB.EQ.1) GO TO 300
                 T = VSIN
                 VSIN = -VCOS
                 VCOS = T
                 GNU = GNU + TWO
   80         CONTINUE
C---------------------------------------------------------------------
C If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1
C---------------------------------------------------------------------
               IF (NB .GT. 2) THEN
                  GNU = ALPHA + ALPHA + TWO
                  DO 90 J=3,NB
                    B(J) = GNU*B(J-1)/X - B(J-2)
                    GNU = GNU + TWO
   90             CONTINUE
               END IF
C---------------------------------------------------------------------
C Use recurrence to generate results.  First initialize the
C calculation of P*S.
C---------------------------------------------------------------------
            ELSE
               NBMX = NB - MAGX
               N = MAGX + 1
               EN = CONV(N+N) + (ALPHA+ALPHA)
               PLAST = ONE
               P = EN/X
C---------------------------------------------------------------------
C Calculate general significance test.
C---------------------------------------------------------------------
               TEST = ENSIG + ENSIG
               IF (NBMX .GE. 3) THEN
C---------------------------------------------------------------------
C Calculate P*S until N = NB-1.  Check for possible overflow.
C---------------------------------------------------------------------
                  TOVER = ENTEN/ENSIG
                  NSTART = MAGX + 2
                  NEND = NB - 1
                  EN = CONV(NSTART+NSTART) - TWO + (ALPHA+ALPHA)
                  DO 130 K=NSTART,NEND
                     N = K
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN*PLAST/X - POLD
                     IF (P.GT.TOVER) THEN
C---------------------------------------------------------------------
C To avoid overflow, divide P*S by TOVER.  Calculate P*S until
C ABS(P) .GT. 1.
C---------------------------------------------------------------------
                        TOVER = ENTEN
                        P = P/TOVER
                        PLAST = PLAST/TOVER
                        PSAVE = P
                        PSAVEL = PLAST
                        NSTART = N + 1
  100                   N = N + 1
                           EN = EN + TWO
                           POLD = PLAST
                           PLAST = P
                           P = EN*PLAST/X - POLD
                        IF (P.LE.ONE) GO TO 100
                        TEMPB = EN/X
C---------------------------------------------------------------------
C Calculate backward test and find NCALC, the highest N such that
C the test is passed.
C---------------------------------------------------------------------
                        TEST = POLD*PLAST*(HALF-HALF/(TEMPB*TEMPB))
                        TEST = TEST/ENSIG
                        P = PLAST*TOVER
                        N = N - 1
                        EN = EN - TWO
                        NEND = MIN(NB,N)
                        DO 110 L=NSTART,NEND
                           POLD = PSAVEL
                           PSAVEL = PSAVE
                           PSAVE = EN*PSAVEL/X - POLD
                           IF (PSAVE*PSAVEL.GT.TEST) THEN
                              NCALC = L - 1
                              GO TO 190
                           END IF
  110                   CONTINUE
                        NCALC = NEND
                        GO TO 190
                     END IF
  130             CONTINUE
                  N = NEND
                  EN = CONV(N+N) + (ALPHA+ALPHA)
C---------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C---------------------------------------------------------------------
                  TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
               END IF
C---------------------------------------------------------------------
C Calculate P*S until significance test passes.
C---------------------------------------------------------------------
  140          N = N + 1
                  EN = EN + TWO
                  POLD = PLAST
                  PLAST = P
                  P = EN*PLAST/X - POLD
               IF (P.LT.TEST) GO TO 140
C---------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C---------------------------------------------------------------------
  190          N = N + 1
               EN = EN + TWO
               TEMPB = ZERO
               TEMPA = ONE/P
               M = 2*N - 4*(N/2)
               SUM = ZERO
               EM = CONV(N/2)
               ALPEM = (EM-ONE) + ALPHA
               ALP2EM = (EM+EM) + ALPHA
               IF (M .NE. 0) SUM = TEMPA*ALPEM*ALP2EM/EM
               NEND = N - NB
               IF (NEND .GT. 0) THEN
C---------------------------------------------------------------------
C Recur backward via difference equation, calculating (but not
C storing) B(N), until N = NB.
C---------------------------------------------------------------------
                  DO 200 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     TEMPC = TEMPB
                     TEMPB = TEMPA
                     TEMPA = (EN*TEMPB)/X - TEMPC
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        IF (N.EQ.1) GO TO 210
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+TEMPA*ALP2EM)*ALPEM/EM
                     END IF
  200             CONTINUE
               END IF
C---------------------------------------------------------------------
C Store B(NB).
C---------------------------------------------------------------------
  210          B(N) = TEMPA
               IF (NEND .GE. 0) THEN
                  IF (NB .LE. 1) THEN
                        ALP2EM = ALPHA
                        IF ((ALPHA+ONE).EQ.ONE) ALP2EM = ONE
                        SUM = SUM + B(1)*ALP2EM
                        GO TO 250
                     ELSE
C---------------------------------------------------------------------
C Calculate and store B(NB-1).
C---------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*TEMPA)/X - TEMPB
                        IF (N.EQ.1) GO TO 240
                        M = 2 - M
                        IF (M .NE. 0) THEN
                           EM = EM - ONE
                           ALP2EM = (EM+EM) + ALPHA
                           ALPEM = (EM-ONE) + ALPHA
                           IF (ALPEM.EQ.ZERO) ALPEM = ONE
                           SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                        END IF
                  END IF
               END IF
               NEND = N - 2
               IF (NEND .NE. 0) THEN
C---------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C---------------------------------------------------------------------
                  DO 230 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     B(N) = (EN*B(N+1))/X - B(N+2)
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                     END IF
  230             CONTINUE
               END IF
C---------------------------------------------------------------------
C Calculate B(1).
C---------------------------------------------------------------------
               B(1) = TWO*(ALPHA+ONE)*B(2)/X - B(3)
  240          EM = EM - ONE
               ALP2EM = (EM+EM) + ALPHA
               IF (ALP2EM.EQ.ZERO) ALP2EM = ONE
               SUM = SUM + B(1)*ALP2EM
C---------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C---------------------------------------------------------------------
  250          IF ((ALPHA+ONE).NE.ONE)
     1              SUM = SUM*FUNC(ALPHA)*(X*HALF)**(-ALPHA)
               TEMPA = ENMTEN
               IF (SUM.GT.ONE) TEMPA = TEMPA*SUM
               DO 260 N=1,NB
                 IF (ABS(B(N)).LT.TEMPA) B(N) = ZERO
                 B(N) = B(N)/SUM
  260          CONTINUE
            END IF
C---------------------------------------------------------------------
C Error return -- X, NB, or ALPHA is out of range.
C---------------------------------------------------------------------
         ELSE
            B(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
C---------------------------------------------------------------------
C Exit
C---------------------------------------------------------------------
  300 RETURN
C ---------- Last line of RJBESL ----------
      END

CS    REAL FUNCTION GAMMA(X)
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C----------------------------------------------------------------------
C
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C
C  Intrinsic functions required are:
C
C     INT, DBLE, EXP, LOG, REAL, SIN
C
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: October 12, 1989
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
CS    REAL 
      DOUBLE PRECISION 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
CS   1     SQRTPI/0.9189385332046727417803297E0/,
CS   2     PI/3.1415926535897932384626434E0/
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
CS   1     XINF/3.4E38/
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C----------------------------------------------------------------------
CS    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
CS   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
CS   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
CS   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
CS    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
CS   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
CS   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
CS   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
C  Statement functions for conversion between integer and float
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C----------------------------------------------------------------------
C  Argument is negative
C----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Argument is positive
C----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
C----------------------------------------------------------------------
C  Argument .LT. EPS
C----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C----------------------------------------------------------------------
C  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
C----------------------------------------------------------------------
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C----------------------------------------------------------------------
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
CS900 GAMMA = RES
  900 DGAMMA = RES
      RETURN
C ---------- Last line of GAMMA ----------
      END


      SUBROUTINE RYBESL(X,ALPHA,NB,BY,NCALC)
C----------------------------------------------------------------------
C
C  This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision non-negative real argument for which
C         Y's are to be calculated.
C ALPHA - Working precision fractional part of order for which
C         Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the 
C         last is of order (NB - 1 + ALPHA).
C BY    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BY
C         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),
C         If (0 .LT. NCALC .LT. NB), BY(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BY, the user should check that 
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   p      = Number of significant base-beta digits in the
C            significand of a floating-point number
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   EPS    = beta ** (-p)
C   DEL    = Machine number below which sin(x)/x = 1; approximately
C            SQRT(EPS).
C   XMIN   = Smallest acceptable argument for RBESY; approximately
C            max(2*beta**minexp,2/XINF), rounded up
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   THRESH = Lower bound for use of the asymptotic form; approximately
C            AINT(-LOG10(EPS/2.0))+1.0
C   XLARGE = Upper bound on X; approximately 1/DEL, because the sine
C            and cosine functions have lost about half of their 
C            precision at that point.
C
C
C     Approximate values for some important machines are:
C
C                        beta    p     minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15
C  Cyber 180/185 
C    under NOS   (S.P.)    2    48      -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16
C  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17
C  VAX           (S.P.)    2    24      -128         127    5.96E-8
C  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17
C  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16
C
C
C                         DEL      XMIN      XINF     THRESH  XLARGE
C
C CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7
C Cyber 180/855
C   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8
C IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8
C VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4
C VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9
C VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all Y's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. -1:  An argument is out of range. For example,
C       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
C       XMAX.  In this case, BY(1) = 0.0, the remainder of the
C       BY-vector is not calculated, and NCALC is set to
C       MIN0(NB,0)-2  so that NCALC .NE. NB.
C  NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function
C       values are set to 0.0.
C  1 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BY(I) contains correct function
C       values for I .LE. NCALC, and and the remaining NB-NCALC
C       array elements contain 0.0.
C
C
C Intrinsic functions required are:
C
C     DBLE, EXP, INT, MAX, MIN, REAL, SQRT
C
C
C Acknowledgement
C
C  This program draws heavily on Temme's Algol program for Y(a,x)
C  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
C  scheme is used for  x < THRESH, and Campbell's scheme is used
C  in the asymptotic region.  Segments of code from both sources
C  have been translated into Fortran 77, merged, and heavily modified.
C  Modifications include parameterization of machine dependencies,
C  use of a new approximation for ln(gamma(x)), and built-in
C  protection against over/underflow.
C
C References: "Bessel functions J_nu(x) and Y_nu(x) of real
C              order and real argument," Campbell, J. B.,
C              Comp. Phy. Comm. 18, 1979, pp. 133-142.
C
C             "On the numerical evaluation of the ordinary
C              Bessel function of the second kind," Temme,
C              N. M., J. Comput. Phys. 21, 1976, pp. 343-350.
C
C  Latest modification: March 19, 1990
C
C  Modified by: W. J. Cody
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C----------------------------------------------------------------------
      INTEGER I,K,NA,NB,NCALC
CS    REAL
      DOUBLE PRECISION
     1  ALFA,ALPHA,AYE,B,BY,C,CH,COSMU,D,DEL,DEN,DDIV,DIV,DMU,D1,D2,
     2  E,EIGHT,EN,ENU,EN1,EPS,EVEN,EX,F,FIVPI,G,GAMMA,H,HALF,ODD,
     3  ONBPI,ONE,ONE5,P,PA,PA1,PI,PIBY2,PIM5,Q,QA,QA1,Q0,R,S,SINMU,
     4  SQ2BPI,TEN9,TERM,THREE,THRESH,TWO,TWOBYX,X,XINF,XLARGE,XMIN,
     5  XNA,X2,YA,YA1,ZERO
      DIMENSION BY(NB),CH(21)
C----------------------------------------------------------------------
C  Mathematical constants
C    FIVPI = 5*PI
C    PIM5 = 5*PI - 15
C    ONBPI = 1/PI
C    PIBY2 = PI/2
C    SQ2BPI = SQUARE ROOT OF 2/PI
C----------------------------------------------------------------------
CS    DATA ZERO,HALF,ONE,TWO,THREE/0.0E0,0.5E0,1.0E0,2.0E0,3.0E0/
CS    DATA EIGHT,ONE5,TEN9/8.0E0,15.0E0,1.9E1/
CS    DATA FIVPI,PIBY2/1.5707963267948966192E1,1.5707963267948966192E0/
CS    DATA PI,SQ2BPI/3.1415926535897932385E0,7.9788456080286535588E-1/
CS    DATA PIM5,ONBPI/7.0796326794896619231E-1,3.1830988618379067154E-1/
      DATA ZERO,HALF,ONE,TWO,THREE/0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
      DATA EIGHT,ONE5,TEN9/8.0D0,15.0D0,1.9D1/
      DATA FIVPI,PIBY2/1.5707963267948966192D1,1.5707963267948966192D0/
      DATA PI,SQ2BPI/3.1415926535897932385D0,7.9788456080286535588D-1/
      DATA PIM5,ONBPI/7.0796326794896619231D-1,3.1830988618379067154D-1/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA DEL,XMIN,XINF,EPS/1.0E-4,2.36E-38,3.40E38,5.96E-8/
CS    DATA THRESH,XLARGE/8.0E0,1.0E4/
      DATA DEL,XMIN,XINF,EPS/1.0D-8,4.46D-308,1.79D308,1.11D-16/
      DATA THRESH,XLARGE/16.0D0,1.0D8/
C----------------------------------------------------------------------
C  Coefficients for Chebyshev polynomial expansion of 
C         1/gamma(1-x), abs(x) .le. .5
C----------------------------------------------------------------------
CS    DATA CH/-0.67735241822398840964E-23,-0.61455180116049879894E-22,
CS   1         0.29017595056104745456E-20, 0.13639417919073099464E-18,
CS   2         0.23826220476859635824E-17,-0.90642907957550702534E-17,
CS   3        -0.14943667065169001769E-14,-0.33919078305362211264E-13,
CS   4        -0.17023776642512729175E-12, 0.91609750938768647911E-11,
CS   5         0.24230957900482704055E-09, 0.17451364971382984243E-08,
CS   6        -0.33126119768180852711E-07,-0.86592079961391259661E-06,
CS   7        -0.49717367041957398581E-05, 0.76309597585908126618E-04,
CS   8         0.12719271366545622927E-02, 0.17063050710955562222E-02,
CS   9        -0.76852840844786673690E-01,-0.28387654227602353814E+00,
CS   A         0.92187029365045265648E+00/
      DATA CH/-0.67735241822398840964D-23,-0.61455180116049879894D-22,
     1         0.29017595056104745456D-20, 0.13639417919073099464D-18,
     2         0.23826220476859635824D-17,-0.90642907957550702534D-17,
     3        -0.14943667065169001769D-14,-0.33919078305362211264D-13,
     4        -0.17023776642512729175D-12, 0.91609750938768647911D-11,
     5         0.24230957900482704055D-09, 0.17451364971382984243D-08,
     6        -0.33126119768180852711D-07,-0.86592079961391259661D-06,
     7        -0.49717367041957398581D-05, 0.76309597585908126618D-04,
     8         0.12719271366545622927D-02, 0.17063050710955562222D-02,
     9        -0.76852840844786673690D-01,-0.28387654227602353814D+00,
     A         0.92187029365045265648D+00/
C----------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      IF ((NB .GT. 0) .AND. (X .GE. XMIN) .AND. (EX .LT. XLARGE)
     1       .AND. (ENU .GE. ZERO) .AND. (ENU .LT. ONE))  THEN
            XNA = AINT(ENU+HALF)
            NA = INT(XNA)
            IF (NA .EQ. 1) ENU = ENU - XNA
            IF (ENU .EQ. -HALF) THEN
                  P = SQ2BPI/SQRT(EX)
                  YA = P * SIN(EX)
                  YA1 = -P * COS(EX)
               ELSE IF (EX .LT. THREE) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for small X
C----------------------------------------------------------------------
                  B = EX * HALF
                  D = -LOG(B)
                  F = ENU * D
                  E = B**(-ENU)
                  IF (ABS(ENU) .LT. DEL) THEN
                        C = ONBPI
                     ELSE
                        C = ENU / SIN(ENU*PI)
                  END IF
C----------------------------------------------------------------------
C  Computation of sinh(f)/f
C----------------------------------------------------------------------
                  IF (ABS(F) .LT. ONE) THEN
                        X2 = F*F
                        EN = TEN9
                        S = ONE
                        DO 80 I = 1, 9
                           S = S*X2/EN/(EN-ONE)+ONE
                           EN = EN - TWO
   80                   CONTINUE
                     ELSE 
                        S = (E - ONE/E) * HALF / F
                  END IF
C----------------------------------------------------------------------
C  Computation of 1/gamma(1-a) using Chebyshev polynomials
C----------------------------------------------------------------------
                  X2 = ENU*ENU*EIGHT
                  AYE = CH(1)
                  EVEN = ZERO
                  ALFA = CH(2)
                  ODD = ZERO
                  DO 40 I = 3, 19, 2
                     EVEN = -(AYE+AYE+EVEN)
                     AYE = -EVEN*X2 - AYE + CH(I)
                     ODD = -(ALFA+ALFA+ODD)
                     ALFA = -ODD*X2 - ALFA + CH(I+1)
   40             CONTINUE
                  EVEN = (EVEN*HALF+AYE)*X2 - AYE + CH(21)
                  ODD = (ODD+ALFA)*TWO
                  GAMMA = ODD*ENU + EVEN
C----------------------------------------------------------------------
C  End of computation of 1/gamma(1-a)
C----------------------------------------------------------------------
                  G = E * GAMMA
                  E = (E + ONE/E) * HALF
                  F = TWO*C*(ODD*E+EVEN*S*D)
                  E = ENU*ENU
                  P = G*C
                  Q = ONBPI / G
                  C = ENU*PIBY2
                  IF (ABS(C) .LT. DEL) THEN
                        R = ONE
                     ELSE 
                        R = SIN(C)/C
                  END IF
                  R = PI*C*R*R
                  C = ONE
                  D = - B*B
                  H = ZERO
                  YA = F + R*Q
                  YA1 = P
                  EN = ZERO
  100             EN = EN + ONE
                  IF (ABS(G/(ONE+ABS(YA)))
     1                      + ABS(H/(ONE+ABS(YA1))) .GT. EPS) THEN
                        F = (F*EN+P+Q)/(EN*EN-E)
                        C = C * D/EN
                        P = P/(EN-ENU)
                        Q = Q/(EN+ENU)
                        G = C*(F+R*Q)
                        H = C*P - EN*G
                        YA = YA + G
                        YA1 = YA1+H
                        GO TO 100
                  END IF
                  YA = -YA
                  YA1 = -YA1/B
               ELSE IF (EX .LT. THRESH) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for moderate X
C----------------------------------------------------------------------
                  C = (HALF-ENU)*(HALF+ENU)
                  B = EX + EX
                  E = (EX*ONBPI*COS(ENU*PI)/EPS)
                  E = E*E
                  P = ONE
                  Q = -EX
                  R = ONE + EX*EX
                  S = R
                  EN = TWO
  200             IF (R*EN*EN .LT. E) THEN
                        EN1 = EN+ONE
                        D = (EN-ONE+C/EN)/S
                        P = (EN+EN-P*D)/EN1
                        Q = (-B+Q*D)/EN1
                        S = P*P + Q*Q
                        R = R*S
                        EN = EN1
                        GO TO 200
                  END IF
                  F = P/S
                  P = F
                  G = -Q/S
                  Q = G
  220             EN = EN - ONE  
                  IF (EN .GT. ZERO) THEN
                        R = EN1*(TWO-P)-TWO
                        S = B + EN1*Q
                        D = (EN-ONE+C/EN)/(R*R+S*S)
                        P = D*R
                        Q = D*S
                        E = F + ONE
                        F = P*E - G*Q
                        G = Q*E + P*G
                        EN1 = EN
                        GO TO 220
                  END IF
                  F = ONE + F
                  D = F*F + G*G
                  PA = F/D
                  QA = -G/D
                  D = ENU + HALF -P
                  Q = Q + EX
                  PA1 = (PA*Q-QA*D)/EX
                  QA1 = (QA*Q+PA*D)/EX
                  B = EX - PIBY2*(ENU+HALF)
                  C = COS(B)
                  S = SIN(B)
                  D = SQ2BPI/SQRT(EX)
                  YA = D*(PA*S+QA*C)
                  YA1 = D*(QA1*S-PA1*C)
               ELSE
C----------------------------------------------------------------------
C  Use Campbell's asymptotic scheme.
C----------------------------------------------------------------------
                  NA = 0
                  D1 = AINT(EX/FIVPI)
                  I = INT(D1)
                  DMU = ((EX-ONE5*D1)-D1*PIM5)-(ALPHA+HALF)*PIBY2
                  IF (I-2*(I/2) .EQ. 0) THEN
                        COSMU = COS(DMU)
                        SINMU = SIN(DMU)
                     ELSE
                        COSMU = -COS(DMU)
                        SINMU = -SIN(DMU)
                  END IF
                  DDIV = EIGHT * EX
                  DMU = ALPHA
                  DEN = SQRT(EX)
                  DO 350 K = 1, 2
                     P = COSMU
                     COSMU = SINMU
                     SINMU = -P
                     D1 = (TWO*DMU-ONE)*(TWO*DMU+ONE)
                     D2 = ZERO
                     DIV = DDIV
                     P = ZERO
                     Q = ZERO
                     Q0 = D1/DIV
                     TERM = Q0
                     DO 310 I = 2, 20
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = -TERM*D1/DIV
                        P = P + TERM
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = TERM*D1/DIV
                        Q = Q + TERM
                        IF (ABS(TERM) .LE. EPS) GO TO 320
  310                CONTINUE
  320                P = P + ONE
                     Q = Q + Q0
                     IF (K .EQ. 1) THEN
                           YA = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                        ELSE
                           YA1 = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                     END IF
                     DMU = DMU + ONE
  350             CONTINUE
            END IF
            IF (NA .EQ. 1) THEN
               H = TWO*(ENU+ONE)/EX
               IF (H .GT. ONE) THEN
                  IF (ABS(YA1) .GT. XINF/H) THEN
                     H = ZERO
                     YA = ZERO
                  END IF
               END IF
               H = H*YA1 - YA
               YA = YA1
               YA1 = H
            END IF
C----------------------------------------------------------------------
C  Now have first one or two Y's
C----------------------------------------------------------------------
            BY(1) = YA
            BY(2) = YA1
            IF (YA1 .EQ. ZERO) THEN
                  NCALC = 1
               ELSE
                  AYE = ONE + ALPHA
                  TWOBYX = TWO/EX
                  NCALC = 2
                  DO 400 I = 3, NB
                     IF (TWOBYX .LT. ONE) THEN
                           IF (ABS(BY(I-1))*TWOBYX .GE. XINF/AYE)
     1                                                     GO TO 450
                        ELSE
                           IF (ABS(BY(I-1)) .GE. XINF/AYE/TWOBYX )
     1                                                     GO TO 450
                     END IF
                     BY(I) = TWOBYX*AYE*BY(I-1) - BY(I-2) 
                     AYE = AYE + ONE
                     NCALC = NCALC + 1
  400             CONTINUE
            END IF
  450       DO 460 I = NCALC+1, NB
               BY(I) = ZERO
  460       CONTINUE
         ELSE
            BY(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
  900 RETURN
C---------- Last line of RYBESL ----------
      END



      subroutine potential(Numchannels,NumDataPoints,NumTotStates,
     >    PFile,QFile,xdata,VQmat,Pmat)

      integer Numchannels,NumDataPoints,NumTotStates

      double precision VQmat(NumDataPoints,Numchannels,Numchannels)
      double precision xdata(NumDataPoints)
      double precision Pmat(NumDataPoints,NumChannels,NumChannels)
      character*64 PFile,QFile

      open(200,file=PFile,status='old')
      open(300,file=QFile,status='old')


      Pmat = 0.d0
      VQmat = 0.d0

      do i=1,NumDataPoints
       read(200,*) xdata(i)
       read(300,*)
       do j=1,Numchannels
        read(200,3001)(Pmat(i,j,k),k=1,Numchannels)
        read(300,3001)(VQmat(i,j,k),k=1,Numchannels)
       enddo
       do j=Numchannels+1,NumTotStates
        read(200,3001)
        read(300,3001)
       enddo
      enddo


cc      XXX
c       do i = 1,NumDataPoints
c       do j = 1,NumChannels
c       do k = j+1,NumChannels
c       j=1; k=3;
c       VQmat(i,j,k) = 0.d0
c       VQmat(i,k,j) = 0.d0
c       Pmat(i,j,k) = 0.d0
c       Pmat(i,k,j) = 0.d0
cc       VQmat(i,j,k) = VQmat(i,k,j) 
cc        Pmat(i,j,k) = -Pmat(i,k,j) 
c       enddo
c       enddo
c       enddo

cc     Filtering noise (jpdincao)
c      do i=1,NumDataPoints
c       do j=1,Numchannels
c       do k=1,Numchannels
c           if (dabs(Pmat(i,j,k)).le.1.d-14) Pmat(i,j,k) = 0.d0 !sign(1.d-14,Pmat(i,j,k))
c           if (dabs(VQmat(i,j,k)).le.1.d-14) VQmat(i,j,k) = 0.d0 !sign(1.d-14,VQmat(i,j,k))
c       enddo
c       enddo
c      enddo


      close(200)
      close(300)
c3001 format(100(e20.12,1x))
 3001 format(100(e18.12,1x))
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function u(n,x)

      integer n
      double precision x

      if (n .eq. 1) u = x**2-1.25d0*x**3-0.5d0*x**4+0.75d0*x**5
      if (n .eq. 2) u = (x**2-x**3-x**4+x**5)/4.0d0
      if (n .eq. 3) u = 1-2.0d0*x**2+x**4
      if (n .eq. 4) u = x-2.0d0*x**3+x**5
      if (n .eq. 5) u = x**2+1.25d0*x**3-0.5d0*x**4-0.75d0*x**5
      if (n .eq. 6) u = (-x**2-x**3+x**4+x**5)/4.0d0

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BasicOverlapMatrix(BasicS)

      double precision BasicS(6,6)

      BasicS(1,1) = 1046.0d0/3465.0d0
      BasicS(1,2) = 38.0d0/1155.0d0
      BasicS(1,3) = 8.0d0/63.0d0
      BasicS(1,4) = -32.0d0/693.0d0
      BasicS(1,5) = 131.0d0/3465.0d0
      BasicS(1,6) = -29.0d0/3465.0d0

      BasicS(2,2) = 16.0d0/3465.0d0
      BasicS(2,3) = 8.0d0/315.0d0
      BasicS(2,4) = -8.0d0/1155.0d0
      BasicS(2,5) = 29.0d0/3465.0d0
      BasicS(2,6) = -2.0d0/1155.0d0

      BasicS(3,3) = 256.0d0/315.0d0
      BasicS(3,4) = 0.0d0
      BasicS(3,5) = 8.0d0/63.0d0
      BasicS(3,6) = -8.0d0/315.0d0

      BasicS(4,4) = 256.0d0/3465.0d0
      BasicS(4,5) = 32.0d0/693.0d0
      BasicS(4,6) = -8.0d0/1155.0d0

      BasicS(5,5) = 1046.0d0/3465.0d0
      BasicS(5,6) = -38.0d0/1155.0d0

      BasicS(6,6) = 16.0d0/3465.0d0

      do i = 2,6
       do j = 1,i-1
        BasicS(i,j) = BasicS(j,i)
       enddo
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BasicKineticMatrix(BasicT)

c -1/2 factor is included

      double precision BasicT(6,6)

      BasicT(1,1) = 139.0d0/210.0d0
      BasicT(1,2) = 223.0d0/420.0d0
      BasicT(1,3) = -64.0d0/105.0d0
      BasicT(1,4) = 4.0d0/21.0d0
      BasicT(1,5) = -11.0d0/210.0d0
      BasicT(1,6) = -1.0d0/140.0d0

      BasicT(2,2) = 2.0d0/45.0d0
      BasicT(2,3) = -4.0d0/105.0d0
      BasicT(2,4) = -4.0d0/315.0d0
      BasicT(2,5) = 1.0d0/140.0d0
      BasicT(2,6) = -1.0d0/126.0d0

      BasicT(3,3) = 128.0d0/105.0d0
      BasicT(3,4) = 0.0d0
      BasicT(3,5) = -64.0d0/105.0d0
      BasicT(3,6) = 4.0d0/105.0d0

      BasicT(4,4) = 128.0d0/315.0d0
      BasicT(4,5) = -4.0d0/21.0d0
      BasicT(4,6) = -4.0d0/315.0d0

      BasicT(5,5) = 139.0d0/210.0d0
      BasicT(5,6) = -223.0d0/420.0d0

      BasicT(6,6) = 2.0d0/45.0d0

      do i = 2,6
       do j = 1,i-1
        BasicT(i,j) = BasicT(j,i)
       enddo
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      subroutine GenMEMap(Numchannels,NumSectors,MEMap)

      integer NumSectors,Numchannels
      integer MEMap(Numchannels*Numsectors,6)

c form matrix element map

      icnt = 1
      do k=1,Numsectors
       do j=1,6
        do i=1,NumChannels

         index = Numchannels*(k-1) + i

         if (k .eq. 1 .and. j .eq. 1)then
          MEMap(index,j) = 0
         elseif (k .eq. Numsectors .and. j .eq. 5)then
          MEMap(index,j) = 0
         else
          MEMap(index,j) = icnt
          icnt = icnt + 1
         endif

         enddo
        enddo
        icnt = icnt - 2*Numchannels
       enddo

1901  format(6(i4,1x))

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcClosed(Numchannels,NumSectors,MatrixDim,kl,iband,Map,Mass,energy,
     >   xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,PFit,VQFit,
     >   WhereStart,PExponents,VQExponents,Gcc)

      integer NumSectors,Numchannels,MatrixDim,kl,iband,LegPoints
      integer MEMap(NumSectors*NumChannels,6)
      double precision Mass,xSector(*),xLeg(*),wLeg(*),energy
      double precision Gcc(iband,MatrixDim)

      double precision Scale(6),IntScale,BasicT(6,6)
      double precision BasicVQ(6,6),BasicS(6,6),BasicP(6,6)
      double precision ScaledZero,x(LegPoints),xL(LegPoints),xR(LegPoints)
      double precision u(6,LegPoints),VQInterp(LegPoints)
      double precision up(6,LegPoints),PInterp(LegPoints)
      double precision dPInterp(LegPoints),PLInterp(LegPoints)
      double precision PRInterp(LegPoints)
      double precision VQmat(NumDataPoints,Numchannels,Numchannels)
      double precision Pmat(NumDataPoints,Numchannels,Numchannels)
      double precision xdata(NumDataPoints),ydata(NumDataPoints)

      double precision PFit(Numchannels,Numchannels,3),VQFit(Numchannels,Numchannels,3)
      double precision WhereStart,vvj,ppj
      double precision PExponents(Numchannels,Numchannels,3)
      double precision VQExponents(Numchannels,Numchannels,3)

      call BasicKineticMatrix(BasicT)

      call CalcLocalBasis(LegPoints,xLeg,u,up)

      call BasicOverlapMatrix(BasicS)

      ku = kl
      kb = kl + ku + 1

      Gcc = 0.0d0
      Scale = 1.0d0

      knc = 1
      kjc = 1

c index one sector at a time
      do n = 1,NumSectors

       IntScale = (xSector(n+1)-xSector(n))/2.0d0
       ScaledZero = (xSector(n+1)+xSector(n))/2.0d0

       do nc = 1,Numchannels

c integrate potentials in each block

        do jc = 1,Numchannels

        do ik=1,NumDataPoints
         ydata(ik)=VQmat(ik,nc,jc)
        enddo

        do j = 1,LegPoints
         x(j) = IntScale*xLeg(j)+ScaledZero

        IF(x(j) .ge. WhereStart) THEN

	   vvj=VQFit(nc,jc,1)/x(j)**VQExponents(nc,jc,1)
	   vvj=vvj+VQFit(nc,jc,2)/x(j)**VQExponents(nc,jc,2)
	   vvj=vvj+VQFit(nc,jc,3)/x(j)**VQExponents(nc,jc,3)
	   VQInterp(j)=vvj

	   ppj=PFit(nc,jc,1)/x(j)**PExponents(nc,jc,1)
	   ppj=ppj+PFit(nc,jc,2)/x(j)**PExponents(nc,jc,2)
	   ppj=ppj+PFit(nc,jc,3)/x(j)**PExponents(nc,jc,3)
	   PInterp(j)=ppj

         endif
        enddo

	IF(x(1) .LT. WhereStart) THEN
         call Akima(NumDataPoints,xdata,ydata,LegPoints,x,VQInterp)
	ENDIF

cc       Filtering noise (jpdincao)
c        do j = 1,LegPoints
c           if (dabs(VQInterp(j)).le.1.d-14) VQInterp(j) = 0.d0 !sign(1.d-14,VQInterp(j))
c        enddo

cc       Cutoff Long-range Q-coupling (jpdincao)
c        if (nc.eq.1.or.nc.eq.3) then
c        if (jc.eq.1.or.jc.eq.3) then
c        do j = 1,LegPoints
c           if (x(j).ge.WhereStart) VQInterp(j) = 0.d0 
c        enddo
c        endif
c        endif

       do i = 1,6
        do ip = 1,6
         BasicVQ(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicVQ(i,ip) = BasicVQ(i,ip) + wLeg(j)*u(i,j)*VQInterp(j)*u(ip,j)
         enddo
        enddo
       enddo

       do ik=1,NumDataPoints
        ydata(ik)=Pmat(ik,nc,jc)
       enddo
       
       IF(x(1) .LT. WhereStart) THEN
        call akima(NumDataPoints,xdata,ydata,LegPoints,x,Pinterp)
       ENDIF

cc       Filtering noise (jpdincao)
c        do j = 1,LegPoints
c           if (dabs(PInterp(j)).le.1.d-14) PInterp(j) = 0.d0 !sign(1.d-14,PInterp(j))
c        enddo   

cc       Cutoff Long-range P-coupling (jpdincao)
c        if (nc.eq.1.or.nc.eq.3) then
c        if (jc.eq.1.or.jc.eq.3) then
c        do j = 1,LegPoints
c           if (x(j).ge.WhereStart) PInterp(j) = 0.d0 
c        enddo
c        endif
c        endif

       do i = 1,6
        do ip = 1,6
         BasicP(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicP(i,ip) = BasicP(i,ip) + wLeg(j)*PInterp(j)*
     >        (u(i,j)*up(ip,j) - up(i,j)*u(ip,j))
         enddo
        enddo
       enddo

       if (n .ne. 1) then
        Scale(2) = (xSector(n+1)-xSector(n))/(xSector(n)-xSector(n-1))
       else
        Scale(2) = 1.0d0
       endif

       do i = 1,6
        do ip = 1,6

        if ((MEMap(knc,i) .ne. 0) .AND. (MEMap(kjc,ip) .ne. 0)) then

         imap = MEMap(knc,i)
         jmap = MEMap(kjc,ip)

         if (knc .eq. kjc)then

c diagonal blocks

         Gcc(kb+imap-jmap,jmap) = Gcc(kb+imap-jmap,jmap) + 
     >         2.0d0*mass*Scale(i)*Scale(ip)*
     >        (energy*Intscale*BasicS(i,ip) -
     >        (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >         + BasicP(i,ip)))

          else 

c off-diagonal blocks

         Gcc(kb+imap-jmap,jmap) = Gcc(kb+imap-jmap,jmap) - 
     >        2.0d0*mass*Scale(i)*Scale(ip)*(BasicVQ(i,ip)*IntScale
     >        + BasicP(i,ip))

          endif
         endif
        enddo
       enddo

        kjc = kjc + 1
        enddo

        kjc = Numchannels*(n-1) + 1
        knc = knc + 1
       enddo

       knc = Numchannels*n + 1
       kjc = knc
      enddo

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcLocalBasis(LegPoints,xLeg,u,up)

      integer LegPoints
      double precision xLeg(LegPoints),u(6,LegPoints)
      double precision up(6,LegPoints)

      integer n
      double precision x,x2,x3,x4,x5

      do n = 1,LegPoints

       x = xLeg(n)
       x2 = x*x
       x3 = x2*x
       x4 = x3*x
       x5 = x4*x

       u(1,n) = x2 - 1.25d0*x3 - 0.5d0*x4 + 0.75d0*x5
       up(1,n) = 2.0d0*x - 3.750*x2 - 2*x3 + 3.750d0*x4
c       upp(1,n) = 2.0d0 - 7.5d0*x - 6.0d0*x2 + 15.0d0*x3

       u(2,n) = 0.25d0*(x2 - x3 - x4 + x5)
       up(2,n) = 0.5d0*x - 0.75d0*x2 - x3 + 1.25d0*x4
c       upp(2,n) = 0.5d0 - 1.5d0*x - 3.0d0*x2 + 5.0d0*x3

       u(3,n) = 1.0d0 - 2.0d0*x2 + x4
       up(3,n) =  -4.0d0*x + 4.0d0*x3
c       upp(3,n) =  -4.0d0 + 12.0d0*x2

       u(4,n) = x - 2.0d0*x3 + x5
       up(4,n) = 1.0d0 - 6.0d0*x2 + 5.0d0*x4
c       upp(4,n) =  -12.0d0*x + 20.0d0*x3

       u(5,n) = x2 + 1.25d0*x3 - 0.5d0*x4 - 0.75d0*x5
       up(5,n) = 2.0d0*x + 3.75d0*x2 - 2.0d0*x3 - 3.75d0*x4
c       upp(5,n) = 2.0d0 + 7.5d0*x - 6.0d0*x2 - 15.0d0*x3

       u(6,n) = 0.25d0*( -x2 - x3 + x4 + x5)
       up(6,n) =  -0.5d0*x - 0.75d0*x2 + x3 + 1.25d0*x4
c       upp(6,n) =  -0.5d0 - 1.5d0*x + 3.0d0*x2 + 5.0d0*x3

      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcOpen(NumOpenL,NumOpenR,Numchannels,NumSectors,MatrixDim,MEMap,Mass,
     >  energy,xsector,LegPoints,xLeg,wLeg,NumDataPoints,xdata,VQmat,Pmat,
     >  PFit,VQFit,WhereStart,PExponents,VQExponents,Gco,Goo)

      integer NumSectors,NumOpenL,NumOpenR,Numchannels,LegPoints,NumDataPoints
      integer MEMap(NumSectors*Numchannels,6),MatrixDim
      double precision Mass,xSector(*),energy
      double precision Gco(MatrixDim,NumOpenL+NumOpenR)
      double precision Goo(NumOpenL+NumOpenR,NumOpenL+NumOpenR)
      double precision xLeg(LegPoints),wLeg(LegPoints)

      double precision Scale(6),IntScale,BasicT(6,6)
      double precision BasicP(6,6),BasicVQ(6,6),BasicS(6,6)
      double precision ScaledZero,x(LegPoints),up(6,LegPoints)
      double precision u(6,LegPoints),VQInterp(LegPoints)
      double precision PInterp(LegPoints)
      double precision VQmat(NumDataPoints,NumChannels,NumChannels)
      double precision Pmat(NumDataPoints,NumChannels,NumChannels)
      double precision xdata(NumDataPoints)
      double precision ydata(NumDataPoints)

      double precision PFit(NumChannels,NumChannels,3),VQFit(NumChannels,NumChannels,3)
      double precision WhereStart,vvj,ppj,Dpj	
      double precision PExponents(NumChannels,NumChannels,3)
      double precision VQExponents(NumChannels,NumChannels,3)

      call BasicKineticMatrix(BasicT)

      call CalcLocalBasis(LegPoints,xLeg,u,up)

      call BasicOverlapMatrix(BasicS)

      NumOpen = NumOpenL + NumOpenR

      Gco = 0.0d0
      Goo = 0.0d0

      Scale = 1.0d0

c open channels on left boundary

      ip = 1
      n  = 1

      do nc = 1,NumOpenL
        IntScale = (xSector(n+1)-xSector(n))/2.0d0
        ScaledZero = (xSector(n+1)+xSector(n))/2.0d0

c integrate potentials in each block

        do jc = 1,NumChannels

        do ik=1,NumDataPoints
         ydata(ik) = VQmat(ik,nc,jc)
        enddo

        do j = 1,LegPoints
           x(j) = IntScale*xLeg(j)+ScaledZero

	   vvj=VQFit(nc,jc,1)/x(j)**VQExponents(nc,jc,1)
	   vvj=vvj+VQFit(nc,jc,2)/x(j)**VQExponents(nc,jc,2)
	   vvj=vvj+VQFit(nc,jc,3)/x(j)**VQExponents(nc,jc,3)
	   VQInterp(j)=vvj

	   ppj=PFit(nc,jc,1)/x(j)**PExponents(nc,jc,1)
	   ppj=ppj+PFit(nc,jc,2)/x(j)**PExponents(nc,jc,2)
	   ppj=ppj+PFit(nc,jc,3)/x(j)**PExponents(nc,jc,3)
	   PInterp(j)=ppj

        enddo

	 IF(x(1) .LT. WhereStart) THEN
         call Akima(NumDataPoints,xdata,ydata,LegPoints,x,VQInterp)  ! do the interpolation for this particular matrix element only (nc,jc)? (npmehta)
	 ENDIF

cc       Filtering noise (jpdincao)
c        do j = 1,LegPoints
c           if (dabs(VQInterp(j)).le.1.d-14) VQInterp(j) = 0.d0 !sign(1.d-14,VQInterp(j))
c        enddo   

cc       Cutoff Long-range Q-coupling (jpdincao)
c        if (nc.eq.1.or.nc.eq.3) then
c        if (jc.eq.1.or.jc.eq.3) then
c        do j = 1,LegPoints
c           if (x(j).ge.WhereStart) VQInterp(j) = 0.d0 
c        enddo
c        endif
c        endif

       do i = 1,6
         BasicVQ(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicVQ(i,ip) = BasicVQ(i,ip) + wLeg(j)*u(i,j)*VQInterp(j)*u(ip,j)
         enddo
       enddo

        do ik=1,NumDataPoints
         ydata(ik)=Pmat(ik,nc,jc)
        enddo

	 IF(x(1) .LT. WhereStart) THEN
         call akima(NumDataPoints,xdata,ydata,LegPoints,x,Pinterp)
	 ENDIF


cc       Filtering noise (jpdincao)
c        do j = 1,LegPoints
c           if (dabs(PInterp(j)).le.1.d-14) PInterp(j) = 0.d0 !sign(1.d-14,PInterp(j))
c        enddo   

cc       Cutoff Long-range P-coupling (jpdincao)
c        if (nc.eq.1.or.nc.eq.3) then
c        if (jc.eq.1.or.jc.eq.3) then
c        do j = 1,LegPoints
c           if (x(j).ge.WhereStart) PInterp(j) = 0.d0 
c        enddo
c        endif
c        endif

       do i = 1,6
         BasicP(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicP(i,ip) = BasicP(i,ip) + wLeg(j)*PInterp(j)*
     >        (u(ip,j)*up(i,j) - up(ip,j)*u(i,j))
         enddo
       enddo

       do i = 1,6
        irow = MEMap(jc,i)

c diagonal blocks
           if (nc .eq. jc)then
            if (i .ne. ip)then
             Gco(irow,nc) = 2.0d0*mass*Scale(i)*Scale(ip)*
     >        (energy*Intscale*BasicS(i,ip) -
     >        (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >         + BasicP(i,ip)))
            else
             if (jc .le. NumOpenL)then
             Goo(jc,nc) =  2.0d0*mass*Scale(i)*Scale(ip)*
     >        (energy*Intscale*BasicS(i,ip) -
     >        (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >         + BasicP(i,ip))) 
             endif
            endif

c bloch operator
           if (i .eq. 2)Gco(irow,nc)=Gco(irow,nc) + 1.0d0/Intscale

c off-diagonal blocks
           else
            if (i .ne. ip)then
             Gco(irow,nc) = -2.0d0*mass*Scale(i)*Scale(ip)*
     >             (BasicVQ(i,ip)*IntScale + BasicP(i,ip))
            else
             if (jc .le. NumOpenL)then
             Goo(jc,nc) =  -2.0d0*mass*Scale(i)*Scale(ip)*
     >                (BasicVQ(i,ip)*IntScale + BasicP(i,ip))
             endif
            endif
           endif

          enddo

        enddo
       enddo

c now do open channels on right boundary (Gco)

      ip = 5
      n = Numsectors

      do nc = NumOpenL+1,NumOpen
        IntScale = (xSector(n+1)-xSector(n))/2.0d0
        ScaledZero = (xSector(n+1)+xSector(n))/2.0d0
        knv = nc - NumOpenL

c integrate potentials in each block

        do jc = 1,Numchannels
         kjcs = jc  + Numchannels*(Numsectors-1)
         kjv = jc 

        do ik=1,NumDataPoints
         ydata(ik) = VQmat(ik,kjv,knv)
        enddo

        do j = 1,LegPoints
         x(j) = IntScale*xLeg(j)+ScaledZero

	   vvj=VQFit(kjv,knv,1)/x(j)**VQExponents(kjv,knv,1)
	   vvj=vvj+VQFit(kjv,knv,2)/x(j)**VQExponents(kjv,knv,2)
	   vvj=vvj+VQFit(kjv,knv,3)/x(j)**VQExponents(kjv,knv,3)
	   VQInterp(j)=vvj

	   ppj=PFit(kjv,knv,1)/x(j)**PExponents(kjv,knv,1)
	   ppj=ppj+PFit(kjv,knv,2)/x(j)**PExponents(kjv,knv,2)
	   ppj=ppj+PFit(kjv,knv,3)/x(j)**PExponents(kjv,knv,3)
	   PInterp(j)=ppj
        enddo

	IF(x(1) .LT. WhereStart) THEN
        call Akima(NumDataPoints,xdata,ydata,LegPoints,x,VQInterp)
	ENDIF

       do i = 1,6
         BasicVQ(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicVQ(i,ip) = BasicVQ(i,ip) + wLeg(j)*u(i,j)*VQInterp(j)*u(ip,j)
         enddo
       enddo

        do ik=1,NumDataPoints
         ydata(ik)=Pmat(ik,kjv,knv)
        enddo

	   
	 IF(x(1) .LT. WhereStart) THEN
         call akima(NumDataPoints,xdata,ydata,LegPoints,x,Pinterp)
	 ENDIF

       do i = 1,6
         BasicP(i,ip) = 0.0d0
         do j = 1,LegPoints
          BasicP(i,ip) = BasicP(i,ip) + wLeg(j)*PInterp(j)*
     >    (u(i,j)*up(ip,j) - up(i,j)*u(ip,j))
         enddo
       enddo

        If(n .ne. 1)Scale(2) = (xSector(n+1)-xSector(n))/(xSector(n)-xSector(n-1))

       do i = 1,6
        irow = MEMap(kjcs,i)

c diagonal blocks
           if (nc .eq. jc+NumOpenL)then
            if (i .ne. ip)then
             Gco(irow,nc) = 2.0d0*mass*Scale(i)*Scale(ip)*
     >        (energy*Intscale*BasicS(i,ip) -
     >        (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >         + BasicP(i,ip)))
            else 
             if (jc .le. NumOpenR)then
             Goo(jc+NumOpenL,nc) =  2.0d0*mass*Scale(i)*Scale(ip)*
     >        (energy*Intscale*BasicS(i,ip) -
     >        (BasicT(i,ip)/(Mass*IntScale)+BasicVQ(i,ip)*IntScale
     >         + BasicP(i,ip))) 
             endif
            endif

c bloch operator
           if (i .eq. 6)Gco(irow,nc)=Gco(irow,nc) - 1.0d0/Intscale


c off-diagonal blocks
           else
            if (i .ne. ip)then
             Gco(irow,nc) = -2.0d0*mass*Scale(i)*Scale(ip)*
     >                (BasicVQ(i,ip)*IntScale+ BasicP(i,ip))
            else
             if (jc .le. NumOpenR)then
             Goo(jc+NumOpenL,nc) =  -2.0d0*mass*Scale(i)*Scale(ip)*
     >                  (BasicVQ(i,ip)*IntScale + BasicP(i,ip))
             endif
            endif
           endif

          enddo

        enddo
       enddo

      return
      end

c     ***********************************************************************

      subroutine CalcSmatrix(NumOpenR,leff,Thresholds,TotalAngMomentum,Mass,
     >                       Energy,RMatch,solution,deriv,Smatrix,Tmatrix,CrossSections)
      integer TotalAngMomentum
      integer NumOpenR,ipiv(NumOpenR)
      double precision leff(NumOpenR),Thresholds(NumOpenR),Mass,RMatch
      double precision solution(NumOpenR,NumOpenR)
      double precision deriv(NumOpenR),cmPERau,secPERau,convertCGS
      double precision convertCGS_K3,convertCGS_RX
      double precision CrossSections(NumOpenR,NumOpenR)
      double precision kVector(NumOpenR),f,fp,g,gp,Wronskian,pi,Energy
      double precision Imat(NumOpenR,NumOpenR)
      double precision Jmat(NumOpenR,NumOpenR)
      double precision enorm(NumOpenR),x
      double precision Tmatrix(NumOpenR,NumOpenR),TmatrixAux(NumOpenR,NumOpenR),ImatAux(NumOpenR,NumOpenR)
      double complex Smatrix(NumOpenR,NumOpenR)
      double complex IJmat(NumOpenR,NumOpenR)

      pi = dacos(-1.0d0)
      cmPERau = 5.29178d-09
      secPERau = 2.4189d-17
      convertCGS_K3 = cmPERau**6/secPERau
      convertCGS_RX = cmPERau**3/secPERau

      do i=1,NumOpenR

       kVector(i) = dsqrt(2.0d0*Mass*(Energy-Thresholds(i)))
       enorm(i) = dsqrt(2.0d0*kVector(i)/pi)

c     call BesselBasePair(kVector(i),RMatch,leff(i),f,fp,g,gp)
 
      x = kVector(i)*Rmatch
      call sphbes(leff(i),x,f,g,fp,gp)
      fp = enorm(i)*(f + x*fp)
      f = enorm(i)*Rmatch*f
      gp = enorm(i)*(g + x*gp)
      g = enorm(i)*Rmatch*g

       Wronskian = 1.0d0/(g*fp-gp*f)
       !write(6,*) 'Wronskian check:    ',2.0d0*Wronskian/Pi
36    format(5(e12.6,1x))
 
       !write(6,*) 'Bessel funcs:'
       !write(6,*) i,f,fp
       !write(6,*) i,g,gp
       !write(6,*)

       do j = 1,NumOpenR
        IMat(i,j) = -(g*deriv(j)+gp)*Solution(i,j)*Wronskian
        JMat(i,j) = -(f*deriv(j)+fp)*Solution(i,j)*Wronskian
       enddo
 
      enddo 

c     calculate S-matrix

      IJMat = (0.0d0,0.0d0)
 
       do i = 1,NumOpenR
        do j = 1,NumOpenR
         IJMat(i,j)   = IMat(j,i)-(0.0d0,1.0d0)*JMat(j,i)
         Smatrix(i,j) = IMat(j,i)+(0.0d0,1.0d0)*JMat(j,i)
        enddo
       enddo
 
       call zgesv(NumOpenR,NumOpenR,IJMat,NumOpenR,ipiv,Smatrix,NumOpenR,info)

       do i = 1,NumOpenR
        Smatrix(i,i) = Smatrix(i,i) - (1.0d0,0.0d0)
       enddo


c     calculate T-matrix

       do i = 1,NumOpenR
        do j = 1,NumOpenR
         ImatAux(i,j) = Imat(j,i)
         TmatrixAux(i,j) = Jmat(j,i)
        enddo
       enddo

       call dgesv(NumOpenR,NumOpenR,ImatAux,NumOpenR,ipiv,TmatrixAux,NumOpenR,info)

       do i = 1,NumOpenR
        do j = 1,NumOpenR
         Tmatrix(i,j) = TmatrixAux(j,i)
        enddo
       enddo

 
c i -> final state
c j -> initial state
 
       do i = 1,NumOpenR
        do j = 1,NumOpenR
c       For Recombination 
        CrossSections(i,j) = convertCGS_K3*dfloat(2*TotalAngMomentum+1)*
     .                       192.d0*pi*pi/(Mass*kVector(j)**4)*
     .                       dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
c       For Relaxation
c       CrossSections(i,j) = convertCGS_RX*dfloat(2*TotalAngMomentum+1)*
c     .                        pi/(Mass*kVector(j)**1)*
c     .                        dreal(dconjg(Smatrix(i,j))*Smatrix(i,j))
        enddo
       enddo 

       return
       end


      
C***********************************************************************
C                                                                      *
C     SUBROUTINE AKIMA (L,X,Y,N,U,V)                                   *
C     INTERPOLATION OF A SINGLE-VALUED FUNCTION                        *
C                                                                      *
C     L = NUMBER OF INPUT DATA POINTS                                  *
C          (MUST BE 2 OR GREATER)                                      *
C     X = ARRAY OF DIMENSION L STORING THE X VALUES                    *
C         (ABSCISSAS) OF INPUT DATA POINTS                             *
C          (IN ASCENDING ORDER)                                        *
C     Y = ARRAY OF DIMENSION L STORING THE Y VALUES                    *
C         (ORDINATES) OF INPUT DATA POINTS                             *
C     N = NUMBER OF POINTS AT WHICH INTERPOLATION OF THE               *
C         Y VALUE (ORDINATE) IS DESIRED                                *
C          (MUST BE 1 OR GREATER)                                      *
C     U = ARRAY OF DIMENSION N STRORING THE X VALUES                   *
C         (ABSCISSAS) OF THE DESIRED POINTS                            *
C     V = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Y                *
C         VALUES (ORDINATES) ARE TO BE DISPLAYED                       *
C                                                                      *
C***********************************************************************
C
      SUBROUTINE AKIMA (L,X,Y,N,U,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(L),Y(L),U(N),V(N)
      EQUIVALENCE (P0,X3),(Q0,Y3),(Q1,T3)
      REAL*8 M1,M2,M3,M4,M5
      EQUIVALENCE (UK,DX),(IMN,X2,A1,M1),(IMX,X5,A5,M5),
     1            (J,SW,SA),(Y2,W2,W4,Q2),(Y5,W3,Q3)
 10   L0=L
      LM1=L0-1
      LM2=LM1-1
      LP1=L0+1
      N0=N
      IF (LM2.LT.0) GO TO 90
      IF (N0.LE.0) GO TO 91
      DO 11 I=2,L0
      IF (X(I-1)-X(I)) 11,95,96
 11   CONTINUE
      IPV=0
      DO 80 K=1,N0
      UK=U(K)
 20   IF (LM2.EQ.0) GO TO 27
      IF (UK.GE.X(L0)) GO TO 26
      IF (UK.LT.X(1)) GO TO 25
      IMN=2
      IMX=L0
 21   I=(IMN+IMX)/2
      IF (UK.GE.X(I)) GO TO 23
 22   IMX=I
      GO TO 24
 23   IMN=I+1
 24   IF (IMX.GT.IMN) GO TO 21
      I=IMX
      GO TO 30
 25   I=1
      GO TO 30
 26   I=LP1
      GO TO 30
 27   I=2
 30   IF (I.EQ.IPV) GO TO 70
      IPV=I
 40   J=I
      IF (J.EQ.1) J=2
      IF (J.EQ.LP1) J=L0
      X3=X(J-1)
      Y3=Y(J-1)
      X4=X(J)
      Y4=Y(J)
      A3=X4-X3
      M3=(Y4-Y3)/A3
      IF (LM2.EQ.0) GO TO 43
      IF (J.EQ.2) GO TO 41
      X2=X(J-2)
      Y2=Y(J-2)
      A2=X3-X2
      M2=(Y3-Y2)/A2
      IF (J.EQ.L0) GO TO 42
 41   X5=X(J+1)
      Y5=Y(J+1)
      A4=X5-X4
      M4=(Y5-Y4)/A4
      IF (J.EQ.2) M2=M3+M3-M4
      GO TO 45
 42   M4=M3+M3-M2
      GO TO 45
 43   M2=M3
      M4=M3
 45   IF (J.LE.3) GO TO 46
      A1=X2-X(J-3)
      M1=(Y2-Y(J-3))/A1
      GO TO 47
 46   M1=M2+M2-M3
 47   IF (J.GE.LM1) GO TO 48
      A5=X(J+2)-X5
      M5=(Y(J+2)-Y5)/A5
      GO TO 50
 48   M5=M4+M4-M3
 50   IF (I.EQ.LP1) GO TO 52
      W2=DABS(M4-M3)
      W3=DABS(M2-M1)
      SW=W2+W3
      IF (SW.NE.0.D0) GO TO 51
      W2=.5D0
      W3=.5D0
      SW=1.D0
 51   T3=(W2*M2+W3*M3)/SW
      IF (I.EQ.1) GO TO 54
 52   W3=DABS(M5-M4)
      W4=DABS(M3-M2)
      SW=W3+W4
      IF (SW.NE.0.D0) GO TO 53
      W3=.5D0
      W4=.5D0
      SW=1.D0
 53   T4=(W3*M3+W4*M4)/SW
      IF (I.NE.LP1) GO TO 60
      T3=T4
      SA=A2+A3
      T4=.5D0*(M4+M5-A2*(A2-A3)*(M2-M3)/(SA*SA))
      X3=X4
      Y3=Y4
      A3=A2
      M3=M4
      GO TO 60
 54   T4=T3
      SA=A3+A4
      T3=.5D0*(M1+M2-A4*(A3-A4)*(M3-M4)/(SA*SA))
      X3=X3-A4
      Y3=Y3-M2*A4
      A3=A4
      M3=M2
 60   Q2=(2.D0*(M3-T3)+M3-T4)/A3
      Q3=(-M3-M3+T3+T4)/(A3*A3)
 70   DX=UK-P0
 80   V(K)=Q0+DX*(Q1+DX*(Q2+DX*Q3))
      RETURN
 90   WRITE (6,2090)
      GO TO 99
 91   WRITE (6,2091)
      GO TO 99
 95   WRITE (6,2095)
      GO TO 97
 96   WRITE (6,2096)
 97   WRITE (6,2097) I,X(I)
 99   WRITE (6,2099) L0,N0
      RETURN
 2090 FORMAT ('***   L = 1 OR LESS.')
 2091 FORMAT (' ***   N = 0 OR LESS.')
 2095 FORMAT (' ***   IDENTICAL X VALUES.')
 2096 FORMAT (' ***   X VALUES OUT OF SEQUENCE.')
 2097 FORMAT ('   I =',I7,10X,'X(I) =',D12.3)
 2099 FORMAT ('   L =',I7,10X,'N =',I7/
     >       'ERROR DETECTED IN ROUTINE     AKIMA')
      END
c************************************************************************************
      subroutine BoxMatching(NumOpenL,NumOpenR,psi_in,b,Goo,logderiv,solution,deriv)
      integer NumOpenChannels,NumOpen,NumOpenR,ierr
      integer Matz,NumCoef,NumOpenTot
      double precision b(NumOpenL),logderiv(NumOpenL+NumOpenR)
      double precision psi_in(NumOpenL,NumOpenL)
      double precision Goo(NumOpenL+NumOpenR,NumOpenL+NumOpenR)
      double precision left(2*NumOpenL+NumOpenR,2*NumOpenL+NumOpenR)
      double precision right(2*NumOpenL+NumOpenR,2*NumOpenL+NumOpenR)
      double precision areal(2*NumOpenL+NumOpenR)
      double precision aimag(2*NumOpenL+NumOpenR)
      double precision cpsi(2*NumOpenL+NumOpenR,2*NumOpenL+NumOpenR)
      double precision denom(2*NumOpenL+NumOpenR)
      double precision sum,solution(NumOpenR,NumOpenR)
      double precision coef(2*NumOpenL+NumOpenR,NumOpenR)
      double precision deriv(NumOpenR)

      NumOpenTot = NumOpenL+NumOpenR

      do i=1,NumOpenL
        do j=1,NumOpenL
         left(i,j) = psi_in(i,j)
         right(i,j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenL
        do j=1,NumOpenTot
         left(i,NumOpenL+j) = -Goo(i,j)
         right(i,NumOpenL+j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenL
        do j=1,NumOpenL
         left(i+NumOpenL,j) = -psi_in(i,j)*b(j)
         right(i+NumOpenL,j) = 0.0d0 
        enddo
       enddo
 
       do i=1,NumOpenL
        do j=1,NumOpenTot
         left(i+NumOpenL,NumOpenL+j) = -Goo(i,j)*logderiv(j)
         right(i+NumOpenL,NumOpenL+j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenR
        do j=1,NumOpenTot
         left(i+2*NumOpenL,j) = 0.0d0
         right(i+2*NumOpenL,j) = 0.0d0
        enddo
       enddo
 
       do i=1,NumOpenR
        do j=1,NumOpenTot
         left(i+2*NumOpenL,NumOpenL+j) =
     >       -Goo(i+NumOpenL,j)*logderiv(j)
         right(i+2*NumOpenL,NumOpenL+j) = Goo(i+NumOpenL,j)
        enddo
       enddo 

       NumCoef = NumOpenL + NumOpenTot
       matz = 1
       call rgg(NumCoef,NumCoef,left,right,areal,aimag,denom,matz,cpsi,ierr)
       !write(6,*)
       !write(6,*)'return from rgg, ierr=',ierr
       !write(6,*)

        iflag = 1
        do i=1,NumCoef
         if (abs(denom(i)) .gt. 1.0d-8)then
          deriv(iflag) = -areal(i)/denom(i)
          do j=1,NumCoef
           coef(j,iflag) = cpsi(j,i)
          enddo
          iflag = iflag+1
         endif
        enddo   

        do i=1,NumOpenR
         do j=1,NumOpenR
          solution(i,j) = 0.0d0
          do k=1,NumOpentot
           solution(i,j) = solution(i,j) +
     >           Goo(i+NumOpenL,k)*coef(k+NumOpenL,j)
          enddo
         enddo
        enddo
 
       do i=1,NumOpenR
        sum = 0.0d0
        do j=1,NumOpenR
         sum = sum + solution(j,i)*solution(j,i)
        enddo
        sum = dsqrt(sum)
        do j=1,NumOpenR
         solution(j,i) = solution(j,i)/sum
        enddo
       enddo 

       return 
       end


      subroutine rgg(nm,n,a,b,alfr,alfi,beta,matz,z,ierr)               
                                                                        
      integer n,nm,ierr,matz                                            
      real*8 a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n),epsl           
      logical tf                                                        
c                                                                       
c     this subroutine calls the recommended sequence of                 
c     subroutines from the eigensystem subroutine package (eispack)     
c     to find the eigenvalues and eigenvectors (if desired)             
c     for the real general generalized eigenproblem  ax = (lambda)bx.   
c                                                                       
c     on input:                                                         
c                                                                       
c        nm  must be set to the row dimension of the two-dimensional    
c        array parameters as declared in the calling program            
c        dimension statement;                                           
c                                                                       
c        n  is the order of the matrices  a  and  b;                    
c                                                                       
c        a  contains a real general matrix;                             
c                                                                       
c        b  contains a real general matrix;                             
c                                                                       
c        matz  is an integer variable set equal to zero if              
c        only eigenvalues are desired;  otherwise it is set to          
c        any non-zero integer for both eigenvalues and eigenvectors.    
c                                                                       
c     on output:                                                        
c                                                                       
c        alfr  and  alfi  contain the real and imaginary parts,         
c        respectively, of the numerators of the eigenvalues;            
c                                                                       
c        beta  contains the denominators of the eigenvalues,            
c        which are thus given by the ratios  (alfr+i*alfi)/beta.        
c        complex conjugate pairs of eigenvalues appear consecutively    
c        with the eigenvalue having the positive imaginary part first;  
c                                                                       
c        z  contains the real and imaginary parts of the eigenvectors   
c        if matz is not zero.  if the j-th eigenvalue is real, the      
c        j-th column of  z  contains its eigenvector.  if the j-th      
c        eigenvalue is complex with positive imaginary part, the        
c        j-th and (j+1)-th columns of  z  contain the real and          
c        imaginary parts of its eigenvector.  the conjugate of this     
c        vector is the eigenvector for the conjugate eigenvalue;        
c                                                                       
c        ierr  is an integer output variable set equal to an            
c        error completion code described in section 2b of the           
c        documentation.  the normal completion code is zero.            
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
      epsl=1.0d-12
      if (n .le. nm) go to 10                                           
      ierr = 10 * n                                                     
      go to 50                                                          
                                                                        
   10 if (matz .ne. 0) go to 20                                         
c     :::::::::: find eigenvalues only ::::::::::                       
      tf = .false.                                                      
      call  qzhes(nm,n,a,b,tf,z)                                        
      call  qzit(nm,n,a,b,epsl,tf,z,ierr)                              
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)                         
      go to 50                                                          
c     :::::::::: find both eigenvalues and eigenvectors ::::::::::      
   20 tf = .true.                                                       
      call  qzhes(nm,n,a,b,tf,z)                                        
      call  qzit(nm,n,a,b,epsl,tf,z,ierr)                              
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)                         
      if (ierr .ne. 0) go to 50                                         
      call  qzvec(nm,n,a,b,alfr,alfi,beta,z)                            
   50 return                                                            
      end                                                               
**********************************************
      subroutine qzhes(nm,n,a,b,matz,z)                                 
                                                                        
      integer i,j,k,l,n,lb,l1,nm,nk1,nm1,nm2                            
      real*8 a(nm,n),b(nm,n),z(nm,n)                                    
      real*8 r,s,t,u1,u2,v1,v2,rho                                      
      real*8 dsqrt,dabs,dsign                                           
      logical matz                                                      
c                                                                       
c     this subroutine is the first step of the qz algorithm             
c     for solving generalized matrix eigenvalue problems,               
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.      
c                                                                       
c     this subroutine accepts a pair of real general matrices and       
c     reduces one of them to upper hessenberg form and the other        
c     to upper triangular form using orthogonal transformations.        
c     it is usually followed by  qzit,  qzval  and, possibly,  qzvec.   
c                                                                       
c     on input:                                                         
c                                                                       
c        nm must be set to the row dimension of two-dimensional         
c          array parameters as declared in the calling program          
c          dimension statement;                                         
c                                                                       
c        n is the order of the matrices;                                
c                                                                       
c        a contains a real general matrix;                              
c                                                                       
c        b contains a real general matrix;                              
c                                                                       
c        matz should be set to .true. if the right hand transformations 
c          are to be accumulated for later use in computing             
c          eigenvectors, and to .false. otherwise.                      
c                                                                       
c     on output:                                                        
c                                                                       
c        a has been reduced to upper hessenberg form.  the elements     
c          below the first subdiagonal have been set to zero;           
c                                                                       
c        b has been reduced to upper triangular form.  the elements     
c          below the main diagonal have been set to zero;               
c                                                                       
c        z contains the product of the right hand transformations if    
c          matz has been set to .true.  otherwise, z is not referenced. 
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
c     :::::::::: initialize z ::::::::::                                
      if (.not. matz) go to 10                                          
      do 3 i = 1, n                                                     
         do 2 j = 1, n                                                  
            z(i,j) = 0.0d0                                              
    2    continue                                                       
         z(i,i) = 1.0d0                                                 
    3 continue                                                          
c     :::::::::: reduce b to upper triangular form ::::::::::           
   10 if (n .le. 1) go to 170                                           
      nm1 = n - 1                                                       
      do 100 l = 1, nm1                                                 
         l1 = l + 1                                                     
         s = 0.0d0                                                      
                                                                        
         do 20 i = l1, n                                                
            s = s + dabs(b(i,l))                                        
   20    continue                                                       
                                                                        
         if (s .eq. 0.0d0) go to 100                                    
         s = s + dabs(b(l,l))                                           
         r = 0.0d0                                                      
                                                                        
         do 25 i = l, n                                                 
            b(i,l) = b(i,l) / s                                         
            r = r + b(i,l)**2                                           
   25    continue                                                       
                                                                        
         r = dsign(dsqrt(r),b(l,l))                                     
         b(l,l) = b(l,l) + r                                            
         rho = r * b(l,l)                                               
                                                                        
         do 50 j = l1, n                                                
            t = 0.0d0                                                   
            do 30 i = l, n                                              
               t = t + b(i,l) * b(i,j)                                  
   30       continue                                                    
            t = -t / rho                                                
            do 40 i = l, n                                              
               b(i,j) = b(i,j) + t * b(i,l)                             
   40       continue                                                    
   50    continue                                                       
         do 80 j = 1, n                                                 
            t = 0.0d0                                                   
            do 60 i = l, n                                              
               t = t + b(i,l) * a(i,j)                                  
   60       continue                                                    
            t = -t / rho                                                
            do 70 i = l, n                                              
               a(i,j) = a(i,j) + t * b(i,l)                             
   70       continue                                                    
   80    continue                                                       
         b(l,l) = -s * r                                                
         do 90 i = l1, n                                                
            b(i,l) = 0.0d0                                              
   90    continue                                                       
  100 continue                                                          
c     :::::::::: reduce a to upper hessenberg form, while               
c                keeping b triangular ::::::::::                        
      if (n .eq. 2) go to 170                                           
      nm2 = n - 2                                                       
                                                                        
      do 160 k = 1, nm2                                                 
         nk1 = nm1 - k                                                  
c     :::::::::: for l=n-1 step -1 until k+1 do -- ::::::::::           
         do 150 lb = 1, nk1                                             
            l = n - lb                                                  
            l1 = l + 1                                                  
c     :::::::::: zero a(l+1,k) ::::::::::                               
            s = dabs(a(l,k)) + dabs(a(l1,k))                            
            if (s .eq. 0.0d0) go to 150                                 
            u1 = a(l,k) / s                                             
            u2 = a(l1,k) / s                                            
            r = dsign(dsqrt(u1*u1+u2*u2),u1)                            
            v1 =  -(u1 + r) / r                                         
            v2 = -u2 / r                                                
            u2 = v2 / v1                                                
                                                                        
            do 110 j = k, n                                             
               t = a(l,j) + u2 * a(l1,j)                                
               a(l,j) = a(l,j) + t * v1                                 
               a(l1,j) = a(l1,j) + t * v2                               
  110       continue                                                    
                                                                        
            a(l1,k) = 0.0d0                                             
                                                                        
            do 120 j = l, n                                             
               t = b(l,j) + u2 * b(l1,j)                                
               b(l,j) = b(l,j) + t * v1                                 
               b(l1,j) = b(l1,j) + t * v2                               
  120       continue                                                    
c     :::::::::: zero b(l+1,l) ::::::::::                               
            s = dabs(b(l1,l1)) + dabs(b(l1,l))                          
            if (s .eq. 0.0d0) go to 150                                 
            u1 = b(l1,l1) / s                                           
            u2 = b(l1,l) / s                                            
            r = dsign(dsqrt(u1*u1+u2*u2),u1)                            
            v1 =  -(u1 + r) / r                                         
            v2 = -u2 / r                                                
            u2 = v2 / v1                                                
                                                                        
            do 130 i = 1, l1                                            
               t = b(i,l1) + u2 * b(i,l)                                
               b(i,l1) = b(i,l1) + t * v1                               
               b(i,l) = b(i,l) + t * v2                                 
  130       continue                                                    
                                                                        
            b(l1,l) = 0.0d0                                             
                                                                        
            do 140 i = 1, n                                             
               t = a(i,l1) + u2 * a(i,l)                                
               a(i,l1) = a(i,l1) + t * v1                               
               a(i,l) = a(i,l) + t * v2                                 
  140       continue                                                    
            if (.not. matz) go to 150                                   
            do 145 i = 1, n                                             
               t = z(i,l1) + u2 * z(i,l)                                
               z(i,l1) = z(i,l1) + t * v1                               
               z(i,l) = z(i,l) + t * v2                                 
  145       continue                                                    
  150    continue                                                       
  160 continue                                                          
  170 return                                                            
      end                                                               
*********************************************
      subroutine qzit(nm,n,a,b,eps1,matz,z,ierr)                        
                                                                        
      integer i,j,k,l,n,en,k1,k2,ld,ll,l1,na,nm,ish,its,km1,lm1,        
     x        enm2,ierr,lor1,enorn                                      
      real*8 a(nm,n),b(nm,n),z(nm,n)                                    
      real*8 r,s,t,a1,a2,a3,ep,sh,u1,u2,u3,v1,v2,v3,ani,a11,a12,        
     x       a21,a22,a33,a34,a43,a44,bni,b11,b12,b22,b33,b34,           
     x       b44,epsa,epsb,eps1,anorm,bnorm                             
      real*8 dsqrt,dabs,dsign                                           
      integer max0,min0                                                 
      logical matz,notlas                                               
c                                                                       
c     this subroutine is the second step of the qz algorithm            
c     for solving generalized matrix eigenvalue problems,               
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,      
c     as modified in technical note nasa tn d-7305(1973) by ward.       
c                                                                       
c     this subroutine accepts a pair of real matrices, one of them      
c     in upper hessenberg form and the other in upper triangular form.  
c     it reduces the hessenberg matrix to quasi-triangular form using   
c     orthogonal transformations while maintaining the triangular form  
c     of the other matrix.  it is usually preceded by  qzhes  and       
c     followed by  qzval  and, possibly,  qzvec.                        
c                                                                       
c     on input:                                                         
c                                                                       
c        nm must be set to the row dimension of two-dimensional         
c          array parameters as declared in the calling program          
c          dimension statement;                                         
c                                                                       
c        n is the order of the matrices;                                
c                                                                       
c        a contains a real upper hessenberg matrix;                     
c                                                                       
c        b contains a real upper triangular matrix;                     
c                                                                       
c        eps1 is a tolerance used to determine negligible elements.     
c          eps1 = 0.0 (or negative) may be input, in which case an      
c          element will be neglected only if it is less than roundoff   
c          error times the norm of its matrix.  if the input eps1 is    
c          positive, then an element will be considered negligible      
c          if it is less than eps1 times the norm of its matrix.  a     
c          positive value of eps1 may result in faster execution,       
c          but less accurate results;                                   
c                                                                       
c        matz should be set to .true. if the right hand transformations 
c          are to be accumulated for later use in computing             
c          eigenvectors, and to .false. otherwise;                      
c                                                                       
c        z contains, if matz has been set to .true., the                
c          transformation matrix produced in the reduction              
c          by  qzhes, if performed, or else the identity matrix.        
c          if matz has been set to .false., z is not referenced.        
c                                                                       
c     on output:                                                        
c                                                                       
c        a has been reduced to quasi-triangular form.  the elements     
c          below the first subdiagonal are still zero and no two        
c          consecutive subdiagonal elements are nonzero;                
c                                                                       
c        b is still in upper triangular form, although its elements     
c          have been altered.  the location b(n,1) is used to store     
c          eps1 times the norm of b for later use by  qzval  and  qzvec;
c                                                                       
c        z contains the product of the right hand transformations       
c          (for both steps) if matz has been set to .true.;             
c                                                                       
c        ierr is set to                                                 
c          zero       for normal return,                                
c          j          if neither a(j,j-1) nor a(j-1,j-2) has become     
c                     zero after 50 iterations.                         
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
      ierr = 0                                                          
c     :::::::::: compute epsa,epsb ::::::::::                           
      anorm = 0.0d0                                                     
      bnorm = 0.0d0                                                     
                                                                        
      do 30 i = 1, n                                                    
         ani = 0.0d0                                                    
         if (i .ne. 1) ani = dabs(a(i,i-1))                             
         bni = 0.0d0                                                    
                                                                        
         do 20 j = i, n                                                 
            ani = ani + dabs(a(i,j))                                    
            bni = bni + dabs(b(i,j))                                    
   20    continue                                                       
                                                                        
         if (ani .gt. anorm) anorm = ani                                
         if (bni .gt. bnorm) bnorm = bni                                
   30 continue                                                          
                                                                        
      if (anorm .eq. 0.0d0) anorm = 1.0d0                               
      if (bnorm .eq. 0.0d0) bnorm = 1.0d0                               
      ep = eps1                                                         
      if (ep .gt. 0.0d0) go to 50                                       
c     :::::::::: compute roundoff level if eps1 is zero ::::::::::      
      ep = 1.0d0                                                        
   40 ep = ep / 2.0d0                                                   
      if (1.0d0 + ep .gt. 1.0d0) go to 40                               
   50 epsa = ep * anorm                                                 
      epsb = ep * bnorm                                                 
c     :::::::::: reduce a to quasi-triangular form, while               
c                keeping b triangular ::::::::::                        
      lor1 = 1                                                          
      enorn = n                                                         
      en = n                                                            
c     :::::::::: begin qz step ::::::::::                               
   60 if (en .le. 2) go to 1001                                         
      if (.not. matz) enorn = en                                        
      its = 0                                                           
      na = en - 1                                                       
      enm2 = na - 1                                                     
   70 ish = 2                                                           
c     :::::::::: check for convergence or reducibility.                 
c                for l=en step -1 until 1 do -- ::::::::::              
      do 80 ll = 1, en                                                  
         lm1 = en - ll                                                  
         l = lm1 + 1                                                    
         if (l .eq. 1) go to 95                                         
         if (dabs(a(l,lm1)) .le. epsa) go to 90                         
   80 continue                                                          
                                                                        
   90 a(l,lm1) = 0.0d0                                                  
      if (l .lt. na) go to 95                                           
c     :::::::::: 1-by-1 or 2-by-2 block isolated ::::::::::             
      en = lm1                                                          
      go to 60                                                          
c     :::::::::: check for small top of b ::::::::::                    
   95 ld = l                                                            
  100 l1 = l + 1                                                        
      b11 = b(l,l)                                                      
      if (dabs(b11) .gt. epsb) go to 120                                
      b(l,l) = 0.0d0                                                    
      s = dabs(a(l,l)) + dabs(a(l1,l))                                  
      u1 = a(l,l) / s                                                   
      u2 = a(l1,l) / s                                                  
      r = dsign(dsqrt(u1*u1+u2*u2),u1)                                  
      v1 = -(u1 + r) / r                                                
      v2 = -u2 / r                                                      
      u2 = v2 / v1                                                      
                                                                        
      do 110 j = l, enorn                                               
         t = a(l,j) + u2 * a(l1,j)                                      
         a(l,j) = a(l,j) + t * v1                                       
         a(l1,j) = a(l1,j) + t * v2                                     
         t = b(l,j) + u2 * b(l1,j)                                      
         b(l,j) = b(l,j) + t * v1                                       
         b(l1,j) = b(l1,j) + t * v2                                     
  110 continue                                                          
                                                                        
      if (l .ne. 1) a(l,lm1) = -a(l,lm1)                                
      lm1 = l                                                           
      l = l1                                                            
      go to 90                                                          
  120 a11 = a(l,l) / b11                                                
      a21 = a(l1,l) / b11                                               
      if (ish .eq. 1) go to 140                                         
c     :::::::::: iteration strategy ::::::::::                          
      if (its .eq. 300) go to 1000                                       
      if (its .eq. 10) go to 155                                        
c     :::::::::: determine type of shift ::::::::::                     
      b22 = b(l1,l1)                                                    
      if (dabs(b22) .lt. epsb) b22 = epsb                               
      b33 = b(na,na)                                                    
      if (dabs(b33) .lt. epsb) b33 = epsb                               
      b44 = b(en,en)                                                    
      if (dabs(b44) .lt. epsb) b44 = epsb                               
      a33 = a(na,na) / b33                                              
      a34 = a(na,en) / b44                                              
      a43 = a(en,na) / b33                                              
      a44 = a(en,en) / b44                                              
      b34 = b(na,en) / b44                                              
      t = 0.5d0 * (a43 * b34 - a33 - a44)                               
      r = t * t + a34 * a43 - a33 * a44                                 
      if (r .lt. 0.0d0) go to 150                                       
c     :::::::::: determine single shift zeroth column of a ::::::::::   
      ish = 1                                                           
      r = dsqrt(r)                                                      
      sh = -t + r                                                       
      s = -t - r                                                        
      if (dabs(s-a44) .lt. dabs(sh-a44)) sh = s                         
c     :::::::::: look for two consecutive small                         
c                sub-diagonal elements of a.                            
c                for l=en-2 step -1 until ld do -- ::::::::::           
      do 130 ll = ld, enm2                                              
         l = enm2 + ld - ll                                             
         if (l .eq. ld) go to 140                                       
         lm1 = l - 1                                                    
         l1 = l + 1                                                     
         t = a(l,l)                                                     
         if (dabs(b(l,l)) .gt. epsb) t = t - sh * b(l,l)                
         if (dabs(a(l,lm1)) .le. dabs(t/a(l1,l)) * epsa) go to 100      
  130 continue                                                          
                                                                        
  140 a1 = a11 - sh                                                     
      a2 = a21                                                          
      if (l .ne. ld) a(l,lm1) = -a(l,lm1)                               
      go to 160                                                         
c     :::::::::: determine double shift zeroth column of a ::::::::::   
  150 a12 = a(l,l1) / b22                                               
      a22 = a(l1,l1) / b22                                              
      b12 = b(l,l1) / b22                                               
      a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11)    
     x     / a21 + a12 - a11 * b12                                      
      a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11)          
     x     + a43 * b34                                                  
      a3 = a(l1+1,l1) / b22                                             
      go to 160                                                         
c     :::::::::: ad hoc shift ::::::::::                                
  155 a1 = 0.0d0                                                        
      a2 = 1.0d0                                                        
      a3 = 1.1605d0                                                     
  160 its = its + 1                                                     
      if (.not. matz) lor1 = ld                                         
c     :::::::::: main loop ::::::::::                                   
      do 260 k = l, na                                                  
         notlas = k .ne. na .and. ish .eq. 2                            
         k1 = k + 1                                                     
         k2 = k + 2                                                     
         km1 = max0(k-1,l)                                              
         ll = min0(en,k1+ish)                                           
         if (notlas) go to 190                                          
c     :::::::::: zero a(k+1,k-1) ::::::::::                             
         if (k .eq. l) go to 170                                        
         a1 = a(k,km1)                                                  
         a2 = a(k1,km1)                                                 
  170    s = dabs(a1) + dabs(a2)                                        
         if (s .eq. 0.0d0) go to 70                                     
         u1 = a1 / s                                                    
         u2 = a2 / s                                                    
         r = dsign(dsqrt(u1*u1+u2*u2),u1)                               
         v1 = -(u1 + r) / r                                             
         v2 = -u2 / r                                                   
         u2 = v2 / v1                                                   
                                                                        
         do 180 j = km1, enorn                                          
            t = a(k,j) + u2 * a(k1,j)                                   
            a(k,j) = a(k,j) + t * v1                                    
            a(k1,j) = a(k1,j) + t * v2                                  
            t = b(k,j) + u2 * b(k1,j)                                   
            b(k,j) = b(k,j) + t * v1                                    
            b(k1,j) = b(k1,j) + t * v2                                  
  180    continue                                                       
                                                                        
         if (k .ne. l) a(k1,km1) = 0.0d0                                
         go to 240                                                      
c     :::::::::: zero a(k+1,k-1) and a(k+2,k-1) ::::::::::              
  190    if (k .eq. l) go to 200                                        
         a1 = a(k,km1)                                                  
         a2 = a(k1,km1)                                                 
         a3 = a(k2,km1)                                                 
  200    s = dabs(a1) + dabs(a2) + dabs(a3)                             
         if (s .eq. 0.0d0) go to 260                                    
         u1 = a1 / s                                                    
         u2 = a2 / s                                                    
         u3 = a3 / s                                                    
         r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)                         
         v1 = -(u1 + r) / r                                             
         v2 = -u2 / r                                                   
         v3 = -u3 / r                                                   
         u2 = v2 / v1                                                   
         u3 = v3 / v1                                                   
                                                                        
         do 210 j = km1, enorn                                          
            t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)                    
            a(k,j) = a(k,j) + t * v1                                    
            a(k1,j) = a(k1,j) + t * v2                                  
            a(k2,j) = a(k2,j) + t * v3                                  
            t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)                    
            b(k,j) = b(k,j) + t * v1                                    
            b(k1,j) = b(k1,j) + t * v2                                  
            b(k2,j) = b(k2,j) + t * v3                                  
  210    continue                                                       
                                                                        
         if (k .eq. l) go to 220                                        
         a(k1,km1) = 0.0d0                                              
         a(k2,km1) = 0.0d0                                              
c     :::::::::: zero b(k+2,k+1) and b(k+2,k) ::::::::::                
  220    s = dabs(b(k2,k2)) + dabs(b(k2,k1)) + dabs(b(k2,k))            
         if (s .eq. 0.0d0) go to 240                                    
         u1 = b(k2,k2) / s                                              
         u2 = b(k2,k1) / s                                              
         u3 = b(k2,k) / s                                               
         r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)                         
         v1 = -(u1 + r) / r                                             
         v2 = -u2 / r                                                   
         v3 = -u3 / r                                                   
         u2 = v2 / v1                                                   
         u3 = v3 / v1                                                   
                                                                        
         do 230 i = lor1, ll                                            
            t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)                    
            a(i,k2) = a(i,k2) + t * v1                                  
            a(i,k1) = a(i,k1) + t * v2                                  
            a(i,k) = a(i,k) + t * v3                                    
            t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)                    
            b(i,k2) = b(i,k2) + t * v1                                  
            b(i,k1) = b(i,k1) + t * v2                                  
            b(i,k) = b(i,k) + t * v3                                    
  230    continue                                                       
                                                                        
         b(k2,k) = 0.0d0                                                
         b(k2,k1) = 0.0d0                                               
         if (.not. matz) go to 240                                      
                                                                        
         do 235 i = 1, n                                                
            t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)                    
            z(i,k2) = z(i,k2) + t * v1                                  
            z(i,k1) = z(i,k1) + t * v2                                  
            z(i,k) = z(i,k) + t * v3                                    
  235    continue                                                       
c     :::::::::: zero b(k+1,k) ::::::::::                               
  240    s = dabs(b(k1,k1)) + dabs(b(k1,k))                             
         if (s .eq. 0.0d0) go to 260                                    
         u1 = b(k1,k1) / s                                              
         u2 = b(k1,k) / s                                               
         r = dsign(dsqrt(u1*u1+u2*u2),u1)                               
         v1 = -(u1 + r) / r                                             
         v2 = -u2 / r                                                   
         u2 = v2 / v1                                                   
                                                                        
         do 250 i = lor1, ll                                            
            t = a(i,k1) + u2 * a(i,k)                                   
            a(i,k1) = a(i,k1) + t * v1                                  
            a(i,k) = a(i,k) + t * v2                                    
            t = b(i,k1) + u2 * b(i,k)                                   
            b(i,k1) = b(i,k1) + t * v1                                  
            b(i,k) = b(i,k) + t * v2                                    
  250    continue                                                       
                                                                        
         b(k1,k) = 0.0d0                                                
         if (.not. matz) go to 260                                      
                                                                        
         do 255 i = 1, n                                                
            t = z(i,k1) + u2 * z(i,k)                                   
            z(i,k1) = z(i,k1) + t * v1                                  
            z(i,k) = z(i,k) + t * v2                                    
  255    continue                                                       
                                                                        
  260 continue                                                          
c     :::::::::: end qz step ::::::::::                                 
      go to 70                                                          
c     :::::::::: set error -- neither bottom subdiagonal element        
c                has become negligible after 50 iterations ::::::::::   
 1000 ierr = en                                                         
c     :::::::::: save epsb for use by qzval and qzvec ::::::::::        
 1001 if (n .gt. 1) b(n,1) = epsb                                       
      return                                                            
      end                                                               
****************************************
      subroutine qzval(nm,n,a,b,alfr,alfi,beta,matz,z)                  
                                                                        
      integer i,j,n,en,na,nm,nn,isw                                     
      real*8 a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)            
      real*8 c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1,u2,        
     x       v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr,         
     x       ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb         
      real*8 dsqrt,dabs,dsign                                           
      logical matz                                                      
c                                                                       
c     this subroutine is the third step of the qz algorithm             
c     for solving generalized matrix eigenvalue problems,               
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.      
c                                                                       
c     this subroutine accepts a pair of real matrices, one of them      
c     in quasi-triangular form and the other in upper triangular form.  
c     it reduces the quasi-triangular matrix further, so that any       
c     remaining 2-by-2 blocks correspond to pairs of complex            
c     eigenvalues, and returns quantities whose ratios give the         
c     generalized eigenvalues.  it is usually preceded by  qzhes        
c     and  qzit  and may be followed by  qzvec.                         
c                                                                       
c     on input:                                                         
c                                                                       
c        nm must be set to the row dimension of two-dimensional         
c          array parameters as declared in the calling program          
c          dimension statement;                                         
c                                                                       
c        n is the order of the matrices;                                
c                                                                       
c        a contains a real upper quasi-triangular matrix;               
c                                                                       
c        b contains a real upper triangular matrix.  in addition,       
c          location b(n,1) contains the tolerance quantity (epsb)       
c          computed and saved in  qzit;                                 
c                                                                       
c        matz should be set to .true. if the right hand transformations 
c          are to be accumulated for later use in computing             
c          eigenvectors, and to .false. otherwise;                      
c                                                                       
c        z contains, if matz has been set to .true., the                
c          transformation matrix produced in the reductions by qzhes    
c          and qzit, if performed, or else the identity matrix.         
c          if matz has been set to .false., z is not referenced.        
c                                                                       
c     on output:                                                        
c                                                                       
c        a has been reduced further to a quasi-triangular matrix        
c          in which all nonzero subdiagonal elements correspond to      
c          pairs of complex eigenvalues;                                
c                                                                       
c        b is still in upper triangular form, although its elements     
c          have been altered.  b(n,1) is unaltered;                     
c                                                                       
c        alfr and alfi contain the real and imaginary parts of the      
c          diagonal elements of the triangular matrix that would be     
c          obtained if a were reduced completely to triangular form     
c          by unitary transformations.  non-zero values of alfi occur   
c          in pairs, the first member positive and the second negative; 
c                                                                       
c        beta contains the diagonal elements of the corresponding b,    
c          normalized to be real and non-negative.  the generalized     
c          eigenvalues are then the ratios ((alfr+i*alfi)/beta);        
c                                                                       
c        z contains the product of the right hand transformations       
c          (for all three steps) if matz has been set to .true.         
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
      epsb = b(n,1)                                                     
      isw = 1                                                           
c     :::::::::: find eigenvalues of quasi-triangular matrices.         
c                for en=n step -1 until 1 do -- ::::::::::              
      do 510 nn = 1, n                                                  
         en = n + 1 - nn                                                
         na = en - 1                                                    
         if (isw .eq. 2) go to 505                                      
         if (en .eq. 1) go to 410                                       
         if (a(en,na) .ne. 0.0d0) go to 420                             
c     :::::::::: 1-by-1 block, one real root ::::::::::                 
  410    alfr(en) = a(en,en)                                            
         if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)                  
         beta(en) = dabs(b(en,en))                                      
         alfi(en) = 0.0d0                                               
         go to 510                                                      
c     :::::::::: 2-by-2 block ::::::::::                                
  420    if (dabs(b(na,na)) .le. epsb) go to 455                        
         if (dabs(b(en,en)) .gt. epsb) go to 430                        
         a1 = a(en,en)                                                  
         a2 = a(en,na)                                                  
         bn = 0.0d0                                                     
         go to 435                                                      
  430    an = dabs(a(na,na)) + dabs(a(na,en)) + dabs(a(en,na))          
     x      + dabs(a(en,en))                                            
         bn = dabs(b(na,na)) + dabs(b(na,en)) + dabs(b(en,en))          
         a11 = a(na,na) / an                                            
         a12 = a(na,en) / an                                            
         a21 = a(en,na) / an                                            
         a22 = a(en,en) / an                                            
         b11 = b(na,na) / bn                                            
         b12 = b(na,en) / bn                                            
         b22 = b(en,en) / bn                                            
         e = a11 / b11                                                  
         ei = a22 / b22                                                 
         s = a21 / (b11 * b22)                                          
         t = (a22 - e * b22) / b22                                      
         if (dabs(e) .le. dabs(ei)) go to 431                           
         e = ei                                                         
         t = (a11 - e * b11) / b11                                      
  431    c = 0.5d0 * (t - s * b12)                                      
         d = c * c + s * (a12 - e * b12)                                
         if (d .lt. 0.0d0) go to 480                                    
c     :::::::::: two real roots.                                        
c                zero both a(en,na) and b(en,na) ::::::::::             
         e = e + (c + dsign(dsqrt(d),c))                                
         a11 = a11 - e * b11                                            
         a12 = a12 - e * b12                                            
         a22 = a22 - e * b22                                            
         if (dabs(a11) + dabs(a12) .lt.                                 
     x       dabs(a21) + dabs(a22)) go to 432                           
         a1 = a12                                                       
         a2 = a11                                                       
         go to 435                                                      
  432    a1 = a22                                                       
         a2 = a21                                                       
c     :::::::::: choose and apply real z ::::::::::                     
  435    s = dabs(a1) + dabs(a2)                                        
         u1 = a1 / s                                                    
         u2 = a2 / s                                                    
         r = dsign(dsqrt(u1*u1+u2*u2),u1)                               
         v1 = -(u1 + r) / r                                             
         v2 = -u2 / r                                                   
         u2 = v2 / v1                                                   
                                                                        
         do 440 i = 1, en                                               
            t = a(i,en) + u2 * a(i,na)                                  
            a(i,en) = a(i,en) + t * v1                                  
            a(i,na) = a(i,na) + t * v2                                  
            t = b(i,en) + u2 * b(i,na)                                  
            b(i,en) = b(i,en) + t * v1                                  
            b(i,na) = b(i,na) + t * v2                                  
  440    continue                                                       
                                                                        
         if (.not. matz) go to 450                                      
                                                                        
         do 445 i = 1, n                                                
            t = z(i,en) + u2 * z(i,na)                                  
            z(i,en) = z(i,en) + t * v1                                  
            z(i,na) = z(i,na) + t * v2                                  
  445    continue                                                       
                                                                        
  450    if (bn .eq. 0.0d0) go to 475                                   
         if (an .lt. dabs(e) * bn) go to 455                            
         a1 = b(na,na)                                                  
         a2 = b(en,na)                                                  
         go to 460                                                      
  455    a1 = a(na,na)                                                  
         a2 = a(en,na)                                                  
c     :::::::::: choose and apply real q ::::::::::                     
  460    s = dabs(a1) + dabs(a2)                                        
         if (s .eq. 0.0d0) go to 475                                    
         u1 = a1 / s                                                    
         u2 = a2 / s                                                    
         r = dsign(dsqrt(u1*u1+u2*u2),u1)                               
         v1 = -(u1 + r) / r                                             
         v2 = -u2 / r                                                   
         u2 = v2 / v1                                                   
                                                                        
         do 470 j = na, n                                               
            t = a(na,j) + u2 * a(en,j)                                  
            a(na,j) = a(na,j) + t * v1                                  
            a(en,j) = a(en,j) + t * v2                                  
            t = b(na,j) + u2 * b(en,j)                                  
            b(na,j) = b(na,j) + t * v1                                  
            b(en,j) = b(en,j) + t * v2                                  
  470    continue                                                       
                                                                        
  475    a(en,na) = 0.0d0                                               
         b(en,na) = 0.0d0                                               
         alfr(na) = a(na,na)                                            
         alfr(en) = a(en,en)                                            
         if (b(na,na) .lt. 0.0d0) alfr(na) = -alfr(na)                  
         if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)                  
         beta(na) = dabs(b(na,na))                                      
         beta(en) = dabs(b(en,en))                                      
         alfi(en) = 0.0d0                                               
         alfi(na) = 0.0d0                                               
         go to 505                                                      
c     :::::::::: two complex roots ::::::::::                           
  480    e = e + c                                                      
         ei = dsqrt(-d)                                                 
         a11r = a11 - e * b11                                           
         a11i = ei * b11                                                
         a12r = a12 - e * b12                                           
         a12i = ei * b12                                                
         a22r = a22 - e * b22                                           
         a22i = ei * b22                                                
         if (dabs(a11r) + dabs(a11i) + dabs(a12r) + dabs(a12i) .lt.     
     x       dabs(a21) + dabs(a22r) + dabs(a22i)) go to 482             
         a1 = a12r                                                      
         a1i = a12i                                                     
         a2 = -a11r                                                     
         a2i = -a11i                                                    
         go to 485                                                      
  482    a1 = a22r                                                      
         a1i = a22i                                                     
         a2 = -a21                                                      
         a2i = 0.0d0                                                    
c     :::::::::: choose complex z ::::::::::                            
  485    cz = dsqrt(a1*a1+a1i*a1i)                                      
         if (cz .eq. 0.0d0) go to 487                                   
         szr = (a1 * a2 + a1i * a2i) / cz                               
         szi = (a1 * a2i - a1i * a2) / cz                               
         r = dsqrt(cz*cz+szr*szr+szi*szi)                               
         cz = cz / r                                                    
         szr = szr / r                                                  
         szi = szi / r                                                  
         go to 490                                                      
  487    szr = 1.0d0                                                    
         szi = 0.0d0                                                    
  490    if (an .lt. (dabs(e) + ei) * bn) go to 492                     
         a1 = cz * b11 + szr * b12                                      
         a1i = szi * b12                                                
         a2 = szr * b22                                                 
         a2i = szi * b22                                                
         go to 495                                                      
  492    a1 = cz * a11 + szr * a12                                      
         a1i = szi * a12                                                
         a2 = cz * a21 + szr * a22                                      
         a2i = szi * a22                                                
c     :::::::::: choose complex q ::::::::::                            
  495    cq = dsqrt(a1*a1+a1i*a1i)                                      
         if (cq .eq. 0.0d0) go to 497                                   
         sqr = (a1 * a2 + a1i * a2i) / cq                               
         sqi = (a1 * a2i - a1i * a2) / cq                               
         r = dsqrt(cq*cq+sqr*sqr+sqi*sqi)                               
         cq = cq / r                                                    
         sqr = sqr / r                                                  
         sqi = sqi / r                                                  
         go to 500                                                      
  497    sqr = 1.0d0                                                    
         sqi = 0.0d0                                                    
c     :::::::::: compute diagonal elements that would result            
c                if transformations were applied ::::::::::             
  500    ssr = sqr * szr + sqi * szi                                    
         ssi = sqr * szi - sqi * szr                                    
         i = 1                                                          
         tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21           
     x      + ssr * a22                                                 
         ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22               
         dr = cq * cz * b11 + cq * szr * b12 + ssr * b22                
         di = cq * szi * b12 + ssi * b22                                
         go to 503                                                      
  502    i = 2                                                          
         tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21               
     x      + cq * cz * a22                                             
         ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21              
         dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22                
         di = -ssi * b11 - sqi * cz * b12                               
  503    t = ti * dr - tr * di                                          
         j = na                                                         
         if (t .lt. 0.0d0) j = en                                       
         r = dsqrt(dr*dr+di*di)                                         
         beta(j) = bn * r                                               
         alfr(j) = an * (tr * dr + ti * di) / r                         
         alfi(j) = an * t / r                                           
         if (i .eq. 1) go to 502                                        
  505    isw = 3 - isw                                                  
  510 continue                                                          
                                                                        
      return                                                            
      end                                                               
*********************************************
      subroutine qzvec(nm,n,a,b,alfr,alfi,beta,z)                       
c                                                                       
      integer i,j,k,m,n,en,ii,jj,na,nm,nn,isw,enm2                      
      real*8 a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)            
      real*8 d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2,w1,x1,zz,z1,    
     x       alfm,almi,almr,betm,epsb                                   
      real*8 dsqrt,dabs                                                 
c                                                                       
c     this subroutine is the optional fourth step of the qz algorithm   
c     for solving generalized matrix eigenvalue problems,               
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.      
c                                                                       
c     this subroutine accepts a pair of real matrices, one of them in   
c     quasi-triangular form (in which each 2-by-2 block corresponds to  
c     a pair of complex eigenvalues) and the other in upper triangular  
c     form.  it computes the eigenvectors of the triangular problem and 
c     transforms the results back to the original coordinate system.    
c     it is usually preceded by  qzhes,  qzit, and  qzval.              
c                                                                       
c     on input:                                                         
c                                                                       
c        nm must be set to the row dimension of two-dimensional         
c          array parameters as declared in the calling program          
c          dimension statement;                                         
c                                                                       
c        n is the order of the matrices;                                
c                                                                       
c        a contains a real upper quasi-triangular matrix;               
c                                                                       
c        b contains a real upper triangular matrix.  in addition,       
c          location b(n,1) contains the tolerance quantity (epsb)       
c          computed and saved in  qzit;                                 
c                                                                       
c        alfr, alfi, and beta  are vectors with components whose        
c          ratios ((alfr+i*alfi)/beta) are the generalized              
c          eigenvalues.  they are usually obtained from  qzval;         
c                                                                       
c        z contains the transformation matrix produced in the           
c          reductions by  qzhes,  qzit, and  qzval, if performed.       
c          if the eigenvectors of the triangular problem are            
c          desired, z must contain the identity matrix.                 
c                                                                       
c     on output:                                                        
c                                                                       
c        a is unaltered.  its subdiagonal elements provide information  
c           about the storage of the complex eigenvectors;              
c                                                                       
c        b has been destroyed;                                          
c                                                                       
c        alfr, alfi, and beta are unaltered;                            
c                                                                       
c        z contains the real and imaginary parts of the eigenvectors.   
c          if alfi(i) .eq. 0.0, the i-th eigenvalue is real and         
c            the i-th column of z contains its eigenvector.             
c          if alfi(i) .ne. 0.0, the i-th eigenvalue is complex.         
c            if alfi(i) .gt. 0.0, the eigenvalue is the first of        
c              a complex pair and the i-th and (i+1)-th columns         
c              of z contain its eigenvector.                            
c            if alfi(i) .lt. 0.0, the eigenvalue is the second of       
c              a complex pair and the (i-1)-th and i-th columns         
c              of z contain the conjugate of its eigenvector.           
c          each eigenvector is normalized so that the modulus           
c          of its largest component is 1.0 .                            
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
      epsb = b(n,1)                                                     
      isw = 1                                                           
c     :::::::::: for en=n step -1 until 1 do -- ::::::::::              
      do 800 nn = 1, n                                                  
         en = n + 1 - nn                                                
         na = en - 1                                                    
         if (isw .eq. 2) go to 795                                      
         if (alfi(en) .ne. 0.0d0) go to 710                             
c     :::::::::: real vector ::::::::::                                 
         m = en                                                         
         b(en,en) = 1.0d0                                               
         if (na .eq. 0) go to 800                                       
         alfm = alfr(m)                                                 
         betm = beta(m)                                                 
          if(dabs(betm).lt.1.d-12) goto 800
c        !write(6,7711) m,alfm,betm
 7711    format(1x,i5,2(1x,1pd22.15))
c     :::::::::: for i=en-1 step -1 until 1 do -- ::::::::::            
         do 700 ii = 1, na                                              
            i = en - ii                                                 
            w = betm * a(i,i) - alfm * b(i,i)                           
            r = 0.0d0                                                   
                                                                        
            do 610 j = m, en                                            
            r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)           
cccccc             !write(6,7722) j,m,b(i,j),b(j,en),r
 7722       format(1x,'loop 610',2x,2i5,3(1x,1pd12.5))
  610       continue
                                                                        
            if (i .eq. 1 .or. isw .eq. 2) go to 630                     
            if (betm * a(i,i-1) .eq. 0.0d0) go to 630                   
            zz = w                                                      
            s = r                                                       
            go to 690                                                   
  630       m = i                                                       
            if (isw .eq. 2) go to 640                                   
c     :::::::::: real 1-by-1 block ::::::::::                           
            t = w                                                       
            if (w .eq. 0.0d0) t = epsb                                  
            b(i,en) = -r / t                                            
            go to 700                                                   
c     :::::::::: real 2-by-2 block ::::::::::                           
  640       x = betm * a(i,i+1) - alfm * b(i,i+1)                       
            y = betm * a(i+1,i)                                         
            q = w * zz - x * y                                          
            t = (x * s - zz * r) / q                                    
            b(i,en) = t                                                 
            if (dabs(x) .le. dabs(zz)) go to 650                        
            b(i+1,en) = (-r - w * t) / x                                
ccccc       !write(6,7733) r,w,t,x
 7733       format(1x,'640-650',4(1x,1pd12.5))
            go to 690                                                   
  650       b(i+1,en) = (-s - y * t) / zz                               
  690       isw = 3 - isw                                               
  700    continue                                                       
c     :::::::::: end real vector ::::::::::                             
         go to 800                                                      
c     :::::::::: complex vector ::::::::::                              
  710    m = na                                                         
         almr = alfr(m)                                                 
         almi = alfi(m)                                                 
         betm = beta(m)                                                 
c     :::::::::: last vector component chosen imaginary so that         
c                eigenvector matrix is triangular ::::::::::            
         y = betm * a(en,na)                                            
         b(na,na) = -almi * b(en,en) / y                                
         b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y             
         b(en,na) = 0.0d0                                               
         b(en,en) = 1.0d0                                               
         enm2 = na - 1                                                  
         if (enm2 .eq. 0) go to 795                                     
c     :::::::::: for i=en-2 step -1 until 1 do -- ::::::::::            
         do 790 ii = 1, enm2                                            
            i = na - ii                                                 
            w = betm * a(i,i) - almr * b(i,i)                           
            w1 = -almi * b(i,i)                                         
            ra = 0.0d0                                                  
            sa = 0.0d0                                                  
                                                                        
            do 760 j = m, en                                            
               x = betm * a(i,j) - almr * b(i,j)                        
               x1 = -almi * b(i,j)                                      
               ra = ra + x * b(j,na) - x1 * b(j,en)                     
               sa = sa + x * b(j,en) + x1 * b(j,na)                     
  760       continue                                                    
                                                                        
            if (i .eq. 1 .or. isw .eq. 2) go to 770                     
            if (betm * a(i,i-1) .eq. 0.0d0) go to 770                   
            zz = w                                                      
            z1 = w1                                                     
            r = ra                                                      
            s = sa                                                      
            isw = 2                                                     
            go to 790                                                   
  770       m = i                                                       
            if (isw .eq. 2) go to 780                                   
c     :::::::::: complex 1-by-1 block ::::::::::                        
            tr = -ra                                                    
            ti = -sa                                                    
  773       dr = w                                                      
            di = w1                                                     
c     :::::::::: complex divide (t1,t2) = (tr,ti) / (dr,di) ::::::::::  
  775       if (dabs(di) .gt. dabs(dr)) go to 777                       
            rr = di / dr                                                
            d = dr + di * rr                                            
            t1 = (tr + ti * rr) / d                                     
            t2 = (ti - tr * rr) / d                                     
            go to (787,782), isw                                        
  777       rr = dr / di                                                
            d = dr * rr + di                                            
            t1 = (tr * rr + ti) / d                                     
            t2 = (ti * rr - tr) / d                                     
            go to (787,782), isw                                        
c     :::::::::: complex 2-by-2 block ::::::::::                        
  780       x = betm * a(i,i+1) - almr * b(i,i+1)                       
            x1 = -almi * b(i,i+1)                                       
            y = betm * a(i+1,i)                                         
            tr = y * ra - w * r + w1 * s                                
            ti = y * sa - w * s - w1 * r                                
            dr = w * zz - w1 * z1 - x * y                               
            di = w * z1 + w1 * zz - x1 * y                              
            if (dr .eq. 0.0d0 .and. di .eq. 0.0d0) dr = epsb            
            go to 775                                                   
  782       b(i+1,na) = t1                                              
            b(i+1,en) = t2                                              
            isw = 1                                                     
            if (dabs(y) .gt. dabs(w) + dabs(w1)) go to 785              
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)                   
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)                   
            go to 773                                                   
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y             
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y             
  787       b(i,na) = t1                                                
            b(i,en) = t2                                                
  790    continue                                                       
c     :::::::::: end complex vector ::::::::::                          
  795    isw = 3 - isw                                                  
  800 continue                                                          
c     :::::::::: end back substitution.                                 
c                transform to original coordinate system.               
c                for j=n step -1 until 1 do -- ::::::::::               
      do 880 jj = 1, n                                                  
         j = n + 1 - jj                                                 
                                                                        
         do 880 i = 1, n                                                
            zz = 0.0d0                                                  
                                                                        
            do 860 k = 1, j                                             
  860       zz = zz + z(i,k) * b(k,j)                                   
                                                                        
            z(i,j) = zz                                                 
  880 continue                                                          
c     :::::::::: normalize so that modulus of largest                   
c                component of each vector is 1.                         
c                (isw is 1 initially from before) ::::::::::            
      do 950 j = 1, n                                                   
         d = 0.0d0                                                      
         if (isw .eq. 2) go to 920                                      
         if (alfi(j) .ne. 0.0d0) go to 945                              
                                                                        
         do 890 i = 1, n                                                
            if (dabs(z(i,j)) .gt. d) d = dabs(z(i,j))                   
  890    continue                                                       
                                                                        
         do 900 i = 1, n                                                
  900    z(i,j) = z(i,j) / d                                            
                                                                        
         go to 950                                                      
                                                                        
  920    do 930 i = 1, n                                                
            r = dabs(z(i,j-1)) + dabs(z(i,j))                           
            if (r .ne. 0.0d0) r = r * dsqrt((z(i,j-1)/r)**2             
     x                                     +(z(i,j)/r)**2)              
            if (r .gt. d) d = r                                         
  930    continue                                                       
                                                                        
         do 940 i = 1, n                                                
            z(i,j-1) = z(i,j-1) / d                                     
            z(i,j) = z(i,j) / d                                         
  940    continue                                                       
                                                                        
  945    isw = 3 - isw                                                  
  950 continue                                                          
                                                                        
      return                                                            
      end                                                               
c*************************************************************
      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      implicit real*8(a-h,o-z)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
CU    USES chebev
      double precision xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     *-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/

      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      implicit real*8 (a-h,o-z)
      INTEGER MAXIT
      double precision rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=200000,XMIN=2.d0,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     *f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,
     *r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
     *temp,w,x2,xi,xi2,xmu,xmu2

      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      write(*,*)
      write(*,*)'x =',x
      pause 'x too large in bessjy; try asymptotic expansion'
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        pause 'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1

 123  return
      END

      Double Precision FUNCTION chebev(a,b,c,m,x)
      implicit real*8(a-h,o-z)
      INTEGER m
      double precision a,b,x,c(m)
      INTEGER j
      double precision d,dd,sv,y,y2

      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END

      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      implicit real*8(a-h,o-z)
      double precision sj,sjp,sy,syp,x,n
CU    USES bessjy
      double precision factor,order,rj,rjp,ry,ryp,RTPIO2

      pi = acos(-1.0d0)
      rtpio2 = dsqrt(0.5d0*pi)
      if(n.lt.0.d0.or.x.le.0.d0)pause 'bad arguments in sphbes' 
      order=n+0.5
      if (x.le.1.86d5) then
      call bessjy(x,order,rj,ry,rjp,ryp)
      else
      call BesselAsymNew(order,x,rj,ry,rjp,ryp)
      endif
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.d0*x)
      syp=factor*ryp-sy/(2.d0*x)
      return
      END




c**************************************************************

      subroutine BesselAsymNew(v,x,sj,sy,sjp,syp)
      integer k,b,kstop
      double precision v,x,sj,sy,sjp,syp,u
      double precision P,Q,R,S,pi,factor,chi
      double precision denomP,denomQ,denomP2,denomQ2
      double precision Pfunction,Qfunction,sign
      double precision coefR,coefS,factrl

c     asymptotic expansions from Abramowitz and Stegun, pg. 364
c     following assumes v = n+1/2 where n=0,1,2,......

      pi = dacos(-1.0d0)
      factor = dsqrt(2.0d0/(pi*x))
      chi = x - pi*(0.5d0*v + 0.25d0)
      u = 4.0d0*v*v

      P = 0.0d0
      R = 0.0d0
      do k=0,int(v/2.d0-0.25d0)
         sign = (-1.0d0)**k
         denomP = (2.0d0*x)**(2*k)
         denomP2 = factrl(2*k)
         coefR = (u-1.0d0+16.0d0*dfloat(k*k))/
     .       (u-(4.0d0*dfloat(k)-1.0d0)*(4.0d0*dfloat(k)-1.0d0))
         Pfunction = 1.0d0
         do b=0,4*k-1
          Pfunction = Pfunction*(v-0.5d0+float(2*k-b))
         enddo
      P = P + sign*Pfunction/(denomP*denomP2)
      R = R + sign*coefR*Pfunction/(denomP*denomP2)
      enddo

      Q = 0.0d0
      S = 0.0d0
      do k=0,int(v/2.d0-0.75d0)
         sign = (-1.0d0)**k
         denomQ = denomP*(2.0d0*x)
         denomQ2 = factrl(2*k+1)
         coefS = (u+4.0d0*(2.0d0*dfloat(k)+1.0d0)*(2.0d0*dfloat(k)+1.0d0)-
     .           1.0d0)/(u-(4.0d0*dfloat(k)+1.0d0)*(4.0d0*dfloat(k)+1.0d0))
         Qfunction = 1.0d0
         do b=0,4*k+1
         Qfunction = Qfunction*(v+0.5d0+float(2*k-b))
         enddo       
      Q = Q + sign*Qfunction/(denomQ*denomQ2)
      S = S + sign*coefS*Qfunction/(denomQ*denomQ2)
      enddo


      if (v.eq.0.5d0) then
         R = 1.d0
         S = 1.d0/(2.d0*x)
      endif

      sj = factor*(P*dcos(chi) - Q*dsin(chi)) 
      sy = factor*(P*dsin(chi) + Q*dcos(chi))
      sjp=-factor*(R*dsin(chi) - S*dcos(chi))
      syp= factor*(R*dcos(chi) - S*dsin(chi))

      return
      end


      subroutine BesselAsym(v,x,sj,sy,sjp,syp)
      integer k,b
      double precision v,x,sj,sy,sjp,syp,u
      double precision P,Q,R,S,pi,factor,chi
      double precision denomP,denomQ,denomP2,denomQ2
      double precision Pfunction,Qfunction,sign
      double precision coefR,coefS,factrl

c asymptotic expansions from Abramowitz and Stegun, pg. 364
c following assumes v = n+1/2 where n=0,1,2,......

      pi = dacos(-1.0d0)
      factor = dsqrt(2.0d0/(pi*x))
      chi = x - pi*(0.5d0*v + 0.25d0)
      u = 4.0d0*v*v

      P = 0.0d0
      Q = 0.0d0
      R = 0.0d0
      S = 0.0d0
      do k=0,int(v-0.5)
       !write(6,*)'Doing Loop'
       sign = (-1.0d0)**k
       denomP = (2.0d0*x)**(2*k)
       denomQ = denomP*(2.0d0*x)
       denomP2 = factrl(2*k)
       denomQ2 = factrl(2*k+1)
       coefR = (u-1.0d0+16.0d0*float(k*k))/
     > (u-(4.0d0*float(k)-1.0d0)*(4.0d0*float(k)-1.0d0))
       coefS = (u+4.0d0*(2.0d0*float(k)+1.0d0)*(2.0d0*float(k)+1.0d0)-
     >          1.0d0)/(u-(4.0d0*float(k)+1.0d0)*(4.0d0*float(k)+1.0d0))
       !write(6,*)'CoefS',coefR,coefS,k

       Pfunction = 1.0d0
       do b=0,4*k-1
        Pfunction = Pfunction*(v-0.5d0+float(2*k-b))
       enddo
       Qfunction = 1.0d0
       do b=0,4*k+1
        Qfunction = Qfunction*(v+0.5d0+float(2*k-b))
       enddo       
     
       P = P + sign*Pfunction/(denomP*denomP2)
       Q = Q + sign*Qfunction/(denomQ*denomQ2)
       R = R + sign*coefR*Pfunction/(denomP*denomP2)
       S = S + sign*coefS*Qfunction/(denomQ*denomQ2)
      enddo

      !write(6,*)'R,S',R,S*2.0d0*x
      sj = factor*(P*dcos(chi) - Q*dsin(chi)) 
      sy = factor*(P*dsin(chi) + Q*dcos(chi))
      sjp=-factor*(R*dsin(chi) - S*dcos(chi))
      syp= factor*(R*dcos(chi) - S*dsin(chi))
      
      return
      end

      double precision FUNCTION factrl(n)
      INTEGER n
CU    USES gammln
      INTEGER j,ntop
      double precision a(33),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      if (n.lt.0) then
        pause 'negative factorial in factrl'
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(float(n)+1.0d0))
      endif
      return
      END

      double precision FUNCTION gammln(xx)
      double precision xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
