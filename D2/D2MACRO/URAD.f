      program URADmain
      LOGICAL DEBUG
      COMMON /DEBUG/ DEBUG

      parameter (NPDIM=400)
      real Epdata(NPDIM),Mfdata(NPDIM),
     &     srdata(NPDIM),sndata(NPDIM)
      character*32 file

*     -- init /DEBUG/
      DEBUG = .TRUE.
*     -- init /MoTsai/
      call iMoTsai
*     -- init /URADMT/
      call iURADMT

      call ClearMT

      np = 200
      Es = .411
      file = 'Es411.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)
      Es = .510
      file = 'Es510.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)
      Es = .590
      file = 'Es590.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)
      Es = .670
      file = 'Es670.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)

      del = .001
      print*, 'del =', del, '  Enter del (/=unchanged)'
      read*, del
      call URAD(del)

      call WAIT('Write to file(s)')
      call WrURAD
      END

*                                         @METAGS .URAD
*                                         12-02-98 10:51am
*--------------- URAD ---------------
*
      subroutine URAD(del)
      IMPLICIT NONE
*     .. scattering angle thetaMT supplies common /MoTsai/
      real del
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      LOGICAL DEBUG
      COMMON /DEBUG/ DEBUG
*      use MoTsai,TradMT,URADMT, DEBUG
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2 =.511e-3*.511e-3)
      real eps
      parameter (eps=0.00001)
      real Kmin,Mfmin
      real lnEsdel,lnEpdel,lnqme2, IdEs,IdEp
      real fEs,fEp,Spence
      external fEs,fEp
      real Es0,Ep0,b,bw,Es,Ep,sp,tr,fs,fp,Epmax,Esmin,termEs,termEp
      real deltat,deltar,rer,aer,errest,Epl,Eph,Esl,Esh
      real costheta,flag,sum
      integer NS,NSminMf,np,nofun

      real     fEstst,fEptst
      external fEstst,fEptst
      logical testEs,testEp
      testEs   = .FALSE.
      testEp   = .FALSE.
*      testEs   = .TRUE.
*      testEp   = .TRUE.

*     .. save old kinematics of common /MoTsai/
      Es0 = EsMT
      Ep0 = EpMT

      print*, 'URAD: Esdat(NS):', Esdat
      print*, 'URAD: nptot(NS):', nptot
      print*, 'theta =', thetaMT
      print*, 'Mtar =', Mtar
      print*, 'mthr =', mthr
      print*, 'Trad   =', Trad
      print*, 'tiwrad =', tiwrad
      print*, 'tfwrad =', tfwrad
      print*, ' '

      b  = bMT
      bw = bwMT
      costheta = cos(thetaMT*pi/180.)

*     .. minimal energy of equivalent gamma-quantum
      Kmin = ((Mtar+mthr)**2 - Mtar**2)/(2.*Mtar)

*     .. estimate the two highest points using only the first term of (IV.2)
*     .. Note: for the second point we can use some interpolation
      DO NS=1,NStot
      Es = Esdat(NS)
      print*, 'Es =', Es, '   np, Ep, srdat, sudat, sndat:'
      do np=nptot(NS), nptot(NS)-1, -1
         Ep = Epdat(np,NS)
         sp = Es*Ep*(1.-costheta)
         lnqme2 = log(2.*sp/me2)

*        .. equivalent radiator
         tr = (alp/pi)*(lnqme2-1.) / b
*        .. fs,fp terms
         fs = b*tr + bw*tiwrad + b*Trad/2.
         fp = b*tr + bw*tfwrad + b*Trad/2.

         lnEsdel = log(Es/del)
         lnEpdel = log(Ep/del)
         deltat = -( (bw*tiwrad+b*Trad/2.)*lnEsdel +
     +               (bw*tfwrad+b*Trad/2.)*lnEpdel )

         deltar = -(alp/pi)*( 28./9. - (13./6)*lnqme2 +
     +                        (lnEsdel+lnEpdel)*(lnqme2-1.) -
     -                        Spence(-(Es-Ep)/Ep) - Spence((Es-Ep)/Es) )

         sudat(np,NS) = srdat(np,NS)*exp(-deltat-deltar)
         print*, np, Ep, srdat(np,NS), sudat(np,NS), sndat(np,NS)
      enddo
*     .. Note: really in this initial procedure we unfolded two points
      npcur(NS) = nptot(NS)
      ENDDO

*     .. Point ncdata is restored, ncdata-1 is estimated here.

      call DWAIT('Restore the whole spectrum in a common way')
10000 continue
*     .. look for the any spectrum with non unfolded point
      do NS=1,NStot
         if (npcur(NS).GT.1) goto 11000
      enddo

*     .. If we just here then all specra are unfolded
      print*, 'Work is complete'
      goto 99999

11000 NSminMf=NS
*     .. find the number of spectrum have unfolded point with min Mf
      Mfmin = Mfdat(npcur(NSminMf)-1,NSminMf)
      do NS=1,NStot
*        .. the last unfolded point
         if (npcur(NS).EQ.1) goto 10100
         if (Mfdat(npcur(NS)-1,NS).LT.Mfmin) NSminMf=NS
10100 enddo

*     .. current spectrum No.
      NS = NSminMf
*     .. store in common URADMT the current spectrum No.
      NScur = NS
*     .. No. of current unfolding point
      np = npcur(NScur)-1

      Es = Esdat(NS)
      Epmax = Epma(NS)
      Ep = Epdat(np,NS)
*     .. store in common /MoTsai/ Es,Ep of current point
      EsMT = Es
      EpMT = Ep

*      print*, 'Restoring point: NS,np,Ep,Epmax', NS,np,Ep,Epmax

      sp = Es*Ep*(1.-costheta)
      lnqme2 = log(2.*sp/me2)

*     .. equivalent radiator
      tr = (alp/pi)*(lnqme2-1.) / b
*     .. fs,fp terms
      fs = b*tr + bw*tiwrad + b*Trad/2.
      fp = b*tr + bw*tfwrad + b*Trad/2.
*     .. store in common /MoTsai/ for following usage by fEs,fEp
      spMT = sp
      trMT = tr
      fsMT = fs
      fpMT = fp

      rer = 1.e-6
      aer = 0.
*     .. integral over Ep
      termEp = 0.
      Epl = Ep+del
      Eph = Epmax
*     .. see Remarks on Programming at p.235
      Eph = Eph - 1.*eps
      if(DEBUG)print*, 'for IdEp: Ep,Epl,Eph:', Ep,Epl,Eph
      if (Eph-Epl .LE. eps) then
         if(DEBUG)print*, '--- IdEp: Epl > Eph. Skip integration'
         goto 1000
      endif

      IF (testEp) THEN
*        .. for test use non-radiative cross section for Ep in int. dEp
         call quanc8(fEptst,Epl,Eph,aer,rer, IdEp,errest,nofun,flag)
         termEp = (del/Es)**(fs/2.)*IdEp
         goto 1000
      ENDIF

      call quanc8(fEp,Epl,Eph,aer,rer, IdEp,errest,nofun,flag)
      termEp = (del/Es)**(fs/2.)*IdEp
      if(DEBUG)print*, 'IdEp,termEp:', IdEp,termEp

 1000 continue

*     .. integral over Es
      termEs = 0.
      Esmin = (Ep + Kmin)/(1.-(Ep/Mtar)*(1-costheta))
      Esl = Esmin
*     .. see Remarks on Programming at p.235
      Esl = Esl + 1.*eps
      Esh = Es-del

      if(DEBUG)print*, 'for IdEs: Es,Esl,Esh:', Es,Esl,Esh
      if (Esh-Esl .LE. eps) then
         if(DEBUG)print*, '--- IdEs: Esl > Esh. Skip integration'
         goto 2000
      endif

      IF (testEs) THEN
*        .. for test use non-radiative cross section for Es in int. dEs
         call quanc8(fEstst,Esl,Esh,aer,rer, IdEs,errest,nofun,flag)
         termEs = (del/Ep)**(fp/2.)*IdEs
         goto 2000
      ENDIF

      call quanc8(fEs,Esl,Esh,aer,rer, IdEs,errest,nofun,flag)
      termEs = (del/Ep)**(fp/2.)*IdEs
      if(DEBUG)print*, 'IdEs,termEs:', IdEs,termEs

 2000 continue
      lnEsdel = log(Es/del)
      lnEpdel = log(Ep/del)
      deltat = -( (bw*tiwrad+b*Trad/2.)*lnEsdel +
     +            (bw*tfwrad+b*Trad/2.)*lnEpdel )

      deltar = -(alp/pi)*( 28./9. - (13./6)*lnqme2 +
     +                     (lnEsdel+lnEpdel)*(lnqme2-1.) -
     -                     Spence(-(Es-Ep)/Ep) - Spence((Es-Ep)/Es) )
      sum = (srdat(np,NS) - termEs - termEp)*exp(-deltat-deltar)

      if (sum.LE.0.) then
         print*, 'sum < 0.'
         print*,'termEp,termEs,srdat,sum,sndat for NS,np,Ep', NS,np,Ep
         print*,termEp,termEs,srdat(np,NS),sum,sndat(np,NS)
*         DEBUG=.TRUE.
         call WAIT('URAD')
         sum = 0.
      endif

      sudat(np,NS) = sum
*     .. point np is unfolded now. Store its number in npcur(NS)
      npcur(NS) = np

      if(DEBUG) then
         print*,'termEp,termEs,srdat,sudat,sndat for NS,np,Ep', NS,np,Ep
         print*,termEp,termEs,srdat(np,NS),sudat(np,NS),sndat(np,NS)
      else
      print*, 'NS,np,Ep,sudat,sndat',NS,np,Ep,sudat(np,NS),sndat(np,NS)
      endif
      call DWAIT('URAD: point is complete')

      goto 10000

99999 continue
*     .. restore old kinematics
      EsMT = Es0
      EpMT = Ep0
      END
*                                         @METAGS .fEs
*                                         12-01-98 01:06pm
*--------------- fEs ---------------
*
      real function fEs(Es1)
*     .. The integrand over Es1. Es,Ep,theta supplies commom /MoTsai/
      IMPLICIT NONE
      real Es1
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      LOGICAL DEBUG
      COMMON /DEBUG/ DEBUG
      real vEs(NSDIM),vsec(NSDIM)
*      use MoTsai,TradMT,URADMT
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2=.511e-3*.511e-3)
      real eps
      parameter (eps=0.00001)
      real sMf,slin
      real Is
      integer NS
      real Es,Ep,theta,b,bw,xs,sp,fs,ts,sec
      real Mf,Mf2
      integer ipoint

      fEs = 0.

***      if (Ep1.GE.Epma(NScur)) RETURN

*     .. Es,Ep,theta from common /MoTsai/
      Es = EsMT
      Ep = EpMT
      theta = thetaMT
      b  = bMT
      bw = bwMT

      xs = Es1/Es
*     .. get already calculated sp,fp from common /MoTsai/
      sp = spMT
      fs = fsMT

*     .. interpolate along the Mf already unfolded cross section
      Mf2=Mtar**2 + 2.*Mtar*(Es1-Ep) - 2.*Es1*Ep*(1.-cos(theta*pi/180.))
      if (Mf2.LT.0.) then
         print*, 'ERROR fEs: Mf2 < 0.'
         STOP
      endif
      Mf = sqrt(Mf2)
*      print*, 'fEs: Es1=', Es1, ' Mf=', Mf

      ipoint = 0
      do NS=1,NStot
         if (Mf.GT.Mfdat(1,NS)) goto 100
         sec = sMf(Mf,NS)
*        .. norm to Mott
         sec = sec*Esdat(NS)**2
         if (sec .LT. 0.) then
            print*, 'fEs: sec < 0: Out of range'
            call WAIT(' ')
         endif
*         if(DEBUG)print*, 'Mf,NS,sec', Mf,NS,sec
*         if (sec .GT. 0.) then
            ipoint = ipoint+1
            vEs(ipoint)  = Esdat(NS)
            vsec(ipoint) = sec
*         endif
  100 enddo

      if (ipoint.GT.0) then
         sec = slin(Es1,ipoint,vEs,vsec)
*        .. return from the Mott norm
         sec = sec/Es1**2
      else
         DEBUG = .TRUE.
         if(DEBUG)print*,'ipoint=0, Es1,Mf', Es1,Mf
         call DWAIT('ERROR fEs: ipoint=0.')
         fEs = 0.
         RETURN
      endif

      ts = (alp/pi)*(.5*(1.+xs**2)*log(2.*sp/me2) - xs)
*     .. probability to emit photon in the final state
      Is = (ts + (bw*tiwrad+b*Trad/2.)*(xs + (3./4.)*(1-xs)**2))*
     &     log(1./xs)**fs / (Es-Es1)

      fEs = Is*sec
      END
*                                         @METAGS .sMf
*                                         12-30-98 12:27pm
*--------------- sMf ---------------
*
      real function sMf(Mf,NS)
      IMPLICIT NONE
      integer NS
      real Mf
      real EsMT,EpMT,thetaMT,Mtar,mthr,spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr,spMT,trMT,fsMT,fpMT
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      LOGICAL DEBUG
      COMMON /DEBUG/ DEBUG
*      use MoTsai,URADMT
      real eps
      parameter (eps=0.00001)
      real sec,dS,dMf,Mfmin,Mfmax
      integer np,nlow,nhigh
      sMf = 0.
      if (Mf.LE.(Mtar+mthr)) then
         if(DEBUG)print*, 'sMf: Mf <= Mtar+mthr', Mf, ' <= ', Mtar+mthr
         call WAIT(' ')
         RETURN
      endif
      if (Mf.LE.Mfdat(npcur(NS),NS)) then
*        .. Note: Mf decrease along the spectrum
         do np=npcur(NS),nptot(NS)
            nlow = np
            if (Mfdat(nlow,NS).LT.Mf) goto 1000
         enddo
      else
*        .. extrapolate using points npcur(NS)-1, npcur(NS)
         nlow = npcur(NS)+1
      endif
 1000 nhigh = nlow-1
*     .. interpolation(extrapolation)
      dMf = Mfdat(nhigh,NS)-Mfdat(nlow,NS)
      if (dMf.GT.eps) then
         dS  = sudat(nhigh,NS)-sudat(nlow,NS)
         sec = sudat(nlow,NS) + dS*(Mf-Mfdat(nlow,NS))/dMf
      else
         sec = (sudat(nhigh,NS) + sudat(nlow,NS))/2.
      endif
      if (sec .LT. 0.) then
*         print*, 'NS=', NS, ' for Mf=', Mf, ' sec=', sec
         sec=0.
      endif
*     .. realible Mf range for extrapolation
*      dMf = Mfdat(1,NS)-Mfdat(nptot(NS),NS)
*      Mfmin = Mfdat(nptot(NS),NS) - 0.10*dMf
*      Mfmax = Mfdat(1,NS) + 0.10*dMf
      dMf = Mfdat(nhigh,NS)-Mfdat(nlow,NS)
      if (dMf.LT.0.) then
         print*, 'ERROR sMf: dMf < 0:'
         print*, 'NS,nlow,nhigh,Mfdat(nlow,NS),Mfdat(nhigh,NS)'
         print*, NS,nlow,nhigh,Mfdat(nlow,NS),Mfdat(nhigh,NS)
         STOP
      endif
      Mfmin = Mfdat(nlow,NS)  - 4.*dMf
      Mfmax = Mfdat(nhigh,NS) + 4.*dMf
**      if ((Mf.GE.Mfmin) .AND. (Mf.LE.Mfmax)) then
*      if (Mf.LE.Mfmax) then
*         sMf = sec
*      else
**        .. invert sign of cross section for "unrealible" extrapolation
*         sMf = -sec
*         if (DEBUG) then
*            print*, 'sMf: out of range. NS,nlow,nhigh:'
*            print*, NS,nlow,nhigh
*            print*, 'Mf,Mfdat(nlow,NS),Mfdat(nhigh,NS),Mfmin,Mfmax'
*            print*, Mf,Mfdat(nlow,NS),Mfdat(nhigh,NS),Mfmin,Mfmax
*            call DWAIT(' ')
*         endif
*      endif
      sMf = sec
      END

*                                         @METAGS .slin
*                                         12-30-98 05:48pm
*--------------- slin ---------------
*
      real function slin(Es1,np,Es,Sn)
*     .. linear interpolation between two closest incident energies
      IMPLICIT NONE
      integer np
      real Es1,Es(np),Sn(np)
      real dE,dS,sec
      integer i,nlow,nhigh
      slin = Sn(1)
      if (np.LT.2) RETURN
      if (Es1.GE.Es(1)) then
         do i=1,np
            nhigh = i
            if (Es(nhigh).GT.Es1) goto 1000
         enddo
      else
         nhigh=2
      endif
 1000 nlow = nhigh-1
      dE = Es(nhigh)-Es(nlow)
      dS = Sn(nhigh)-Sn(nlow)
      sec = Sn(nlow) + (Es1-Es(nlow))*dS/dE
      if (sec .LT. 0.) sec=0.
      slin = sec
      END

*                                         @METAGS .sinter
*                                         12-30-98 05:48pm
*--------------- sinter ---------------
*
      real function sinter(E0,np,Es,sec)
      IMPLICIT NONE
      integer np
      real E0,Es(np),sec(np)
      real a,b
      real x,sumx,sumx2,lny,sumlny,sumxlny,arg
      integer i
      sinter = 0.
      if (np.LT.2) RETURN
*     .. fit exp(a+b*x)
      sumlny  = 0.
      sumxlny = 0.
      sumx  = 0.
      sumx2 = 0.
      do i=1,np
         x = Es(i)
         sumx  = sumx  + x
         sumx2 = sumx2 + x**2
         lny = log(sec(i))
         sumlny  = sumlny  + lny
         sumxlny = sumxlny  + x*lny
      enddo
      b = (sumxlny - sumlny*sumx/np)/(sumx2-sumx**2/np)
      a = (sumlny - b*sumx)/np
*     .. get interpolated value
      arg = a + b*E0
      if (ABS(arg).GT.50.) RETURN
      sinter = exp(arg)
      END

*                                         @METAGS .fEp
*                                         12-01-98 01:06pm
*--------------- fEp ---------------
*
      real function fEp(Ep1)
*     .. The integrand over Ep1. Es,Ep,theta supplies commom /MoTsai/
      IMPLICIT NONE
      real Ep1
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
*      use MoTsai,TradMT,URADMT
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2=.511e-3*.511e-3)
      real eps
      parameter (eps=0.00001)
      real Ip
      integer NS,np,nlow,nhigh
      real Es,Ep,theta,b,bw,xp,sp,fp,tp,sec,dS,dE

*     .. current spectrum No.
      NS = NScur

      fEp = 0.
      if (Ep1.GE.Epma(NScur)) RETURN

*     .. Es,Ep,theta from common /MoTsai/
      Es = EsMT
      Ep = EpMT
      theta = thetaMT
      b  = bMT
      bw = bwMT

      xp = Ep/Ep1
*     .. get already calculated sp,fp from common /MoTsai/
      sp = spMT
      fp = fpMT

      tp = (alp/pi)*(.5*(1.+xp**2)*log(2.*sp/me2) - xp)
*     .. probability to emit photon in the final state
      Ip = (tp + (bw*tfwrad+b*Trad/2.)*(xp + (3./4.)*(1-xp)**2))*
     &     log(1./xp)**fp / (Ep1-Ep)

*     .. interpolate already unfolded cross section
      if (Epdat(npcur(NS),NS).LE.Ep1) then
*        .. search forward
         do np=npcur(NS),nptot(NS)
            nhigh = np
            if (Epdat(nhigh,NS).GT.Ep1) goto 1000
         enddo
      else
*        .. extrapolate using points npcur(NS), npcur(NS)+1
         nhigh = npcur(NS)+1
      endif

 1000 nlow = nhigh-1
*     .. interpolation(extrapolation)
      dE = Epdat(nhigh,NS)-Epdat(nlow,NS)
      if (dE.GT.eps) then
         dS = sudat(nhigh,NS)-sudat(nlow,NS)
         sec = sudat(nlow,NS) + (Ep1-Epdat(nlow,NS))*dS/dE
      else
         sec = (sudat(nhigh,NS) + sudat(nlow,NS))/2.
      endif 
      if (sec .LT. 0.) sec=0.
      fEp = Ip*sec
      END

*                                         @METAGS .Spence
*                                         12-02-98 12:07pm
*--------------- Spence ---------------
*
      FUNCTION Spence(X) 
*     .. Spence function
      parameter (pi2=3.141593*3.141593)
      IF (ABS(X).LE.1.) GO TO 1 
      IF (X.GT.1.) GO TO 2

      Y = -pi2/6. - .5*LOG(-X)**2 
      Y1=-1.
      X=1./X
      GOTO 3 

1     Y=0.
      Y1=1. 
      GOTO 3

2     Y =  pi2/3. - .5*LOG(X)**2 
      Y1=-1.
      X=1./X

3     CONTINUE
      N=1
4     Y1=Y1*X 
      Y0=Y1/real(N**2)
      Y=Y+Y0
      N=N+1
      IF(ABS(Y0).GT.ABS(.00001*Y)) GOTO 4 
      Spence=Y 
      END 

*                                         @METAGS .fEptst
*                                         12-01-98 01:06pm
*--------------- fEptst ---------------
*
      real function fEptst(Ep1)
*
*     Uses non-radiative cross sections for test
*
*     .. The integrand over Ep1. Es,Ep,theta supplies commom /MoTsai/
      IMPLICIT NONE
      real Ep1
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
*      use MoTsai,TradMT,URADMT
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2=.511e-3*.511e-3)
      real eps
      parameter (eps=0.00001)
      real Ip
*      integer NS,np,nlow,nhigh
      real Es,Ep,theta,b,bw,xp,sp,fp,tp,sec
*      real dS,dE
      real D2TOT

*     .. Es,Ep,theta from common /MoTsai/
      Es = EsMT
      Ep = EpMT
      theta = thetaMT
      b  = bMT
      bw = bwMT

      xp = Ep/Ep1
*     .. get already calculated sp,fp from common /MoTsai/
      sp = spMT
      fp = fpMT

      sec = D2TOT(Es,Ep1,theta)

      tp = (alp/pi)*(.5*(1.+xp**2)*log(2.*sp/me2) - xp)
*     .. probability to emit photon in the final state
      Ip = (tp + (bw*tfwrad+b*Trad/2.)*(xp + (3./4.)*(1-xp)**2))*
     &     log(1./xp)**fp / (Ep1-Ep)

      fEptst = Ip*sec
      END

*                                         @METAGS .fEstst
*                                         12-01-98 01:06pm
*--------------- fEstst ---------------
*
      real function fEstst(Es1)
*
*     Uses non-radiative cross sections for test
*
*     .. The integrand over Es1. Es,Ep,theta supplies commom /MoTsai/
      IMPLICIT NONE
      real Es1
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      LOGICAL DEBUG
      COMMON /DEBUG/ DEBUG
*      real vEs(NSDIM),vsec(NSDIM)
*      use MoTsai,TradMT,URADMT
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2=.511e-3*.511e-3)
      real eps
      parameter (eps=0.00001)
      real D2TOT
*      real sMf,slin
      real Is
*      integer NS
      real Es,Ep,theta,b,bw,xs,sp,fs,ts,sec
*      real Mf,Mf2
*      integer ipoint

      fEstst = 0.

***      if (Ep1.GE.Epma(NScur)) RETURN

*     .. Es,Ep,theta from common /MoTsai/
      Es = EsMT
      Ep = EpMT
      theta = thetaMT
      b  = bMT
      bw = bwMT

      xs = Es1/Es
*     .. get already calculated sp,fp from common /MoTsai/
      sp = spMT
      fs = fsMT

      sec = D2TOT(Es1,Ep,theta)

      ts = (alp/pi)*(.5*(1.+xs**2)*log(2.*sp/me2) - xs)
*     .. probability to emit photon in the final state
      Is = (ts + (bw*tiwrad+b*Trad/2.)*(xs + (3./4.)*(1-xs)**2))*
     &     log(1./xs)**fs / (Es-Es1)

      fEstst = Is*sec
      END
*                                         @METAGS .D2TOT
*                                         11-19-98 04:53pm
*--------------- D2TOT ---------------
*
      function D2TOT(E0,E,THETA)
*     .. D2TOT in nb
*
*   Programm calculation cross section for deutron,namely:
*   cross section of delta exitation(in approximation free nucleons),
*   cross section of quasielastic peak(Durand)
*   and they sum.
*            30-jun-89, A.Zatserkljanyi.
*
*      REAL MOTT
      real PROTON,NEUTR
*      real MATR,K,KMEV
      REAL M,MD
*      REAL MPION,KPION,KDELTA
*
      DATA PROTON,NEUTR/1.,1./
      DATA M/.93828/,MD/1.87563/
*      DATA MPION/.1396/,KPION/.150/

*      print*, 'D2TOT: E0,E,THETA:', E0,E,THETA
*      call WAIT('D2TOT begin')
      
      D2TOT = 0.

      PI=acos(-1.)

      WR=1.236 !| will be inputed in future,disarable
      KDELTA=.345
*
      TETA=THETA*PI/180.
      E1 = E0
      E2 = E

*     .. elastic peak
      Eel = E1/(1.+(2.*E1/MD)*(sin(TETA/2.))**2)
      if (E2.GE.Eel) then
*         print*, 'D2TOT: E2 > Eel: ', E2*1000., ' >', Eel*1000.
         RETURN
      endif
*
      CALL ROSTOT(E1,TETA,PROTON,NEUTR, EH,ROSP,ROSN,ROSPN)
*
      SIN2=SIN(TETA/2.)**2
      ETAH=1.+2.*E1*SIN2/M
      ETAD=1.+2.*E1*SIN2/MD
*      EH=E1/ETAH   ! Calculated by subroutine ROSTOT
      ED=E1/ETAD
*      EDELTA=(E1-KDELTA)/ETAH
*      EPION=(E1-KPION)/ETAH
*
*
*   BEGINING OF CALCULATION
*
      SDUR = 0.
      SDUR = DURAND(E1,E2,TETA)
      SDUR=SDUR*ROSPN ! cm2/GeV/sr
      SDUR=SDUR*(E1/E2) ! For omitting of recoil

*
*      CALL DELTATIT(E1,E2,TETA,WR,  K,PRES,PNONR,PDELTA)
**
**   We are assumed,that contributions in delta region
**   of proton and neutron are equal.
**
*      D2RES =2.*PRES
*      D2NONR=2.*PNONR
*      D2DELT=2.*PDELTA
*      QM2=4.*E1*E2*SIN2
*      Q=SQRT(QM2+(E1-E2)**2)
*      E2MEV=E2*1000.
*      ELOSS=(E1-E2)*1000. ! ELOSS in MeV !!!
*      KMEV=K*1000.

*     .. delta from Mo & Tsai
      d2delt = 2.*delta(E1,E2,THETA)
*
      SUM=SDUR+D2DELT      
*      print*, 'D2TOT: SUM =', SUM
      D2TOT=SUM
      END
*                                         @METAGS .DURAND
*                                         11-23-98 09:22pm
*--------------- DURAND ---------------
*
      function DURAND(E1,E2,THETA)
*     
*   This subroutine calculate factor SDUR for Durand cross section.
*   Durand cross section=SDUR*(sum Rosenbluths for proton and neutron)       
*   NB! In Rosenbluth recoiling factor E2/E1 must be omitted!
*   (In general,in Durand recoiling factor E2/E1 is absend!)
*   Inputing parameters(all REAL*4):
*     E1-energy of incident electron,GeV;
*     E2-energy of scattering electron,GeV;
*     THETA-scattering angle,rad
*   Outputing parametr(REAL*4) SDUR,1/GeV  
*              30-jun-89, A.Zatserkljanyi
*
      REAL M,M2,N2,MATR
*
      DATA PI/3.14159/
      DATA M/.9389/,EBIND/.002226/,A2/.0020905/,B2/.075116/,N2/.1533/ 
*
      DURAND=0.
*
      M2=M**2
      QM2=4.*E1*E2*SIN(THETA/2.)**2
      PCM2=M*(E1-E2-EBIND)-QM2/4.
      IF(PCM2.LE.0.) GOTO 10000  ! Threshold of electrodesintegration
      ECM2=PCM2+M2
      QCM2=QM2+(QM2/4.-PCM2-A2)**2/ECM2
*
      X=(A2+PCM2+QM2/4.)/SQRT(PCM2*QCM2)
      Y=(B2+PCM2+QM2/4.)/SQRT(PCM2*QCM2)
      X2=X**2
      Y2=Y**2
*
      MATR=N2/PCM2/QCM2*
     *( 1./(X2-1.)+1./(Y2-1.)-ALOG((X+1.)/(X-1.)*(Y-1.)/(Y+1.))/(Y-X) ) 
*
      SDUR=M2/PI*SQRT(PCM2/ECM2)*MATR
      DURAND = SDUR
*
10000 RETURN
      END

*                                         @METAGS .ROSTOT
*                                         11-23-98 09:22pm
*--------------- ROSTOT ---------------
*
      SUBROUTINE ROSTOT(E1,TETA,PROTON,NEUTR, EH,ROSP,ROSN,ROSPN)
*
*   This subroutine calculate Rosenbluth cross sections
*   for proton,neutron and tney sum in nb/sr.
*   Inputing parameters(all REAL*4):
*       E1-energy of incidet electron,GeV;
*       TETA-scatering angle,rad;
*    if PROTON=1.,calculate Rosenbluth for proton;
*    if NEUTR=1.,calculate rosenbluth for neutron;
*   Outputing parametres(all REAL*4):
*     EH- Position of elastic peak for proton
*       ROSP-Rosenbluth for proton(if at least PROTON=1.);
*       ROSN-Rosenbluth for neutron(if at least NEUTR=1.);
*       ROSPN-sum Rosenbluths for p and n (if PROTON=1.and NEUTR=1.);    
*               30-jun-89, A.Zatserkljanyi
*
      REAL MOTT,PROTON,NEUTR,M,MAGP2,MAGN2
*
*      DATA M/.93828/,MAGP2/7.800/,MAGN2/3.658/,CONST/5.1818E-33/
      DATA M/.93828/,MAGP2/7.800/,MAGN2/3.658/,CONST/5.1818/
*
      ROSP=0.
      ROSN=0.
      ROSPN=0.
*
      SIN2=SIN(TETA/2.)**2
      COS2=COS(TETA/2.)**2
      TAN2=SIN2/COS2
*
      MOTT=CONST/(E1**2)*COS2/SIN2**2
*
      EH=E1/(1.+2.*E1*SIN2/M)
      QM2=4.*E1*EH*SIN2
      TAU=QM2/4./M**2
      G2=(1./(1+QM2/0.71)**2)**2
*
      MAGP2=MAGP2*PROTON
      MAGN2=MAGN2*NEUTR
      GPE2=G2
      GPM2=MAGP2*G2
      GNE2=0.
      GNM2=MAGN2*G2
*
      SIGNS=MOTT*EH/E1
      ROSP=SIGNS*((GPE2+TAU*GPM2)/(1.+TAU)+2.*TAU*GPM2*TAN2)
      ROSN=SIGNS*((GNE2+TAU*GNM2)/(1.+TAU)+2.*TAU*GNM2*TAN2)
      ROSPN=ROSP+ROSN
*
      RETURN
      END

*                                         @METAGS .delta
*                                         11-25-98 06:06pm
*--------------- delta ---------------
*
      function delta(Es,Ep,theta)
      real Mp
      parameter (Mp=.93827)
      parameter (GeV2nb=.38938e6)
      real Fdelta,Gdelta
      external Fdelta,Gdelta
      delta = GeV2nb*dsdp(Es,Ep,theta,Mp,Fdelta,Gdelta)
      END

*                                         @METAGS .dsdp
*                                         11-25-98 03:58pm
*--------------- dsdp ---------------
*
      function dsdp(Es,Ep,theta,M,F,G)
      real M,Mf2
      parameter (alp=1./137.036)
      real F,G
      dsdp = 0.
      if (Ep.GE.Es) RETURN
      pi = acos(-1.)
      thr2 = theta*pi/360.
      sin2 = (sin(thr2))**2
      cos2 = 1.-sin2
      tan2 = sin2/cos2
      q2 = -4.*Es*Ep*sin2
      Mf2 = M**2 + 2.*M*(Es-Ep) + q2
*      bracket = F(q2,Mf2) + (2./M**2)*tan2*G(q2,Mf2)
      a = F(q2,Mf2)
      b = (2./M**2)*tan2*G(q2,Mf2)
      bracket = a + b
*      print*, '--- bracket, a,b:', bracket, a,b
      dsdp = 2.*(alp*Ep/q2)**2*M*cos2*bracket
      END

*                                         @METAGS .funEs
*                                         12-01-98 01:06pm
*--------------- funEs ---------------
*
      function funEs(Es1)
*     .. The integrand over Es1. Es,Ep,theta supplies commom /MoTsai/
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
*      use MoTsai,TradMT
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2=.511e-3*.511e-3)
      real Is

      funEs = 0.
*     .. Es,Ep,theta from common /MoTsai/
      Es = EsMT
      Ep = EpMT
      theta = thetaMT

      b  = bMT
      bw = bwMT

      xs = Es1/Es
*     .. get already calculated sp,fs from common /MoTsai/
      sp = spMT
      fs = fsMT

      ts = (alp/pi)*(.5*(1.+xs**2)*log(2.*sp/me2) - xs)
*     .. probability to emit photon in the initial state
      Is = (ts + (bw*tiwrad+b*Trad/2.)*(xs + (3./4.)*(1-xs)**2))*
     &     log(1./xs)**fs / (Es-Es1)
      funEs = Is*D2TOT(Es1,Ep,theta)
      END
*                                         @METAGS .funEp
*                                         12-01-98 01:06pm
*--------------- funEp ---------------
*
      function funEp(Ep1)
*     .. The integrand over Ep1. Es,Ep,theta supplies commom /MoTsai/
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
*      use MoTsai,TradMT
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2=.511e-3*.511e-3)
      real Ip

      funEp = 0.
*     .. Es,Ep,theta from common /MoTsai/
      Es = EsMT
      Ep = EpMT
      theta = thetaMT

      b  = bMT
      bw = bwMT

      xp = Ep/Ep1
*     .. get already calculated sp,fp from common /MoTsai/
      sp = spMT
      fp = fpMT

      tp = (alp/pi)*(.5*(1.+xp**2)*log(2.*sp/me2) - xp)
*     .. probability to emit photon in the final state
      Ip = (tp + (bw*tfwrad+b*Trad/2.)*(xp + (3./4.)*(1-xp)**2))*
     &     log(1./xp)**fp / (Ep1-Ep)
      funEp = Ip*D2TOT(Es,Ep1,theta)
      END

*                                         @METAGS .RadMT
*                                         12-02-98 10:51am
*--------------- RadMT ---------------
*
      function RadMT(Es,Ep,del)
*     .. scattering angle thetaMT supplies common /MoTsai/
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
*      use MoTsai,TradMT
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2 =.511e-3*.511e-3)
      real lnEsdel,lnEpdel,lnqme2, IdEs,IdEp
      real Kmin
      external funEs,funEp

      RadMT = 0.
*     .. store given values Es,Ep in common /MoTsai/
      EsMT = Es
      EpMT = Ep

      b  = bMT
      bw = bwMT

      costheta = cos(thetaMT*pi/180.)
      sp = Es*Ep*(1.-costheta)
      lnqme2 = log(2.*sp/me2)

*     .. equivalent radiator
      tr = (alp/pi)*(lnqme2-1.) / b
*     .. fs,fp terms
      fs = b*tr + bw*tiwrad + b*Trad/2.
      fp = b*tr + bw*tfwrad + b*Trad/2.
*     .. store in common /MoTsai/ for following usage by funEs,funEp
      spMT = sp
      trMT = tr
      fsMT = fs
      fpMT = fp

*     .. minimal energy of equivalent gamma-quantum
      Kmin = ((Mtar+mthr)**2 - Mtar**2)/(2.*Mtar)

      relerr = 1.e-6
      abserr = 0.
*     .. integral over Es
      termEs = 0.
      Esmin = (Ep + Kmin)/(1.-(Ep/Mtar)*(1-costheta))
      Es1 = Esmin
      Es2 = Es-del
      if (Es1.GE.Es2) goto 1000
      call quanc8(funEs,Es1,Es2, abserr,relerr, IdEs,errest,nofun,flag)
*      print*, 'IdEs =', IdEs
      termEs = (del/Ep)**(fp/2.)*IdEs
*      print*, 'termEs =', termEs
 1000 continue
*     .. integral over Ep
      Epmax = (Es - Kmin)/(1.+(Es/Mtar)*(1-costheta))
      Ep1 = Ep+del
      Ep2 = Epmax
      if (Ep1.GE.Ep2) goto 2000
      call quanc8(funEp,Ep1,Ep2, abserr,relerr, IdEp,errest,nofun,flag)
*      print*, 'IdEp =', IdEp
      termEp = (del/Es)**(fs/2.)*IdEp
*      print*, 'termEp =', termEp

 2000 continue
      lnEsdel = log(Es/del)
      lnEpdel = log(Ep/del)
      deltat = -( (bw*tiwrad+b*Trad/2.)*lnEsdel +
     +            (bw*tfwrad+b*Trad/2.)*lnEpdel )

      deltar = -(alp/pi) * ( 28./9. - (13./6)*lnqme2 +
     +                       (lnEsdel+lnEpdel)*(lnqme2-1.) -
     -                       Spence(-(Es-Ep)/Ep) - Spence( (Es-Ep)/Es) )

      sum = D2TOT(Es,Ep,thetaMT)*exp(deltat+deltar) + termEs + termEp
*      print*, 'sum =', sum
      RadMT = sum
      END

*                                         @METAGS .Fdelta
*                                         11-25-98 03:47pm
*--------------- Fdelta ---------------
*
      real function Fdelta(q2,Mf2)
*     .. (III.5)-(III.9) for (B1)
*     .. Mo & Tsai metric: q2 < 0. Main dimension - GeV.
      IMPLICIT NONE
      real q2,Mf2
*     .. x denotes variables in the rest system of final hadron state (Nx)
      real px,px2, Gamma, QQ2,QQx2, q, C32, Eix, pi, Breit
      real G1
      real G2
*      real pxR
      real    Mf,Mp, M33, mpi
      real       Mp2,M332,mpi2
      parameter (Mp =.93827,        M33 =1.236,       mpi =.140)
      parameter (Mp2=.93827*.93827, M332=1.236*1.236, mpi2=.140*.140)

      Fdelta = 0.
      if (Mf2.LE.(Mp+mpi)**2) RETURN
      Mf = sqrt(Mf2)
*     .. px is momentum of decaing pion in the rest system of Nx
      px2 = (Mf2-Mp2+mpi2)**2/(4.*Mf2) - mpi2
      px = sqrt(px2)
      Gamma = .1293*(.85*px/mpi)**3 / (1.+(.85*px/mpi)**2)
**--
**        another form of Gamma(Mf2) from
**        A.J.Dufner and Y.S.Tsai, Phys.Rev.,168,1801(1968)
**        The metric of Dufner & Tsai the same:
**        q2 < 0, e**2/4*pi = alpha, r0*me = alpha.
**        Notice: both forms are agree very well with examples from
**        Dufner & Tsai and contradict with Mo & Tsai examples
**        of "Non-Radiative" curve at low Ep slope of 1236.
**        May be they included the higher resonances 1525 and 1688 MeV?
**
**        Here pxR is the value of px at the resonance;
**        i.e., they let Mf=M33=1.236 GeV
*      pxR = sqrt( (M332-Mp2+mpi2)**2/(4.*M332) - mpi2 )
*      Gamma = .12*(px/pxR)**3
**--
*     .. QQ is 3-momentum transferred in lab, QQ2=(Es-Ep)**2-q2
      QQ2  = (Mf2-q2-Mp2)**2/(2.*Mp)**2 - q2
*     .. QQx is 3-momentum transferred in the rest system of Nx
      QQx2 = QQ2*Mp2/Mf2
      q = sqrt(-q2)
      C32 = (1./Mp2) * 2.05**2 * exp(-6.3*q) * (1.+9.0*q)
*     .. Eix is energy of initial proton in the rest system of Nx
      Eix = (Mf2+Mp2-q2)/(2.*Mf)
      pi = acos(-1.)
      Breit = (1./pi) * Gamma*M33 / ((Mf2-M332)**2 + (Gamma*M33)**2)
      G1 = Breit*Mf*QQx2*2.*C32*(Eix+Mp)/(3.*Mp)
      G2 = (-q2/QQ2)*G1
      Fdelta = (2./Mp)*G2
      END
*                                         @METAGS .Gdelta
*                                         11-25-98 03:47pm
*--------------- Gdelta ---------------
*
      real function Gdelta(q2,Mf2)
*     .. (III.5)-(III.9) for (B1)
*     .. Mo & Tsai metric: q2 < 0. Main dimension - GeV.
      IMPLICIT NONE
      real q2,Mf2
*     .. x denotes variables in the rest system of final hadron state (Nx)
      real px,px2, Gamma, QQ2,QQx2, q, C32, Eix, pi, Breit
      real G1
*      real G2
*      real pxR
      real    Mf,Mp, M33, mpi
      real       Mp2,M332,mpi2
      parameter (Mp =.93827,        M33 =1.236,       mpi =.140)
      parameter (Mp2=.93827*.93827, M332=1.236*1.236, mpi2=.140*.140)

      Gdelta = 0.
      if (Mf2.LE.(Mp+mpi)**2) RETURN
      Mf = sqrt(Mf2)
*     .. px is momentum of decaing pion in the rest system of Nx
      px2 = (Mf2-Mp2+mpi2)**2/(4.*Mf2) - mpi2
      px = sqrt(px2)
      Gamma = .1293*(.85*px/mpi)**3 / (1.+(.85*px/mpi)**2)
**--
**        another form of Gamma(Mf2) from
**        A.J.Dufner and Y.S.Tsai, Phys.Rev.,168,1801(1968)
**        The metric of Dufner & Tsai the same:
**        q2 < 0, e**2/4*pi = alpha, r0*me = alpha.
**        Notice: both forms are agree very well with examples from
**        Dufner & Tsai and contradict with Mo & Tsai examples
**        of "Non-Radiative" curve at low Ep slope of 1236.
**        May be they included the higher resonances 1525 and 1688 MeV?
**
**        Here pxR is the value of px at the resonance;
**        i.e., they let Mf=M33=1.236 GeV
*      pxR = sqrt( (M332-Mp2+mpi2)**2/(4.*M332) - mpi2 )
*      Gamma = .12*(px/pxR)**3
**--
*     .. QQ is 3-momentum transferred in lab, QQ2=(Es-Ep)**2-q2
      QQ2  = (Mf2-q2-Mp2)**2/(2.*Mp)**2 - q2
*     .. QQx is 3-momentum transferred in the rest system of Nx
      QQx2 = QQ2*Mp2/Mf2
      q = sqrt(-q2)
      C32 = (1./Mp2) * 2.05**2 * exp(-6.3*q) * (1.+9.0*q)
*     .. Eix is energy of initial proton in the rest system of Nx
      Eix = (Mf2+Mp2-q2)/(2.*Mf)
      pi = acos(-1.)
      Breit = (1./pi) * Gamma*M33 / ((Mf2-M332)**2 + (Gamma*M33)**2)
      G1 = Breit*Mf*QQx2*2.*C32*(Eix+Mp)/(3.*Mp)
*      G2 = (-q2/QQ2)*G1
      Gdelta = 2.*Mp*G1
      END

*                                         @METAGS .ReadRad
*                                         12-24-98 04:13pm
*--------------- ReadRad ---------------
*
      SUBROUTINE ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      character*(*) file
      integer np
      real Epdata(np),Mfdata(np),srdata(np),sndata(np)
      lun = LUNFREE(1)
      open (lun,FILE=file,STATUS='OLD')
      do i=1,np
         read (lun,*,END=1000) Epdata(i),Mfdata(i),sndata(i),srdata(i)
*         if (i.GE.192) then
*            print*, i, Epdata(i),sndata(i),srdata(i)
*            call WAIT('i, Epdata(i),,sndata(i),srdata(i)')
*         endif
      enddo
 1000 close(lun)
      END

*                                         @METAGS .WrURAD
*                                         12-24-98 04:13pm
*--------------- WrURAD ---------------
*
      SUBROUTINE WrURAD()
      IMPLICIT NONE
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      character*32 file
      integer NS,NP
      integer lun, LUNFREE
      integer len, LENOCC
      do NS=1,NStot
         file = ' '
         write(file,10) NS
*GNU!!!   10    format(1X, 'urad', I1, '.dat')
   10    format('urad', I1, '.dat')
         lun = LUNFREE(1)
         len = LENOCC(file)
         print*, 'Writing file ', file(1:len), ' for NS =', NS
         open (lun, FILE=file, STATUS='OLD', ERR=1000)
         close(lun, STATUS='DELETE')
 1000    continue
         open (lun, FILE=file, STATUS='UNKNOWN', ERR=10000)
         do NP=1,nptot(NS)
            write(lun,*, ERR=20000) Epdat(NP,NS),sudat(NP,NS)
         enddo
         close(lun)
      enddo
      RETURN
10000 print*, 'Open error file ', file(1:len)
      RETURN
20000 close(lun)
      print*, 'Write error file', file(1:len), '   NP =', NP
      END

*                                         @METAGS .FillMT
*                                         12-24-98 04:13pm
*--------------- FillMT ---------------
*
      SUBROUTINE FillMT(np,Es,Epdata,srdata,sndata)
      IMPLICIT NONE
      integer np
      real Es,Epdata(np),srdata(np),sndata(np)

      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      integer NSmax,NPmax
      COMMON /URADNSNP/ NSmax,NPmax

      real pi
      parameter (pi=3.141593)
      real Kmin,Mf2
      real Epmax,E2max,Ep,theta,costheta
      integer i

      if (NStot.EQ.NSmax) then
         print*, 'COMMON /URADMT/ is full. No space to arrange data'
         STOP
      endif
      if (np.GT.NPmax) then
         print*, 'Dimension of common /URADMT/ is overflowed'
         print*, 'The maximum number of data point is', NPmax
         STOP
      endif

      NStot = NStot+1
      Esdat(NStot) = Es
      npcur(NStot) = 0

      theta = thetaMT
      costheta = cos(theta*pi/180.)

*     .. minimum energy of equivalent gamma-quantum
      Kmin = ((Mtar+mthr)**2 - Mtar**2)/(2.*Mtar)
*     .. threshold energy
      Epmax = (Es - Kmin)/(1.+(Es/Mtar)*(1-costheta))

      E2max = Epdata(1)
      do i=1,np
         Ep = Epdata(i)
         if ((Ep.EQ.0.) .OR. (Ep.GT.Epmax)) goto 1000
         nptot(NStot) = i
         E2max=Ep
         Epdat(i,NStot) = Ep
         Mf2 = Mtar**2+2.*Mtar*(Es-Ep)-2.*Es*Ep*(1.-costheta)
         Mfdat(i,NStot) = sqrt(Mf2)
         srdat(i,NStot) = srdata(i)
         sudat(i,NStot) = 0.
         sndat(i,NStot) = sndata(i)
      enddo
 1000 continue
      Epma(NStot)  = MIN(E2max,Epmax)

      do i=nptot(NStot)+1,NPmax
         Epdat(i,NStot) = 0.
         Mfdat(i,NStot) = 0.
         srdat(i,NStot) = 0.
         sudat(i,NStot) = 0.
         sndat(i,NStot) = 0.
      enddo
      END

*                                         @METAGS .Getsudat
*                                         12-24-98 04:13pm
*--------------- Getsudat ---------------
*
      SUBROUTINE Getsudat(NS,NP,Epdata,sudata)
      integer NS,NP
      real Epdata(NP),sudata(NP)
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      real Mfdat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      COMMON /URADNSNP/ NSmax,NPmax
*      use URADMT
      if (NS.GT.NStot) then
         print*, 'The maximum number of spectrum is', NStot
         STOP
      endif
      if (NP.LT.nptot(NS)) then
         print*, 'The maximum number of restored points is', nptot(NS)
         STOP
      endif
      do i=1,nptot(NS)
         sudata(i) = sudat(i,NS)
         Epdata(i) = Epdat(i,NS)
      enddo
      END

*                                         @METAGS .ClearMT
*                                         12-24-98 04:21pm
*--------------- ClearMT ---------------
*
      SUBROUTINE ClearMT()
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      real Mfdat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      COMMON /URADNSNP/ NSmax,NPmax
      NStot = 0
      END

*                                         @METAGS .iMoTsai
*                                         12-28-98 09:43pm
*--------------- iMoTsai ---------------
*
      SUBROUTINE iMoTsai()
      IMPLICIT NONE
      real EsMT,EpMT,thetaMT,Mtar,mthr,spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr,spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      real    Mf,MD
      parameter (MD=1.87563)
      real me2
      parameter (me2 =.511e-3*.511e-3)
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real costheta,xi,recoil,Z,b,bw
*     .. fill by some values
      EsMT = .411
      thetaMT = 143.5
      Mtar = MD
      Mthr = .00222
      Mf = Mtar+Mthr
      costheta = cos(thetaMT*pi/180.)
      recoil = 1. + (EsMT/MD)*(1-costheta)
*     .. scattered energy in 1236
      EpMT = (EsMT - (Mf**2-MD**2)/(2.*MD)) / recoil
      spMT = EsMT*EpMT*(1.-costheta)

*     .. target - D2
      Ztar = 1.
*     .. window - Al
      Zwin = 13.

      Z = Ztar
*     .. formula (A5)
      xi = log(1440.*Z**(-2./3.)) / log(183.*Z**(-1./3.))
*     .. formula (A4)
      b  = (4./3.)*(1. + (1./9.)*((Z+1.)/(Z+xi))/log(183.*Z**(-1./3.)))
      bMT = b
      Z = Zwin
      xi = log(1440.*Z**(-2./3.)) / log(183.*Z**(-1./3.))
      bw = (4./3.)*(1. + (1./9.)*((Z+1.)/(Z+xi))/log(183.*Z**(-1./3.)))
      bwMT = bw

      Trad   = .00663
      tiwrad = .0051
      tfwrad = .00067
*     .. equivalent radiator
      trMT = (alp/pi)*(log(2.*spMT/me2)-1.) / b
*     .. fs,fp terms
      fsMT = b*trMT + bw*tiwrad + b*Trad/2.
      fpMT = b*trMT + bw*tfwrad + b*Trad/2.
      END

*                                         @METAGS .iURADMT
*                                         12-29-98 12:44pm
*--------------- iURADMT ---------------
*
      SUBROUTINE iURADMT()
      IMPLICIT NONE
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      integer NSmax,NPmax
      COMMON /URADNSNP/ NSmax,NPmax
      NSmax = NSDIM
      NPmax = NPDIM
      NStot = 0
      END

*                                         @METAGS .quanc8.MAIN
*                                         11-19-98 04:45pm
*--------------- quanc8.MAIN ---------------
*
**
**     Иллюстpитующая пpогpамма для quanc8
**     Demo program for quanc8
**
*      real function fun(x)
*      real x
*      if (x .eq. 0.) fun = 1.
*      if (x .ne. 0.) fun = sin(x)/x
*      return
*      end
**
**     .. main ..
**
*      external fun
*      real a, b, abserr, relerr, result, errest, flag
*      integer nofun
*
*      a = 0.
*      b = 2.
*      relerr = 1.e-10
*      abserr = 0.
*      call quanc8(fun,a,b, abserr,relerr, result,errest, nofun,flag)
*      write(*,1) result, errest
*      if (flag .ne. 0.) write(*,2) flag
*    1 format (8h result=, f15.10, 10h errest=  , e10.2)
*    2 format (44h Warning.. result may be unreliable.Flag=    , f6.2)
*      stop
*      end
*
*
*
*                                         @METAGS .quanc8
*                                         11-19-98 04:45pm
*--------------- quanc8 ---------------
*
      subroutine quanc8
     &           (fun,a,b, abserr,relerr, result,errest, nofun,flag)
*
      real fun,a,b,abserr,relerr,result,errest,flag
      integer nofun
*
*     Estimate the integral for fun(x) from a to b to given accuracy.
*     Automatic adaptive routine bases on the Newton-Cotes formula of
*     8th order
*
*     Input parameters ..
*
*     fun         function name
*     a           low limit
*     b           upper limit (may be less then a)
*     abserr      absolute error boundary (must be positive)
*     relerr      relative error boundary (must be positive)
*
*     Output parameters ..
*
*     result      approximation to the integral that satisfy (we hope)
*                 less strict error boundary from the two ones
*     errest      astimation of the actual error
*     nofun       the number of function calls have been used for the
*                 calculation of result
*     flag        the indicator of the reliability. If flag equal zero then
*                 result probably satisfies the given error boundary.
*                 If flag = XXX.YYY then XXX = the number of the intervals
*                 without convergence and 0.YYY = part of the main interval
*                 have been remained to process when the routine reaches the
*                 upper value for nofun.
*
      real w0, w1, w2, w3, w4, area, x0, f0, stone, step, cor11, temp
      real qprev, qnow, qdiff, qleft, esterr, tolerr
      real qright(31), f(16), x(16), fsave(8,30), xsave(8,30)
      integer levmin, levmax, levout, nomax, nofin, lev, nim, i, j
*
*     *** Step 1 *** Assignment the initial values to the variables that
*     not depend from the interval. Constants generation.
*
      levmin = 1
      levmax = 30
      levout = 6
      nomax = 5000
      nofin = nomax - 8*(levmax - levout + 2**(levout+1))
*
*     If nofun reaches the nofun then alarm
*
      w0 = 3956./14175.
      w1 = 23552./14175.
      w2 = -3712./14175.
      w3 = 41984./14175.
      w4 = -18160./14175.
*
*     Assign the zero values to the variable sums.
*
      flag = 0.
      result = 0.
      cor11 = 0.
      errest = 0.
      area = 0.
      nofun = 0
      if (a .eq. b) return
*
*     *** Step 2 *** Assignment the initial values to the variables that
*     depend from the interval in accordance to the first interval
*
      lev = 0
      nim = 1
      x0 = a
      x(16) = b
      qprev = 0.
      f0 = fun(x0)
      stone = (b - a) / 16.
      x(8) = (x0 + x(16)) / 2.
      x(4) = (x0 + x(8)) /2.
      x(12) = (x(8) + x(16)) /2.
      x(2) = (x0 + x(4)) /2.
      x(6) = (x(4) + x(8)) /2.
      x(10) = (x(8) + x(12)) /2.
      x(14) = (x(12) + x(16)) /2.
      do 25 j = 2, 16, 2
         f(j) = fun(x(j))
   25 continue
      nofun = 9
*
*     *** Step 3 *** Main calculations.
*     Have been require qprev, x0,x(2),x(4),...,x(16), f0,f(2),f(4),...,f(16)
*     Have been calculate x(1),x(3),...,x(15), f(1),f(3),...,f(15),
*     qleft, qright, qnow, qdiff, area
*
   30 x(1) = (x0 + x(2)) / 2.
      f(1) = fun(x(1))
      do 35 j = 3, 15, 2
         x(j) = (x(j-1) + x(j+1)) / 2.
         f(j) = fun(x(j))
   35 continue
      nofun = nofun+8
      step = (x(16) - x0) / 16.
      qleft = ( w0*(f0+f(8)) + w1*(f(1)+f(7)) + w2*(f(2)+f(6)) +
     +          w3*(f(3)+f(5)) + w4*f(4) ) * step
      qright(lev+1) = ( w0*(f(8)+f(16)) + w1*(f(9)+f(15)) +
     +                  w2*(f(10)+f(14)) + w3*(f(11)+f(13)) +
     +                  w4*f(12) ) * step
      qnow = qleft + qright(lev+1)
      qdiff = qnow - qprev
      area = area + qdiff
*
*     *** Step 4 *** Test of convergence for the interval
*
      esterr = abs(qdiff) / 1023.
      tolerr = amax1(abserr, relerr*abs(area)) * (step/stone)
      if (lev .lt. levmin) goto 50
      if (lev .ge. levmax) goto 62
      if (nofun .gt. nofin) goto 60
      if (esterr .le. tolerr) goto 70
*
*     *** Step 5 *** Convergence is absend.
*     Set the next interval.
*
   50 nim = 2*nim
      lev = lev+1
*
*     Store for the future using the elements fall into the right
*     half of the interval
*
      do 52 i = 1, 8
         fsave(i,lev) = f(i+8)
         xsave(i,lev) = x(i+8)
   52 continue
*
*     Collect for the immediately using the elements fall into the left
*     half of the interval
*
      qprev = qleft
      do 55 i = 1, 8
         j = -i
         f(2*j + 18) = f(j+9)
         x(2*j + 18) = x(j+9)
   55 continue
      goto 30
*
*     *** Step 6 *** "Fire part".
*     The number of the function calls is close to overflow the
*     established limit
*
   60 nofin = 2.*nofin
      levmax = levout
      flag = flag + (b-x0) / (b-a)
      goto 70
*
*     The current limit value of the halving depth is equal to levmax
*
   62 flag = flag+1.
*
*     *** Step 7 *** Convergence for the interval have been occured.
*     Add the next terms to the variable sums.
*
   70 result = result + qnow
      errest = errest + esterr
      cor11 = cor11 + qdiff/1023.
*
*     Set the next interval
*
   72 if (nim .eq. 2*(nim/2)) goto 75
      nim = nim/2
      lev = lev-1
      goto 72
   75 nim = nim+1
      if (lev .le. 0) goto 80
*
*     Collect the elemets have been needed for the next interval
*
      qprev = qright(lev)
      x0 = x(16)
      f0 = f(16)
      do 78 i = 1, 8
         f(2*i) = fsave(i,lev)
         x(2*i) = xsave(i,lev)
   78 continue
      goto 30
*
*     *** Step 8 *** Final operations and exit
*
   80 result = result + cor11
*
*     Supply the value of errest at least the round level
*
      if (errest .eq. 0.) return
   82 temp = abs(result) + errest
      if (temp .ne. abs(result)) return
      errest = 2.*errest
      goto 82
      end

*                                         @METAGS .WAIT.COMIS
*                                         11-15-96 09:45pm
*--------------- WAIT ---------------
*
      SUBROUTINE WAIT(mess)
      character mess*(*), ch*1
      length = LENOCC(mess)
      if (length.GT.0) print 1, mess(1:length)
      print*, '<CR>=Continue, Q=Quit'
      read 1, ch
      if ((ch.EQ.'q') .OR. (ch.EQ.'Q')) STOP
      RETURN
    1 FORMAT(A)
      END
*                                         @METAGS .DWAIT.COMIS
*                                         11-15-96 09:45pm
*--------------- DWAIT ---------------
*
      SUBROUTINE DWAIT(mess)
      LOGICAL DEBUG
      COMMON /DEBUG/ DEBUG
*      use DEBUG
      character mess*(*), ch*32
      character ON*2,OFF*2, SHOW*2, KEY*1
      data ON/'$1'/, OFF/'$0'/, SHOW/'$$'/, KEY/'$'/
      length = LENOCC(mess)
      if (INDEX(mess,KEY).GT.0) then
         if (mess.EQ.ON) then
            DEBUG = .TRUE.
            RETURN
         endif
         if (mess.EQ.OFF) then
            DEBUG = .FALSE.
            RETURN
         endif
         if (mess.EQ.SHOW) then
            print*, 'DEBUG: Debugging is ', DEBUG
            RETURN
         endif
      endif

      if (.NOT.DEBUG) RETURN

      if (length.GT.0) print 1, mess(1:length)
      print*, '<CR>=Continue, Q=Quit, 0=Cancel debugging'
      read 1, ch
      call CLTOU(ch)
      if (INDEX(ch,'0').GT.0) DEBUG=.FALSE.
      if (INDEX(ch,'Q').GT.0) STOP
      RETURN
    1 FORMAT(A)
      END
*                                         @METAGS .MESS.COMIS
*                                         01-22-98 02:11pm
*--------------- MESS ---------------
*
      SUBROUTINE MESS(line)
      character line*(*)
      length = LENOCC(line)
      if (length.GT.0) print 1, line(1:length)
      RETURN
    1 FORMAT(A)
      END
*                                         @METAGS .CLOSE.COMIS
*                                         01-14-98 05:44pm
*--------------- CLOSE ---------------
*
      SUBROUTINE CLOSE(lun)
      close(lun)
      END

*                                         @METAGS .LUNFREE
*                                         01-20-98 02:25pm
*--------------- LUNFREE ---------------
*
      integer function LUNFREE(lunstart)
      logical used
      lun = lunstart-1
      if (lun.LT.0) lun=0
  100 lun = lun+1
      inquire (UNIT=lun, OPENED=used)
      if (used) goto 100
      LUNFREE = lun
      END
*                                         @METAGS .LUNFILE
*                                         01-20-98 02:40pm
*--------------- LUNFILE ---------------
*
      integer function LUNFILE(file)
      character*32 file,lofile
      lofile = file
      call CUTOL(lofile)
      inquire (NUMBER=lun, FILE=lofile)
      LUNFILE = lun
      END
*                                         @METAGS .FILELUN
*                                         01-20-98 02:40pm
*--------------- FILELUN ---------------
*
      SUBROUTINE FILELUN(lun,file)
      character*(*) file, lofile*32
      inquire (UNIT=lun, NAME=lofile)
      call CUTOL(lofile)
      lenght = LEN(file)
      file(1:length) = lofile
      END
*                                         @METAGS .FCLEAR
*                                         01-20-98 01:27pm
*--------------- FCLEAR ---------------
*
      SUBROUTINE FCLEAR(file)
      character*(*) file
      lun = LUNFREE(40)
      open (lun, FILE=file, STATUS='UNKNOWN', ERR=10000)
      endfile lun
      close(lun)
      RETURN
10000 len = LENOCC(file)
      print*, 'ERROR FCLEAR: File ', file(1:len), ' did not found'
      END
*                                         @METAGS .FPARSE
*                                         01-20-98 02:53pm
*--------------- FPARSE ---------------
*
      SUBROUTINE FPARSE(fstr,fname,file,ext)
      character*(*) fstr,fname,file,ext
      character f*32, point
      data point /'.'/
      f = fstr
      lenf = LENOCC(f)
      lene = LENOCC(ext)
      ipoint = INDEX(f,'.')
      if (ipoint.EQ.0) then
         fname = f
         file  = f(1:lenf)//point//ext(1:lene)
      else
         fname = f(1:ipoint-1)
         file  = f
      endif
      call CUTOL(file)
      END

*                                         @METAGS .fEs0
*                                         12-01-98 01:06pm
*--------------- fEs0 ---------------
*
      real function fEs0(Es1)
*     .. The integrand over Es1. Es,Ep,theta supplies commom /MoTsai/
      IMPLICIT NONE
      real Es1
      real EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      COMMON /MoTsai/ EsMT,EpMT,thetaMT,Mtar,mthr, spMT,trMT,fsMT,fpMT
      real Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      COMMON /TradMT/ Trad,tiwrad,tfwrad, Ztar,Zwin, bMT,bwMT
      integer NSDIM,NPDIM
      parameter (NSDIM=9)
      parameter (NPDIM=400)
      integer NStot,NScur,nptot,npcur
      real Esdat,Epma,Epdat,Mfdat,srdat,sudat,sndat
      COMMON/URADMT/
     & NStot,NScur,nptot(NSDIM),npcur(NSDIM),Esdat(NSDIM),Epma(NSDIM),
     & Epdat(NPDIM,NSDIM),Mfdat(NPDIM,NSDIM),
     & srdat(NPDIM,NSDIM),sudat(NPDIM,NSDIM),sndat(NPDIM,NSDIM)
      LOGICAL DEBUG
      COMMON /DEBUG/ DEBUG
      real vEs(NSDIM),vsec(NSDIM)
*      use MoTsai,TradMT,URADMT
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real me,me2
      parameter (me=.511e-3, me2=.511e-3*.511e-3)
      real eps
      parameter (eps=0.00001)
      real sMf,sinter
      real Is
      integer NS
      real Es,Ep,theta,b,bw,xs,sp,fs,ts,sec
      real Mf,Mf2
      integer ipoint

      real Ep1,recoil,dE,dS
      integer np,nlow,nhigh

      fEs0 = 0.
*      print*, 'fEs0: Es1=', Es1

***      if (Ep1.GE.Epma(NScur)) RETURN

*     .. Es,Ep,theta from common /MoTsai/
      Es = EsMT
      Ep = EpMT
      theta = thetaMT
      b  = bMT
      bw = bwMT

      xs = Es1/Es
*     .. get already calculated sp,fp from common /MoTsai/
      sp = spMT
      fs = fsMT

      ts = (alp/pi)*(.5*(1.+xs**2)*log(2.*sp/me2) - xs)
*     .. probability to emit photon in the final state
      Is = (ts + (bw*tfwrad+b*Trad/2.)*(xs + (3./4.)*(1-xs)**2))*
     &     log(1./xs)**fs / (Es-Es1)

*     .. interpolate along the Mf already unfolded cross section
      Mf2=Mtar**2 + 2.*Mtar*(Es1-Ep) - 2.*Es1*Ep*(1.-cos(theta*pi/180.))
      if (Mf2.LT.0.) then
         print*, 'ERROR fEs0: Mf2 < 0.'
         STOP
      endif
      Mf = sqrt(Mf2)

      ipoint = 0
      do NS=1,NStot
         if (Mf.GT.Mfdat(1,NS)) goto 100
         sec = sMf(Mf,NS)
*         if (sec .LT. 0.) then
*            print*, 'fEs0: Out of range'
*            STOP
*         endif
         if (sec .GT. 0.) then
            ipoint = ipoint+1
            vEs(ipoint)  = Esdat(NS)
            vsec(ipoint) = sec
         endif
  100 enddo

      if (ipoint.GT.0) then
         sec = sinter(Es1,ipoint,vEs,vsec)
      else
*         print*, 'ERROR fEs0: ipoint=0.'
         if(DEBUG)print*,'ipoint=0, Es1,Mf', Es1,Mf
*         call DWAIT('ERROR fEs0: ipoint=0.')
*         STOP
*         RETURN

*        .. get cross section at corresponded Ep value
         recoil = 1. + (Es1/Mtar)*(1.-cos(theta*pi/180.))
         Ep1 = (Es1 - (Mf**2-Mtar**2)/(2.*Mtar))/recoil

*        .. interpolate already unfolded cross section
         NS = NScur
         if (Epdat(npcur(NS),NS).LE.Ep1) then
*           .. search forward
            do np=npcur(NS),nptot(NS)
               nhigh = np
               if (Epdat(nhigh,NS).GT.Ep1) goto 1000
            enddo
         else
*           .. extrapolate using points npcur(NS), npcur(NS)+1
            nhigh = npcur(NS)+1
         endif

 1000    nlow = nhigh-1
*        .. interpolation(extrapolation)
         dE = Epdat(nhigh,NS)-Epdat(nlow,NS)
         if (dE.GT.eps) then
            dS = sudat(nhigh,NS)-sudat(nlow,NS)
            sec = sudat(nlow,NS) + dS*(Ep1-Epdat(nlow,NS))/dE
         else
            sec = (sudat(nhigh,NS) + sudat(nlow,NS))/2.
         endif 
         if (sec .LT. 0.) sec=0.
      endif

      fEs0 = Is*sec
      END

