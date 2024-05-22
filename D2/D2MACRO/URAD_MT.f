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
      Es = 7.0
      file = 'Es7.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)
      Es = 10.0
      file = 'Es10.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)
      Es = 13.5
      file = 'Es13.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)
      Es = 16.0
      file = 'Es16.dat'
      call ReadRad(file,np,Epdata,Mfdata,srdata,sndata)
      call FillMT(np,Es,Epdata,srdata,sndata)

      del = .010
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
      real Kmin,Mfmin,Mfmax,Mf2
      real lnEsdel,lnEpdel,lnqme2, IdEs,IdEp,IdEsEp,recoil
      real fEs,fEp,Spence
      external fEs,fEp
      real Es0,Ep0,b,bw,Es,Ep,sp,tr,fs,fp,Epmax,Esmin,termEs,termEp
      real deltat,deltar,relerr,abserr,errest,Epl,Eph,Esl,Esh,Kmax,Es1
      real costheta,flag,sum
      integer i,NS,NSminMf,np,nofun

*     .. save old kinematics of common /MoTsai/
      Es0 = EsMT
      Ep0 = EpMT

      print*, 'URAD: nptot(NS):', nptot

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
*      np=nptot(NS)
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

*     .. Points Nos. ncdata,ncdata-1 are restored here.
*     .. No of last restored point
*      npcur(NS) = nptot(NS)-1
*      npcur(NS) = nptot(NS)

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

      print*, 'NS,np,Ep,Epmax', NS,np,Ep,Epmax

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

      relerr = 1.e-6
      abserr = 0.
*     .. integral over Ep
      termEp = 0.
      Epl = Ep+del
      Eph = Epmax
*     .. see Remarks on Programming at p.235
      Eph = Eph - 10.*eps
      if(DEBUG)print*, 'for IdEp: Ep,Epl,Eph:', Ep,Epl,Eph
      if (Eph-Epl .LE. eps) then
         if(DEBUG)print*, '--- IdEp: Epl > Eph. Skip integration'
         goto 1000
      endif
      call quanc8(fEp,Epl,Eph,abserr,relerr, IdEp,errest,nofun,flag)
      termEp = (del/Es)**(fs/2.)*IdEp
      if(DEBUG)print*, 'IdEp,termEp:', IdEp,termEp

 1000 continue

*     .. integral over Es
      termEs = 0.
      Esmin = (Ep + Kmin)/(1.-(Ep/Mtar)*(1-costheta))
      Esl = Esmin
*     .. see Remarks on Programming at p.235
      Esl = Esl + 10.*eps
      Esh = Es-del

*     .. max Mf available from the rest of spectra
      Mfmax = 0.
      do i=1,NStot
         if (i.EQ.NS) goto 100
         if (Mfdat(i,1).GT.Mfmax) Mfmax=Mfdat(i,1)
  100 enddo

*     .. Es1 corresponding to Mfmax
      Kmax = (Mfmax**2 - Mtar**2)/(2.*Mtar)
      Es1 = (Ep + Kmax)/(1.-(Ep/Mtar)*(1-costheta))

      if(DEBUG)print*, 'for IdEs: Es,Esl,Esh:', Es,Esl,Esh
      if (Esh-Esl .LE. eps) then
         if(DEBUG)print*, '--- IdEs: Esl > Esh. Skip integration'
         goto 2000
      endif

      if (Esh.GT.Es1) then
         IdEs = 0.
         DEBUG = .TRUE.
         print*, 'Split integration over Es: Es1,Esh=', Es1,Esh
         if (Es1-Esl .LE. eps) goto 2002
*        .. integral over available Es region
         call quanc8(fEs,Esl,Es1,abserr,relerr, IdEs,errest,nofun,flag)
*        .. substitute the rest by integration over Ep
 2002    continue
         IdEsEp = 0.
*        .. Eph corresponding Es1
*        .. process missing mass
         Mf2 = Mtar**2 + 2.*Mtar*(Es1-Ep) - 2.*Es1*Ep*(1.-costheta)
*        .. Ep of this process for incident energy Es
         recoil = 1. + (Es1/Mtar)*(1.-costheta)
         Eph = (Es1 - (Mf2-Mtar**2)/(2.*Mtar))/recoil

*        .. Epl corresponding Esh
*        .. process missing mass
         Mf2 = Mtar**2 + 2.*Mtar*(Esh-Ep) - 2.*Esh*Ep*(1.-costheta)
*        .. Ep of this process for incident energy Es
         recoil = 1. + (Esh/Mtar)*(1.-costheta)
         Epl = (Esh - (Mf2-Mtar**2)/(2.*Mtar))/recoil

         if(DEBUG)print*, '--- IdEsEp: Ep,Epl,Eph:', Ep,Epl,Eph
         if (Eph-Epl .LE. eps) then
            if(DEBUG)print*, '--- IdEsEp: Epl > Eph. Skip integration'
            goto 2001
         endif
         call quanc8(fEp,Epl,Eph,abserr,relerr,IdEsEp,errest,nofun,flag)
 2001    continue
         termEs = (del/Ep)**(fp/2.)*(IdEs+IdEsEp)
      else
         call quanc8(fEs,Esl,Esh,abserr,relerr, IdEs,errest,nofun,flag)
         termEs = (del/Ep)**(fp/2.)*IdEs
      endif
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
      sudat(np,NS) = sum
*     .. point np is unfolded now. Store its number in npcur(NS)
      npcur(NS) = np

      if (sum.LE.0.) DEBUG=.TRUE.

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
      real sMf,sinter
      real Is
      integer NS
      real Es,Ep,theta,b,bw,xs,sp,fs,ts,sec
      real Mf,Mf2
      integer ipoint

      fEs = 0.
*      print*, 'fEs: Es1=', Es1

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
         print*, 'ERROR fEs: Mf2 < 0.'
         STOP
      endif
      Mf = sqrt(Mf2)

      ipoint = 0
      do NS=1,NStot
         if (Mf.GT.Mfdat(1,NS)) goto 100
         sec = sMf(Mf,NS)
*         if (sec .LT. 0.) then
*            print*, 'fEs: Out of range'
*            STOP
*         endif
         if (sec .GT. 0.) then
            ipoint = ipoint+1
            vEs(ipoint)  = Esdat(NS)
            vsec(ipoint) = sec
         endif
  100 enddo

      sec = 0.
      if (ipoint.GT.0) then
         sec = sinter(Es1,ipoint,vEs,vsec)
      else
*         print*, 'ERROR fEs: ipoint=0.'
         if(DEBUG)print*,'ipoint=0, Es1,Mf', Es1,Mf
         call DWAIT('ERROR fEs: ipoint=0.')
      endif

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
      if (sec .LT. 0.) sec=0.
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
      Mfmin = Mfdat(nlow,NS)  - 2.*dMf
      Mfmax = Mfdat(nhigh,NS) + 2.*dMf
      if ((Mf.GE.Mfmin) .AND. (Mf.LE.Mfmax)) then
         sMf = sec
      else
*        .. invert sign of cross section for "unrealible" extrapolation
         sMf = -sec
         if (DEBUG) then
            print*, 'sMf: out of range. NS,nlow,nhigh:'
            print*, NS,nlow,nhigh
            print*, 'Mf,Mfdat(nlow,NS),Mfdat(nhigh,NS),Mfmin,Mfmax'
            print*, Mf,Mfdat(nlow,NS),Mfdat(nhigh,NS),Mfmin,Mfmax
         endif
      endif
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
         sec = sudat(nlow,NS) + dS*(Ep1-Epdat(nlow,NS))/dE
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
      real    Mf,Mp,M33
      parameter (Mp=.93827, M33=1.236)
      real me2
      parameter (me2 =.511e-3*.511e-3)
      real alp,pi
      parameter (alp=1./137.036, pi=3.141593)
      real costheta,xi,recoil,Z,b,bw
*     .. fill by some values
      EsMT = 20.
***      thetaMT = 5.
      thetaMT = 6.
      Mtar = Mp
      Mthr = .140
      Mf = M33
      costheta = cos(thetaMT*pi/180.)
      recoil = 1. + (EsMT/Mp)*(1-costheta)
*     .. scattered energy in 1236
      EpMT = (EsMT - (Mf**2-Mp**2)/(2.*Mp)) / recoil
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

      Trad   = .02
      tiwrad = .005
      tfwrad = .005
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

