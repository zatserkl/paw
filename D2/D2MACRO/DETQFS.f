      program DETqfs
*      IMPLICIT NONE
      PARAMETER (NWPAWC = 300000)
      PARAMETER (LRECL  = 1024)
      PARAMETER (id = 1)
      PARAMETER (idh = 2, NXE=100, NYth=40)
      COMMON /PAWC/ IPAW(NWPAWC)

      character CHTOP*8, CHFILE*32
      character*80 chtitl
      REAL F,E1,E2,THin,THout,Xin,Xout
      INTEGER pdf,PChit,PCch,PCcl, Thit,Tch,Tcl, err
*      CHARACTER distr*4
      COMMON/NTDATA/ F,E1,E2,THin,THout,Xin,Xout,
     &               pdf,PChit,PCch,PCcl,Thit,Tch,Tcl, err
*      COMMON/NTDATAC/ distr

      external QFS
      logical EVENT
      character*32 thisf,thisn,ifile,bfile, afile
      character*80 tit
      character*32 spdf
      character*127 CHPATH

      data CHTOP /'NTDIR'/
      data CHFILE/'NTqfs.hbook'/
      data afile /' '/
      data nhitROS,lostROS,nhitSPE,lostSPE,
     &     nhitQFS,lostQFS,nhitUNI,lostUNI /8*0/

      call GETARG(0,thisf)
      call FPARSE(thisf,thisn,thisf,'x')
      lename = LENOCC(thisn)
      ifile = thisn(1:lename)//'.ini'
      bfile = thisn(1:lename)//'_t.for'

      print*, ' '
      print*, 'This file ', thisf
      print*, 'Ini file  ', ifile
      print*, 'HBOOK file with Ntuple id =', id, '   ', CHFILE
      print*, 'Batch file for PAW usage     ', bfile
      print*, ' '

      lun = LUNFREE(1)
      open (lun, FILE=ifile, STATUS='OLD', ERR=10000)
      read (lun,*) nran
      read (lun,*) init
      read (lun,*) F
      read (lun,*) Rchamb
      read (lun,*) dRchamb
      read (lun,*) ALPdeg
      read (lun,*) Rtitov
      read (lun,*) dRtitov
      read (lun,*) Espread
      read (lun,*) Xbeam
      read (lun,*) Rbeam
      read (lun,*) E0
      read (lun,*) theta
      read (lun,*) thacc
      close(lun)

*     .. keep 3 decimal digits
      F = AINT(0.5 + F*1000.)/1000.

      print*, 'Ini file data:'
      print*, 'N randoms/bin =', nran
      print*, 'F =', F, ' MHz'
      print*, 'Rchamb  =', Rchamb, ' cm'
      print*, 'dRchamb =', dRchamb, ' cm'
      print*, 'ALPdeg  =', ALPdeg, ' degree'
      print*, 'Rtitov  =', Rtitov, ' cm'
      print*, 'dRtitov =', dRtitov, ' cm'
      print*, 'Espread =', Espread, '%'
      print*, 'Xbeam   =', Xbeam, ' mm'
      print*, 'Rbeam   =', Rbeam, ' mm'
      print*, 'E0   =', E0, ' MeV'
      print*, 'theta and acceptance:', theta, thacc

c     .. convert
      pi = ACOS(-1.)
      torad = pi/180.
      todeg = 180./pi
      ALPr  = ALPdeg*torad
      thetar  = torad*theta
      thaccr  = torad*thacc
      thmin = theta - thacc/2.
      thmax = theta + thacc/2.
      thminr = torad*thmin
      thmaxr = torad*thmax

      Ec = 17.61*F
*     .. Energy limits used by QFS and UNI distributions
      Emin = 0.90*Ec
      Emax = 1.05*Ec

      print*, ' '
      npar = IARGC()
      if (npar.GT.0) then
         call GETARG(1,spdf)
      else
         spdf = 'Q'
         print*, 'Enter the number of probability density function:'
         print*, ' R:    for Rosenbluth distribution'
         print*, ' S:    for Rosenbluth energy and uniform angle distr.'
         print*, ' Q:    for QFS distribution'
         print*, ' U:    for uniform distribution'
         print*, 'or any combination'
         print*, '<CR>=', spdf
         read '(A)', spdf
      endif
      call CLTOU(spdf)
      
      write(chtitl,1) E0, F
    1 format ('Exposition for E0 =', F6.1, ' MeV, F = ', F6.3, ' MHz')
      call MESS('Ntuple title:')
      call MESS(chtitl)

      CALL HLIMIT(NWPAWC)
*     .. open a new RZ file
      ntlun = LUNFREE(1)
      CALL HROPEN(ntlun,CHTOP,CHFILE,'N',LRECL,ISTAT)
*     .. book Ntuple
      CALL HBNT(id,chtitl,' ')
*     .. define Ntuple
      CALL HBNAME(id,'NTDATA', F, 'F:R, E1:R, E2:R,
     &     THin:R, THout:R, Xin:R, Xout:R, pdf[1,4]:U,
     &     PChit[0,1]:U, PCch[1,96]:U, PCcl[0,1]:U,
     &      Thit[0,1]:U,  Tch[1,11]:U,  Tcl[0,1]:U, err[0,1]:U')
*      CALL HBNAMC(id, 'NTDATA', distr, 'distr:C')
      
*     .. Fill SP commons
      call FILLSP()

      call BAGINT(init)
      print*, ' '

      call FILELUN(ntlun,afile)
      print*, 'With unit', ntlun, ' connected file ', afile
      print*, ' '
      call WAIT('Begin the calculation')
      print*, 'Running ..'

      tlim = 10000.
      call TIMEST(tlim)

      if (INDEX(spdf,'Q').GT.0) then
         print*, 'Create the sample hist'
         write(tit,2) E0, F
    2    format ('Spectrum for E0 =', F6.1, ' MeV, F = ', F6.3, ' MHz')
*        .. set QFS parameters
         Z  = 1.
         A  = 2.
         EPS = 2.22
         EPSD = 2.
         Pf = 40.
         irad = 0
         call SetQFS(E0,Z,A,EPS,EPSD,Pf,irad)
         vmax=QFSMAX(idh,tit,E0,Espread,NXE,Emin,Emax,NYth,thmin,thmax)
         print*, 'Choose vmax for random generator as ', vmax
         print*, 'Enter desirable vmax value or <CR> to continue'
         read '(G15.0)', val
         if (val.NE.0.) then
            vmax = val
            print*, 'Set vmax for random generator for ', vmax
            print*, 'Press <CR> to continue'
            read '(I5)', i
         endif
      endif

*     .. fit to experimental data
      call FOCPAR(Ec, Dexp,Rfoc,amudeg)
*     .. set experimental dispersion
      call SetDmm(Dexp)

      do i=1,nran
         idummy = i
*        .. incident particle
         call BEAM(E0,Espread,Xbeam,Rbeam, E1MeV,X0mm)

*                                         @METAGS ROS_distr
*        Rosenbluth distributed angle
*        ----------------------------
         if (INDEX(spdf,'R').EQ.0) goto 200
         pdf = 1
*         distr = 'FROS'
*        .. My compiler inverts letter sequence in character fields of Ntuple!
*        .. Correction:
*         call INVNTC(distr)
         call ROSRAN(E1MeV,thminr,thmaxr, E2MeV,thranr)
*        .. ray angle with respect to spectrometer angle
*        .. thranr already in rad
         th0 = thranr-thetar
*        .. get output ray
         call RAY(E2MeV,Ec,X0mm,th0, x1,th1)
*        .. store ray
         E1    = E1MeV
         E2    = E2MeV
         Xin   = X0mm
         Xout  = x1
         THin  = th0*todeg
         THout = th1*todeg
         err   = 0

         if (EVENT(x1,th1,Rfoc,Rchamb,dRchamb,Rtitov,dRtitov,ALPr)) then
            nhitROS = nhitROS+1
         else
            lostROS = lostROS+1
         endif

*        .. store Rosenbluth event
         call HFNT(id)

*                                         @METAGS SPE_distr
*        Rosenbluth energy and uniform distributed angle
*        -----------------------------------------------
  200    continue
         if (INDEX(spdf,'S').EQ.0) goto 300
         pdf = 2
*         distr = 'FSPE'
*        .. My compiler inverts letter sequence in character fields of Ntuple!
*        .. Correction:
*         call INVNTC(distr)

         idummy = idummy+1
         th0 = thaccr*(RANBAG(idummy)-.5)
*        .. get output ray
         call RAY(E2MeV,Ec,X0mm,th0, x1,th1)
*        .. store ray
         E1    = E1MeV
         E2    = E2MeV
         Xin   = X0mm
         Xout  = x1
         THin  = th0*todeg
         THout = th1*todeg
         err   = 0

         if (EVENT(x1,th1,Rfoc,Rchamb,dRchamb,Rtitov,dRtitov,ALPr)) then
            nchitSPE = nhitSPE+1
         else
            lostSPE = lostSPE+1
         endif

*        .. store uniform event
         call HFNT(id)

*                                         @METAGS QFS_distr
*        QFS distribution
*        ----------------
  300    continue
         if (INDEX(spdf,'Q').EQ.0) goto 400
         pdf = 3
*         distr = 'FQFS'
*        .. My compiler inverts letter sequence in character fields of Ntuple!
*        .. Correction:
*         call INVNTC(distr)

**        .. random point for constant initial energy E0
*         call RAN2DH(idh,E2Mev,thran)
**        .. ray angle with respect to spectrometer angle
*         th0 = torad*(thran-theta)

*        .. random point for random initial energy E1MeV
*        .. set E1MeV as initial energy
         call SetQFS(E1MeV,Z,A,EPS,EPSD,Pf,irad)
         call RFUN2(QFS,vmax,Emin,Emax,thmin,thmax, E2MeV,thran,iover)

*        .. ray angle with respect to spectrometer angle
         th0 = torad*(thran-theta)
*        .. get output ray
         call RAY(E2MeV,Ec,X0mm,th0, x1,th1)
*        .. store ray
         E1    = E1MeV
         E2    = E2MeV
         Xin   = X0mm
         Xout  = x1
         THin  = th0*todeg
         THout = th1*todeg
         err   = iover

         if (EVENT(x1,th1,Rfoc,Rchamb,dRchamb,Rtitov,dRtitov,ALPr)) then
            nhitQFS = nhitQFS+1
         else
            lostQFS = lostQFS+1
         endif

*        .. store QFS event
         call HFNT(id)

*                                         @METAGS UNI_distr
*        UNI distribution
*        ----------------
  400    continue
         if (INDEX(spdf,'U').EQ.0) goto 100
         pdf = 4
*         distr = 'FUNI'
*        .. My compiler inverts letter sequence in character fields of Ntuple!
*        .. Correction:
*         call INVNTC(distr)

*        .. random scattered energy Emin..Emax
         idummy = idummy+1
         E2MeV = Emin + (Emax-Emin)*RANBAG(idummy)
*        .. random ray angle thminr..thmaxr
         idummy = idummy+1
         th0 = thaccr*(RANBAG(idummy) - 0.5)
*        .. get output ray
         call RAY(E2MeV,Ec,X0mm,th0, x1,th1)
*        .. store ray
         E1    = 0.
         E2    = E2MeV
         Xin   = X0mm
         Xout  = x1
         THin  = th0*todeg
         THout = th1*todeg
         err   = 0

         if (EVENT(x1,th1,Rfoc,Rchamb,dRchamb,Rtitov,dRtitov,ALPr)) then
            nhitUNI = nhitUNI+1
         else
            lostUNI = lostUNI+1
         endif

*        .. store UNI event
         call HFNT(id)

  100    if (MOD(i,1000).EQ.0) print*, 'Processed', i, ' events'
      enddo

      if (INDEX(spdf,'R').NE.0) then
         print*, 'nhitROS =', nhitROS, ',  lostROS =', lostROS
         call ROSREP(nhit,miss)
      endif

      if (INDEX(spdf,'S').NE.0)
     &         print*, 'nhitSPE =', nhitSPE, ',  lostSPE =', lostSPE
      
      if (INDEX(spdf,'Q').NE.0) then
         print*, 'nhitQFS =', nhitQFS, ',  lostQFS =', lostQFS
         call RFUN2REP(nhit,miss,nover)
      endif
      
      if (INDEX(spdf,'U').NE.0)
     &         print*, 'nhitUNI =', nhitUNI, ',  lostUNI =', lostUNI
      
      dummy = RNOREP(nhit,miss)
      call BAGREP()
      
      print*, ' '
      call TIMEX(time)
      print*,'Job time', time, ' seconds'
      print*, ' '

      call LOOK(id)

      CALL HPRNT(id)
*
*-- write batch version of analisys routine to file
*
      lunb = LUNFREE(1)
      OPEN(lunb, FILE=bfile, STATUS='UNKNOWN')
      CALL HUWFUN(lunb, id, thisn(1:lename), 0, 'B')
*
*-- write Ntuple buffer to disk and close RZ file
*
      CALL HROUT(0,ICYCLE,'T')
      print*, 'After HROUT(0,ICYCLE,''T'') ICYCLE is', ICYCLE

 1000 continue
      CHPATH = ' '
      call HCDIR(CHPATH,'R')
      print*, ' '
      call HLDIR('//PAWC','T')
      print '('' Current directory is '', A)', CHPATH(1:LENOCC(CHPATH))
      call HLDIR(' ','T')

      print*, 'Enter ID to write to file. (0=Continue)'
      read '(I8)', idwr
      if (idwr.NE.0) then
         CALL HROUT(idwr,ICYCLE,' ')
         print*, 'ICYCLE =', ICYCLE
         goto 1000
      endif

      CALL HREND(CHTOP)
      STOP
10000 print*, 'Ini file ', ifile(1:LENOCC(ifile)), ' did not found'
      END

*                                         @METAGS QFSMAX
*                                         04-17-98 05:04pm
*--------------- QFSMAX ---------------
*
      function QFSMAX(idh,tit,E0,Espread,NXE,Emin,Emax,NYth,thmin,thmax)
      external QFS
      character*80 tit
      logical HEXIST

*     .. get QFS parameters (store energy E00)
      call GetQFS(E00,Z,A,EPS,EPSD,Pf,irad)
*     .. set QFS parameters
      E01 = E0
      call SetQFS(E01,Z,A,EPS,EPSD,Pf,irad)
      call HBFUN2(idh,tit,NXE,Emin,Emax,NYth,thmin,thmax,QFS)
*     .. simple search of function maximum in E0 region
      vmax0 = HMAX(idh)
      print*, ' '
      print*, 'The sample hist is created'
      print*, 'Scan to estimate the maximum'
      print*, ' '


      idtmp = 1
      do while (HEXIST(idtmp))
         idtmp = idtmp+1
      enddo
      print*, 'QFSMAX: Use for temporary hists ID =', idtmp
      print*, 'Working..'

*     .. may be min..
      E02 = E0*(1. + 3.*0.01*Espread)
      call SetQFS(E02,Z,A,EPS,EPSD,Pf,irad)
      call HBFUN2(idtmp,' ',NXE,Emin,Emax,NYth,thmin,thmax,QFS)
      vmax1 = HMAX(idtmp)
      call HDELET(idtmp)
*     .. may be max..
      E03 = E0*(1. - 3.*0.01*Espread)
      call SetQFS(E03,Z,A,EPS,EPSD,Pf,irad)
      call HBFUN2(idtmp,' ',NXE,Emin,Emax,NYth,thmin,thmax,QFS)
      vmax2 = HMAX(idtmp)
      call HDELET(idtmp)

*     .. restore QFS energy
      call SetQFS(E00,Z,A,EPS,EPSD,Pf,irad)

      vmin = MIN(vmax1,vmax0,vmax2)
      vmax = MAX(vmax1,vmax0,vmax2)
*     .. choose as the upper limit in E0 region
      vmax = vmax + (vmax-vmin)

      print*, 'QFSMAX: The maximum search results:'
      print '(''E0 ='', F6.1,'' MeV: vmax ='', E9.2)', E03, vmax2
      print '(''E0 ='', F6.1,'' MeV: vmax ='', E9.2)', E02, vmax0
      print '(''E0 ='', F6.1,'' MeV: vmax ='', E9.2)', E01, vmax1
      print*, 'As maximum value looks ', vmax
      QFSMAX = vmax
      END

*                                         @METAGS EVENT
*                                         03-26-98 00:11am
*--------------- EVENT ---------------
*
      logical function EVENT(x1,th1,Rfoc,RC,dRC,RT,dRT,ALPr)
      logical Cevent,Tevent
*      CHARACTER distr*4
      REAL F,E1,E2,THin,THout,Xin,Xout
      INTEGER pdf,PChit,PCch,PCcl, Thit,Tch,Tcl, err
*      CHARACTER distr*4
      COMMON/NTDATA/ F,E1,E2,THin,THout,Xin,Xout,
     &               pdf,PChit,PCch,PCcl,Thit,Tch,Tcl, err
*      COMMON/NTDATAC/ distr
      PChit = 0
      PCch = 0
      PCcl = 0
      Thit = 0
      Tch = 0
      Tcl = 0
*---                       @METAGS PC_ClasterResults
*
*        for ALPdeg=0
*     .. dRc = 0.1 cm gets: ROS -  2.5%, UNI -  2.5%
*     .. dRc = 0.2 cm gets: ROS -  5.1%, UNI -  4.9%
*     .. dRc = 0.4 cm gets: ROS - 10.5%, UNI - 10.1%
*     .. dRc = 1.0 cm gets: ROS - 26.2%, UNI - 25.1%
*     .. dRc = 2.0 cm gets: ROS - 51.0%, UNI - 49.6%
*        for ALPdeg=4
*     .. dRc = 2.0 cm gets: ROS - 51.3%, UNI - 49.7%
*---
      RC1 = RC
      if (Cevent(x1,th1,Rfoc,RC1,ALPr, x,n1)) then
         EVENT = .TRUE.
         PChit = 1
         PCch = n1
         PCcl = 0
*        .. test PC clasters
         RC2 = RC+dRC
         if (Cevent(x1,th1,Rfoc,RC2,ALPr, x,n2)) then
            if (n2.NE.n1) PCcl=1
         endif
*        .. hit Titov
*---                       @METAGS Titov_ClasterResults
*
*        .. dRT = 0.8 cm gets: ROS - 3.0%, UNI - 2.5%
*---
         RT1 = RT+dRT/2.
         if (Tevent(x1,th1,Rfoc,RT1, x,n1)) then
            Thit = 1
            Tch  = n1
            Tcl  = 0
*           .. test Titov clasters
            RT2 = RT-dRT/2.
            if (Tevent(x1,th1,Rfoc,RT2, x,n2)) then
               if (n2.NE.n1) Tcl=1
            endif
         endif
      else
         EVENT = .FALSE.
      endif
      END

*                                         @METAGS LOOK
*                                         03-26-98 01:30am
*--------------- LOOK ---------------
*
      SUBROUTINE LOOK(id)
      REAL F,E1,E2,THin,THout,Xin,Xout
      INTEGER pdf,PChit,PCch,PCcl, Thit,Tch,Tcl, err
*      CHARACTER distr*4
      COMMON/NTDATA/ F,E1,E2,THin,THout,Xin,Xout,
     &               pdf,PChit,PCch,PCcl,Thit,Tch,Tcl, err
*      COMMON /NTDATAC/ distr

 1000 I = 1
      print*, 'Enter a number of event to look (<CR>=Exit)'
      read 1, I
    1 format(I8)
      if (I.EQ.0) goto 2000
      CALL HGNT(id, I, IER)
      IF (IER .NE. 0) THEN
         PRINT *, 'Error reading row ', I
      ENDIF
      PRINT *,'Event',I,' E1,E2,PCch,pdf:',E1,E2,PCch,pdf
      goto 1000
 2000 continue
      END

*                                         @METAGS Cevent
*                                         03-26-98 00:11am
*--------------- Cevent ---------------
*
      logical function Cevent(x1,th1,Rfoc,Rcm,ALPr, Xch,Nch)
      logical XCROSS
      call CHAMBER(96,36,2.)
      if(XCROSS(x1,th1,Rfoc,Rcm,ALPr, Xch,Nch)) then
         Cevent = .TRUE.
      else
         Xch = -1000.
         Nch = -1000
         Cevent = .FALSE.
      endif
      END
*                                         @METAGS Tevent
*                                         03-26-98 00:11am
*--------------- Tevent ---------------
*
      logical function Tevent(x1,th1,Rfoc,Rcm, Xch,Nch)
      logical XCROSS
      call CHAMBER(11,4,14.)
      if(XCROSS(x1,th1,Rfoc,Rcm,0., Xch,Nch)) then
         Tevent = .TRUE.
      else
         Xch = -1000.
         Nch = -1000
         Tevent = .FALSE.
      endif
      END

*                                         @METAGS RAN2DH
*                                         04-14-98 04:49pm
*--------------- RAN2DH ---------------
*
      SUBROUTINE RAN2DH(id, x,y)
c
c     Random number according to 2-dim histo. Uses RANBAG
c     Uses Von Neumann acceptance-rejection method
c     Adjusting parameters:
c     LEN - bag size
c
      parameter (LEN=1000)
      integer id
      real x,y
      integer id0
      character*80 CHTITL
      real xbag(LEN),ybag(LEN)
      save xbag,ybag,ncur,id0,idummy
*     .. entry R2DHREP
      real ratio
      integer nhit,miss, nhit1,miss1
      save    nhit,miss
c
      data ncur/LEN/, id0/0/, idummy/0/
      data nhit,miss /0,0/
c
      if (id0.EQ.0) id0=id
      if (id.NE.id0) then
         print*, 'RAN2DH: Sample historgam is changed:', id0,' --> ', id
         print*, 'Recalculate the bag'
         id0  = id
         ncur = LEN
         nhit = 0
         miss = 0
      endif
c
      ncur = ncur+1
      if (ncur.GT.LEN) then
c        .. create new bag
         call HGIVE(id,CHTITL,NX,XMI,XMA,NY,YMI,YMA,NWT,LOC)
*        .. the maximum value
         vmax = HMAX(id)
         dX = XMA-XMI
         dY = YMA-YMI
         do i=1,LEN
  100       continue
*           .. choose some point
            idummy = idummy+1
            xtry = XMI + dX*RANBAG(idummy)
            idummy = idummy+1
            ytry = YMI + dY*RANBAG(idummy)
            hval = HXY(id,xtry,ytry)
*           .. get some random 0..vmax
            idummy = idummy+1
            rval = vmax*RANBAG(idummy)
            if (rval.LE.hval) then
*              .. success
               xbag(i) = xtry
               ybag(i) = ytry
               nhit = nhit+1
            else
*              .. try another point
               miss = miss+1
               goto 100
            endif
         enddo
*        .. clear counter
         ncur=1
      endif
      x = xbag(ncur)
      y = ybag(ncur)
      RETURN
*
*---     @METAGS R2DHREP
*
      ENTRY R2DHREP(nhit1,miss1)
*     .. Entry to report RAN2DH
      ratio = REAL(miss)/REAL(nhit)
      print*, '*** RAN2DH report to generated randoms:'
      print*, 'Generated ', nhit, ' points, missing shots =', miss
      print*, 'There are ', ratio, ' missings per hit'
      print*, '***'
      nhit1 = nhit
      miss1 = miss
      END

*                                         @METAGS RFUN2
*                                         04-14-98 04:49pm
*--------------- RFUN2 ---------------
*
      SUBROUTINE RFUN2(fun,vmax,xmin,xmax,ymin,ymax, x,y,iover)
c
c     Random number according to 2-dim function. Uses RANBAG
c     Uses Von Neumann acceptance-rejection method
c     vmax must be GE the max function value in the x-y range
c
      real fun,vmax,xmin,xmax,ymin,ymax,vmax, x,y
      integer iover
      save idummy
*     .. entry RFUN2REP
      real ratio
      integer nhit,miss,nover, nhit1,miss1,nover1
      save    nhit,miss,nover
c
      data idummy/0/
      data nhit,miss,nover /0,0,0/
c
      iover = 0
      dx = xmax-xmin
      dy = ymax-ymin
  100 continue
*     .. choose some point
      idummy = idummy+1
      xtry = xmin + dx*RANBAG(idummy)
      idummy = idummy+1
      ytry = ymin + dy*RANBAG(idummy)
      fval = fun(xtry,ytry)
      if (fval.GT.vmax) then
         iover = 1
         nover = nover+1
         print*, 'RFUN2: Function value', fval, ' for', xtry,ytry
         print*, 'exceeds the fuction limit vmax =', vmax
         call MESS('RFUN2 error: vmax too low')
      endif
*     .. get some random 0..vmax
      idummy = idummy+1
      rval = vmax*RANBAG(idummy)
      if (rval.LE.fval) then
*        .. success
         x = xtry
         y = ytry
         nhit = nhit+1
      else
*        .. try another point
         miss = miss+1
         goto 100
      endif
      RETURN
*
*---     @METAGS RFUN2REP
*
      ENTRY RFUN2REP(nhit1,miss1,nover1)
*     .. Entry to report RFUN2
      ratio = REAL(miss)/REAL(nhit)
      print*, '*** RFUN2 report to generated randoms:'
      print*, 'Generated ', nhit, ' points, missing shots =', miss
      print*, 'There are ', ratio, ' missings per hit'
      if (nover.GT.0) then
         print*, '*********************'
         print*, '*   RFUN2 WARNING   *'
         print*, '*********************'
      endif
      print*, 'The number of overflows ', nover
      print*, '***'
      nhit1  = nhit
      miss1  = miss
      nover1 = nover
      END
