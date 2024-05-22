c                                      @METAGS main
c============= .MAIN. ================
c
      character fstr*2, fname*13
      logical fexist
c
      common /INIRAN/ iniran
      common /IRUSED/ irused
      dimension NRAND(1),RASST(1),ALPHA(1)
      dimension nsec(1),sdelta(100),w(100,100),nsstat(100)
      dimension nfile(1)
      RASST(1) = -24.
      ALPHA(1) = 0.
      print*, 'RASST(1),ALPHA(1):', RASST(1),ALPHA(1)
      print*, 'Enter RASST,ALPHA'
      read *, RASST(1),ALPHA(1)
      print*, 'RASST(1),ALPHA(1):', RASST(1),ALPHA(1)

c     .. read the SP-02 matrix
      open (1, FILE='sp02mat.kum', STATUS='OLD')
      NF = 1
      JTR = 1
      call M02RE(NF,JTR)
      close(1)
c

      ch = 2.
      NRAND(1) = 1000000
*      minran   =   10000
      minran   = 0
      iniran   = 1

10000 print*, '---'
      print*, 'Current settings:'
      print*, 'Energy acceptance (in channels)'
      print*, ch
      print*, 'Number of randoms and minimal statistic/point'
      print*, NRAND(1), minran
      print*, 'Initial random integer number'
      print*, iniran

      print*, 'Enter new settings ("/" retains the old values)'
      print*, ' '
      print*, 'Enter energy acceptance (in channels)'
      read *, ch
*      print*, 'Enter number of randoms and minimal statistic/point'
*      print*, 'and initial random integer number'
*      read *, NRAND(1), minran, iniran
      print*, 'Enter number of randoms and initial integer'
      read *, NRAND(1), iniran

      call WAIT('Compute the new w-matrix')
      print*, '.. begin'
      call wmatr(NRAND(1),minran,RASST(1),ALPHA(1),ch,
     +           ns,sdelta,nsstat,w)
      print*, '.. end'
      nsec(1) = ns

c     .. write w to disk
      fname = 'Wxx.BIN'
 1000 print*, 'Enter file number (1:99) to write matrix w to the disk'
      print*, '  -1: don''t write to the disk'

      read *, nf
      if (nf.GT.99) goto 1000
      if (nf.LT.0) goto 20000
      if (nf.GT.0) then
         nfile(1) = nf
         write (fstr,1) nf
         if (fstr(1:1).EQ.' ') fstr(1:1)='0'
         fname(2:3) = fstr
      endif
      print*, 'Writing the file ', fname
c     ..         
      fexist = .TRUE.
      do while (fexist)
         inquire (FILE=fname, EXIST=fexist)
         if (fexist) then
            print*,'File ', fname, ' already exist.'
 2000       print*,'Enter new file number (1:99). Overwrite - 0'
            read *, nf
            if (nf.GT.99) goto 2000
            if (nf.EQ.0)  goto 3000
            write (fstr,1) nf
            if (fstr(1:1).EQ.' ') fstr(1:1)='0'
            fname(2:3) = fstr
         endif
      enddo

 3000 continue
      print*, irused, ' randoms have been used'
      print*, '--- Writing file ', fname
      call writew(fname,NRAND(1),RASST(1),ALPHA(1),
     +            ch,nsec(1),sdelta,nsstat,w)
20000 CONTINUE

      call WAIT('Calculate another w-matrix?')
      GOTO 10000

*      STOP
    1 FORMAT(I2.2)
      END

c                                   @METAGS WMATR
c------------ WMATR --------------
c
      subroutine wmatr(NRAND,minran,R,alpha,ch, ns,sdelta,nsstat,w)
*      parameter (minran=100000)
      parameter (d = .2)
*      parameter (Dcm = 1.76)
      parameter (Dcm = 1.81)
      parameter (TETA = 5.)
      parameter (NCENTR = 36)
      integer ns, nsstat(100)
      real sdelta(100), w(100,100)
      integer nnn(96)
      logical CROSS, cross1,cross2
      data nnn/96*0/
c
      X0 = 0.
      pi = acos(-1.)
      TETA0 = TETA*pi/180. * 1000.
      alphar = alpha*pi/180.
c
c     Find the min delta (i.e. min E) that crosses the chamber
c
      dstep = ch * 100.*d*(.01/Dcm)
      print*, 'delta step =', dstep
c
      cross1=.TRUE.
      cross2=.TRUE.
      delta0 = 0.
      do while (cross1.OR.cross2)
         delta0 = delta0-dstep
*        LOGICAL FUNCTION CROSS(DEL,Xcm,TH,Ycm,PH,Rcm,Nc, Xch,Nch)
         cross1 = CROSS(delta0,0.,-TETA0,0.,0.,R,NCENTR, x,n)
         cross2 = CROSS(delta0,0., TETA0,0.,0.,R,NCENTR, x,n)
      enddo
c     .. NB: ray delta0 out the chamber (#channel > 96)
      print*, 'min delta =', delta0
      call WAIT('')
c
c     Fill the matrix w
c
      do k=1,100
         sdelta(k) = 0.
         do j=1,100
            w(k,j)=0.
         enddo
      enddo
c
      ns = 0
      cross1=.TRUE.
      cross2=.TRUE.
      do while (cross1.OR.cross2)
         ns = ns+1
c        .. sdelta point to middle of delta interval
         sdelta(ns) = delta0 + dstep/2.
c--
         icount=0
         icross=0
*         do while((icount.LT.NRAND) .OR. (icross.LT.minran))
         do while(icount.LT.NRAND)
            icount = icount+1
            delta = delta0 + dstep*RANBAG(1*icount)
            th = -TETA0 + (2.*TETA0)*RANBAG(2*icount)
*           .. target 1X1 cm**2
            Xcm = -.5 + RANBAG(3*icount)
            Ycm = -.5 + RANBAG(4*icount)
            PH = 0.
*            LOGICAL FUNCTION CROSS(DEL,Xcm,TH,Ycm,PH,Rcm,Nc, Xch,Nch)
            if (CROSS(delta,Xcm,th,Ycm,PH,R,NCENTR, x,n)) then
               icross = icross+1
               nnn(n)=nnn(n)+1
            else
c              .. ray out the chamber
               continue
            endif
*            if (MOD(icount,NRAND/10).EQ.0) then 
*               print*, 'icount,icross:', icount,icross
*            endif
         enddo
         nsstat(ns) = icount
c--
         rand1 = 1./icount
         do k=1,96
            if (nnn(k).GT.0) then
               w(k,ns) = w(k,ns) + rand1*nnn(k)
               nnn(k) = 0
            endif
         enddo
         delta0 = delta0+dstep
         cross1 = CROSS(delta0,0.,-TETA0,0.,0.,R,NCENTR, x,n)
         cross2 = CROSS(delta0,0., TETA0,0.,0.,R,NCENTR, x,n)

         print*, 'ns,icount,icross:', ns,icount,icross
*         if (mod(ns,5).EQ.0) print*, 'Current ns =', ns

         if (ns.EQ.100) then
            print*, 'w array full. ns=100'
            GOTO 1000
         endif
      enddo
 1000 continue
      print*, 'max delta =', delta0-dstep
      print*, 'Resulting ns =', ns
      END

*                                         @METAGS M02RE
*                                         01-24-97 11:18pm
*--------------- M02RE ---------------
*
      SUBROUTINE M02RE(NF,JTR)
*      SUBROUTINE MTRDRF(NF,JTR)                                        00001000
*
*     DOUBLE PRECISION version
*
* 24-JAN-97. Changed variable names:
*            XF0 -> X0, THF0 -> TH0, YF0 -> Y0, PHF0 -> PH0
*            and common declaration for SPADIM
*            to make identical common declaration for SPADIM
*            in the routines MTRDRF and SPMRAY
*
*            Changed common declaration for SPARIM:
*            statement DIMENSION is removed
*            and array dimensions are declared in COMMON
*
C-----------------------------------------------------------------------00002000
C  HAáHAóEHàE: óTEHàE MATPàñõ Cè àá îOPMATHOâ áAèàCà                    00003000
C  èAPAMETPõ: NF: N îAâãA C MATPàñEâ                                    00004000
C     JTR: =1(2): Ñãü èPüMOâ (OÅPATHOâ) MATPàñõ                         00005000
C  COMMON-Åã: /SPADIM/: èEPEÑAET èAPAMETPõ MATPàñõ: DIM(55)+F0          00006000
C     /SPARIM/: èEPEÑAET èAPAMETPõ MATPàñõ: DIM(66)                     00007000
C-----------------------------------------------------------------------00008000
C                                                                       00009000
C  OÅECèEóEHàE óTEHàü MATPàñ Cè àá îAâãA N NF                           00010000
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      DIMENSION F0X(10),F0T(10),F0Y(10),F0P(10)                         00011000
*      DIMENSION AAX(55),AAT(55),AAY(55),AAP(55)                        00012000
*      COMMON /SPADIM/ AAX,XF0,AAT,THF0,AAY,YF0,AAP,PHF0                00013000
      COMMON /SPADIM/ AAX(55),X0,AAT(55),TH0,AAY(55),Y0H,AAP(55),PH0    ********
*       DIMENSION AX(66),AT(66),AY(66),AP(66)                           00014000
*       COMMON /SPARIM/ AX,AT,AY,AP                                     00015000
       COMMON /SPARIM/ AX(66),AT(66),AY(66),AP(66)                      00015000
*      use SPADIM,SPARIM
C                                                                       00016000
C  CóàTõBAHàE MACCàBOB AAX,..,AAP àá îAâãA N NF                         00017000
      REWIND NF                                                         00018000
      NPRST=5                                                           00019000
      DO 10 IAA=1,NPRST                                                 00020000
   10 READ (NF,20)                                                      00021000
   20 FORMAT(1X)                                                        00022000
      DO 30 IAA=1,55                                                    00023000
      READ (NF,40) AAX(IAA),AAT(IAA),AAY(IAA),AAP(IAA)                  00024000
       AX(IAA)=AAX(IAA)                                                 00025000
       AT(IAA)=AAT(IAA)                                                 00026000
       AY(IAA)=AAY(IAA)                                                 00027000
   30  AP(IAA)=AAP(IAA)                                                 00028000
   40 FORMAT(19X,4(1X,E12.5))                                           00029000
      DO 50 IAA=1,10                                                    00030000
      READ (NF,40) F0X(IAA),F0T(IAA),F0Y(IAA),F0P(IAA)                  00031000
       IAB=55+IAA                                                       00032000
       AX(IAB)=F0X(IAA)                                                 00033000
       AT(IAB)=F0T(IAA)                                                 00034000
       AY(IAB)=F0Y(IAA)                                                 00035000
   50  AP(IAB)=F0P(IAA)                                                 00036000
       X0=F0X(8)                                                        00037000
      TH0=F0T(7)                                                        00038000
       Y0=F0Y(10)                                                       00039000
      PH0=F0P(9)                                                        00040000
      IF(JTR.EQ.2)  X0=F0X(3)                                           00041000
      IF(JTR.EQ.2) TH0=F0T(2)                                           00042000
      IF(JTR.EQ.2)  Y0=F0Y(5)                                           00043000
      IF(JTR.EQ.2) PH0=F0P(4)                                           00044000
C                                                                       00045000
      RETURN                                                            00046000
      END                                                               00047000

*                                         @METAGS SPMRAY
*                                         01-24-97 10:16pm
*--------------- SPMRAY ---------------
*
      SUBROUTINE SPMRAY(DEL,X,TH,Y,PH,XS,THS,YS,PHS)                    00001000
*
*     DOUBLE PRECISION version
*
C-----------------------------------------------------------------------00002000
C  HAáHAóEHàE: èO èAPAMETPAM TPEKA óACTàñõ, èPà JTR=1(2): B CK 'BXOÑ    00003000
C     Cè' ('BõXOÑ Cè') BOCCTAHABãàBAET èAPAMETPõ TPEKA óACTàñõ B        00004000
C     CK 'BõXOÑ Cè' (CK 'BXOÑ Cè') MATPàóHõM CèOCOÅOM                   00005000
C  èAPAMETPõ: DEL: àMèìãúC óACTàñõ                          / % /       00006000
C     X,Y: KOOPÑ.èAPAMETPõ TPEKA óACTàñõ                    / MM /      00007000
C     TH,PH: ìÉãOBõE èAPAMETPõ TPEKA óACTàñõ                / MPAÑ /    00008000
C      /B CàCT.KOOPÑ.: JTR=1(2): 'BXOÑ Cè'('BõXOÑ Cè')/                 00009000
C     XS,YS: KOOPÑ.èAPAMETPõ TPEKA óACTàñõ                  / MM /      00010000
C     THS,PHS: ìÉãOBõE èAPAMETPõ TPEKA óACTàñõ              / MPAÑ /    00011000
C      /B CàCT.KOOPÑ.: JTR=1(2): 'BõXOÑ Cè'('BXOÑ Cè')/                 00012000
C  COMMON-Åã.: /SPADIM/: èEPEÑAET MATPàñì Cè: AX,ATH,AY,APH             00013000
C-----------------------------------------------------------------------00014000
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real DEL,X,TH,Y,PH,XS,THS,YS,PHS
c
      COMMON /SPADIM/ AAX(55),X0,AAT(55),TH0,AAY(55),Y0,AAP(55),PH0     00015000
      DIMENSION XX(5)                                                   00016000
*      use SPADIM
C                                                                       00017000
      NMAX=5                                                            00018000
      L=NMAX                                                            00019000
      MN=14                                                             00020000
      IF(NMAX.EQ.5) MN=20                                               00021000
C  èPEÑBAPàTEãúHOE áAHìãEHàE                                            00022000
      XS8=0.D0                                                          00023000
      THS8=0.D0                                                         00024000
      YS8=0.D0                                                          00025000
      PHS8=0.D0                                                         00026000
C  èOÑÉOTOBKA MACCàBA èEPEMEHHõX XX                                     00027000
      XX(1)=X                                                           00028000
      XX(2)=TH                                                          00029000
      XX(3)=Y                                                           00030000
      XX(4)=PH                                                          00031000
      XX(5)=DEL                                                         00032000
C  èEPEBOÑ KOOPÑàHAT àá MM B CM                                         00033000
      XX(1)=XX(1)/10.D0                                                 00034000
      XX(3)=XX(3)/10.D0                                                 00035000
C  BOCCT.ÉOPàá.èPOEKñàû TPEKA                                           00036000
      G0=0.D0                                                           00037000
      THS8=TH0                                                          00038000
      DO 10 I=1,NMAX                                                    00039000
      G0=G0+AAX(I)*XX(I)                                                00040000
      THS8=THS8+AAT(I)*XX(I)                                            00041000
      DO 10 K=I,NMAX                                                    00042000
      L=L+1                                                             00043000
      G0=G0+AAX(L)*XX(I)*XX(K)                                          00044000
      THS8=THS8+AAT(L)*XX(I)*XX(K)                                      00045000
      DO 10 J=K,NMAX                                                    00046000
      MN=MN+1                                                           00047000
      G0=G0+AAX(MN)*XX(I)*XX(K)*XX(J)                                   00048000
      THS8=THS8+AAT(MN)*XX(I)*XX(K)*XX(J)                               00049000
   10 CONTINUE                                                          00050000
C  èEPEBOÑ KOOPÑàHATõ X àá CM B MM                                      00051000
      XS8=(X0+G0)*10.D0                                                 00052000
C                                                                       00053000
      L=NMAX                                                            00054000
      MN=14                                                             00055000
      IF(NMAX.EQ.5) MN=20                                               00056000
C  BOCCT.BEPTàK.èPOEKñàû TPEKA                                          00057000
      V0=0.D0                                                           00058000
      PHS8=PH0                                                          00059000
      DO 20 I=1,NMAX                                                    00060000
      V0=V0+AAY(I)*XX(I)                                                00061000
      PHS8=PHS8+AAP(I)*XX(I)                                            00062000
      DO 20 K=I,NMAX                                                    00063000
      L=L+1                                                             00064000
      V0=V0+AAY(L)*XX(I)*XX(K)                                          00065000
      PHS8=PHS8+AAP(L)*XX(I)*XX(K)                                      00066000
      DO 20 J=K,NMAX                                                    00067000
      MN=MN+1                                                           00068000
      V0=V0+AAY(MN)*XX(I)*XX(K)*XX(J)                                   00069000
      PHS8=PHS8+AAP(MN)*XX(I)*XX(K)*XX(J)                               00070000
   20 CONTINUE                                                          00071000
C  èEPEBOÑ KOOPÑàHATõ Y àá CM B MM                                      00072000
      YS8=(Y0+V0)*10.D0                                                 00073000
c
      XS  = XS8
      THS = THS8
      YS  = YS8
      PHS = PHS8
C                                                                       00074000
      RETURN                                                            00075000
C     DEBUG INIT,SUBTRACE                                               00076000
      END                                                               00077000

*                                         @METAGS Mprint
*                                         01-25-97 00:24am
*--------------- Mprint ---------------
*
      SUBROUTINE Mprint()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /SPADIM/ AAX(55),X0,AAT(55),TH0,AAY(55),Y0,AAP(55),PH0
*      use SPADIM
      do IAA=1,55
         print 40, AAX(IAA),AAT(IAA),AAY(IAA),AAP(IAA)
      enddo
   40 FORMAT(4(1X,D12.5))
      END

*                                         @METAGS CROSS
*                                         01-25-97 03:28am
*--------------- CROSS ---------------
*
      LOGICAL FUNCTION CROSS(DEL,Xcm,TH,Ycm,PH,Rcm,Nc, Xch,Nch)
c     .. input pars
c     DEL = 100.*(E-Ec)/Ec, %
c     Xcm,Ycm in cm
c     TH,PH in mrad
c     Rcm - the coordinate of the center of chamber, cm
c     Nc - the central channel. Central particle cross the center of Nc
c     .. output pars
c     Xch - coordinate of the cross point, cm
c     Nch - cross channel
c
c     Returns .TRUE. if particle cross the chamber
c
      parameter (dcm=.2)
      CROSS = .FALSE.
      x0  = -X*10.
      y0  = -Y*10.
      ph0 = -PH
      call SPMRAY(DEL,x0,TH,y0,ph0,XS,THS,YS,PHS)
      Xch = -.1*XS - Rcm*TAN(.001*THS)
      Nch = .5 + Nc+(Xch/dcm)
      if ((Nch.GE.1) .AND. (Nch.LE.96)) CROSS=.TRUE.
      END
*                                         @METAGS WAIT
*                                         11-15-96 09:45pm
*--------------- WAIT ---------------
*
      SUBROUTINE WAIT(mess)
      character mess*(*), ch*1
      print*, mess
      print*, '<CR>=Continue, Q=Quit'
      read 1, ch
      if ((ch.EQ.'q') .OR. (ch.EQ.'Q')) STOP
      RETURN
    1 FORMAT(A)
      END

c                                   @METAGS READW
c------------- READW ----------------
c
      subroutine readw(fname,NRAND,R,alpha,ch,ns,sdelta,nsstat,w)
      parameter (lun = 40)
      character*(*) fname
      integer NRAND,ns,nsstat(100)
      real R,alpha,ch,sdelta(100),w(100,100)
c
      open (UNIT=lun, FILE=fname, STATUS='OLD',
     +      ACCESS='DIRECT', RECL=4)
      read (lun,REC=1) rand
      read (lun,REC=2) R
      read (lun,REC=3) alpha
      read (lun,REC=4) ch
      read (lun,REC=5) totns
      ns  = totns
      irec=5
      do j=1,ns
         irec = irec+1
         read (lun,REC=irec) sdelta(j)
      enddo
      do j=1,ns
         irec = irec+1
         read (lun,REC=irec) nsstat(j)
      enddo
      do j=1,100
         do i=1,100
            irec = irec+1
            read (lun,REC=irec) w(i,j)
         enddo
      enddo
      close (UNIT=lun)
      NRAND  = rand
      ns     = totns
      END ! READW
c                                      @METAGS WRITEW
c------------- WRITEW ----------------
c
      subroutine writew(fname,NRAND,R,alpha,ch,ns,sdelta,nsstat,w)
      parameter (lun = 40)
      character*(*) fname
      integer NRAND,ns,nsstat(100)
      real R,alpha,ch,sdelta(100),w(100,100)
c
      rand  = NRAND
      totns = ns
      open (UNIT=lun, FILE=fname, STATUS='UNKNOWN',
     +      ACCESS='DIRECT', RECL=4)
      write (lun,REC=1) rand
      write (lun,REC=2) R
      write (lun,REC=3) alpha
      write (lun,REC=4) ch
      write (lun,REC=5) totns
      ns  = totns
      irec=5
      do j=1,ns
         irec = irec+1
         write (lun,REC=irec) sdelta(j)
      enddo
      do j=1,ns
         irec = irec+1
         write (lun,REC=irec) nsstat(j)
      enddo
      do j=1,100
         do i=1,100
            irec = irec+1
            write (lun,REC=irec) w(i,j)
         enddo
      enddo
      close (UNIT=lun)
      END ! WRITEW

c                                   @METAGS RNOR
c------------- RNOR ---------------
c
      function RNOR(idummy)
*
*     Changed form of RNORMX to use RANBAG.
*
C        Generator of a vector of independent Gaussian-distributed 
C        (pseudo-)random numbers, of mean zero and variance one,
C        making use of a uniform pseudo-random generator (RANMAR).
C        The algorithm for converting uniform numbers to Gaussian
C        is that of "Ratio of Uniforms with Quadratic Bounds."  The
C        method is in principle exact (apart from rounding errors),
C        and is based on the variant published by Joseph Leva in
C        ACM TOMS vol. 18(1992), page 449 for the method and 454 for
C        the Fortran algorithm (ACM No. 712).
C        It requires at least 2 and on average 2.74 uniform deviates
C        per Gaussian (normal) deviate.
C   WARNING -- The uniform generator should not produce exact zeroes,
C   since the pair (0.0, 0.5) provokes a floating point exception.
      common /NORPAR/ amean,sigma
      SAVE  S, T, A, B, R1, R2
      DIMENSION U(2)
      DATA  S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA  R1, R2/ 0.27597, 0.27846/
C         generate pair of uniform deviates
      jdummy = idummy
   50 U(1) = ranbag(jdummy)
      idummy = idummy+1
      U(2) = ranbag(jdummy)
      V = 1.7156 * (U(2) - 0.5)
      X = U(1) - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)
C           accept P if inside inner ellipse
      IF (Q .LT. R1)  GO TO 100
C           reject P if outside outer ellipse
      IF (Q .GT. R2)  GO TO 50
C           reject P if outside acceptance region
      IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2)  GO TO 50
C           ratio of P's coordinates is normal deviate
  100 DEVIAT = V/U(1)
      rnor = amean + sigma*DEVIAT
      RETURN
      END ! rnor
c                                      @METAGS RANBAG
c---------- RANBAG ----------
c
      function ranbag(idummy)
      parameter (LEN=100000)
      parameter (LUX=3)
      real rvec(LEN)
      integer ivec(25)
      logical first
      common /INIRAN/ iniran
      common /IRUSED/ irused
      save rvec,ivec,ncur
      data first/.TRUE./
c
      irused = irused+1
      if (first) then
c        .. initialization
         first=.FALSE.
         call RLUXGO(LUX,iniran,0,0)
c        .. create first bag
         call RANLUX(rvec,LEN)
         call RLUXUT(ivec)
         ncur=0
      endif
      ncur = ncur+1
      if (ncur.GT.LEN) then
*         print*, 'Create next bag'
c        .. create next bag
         call RLUXIN(ivec)         
         call RANLUX(rvec,LEN)
         call RLUXUT(ivec)
         ncur=1
      endif
      ranbag = rvec(ncur)
      END ! ranbag
c
c                                      @METAGS INIBAG
c
      BLOCK DATA INIBAG
      common /INIRAN/ iniran
      common /IRUSED/ irused
      common /NORPAR/ amean,sigma
      data iniran/1/
      data irused/0/
      data amean/0./, sigma/1./
      end
C
C                                   @METAGS RANLUX
C
      SUBROUTINE RANLUX(RVEC,LENV)
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++ 
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
C                               default
C  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT  
         INSEED = JSEED
***         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
***         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
***     +        LUXLEV,'      p =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
  195    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
***      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
*      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
***          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
***     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
***          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
***         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
***     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
***          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')   
     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
***        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
***     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')  
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
     +     K1, K2, ' cannot occur at luxury level', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
