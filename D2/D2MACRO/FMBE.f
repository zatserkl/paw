      program FMBEcommon
*
*     Summary FMBE file name
*
*     Input information
*
      character Energy*4, theta*5
      character*32 fHIS,fMBE, fname
*
*--
*
      integer NWPAWC
      parameter (NWPAWC=30000)
      common /PAWC/ PAW(NWPAWC)

      parameter (MAXDIR=100)
      character*8 chpath(MAXDIR), uppath,lopath
      character*8 chtype,chopt, chtitl*80

      logical fexist
      character*8 MBE(4), MBEtmp
      character MBEs*4, MBEs1*4, MBEs2*5
      data MBEs1/'MBE='/, MBEs2/'MBE ='/

      print*, ' '
      npar = IARGC()
      if (npar.GT.0) then
         call GETARG(1,fHIS)
      else
         fHIS = ' '
         print*, 'Enter the histogram file name'
         read '(A)', fHIS
      endif

      call FPARSE(fHIS,fname,fHIS,'HIS')
      lename = LENOCC(fname)
      fmbe = fname(1:lename)//'.MBE'

      print*, ' '
      print*, 'Histogram file ', fHIS
      print*, 'Output file    ', fMBE
      print*, ' '

      print*, 'For comment'
      Energy = ' '
      print*, 'Enter energy'
      read '(A)', Energy
      theta = ' '
      print*, 'Enter theta'
      read '(A)', theta

      print*, ' '
      print*, 'Energy = ', Energy, ' MeV, theta = ', theta, ' degree'
      print*, 'Histogram file ', fHIS
      call WAIT(' ')

      lunhis = LUNFREE(1)
*--   .. Tell HBOOK how many words are in PAWC
      call HLIMIT(NWPAWC)
      call HROPEN(lunhis,'D2',fHIS,' ',1024,ISTAT)
      if (ISTAT.NE.0) goto 10000

      inquire (FILE=fMBE, EXIST=fexist)
      if (fexist) then
         print*, 'File ', fMBE, ' already exist!'
         STOP
      endif

      lunout = LUNFREE(2)
      open (lunout, FILE=fMBE, STATUS='NEW', ERR=20000)
      write(lunout,1)
      write(lunout,2) Energy,theta
      write(lunout,3)
      write(lunout,4)
      write(lunout,5)
    1 format ('*')
    2 format ('* F and MBE/1e6 for E = ', A, ' MeV, ', A, ' degree')
    3 format ('*')
    4 format ('* F     D2      AL      D2+     AL+')
    5 format ('*')

*     .. store ndir directory names in the array chpath
      call HRDIR(MAXDIR,chpath,ndir)

      do i=1,ndir
         print*, chpath(i)
      enddo
      call WAIT('List of directories in histogram file '//fHIS)

      do i=1,ndir
         print*, chpath(i)
         uppath = chpath(i)
         lopath = chpath(i)
         call CLTOU(uppath)
         call CUTOL(lopath)
         if (uppath.NE.lopath) then
*           .. Non-numeric path
            goto 2000
         endif

         call HCDIR('\\' // chpath(i),' ')
         do nMBE=1,4
            MBE(nMBE) = '0.'
         enddo
         
         idh = 0
  100    continue
         chtype = '1'
         chopt  = ' '
         call HLNEXT(idh,chtype,chtitl,chopt)
         if (idh.NE.0) then
            if (idh.EQ.1000) then
               nMBE = 1
            elseif (idh.EQ.2000) then
               nMBE = 2
            elseif (idh.EQ.3000) then
               nMBE = 3
            elseif (idh.EQ.4000) then
               nMBE = 4
            else
               print*, 'Skip histogram', idh
               goto 100
            endif

            qMBE = 1e9
            length = LEN(chtitl)
*           .. match 'MBE='
            MBEs = MBEs1
            nshift = LEN(MBEs1)
            n = INDEX(chtitl,MBEs1)
            if (n.GT.0) then
               read (chtitl(n+LEN(MBEs1):length),*) qMBE
               goto 1000
            endif

*           .. match 'MBE ='
            MBEs = MBEs2
            nshift = LEN(MBEs2)
            n = INDEX(chtitl,MBEs2)
            if (n.GT.0) then
               read (chtitl(n+LEN(MBEs1):length),*) qMBE
               goto 1000
            endif

            print*, 'ID =', idh
            print*, chtitl
            call WAIT('MBE field did not found')
 
 1000       rMBE = qMBE/1e6
            write (MBEtmp,10) rMBE
   10       format (F8.3)
            n=1
            do while (MBEtmp(n:n).EQ.' ')
               n = n+1
            enddo
            MBE(nMBE) = MBEtmp(n:LEN(MBE))
            goto 100
         endif

         print 4
         print 20, chpath(i), MBE
         write (lunout,20) chpath(i), MBE
   20    format(5A)
 2000    call WAIT(chpath(i)//' have been processed')
      enddo
      close(lunout)

*--   .. Flush remaining buffers to disk
      call HREND('D2')
      stop
10000 print*, 'Could not open file ', fHIS
      stop
20000 print*, 'Could not create NEW file ', fMBE
      end

*                                         @METAGS WAIT
*                                         11-15-96 09:45pm
*--------------- WAIT ---------------
*
      SUBROUTINE WAIT(mess)
      character mess*(*), ch*1
      len = LENOCC(mess)
      if (len.GT.0) print*, mess(1:len)
      print*, '<CR>=Continue, Q=Quit'
      read 1, ch
      if ((ch.EQ.'q') .OR. (ch.EQ.'Q')) STOP
      RETURN
    1 FORMAT(A)
      END

*                                         @METAGS FPARSE
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

*                                         @METAGS LUNFREE
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

*                                         @METAGS FILELUN
*                                         04-07-98 01:01pm
*--------------- FILELUN ---------------
*
      SUBROUTINE FILELUN(lun,file)
      character*(*) file
      inquire (UNIT=lun, NAME=file)
      call CUTOL(file)
      END
