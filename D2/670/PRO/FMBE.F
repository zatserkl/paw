      program FMBE
*
*     Summary FMBE file name
*
      character fname*12
      parameter (fname='FMBE.dat')
*
*     Input information
*
      character Energy*4, theta*5
      parameter (Energy='550', theta='143.5')
*     .. Histogram file
      character hfile*12
      parameter (hfile='h550.his')
*
*--
*
      integer NWPAWC
      parameter (NWPAWC=30000)
      common /PAWC/ PAW(NWPAWC)

      parameter (MAXDIR=100)
      character*8 chpath(MAXDIR)
      character*8 chtype,chopt, chtitl*60
      
      logical fexist
      character*8 MBE(4), MBEtmp
      character MBEs*4, MBEs1*4, MBEs2*5
      data MBEs1/'MBE='/, MBEs2/'MBE ='/

*--   .. Tell HBOOK how many words are in PAWC
      call HLIMIT(NWPAWC)
      call HROPEN(1,'D2',hfile,' ',1024,ISTAT)
      if (ISTAT.NE.0) goto 10000

      print*, ' '
      print*, 'Energy = ', Energy, ' MeV, theta = ', theta, ' degree'
      print*, 'Histogram file ', hfile
      call WAIT('To change correct the source file and recompile')

      inquire (FILE=fname, EXIST=fexist)
      if (fexist) then
         print*, 'File ', fname, ' already exist!'
         STOP
      endif

      open (2, FILE=fname, STATUS='NEW', ERR=20000)
      write(2,1)
      write(2,2) Energy,theta
      write(2,3)
      write(2,4)
      write(2,5)
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
      call WAIT('List of directories in histogram file '//hfile)

      do i=1,ndir
         print*, chpath(i)
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
         write (2,20) chpath(i), MBE
   20    format(5A)
         call WAIT(chpath(i)//' have been processed')
      enddo
      close(2)

*--   .. Flush remaining buffers to disk
      call HREND('D2')
      stop
10000 print*, 'Could not open file ', hfile
      stop
20000 print*, 'Could not create NEW file ', fname
      end

*                                         @METAGS WAIT
*                                         11-15-96 09:45pm
*--------------- WAIT ---------------
*
      SUBROUTINE WAIT(mess)
      character mess*(*), ch*1
      if ((LEN(mess).GT.0) .AND. (mess.NE.' ')) print*, mess
      print*, '<CR>=Continue, Q=Quit'
      read 1, ch
      if ((ch.EQ.'q') .OR. (ch.EQ.'Q')) STOP
      RETURN
    1 FORMAT(A)
      END
