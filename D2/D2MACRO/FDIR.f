      program FDIR
      integer NWPAWC
      parameter (NWPAWC=30000)
      common /PAWC/ PAW(NWPAWC)

      parameter (MAXDIR=100)
      character*8 chpath(MAXDIR)

      character fname*32,file*32, exthis*3
      data exthis/'his'/

      npar = IARGC()
      if (npar.GT.0) then
         call GETARG(1,file)
      else
         print*, 'Enter histogram file name'
         read '(A)', file
      endif
      call FPARSE(file,fname,file,exthis)
      lenf = LENOCC(file)
      print*, 'Histogram file ', file(1:lenf)

*--   .. Tell HBOOK how many words are in PAWC
      call HLIMIT(NWPAWC)
      ntlun = LUNFREE(1)
      call HROPEN(ntlun,'D2',file(1:lenf),' ',1024,ISTAT)
      if (ISTAT.NE.0) goto 10000

*     .. store ndir directory names in the array chpath
      call HRDIR(MAXDIR,chpath,ndir)

      print*, 'List of directories in histogram file ', file
      do i=1,ndir
         print*, chpath(i)
      enddo

      print*, ' '
      print*, ' '
      print*, 'Complete directory tree'
      print*, '-----------------------'
      call HLDIR(' ','T')

*--   .. Flush remaining buffers to disk
      call HREND('D2')
      stop
10000 print*, 'Could not open file ', file(1:lenf)
      end
