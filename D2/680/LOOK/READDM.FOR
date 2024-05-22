      print*, 'Enter nz,nf'
      read*, nz,nf
      call READDM(nz,nf)
      end

      SUBROUTINE READDM(nz,nf)
      character file*32, fname*12
      data fname /'dm**__.dat'/

      write (fname(3:4),3) nf
      file = '\D2\DAT\'//fname
    3 FORMAT (I2.2)

      lun = LUNFREE(1)
      open (lun, FILE=file, STATUS='OLD',
     &      ACCESS='DIRECT', RECL=12, ERR=10000)
      read (lun,REC=nz, ERR=20000) Dmm,F,QMBE
      call FREELUN(lun,file)

      print*, 'Dmm,F,QMBE:', Dmm,F,QMBE
      RETURN
20000 call FREELUN(lun,file)
      print*, 'READ ERROR file ', file
      call WAIT('*** ERROR READDM ***')
      RETURN
10000 call FREELUN(lun,file)
      print*, 'OPEN ERROR file ', file
      call WAIT('*** ERROR READDM ***')
      RETURN
      END

*                                         @METAGS .LUNFREE
*                                         01-20-98 02:25pm
*--------------- LUNFREE ---------------
*
      integer function LUNFREE(lunstart)
      integer lunstart,lun
      logical used
      lun = lunstart
  100 inquire (UNIT=lun, OPENED=used, ERR=10000)
      if (used) then
         lun = lun+1
         if (lun.GT.99) lun=1
         goto 100
      endif
      LUNFREE = lun
      RETURN
10000 LUNFREE = 40
      call WAIT('ERROR LUNFREE')
      END

*                                         @METAGS .FREELUN
*                                         11-11-98 12:33pm
*--------------- FREELUN ---------------
*
      SUBROUTINE FREELUN(lun,file)
*
*     Closes DIRECT ACCESS file and/or correctly frees the unit
*     associated with DIRECT ACCESS file
*
      integer lun
      character*(*) file

*     .. Note: COMIS does not understand parameters in CLOSE others then UNIT
      close(UNIT=lun)
*
*     COMIS does not love DIRECT ACCESS files.
*     Unit lun is free now, but something like following
*     INQUIRE(UNIT=lun, ...) will crash the system.
*     Patch: associate the unit lun with SEQUENTIAL file and free it again.
*
      open (UNIT=lun, FILE=file, STATUS='OLD', ERR=10000)
      close(UNIT=lun)
10000 continue
      END

*                                         @METAGS WAIT.COMIS
*                                         11-15-96 09:45pm
*--------------- WAIT ---------------
*
      SUBROUTINE WAIT(mess)
      character mess*(*), ch*1
      length = LEN(mess)
      if (length.GT.0) print 1, mess(1:length)
      print*, '<CR>=Continue, Q=Quit'
      read 1, ch
      if ((ch.EQ.'q') .OR. (ch.EQ.'Q')) STOP
      RETURN
    1 FORMAT(A)
      END
