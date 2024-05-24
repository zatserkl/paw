*                                         @METAGS MESSLUN
*                                         10-08-98 12:32pm
*--------------- MESSLUN ---------------
*
      SUBROUTINE MESSLUN(lun)
      character*32 file
      file = ' '
      call FILELUN(lun,file)
      length = LENOCC(file)
      print*, 'file name length =', length
      if (length.GT.0) then
         print*,'With lun =', lun, ' is connected file ', file(1:length)
*         do i=1,length
*            ic = ICHAR(file(i:i))
*            print*, file(i:i), ' ic =', ic
*         enddo
      else
         print*,'With lun =', lun, ' is no file connected'
      endif
      END
