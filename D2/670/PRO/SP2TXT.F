c                                      @METAGS .MAIN.
c============ .MAIN. ==============    
c
      integer*2 idata(400)
      integer*2 ps32(32), vp(16), ncount(96), ndelay(96), icomm(70)
      equivalence (idata,ps32), (idata(65),vp),
     +   (idata(81),ncount), (idata(177),icomm), (idata(257),ndelay)
      character comm*70, target*3
      equivalence (icomm,comm), (target,vp)
      character fname*12, fnumb*2, nzap*2
      character mess*8
      equivalence (fname(3:3), fnumb)
      data fname /'sp**__.dat'/
c
      print*, 'Enter # file'
      read *, i
      write (fnumb,10) i
      print*, fname
      open (1, file='\\D2\\DAT\\'//fname, status='OLD',
     *      access='DIRECT', RECL=800, ERR=10000)
c
 1000 print*, '*****'
      print*, 'Data file ', fname, '. Enter #bin. 0 - exit'
      read *, irec
      if (irec.EQ.0) then
         close (1)
         STOP
      endif
c
      read (1,REC=irec, ERR=20000) idata
c
      write (nzap,10) irec
*      open (2, file='SPECTRA\\'//nzap//'.'//fnumb, status='UNKNOWN')
*      open (3, file='SPECTRA\\'//nzap//'i.'//fnumb, status='UNKNOWN')
      open (2, file=nzap//'.'//fnumb, status='UNKNOWN')
      open (3, file=nzap//'i.'//fnumb, status='UNKNOWN')
c
      write (2,20) (ncount(i), ndelay(i), i=1,96)
      close (2)
c
      print*, '---'
      MBE = makei4(ps32(19),ps32(20))
      print*, 'MBE =', MBE
      ips1 = makei4(ps32(21),ps32(22))
      print*, 'ips1 =', ips1
      itime = makei4(ps32(23),ps32(24))
      print*, 'itime =', itime
      iPGS = makei4(ps32(25),ps32(26))
      print*, 'iPGS =', iPGS
      LIVN = makei4(ps32(27),ps32(28))
      print*, 'LIVN =', LIVN
      NS = makei4(ps32(29),ps32(30))
      print*, 'NS =', NS
      print*, 'COMMENT:', comm
      print*, 'VP data:'
      print*, 'target ', target
      print*, 'Nnucl/10**20 =', VP(4)
      print*, 'Frequency =', VP(5), ' kHz'
      print*, '# integrater', vp(6)
      MBEvp = makei4(VP(7),VP(8))
      print*, 'MBEvp =', MBEvp
      MBEvp1 = vp(7)*2**15 + vp(8)
      print*, 'MBEvp1 = vp(7)*2**15 + vp(8):', MBEvp1

      print*, '# coll', vp(9)
      print*, 'Angle =', VP(10)
      print*, 'Teta to beam', vp(11)
      print*, 'E0 =', vp(12)
c
      mess = 'MBE'
      write(3,30) mess, MBE
      mess = 'MBEvp1'
      write(3,30) mess, MBEvp1
      mess = 'TIME'
      write(3,30) mess, itime
      mess = 'GAS'
      write(3,30) mess, iPGS
      mess = 'LIVN'
      write(3,30) mess, LIVN
      mess = 'NS'
      write(3,30) mess, NS
      mess = 'COMMENT:'
      write(3,40) mess, comm
      close(3)
c
      goto 1000
c
      STOP
20000 print*, 'READ ERROR'
      goto 1000
10000 print*, 'OPEN ERROR'
      STOP
   10 FORMAT (I2.2)
   20 FORMAT (2I8)
   30 FORMAT (A8,I8)
   40 FORMAT (A8,A70)
      end
c                                      @METAGS MAKEI4
c------------ MAKEI4 --------------    
c
      integer function makei4(i2old,i2low)
      integer*2 i2old,i2low, i2(2)
      integer*4 i4
      equivalence (i2,i4)
c     .. exchange old <--> low for I4, because I4 = (I2low,I2old)
      i2(1) = i2low
      i2(2) = i2old
      makei4 = i4
      end
