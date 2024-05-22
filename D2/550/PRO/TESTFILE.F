      program test
      integer NWPAWC
      parameter (NWPAWC=15000)
      common /PAWC/ PAW(NWPAWC)
      character*8 chtags(5)
      dimension event(5)
      equivalence (event(1),X), (event(2),Y), (event(3),Z)
      equivalence (event(4),energy), (event(5),Eloss)
      data chtags/'X', 'Y', 'Z', 'Energy', 'Eloss'/
*--   .. Tell HBOOK how many words are in PAWC
      call HLIMIT(NWPAWC)
      call HROPEN(1,'EXAMPLE','EXAMPLE.DAT','N',1024,ISTAT)
      if (ISTAT.NE.0) goto 99
*--   .. Make sub-directory on disk (as HROUT does not do this for us)
      call HMDIR('US','S')
      call HCDIR('//PAWC',' ')
*--   .. Make sub-directory in memory
      call HMDIR('US','S')
      call HCDIR('//PAWC/US',' ')
*--   .. Book Ntuple + 1d histogram
      call HBOOKN(10,'A simple Ntuple',5,'EXAMPLE',5000,chtags)
      call HBOOK1(100,'Energy distribution',100,0.,100.,0.)
*--   .. Fill the Ntuple and histogram
      do 10 i=1,1000
         call RANNOR(X,Y)
         Z = sqrt(X*X + Y*Y)
         energy = 50. + 10.*X
         Eloss = 10.*abs(Y)
         call HFN(10,event)
         call HFILL(100,energy,0.,1.)
   10 continue
*--   .. Juggle top directories (order of this calls is important!)
      call HCDIR('//PAWC',' ')
      call HCDIR('//EXAMPLE',' ')
*--   .. Write out everything to disk
      call HROUT(0,ICYCLE,'T')
*--   .. Flush remaining buffers to disk
      call HREND('EXAMPLE')

   99 continue
      end
