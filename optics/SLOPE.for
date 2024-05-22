      REAL FUNCTION SLOPE(XDUMMY)
      include 'include.inc'
*

*
*--   Enter user code here
*
      id = 1
      print*, 'Look at OPTIC'
      ievent = 1
      call HGNTB(id,'OPTIC',ievent,ier)
      if (ier.NE.0) then
         print*, 'Error reading OPTIC' 
      else
         print*, 'ip,ireg,npart:'
         print*, ip,ireg,npart
         print*, 'raymir:'
         print*, raymir
      endif

      print*, 'Look at PARAM'
      print*, 'Let assign e.g. Rbox=12345. and read Rbox from Ntuple'
      Rbox = 12345.
      ievent = 1
      call HGNTB(id,'PARAM',ievent,ier)
      if (ier.NE.0) then
         print*, 'Error reading PARAM' 
      else
         print*, 'Dbox, Rbox, Dbox, f, Rpmt, Dcol, Rcol2, dhcol:'
         print*, Dbox, Rbox, Dbox, f, Rpmt, Dcol, Rcol2, dhcol
      endif

      SLOPE = 1.
*
      END

