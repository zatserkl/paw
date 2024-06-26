      REAL FUNCTION fun(XDUMMY)
*********************************************************
*                                                       *
* This file was generated by HUWFUN.                    *
*                                                       *
*********************************************************
*
*     N-tuple Id:     1    
*     N-tuple Title:  Spherical mirror + plane reflector
*     Creation:       11/09/98 02.25.08
*
*********************************************************
*
      COMMON /PAWIDN/ IDNEVT,VIDN1,VIDN2,VIDN3,VIDN(10)
*
*--   Ntuple Variable Declarations
*
      REAL part(6),ray0(6),tst(6),ray(6),raycap(6),raypmt(6)
      INTEGER iref,ip,ireg,intref
*
      COMMON /PAWCR4/ iref,ip,part,ray0,tst,ray,ireg,intref,raycap
     + ,raypmt
*
      character*80 tit1,tit2
      character*8 var(4)
*
*--   Enter user code here
*
*      vector phot(1000,1000,1000)
*
*     Read only the three desired columns
*
      var(1) = 'ip'
      var(2) = 'part'
      var(3) = 'ireg'
      call HGNTV(1,var,3,1,ier)
      if (ier.NE.0) then
         print*, 'Error in HGNTV'
      endif
      
      idh1 = 100
      idh2 = 200
      print*, 'The number of entries is', IDNEVT
      print*, 'VIDN1,VIDN2,VIDN3:', VIDN1,VIDN2,VIDN3

      tit1 = 'Photon vs position'
      call HBOOK2(idh1,tit1,60.,-30.,30.,60.,-30.,30.,0.)
      tit2 = 'Registrated photon vs position'
      call HBOOK2(idh2,tit2,60.,-30.,30.,60.,-30.,30.,0.)

      np = 0
      nemit = 0
      nreg  = 0
      x = 100.
      y = 100.
      do 100 nevent=1,IDNEVT
         if (nevent.NE.1) call HGNTF(1,nevent,ier)
         if (ier.NE.0) then
            print*, 'Error reading row', nevent
         endif
         if (ip.EQ.1) goto 100
         if ((part(5).EQ.x).AND.(part(6).EQ.y)) then
            nemit = nemit+1
            if (ireg.EQ.1) nreg = nreg+1
         else
            if (np.GT.0) then
*              .. fill particle
               call HFILL(idh1,x,y,REAL(nemit))
               call HFILL(idh2,x,y,REAL(nreg))
            endif
            np = np+1
            x = part(5)
            y = part(6)
            nemit = 0
            nreg  = 0
         endif
*         print*, 'in cycle np =', np
  100 enddo
*     .. fill last particle
      call HFILL(idh1,x,y,REAL(nemit))
      call HFILL(idh2,x,y,REAL(nreg))

      print*, 'np =', np
      fun = 1.
*
      END

