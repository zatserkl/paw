* 3000        50000   nran    *** WARNING: FIRSTS LINES USED AS INI DATA! 
* 1        2       init INT
* 30.      Rbox, cm
* 500.     Dbox, cm
* 2.5      Rpmt, cm
* 34.9     rindex: index of refraction (n-1)*1e6
* 90.      PMconst
* 0.20     PMconv 
* 10.      angmax, mrad
*
      program right
*      IMPLICIT NONE
      PARAMETER (NWPAWC = 300000)
      PARAMETER (LRECL  = 1024)
      PARAMETER (id = 1)
      COMMON /PAWC/ IPAW(NWPAWC)

      parameter (MAXRAY=100)
      real rayhis(6,MAXRAY+1)
      
      real sysmir(5),sysref(5),syscol(5),syspmt(5)
      character*127 CHPATH

      character CHTOP*8, CHFILE*32
      character*80 chtitl
      
      character c
*
*     Store in Ntuple:
*     iref     - flag hit to reflector
*     ireg     - flag hit to PM
*     part(6)  - vector of incident particle
*     ray0(6)  - vector of initial photon
*     ray(6)   - vector of photon after the reflector
*     xypmt(2) - x,y coordinates of intersection with PM in PM system
*
      COMMON /OPTIC/ iref,ireg,part(6),ray0(6),ray(6),xypmt(2)
      real raylab(6),rayloc(6),rayref(6)

      character*32 thisf,thisn,ifile,bfile, afile

      data CHTOP /'NTDIR'/
      data CHFILE/' '/
      data afile /'unknown'/

      call GETARG(0,thisf)
      call FPARSE(thisf,thisn,thisf,'x')
      lename = LENOCC(thisn)
      CHFILE = thisn(1:lename)//'.hbook'
*      ifile  = thisn(1:lename)//'.ini'
      ifile  = thisn(1:lename)//'.f'
      bfile  = thisn(1:lename)//'_t.for'
*
*     Ntuple
*      
      write(chtitl,10)
   10 format ('Spherical mirror + plane reflector')
      call MESS('Ntuple title:')
      call MESS(chtitl)

      CALL HLIMIT(NWPAWC)
*     .. open a new RZ file
      ntlun = LUNFREE(1)
      CALL HROPEN(ntlun,CHTOP,CHFILE,'N',LRECL,ISTAT)
*     .. book Ntuple
      CALL HBNT(id,chtitl,' ')
*     .. define Ntuple
      CALL HBNAME(id,'OPTIC', iref,
     &'iref[0,1]:U,ireg[0,1]:U,part(6):R,ray0(6):R,ray(6):R,xypmt(2):R')

      print*, ' '
      print*, 'This file ', thisf
      print*, 'Ini file  ', ifile
      print*, 'HBOOK file with Ntuple id =', id, '   ', CHFILE
      print*, 'Batch file for PAW usage     ', bfile
      print*, ' '

      lun = LUNFREE(1)
      open (lun, FILE=ifile, STATUS='OLD', ERR=10000)
      read (lun,*) c,nran
      read (lun,*) c,init
      read (lun,*) c,Rbox
      read (lun,*) c,Dbox
      read (lun,*) c,Rpmt
      read (lun,*) c,rindex
      read (lun,*) c,PMconst
      read (lun,*) c,PMconv
      read (lun,*) c,angmax

      close(lun)

      print*, 'Ini file data:'
      print*, 'N randoms/bin =', nran
      print*, 'Rbox  =', Rbox, ' cm'
      print*, 'Dbox =', Dbox, ' cm'
      print*, 'Rpmt =', Rpmt, ' cm'
      print*, 'rindex =', rindex, ' (n-1)*1e6'
      print*, 'PMconst =', PMconst, ' 1/cm'
      print*, 'PMconv =', PMconv
      print*, 'angmax  =', angmax, ' mrad'

c     .. convert
      pi = ACOS(-1.)
      torad = pi/180.
      todeg = 180./pi
      
      angrad = angmax/1000.
      
*
*     Geometry
*      
*
*     x|                dzpmt
*      |                 +-+
*      |                 | |
*      |                 | |
*      |  /              +-+ dxpmt
*      | /               
*      |/    
*      |                   /
*     0+------------------/----------
*     /|                 /dzpla     z
*    / |\                   
*   /  | \    <---------*---
* y/   |  \  ray(1:3) ray(4:6)
*

*     .. mirror with R = 2f, where f=Dbox
      R = 2.*Dbox + 1.
*     .. mirror frame
      iaxis = 3
      thmir = 0.
      zmir = 0. - 1.
      xmir = 0.
      call mkframe(iaxis,thmir,zmir,xmir, sysmir)

*     .. Reflector
      Rref = 10.
*     .. reflector frame
      iaxis = 3
      thref = 135.
      zref = Dbox
      xref = 0.
      call mkframe(iaxis,thref,zref,xref, sysref)

*     .. paraboloid collector
      Rpar1 = Rpmt
      Rpar2 = 12.
      call PARABpar(Rpar1,Rpar2,fcol,hcol)
      h1par = fcol
      h2par = hcol
*     .. par frame
      iaxis = 3
      thpar = -90.
      zpar = Dbox
      xpar = Rbox + h2par
      call mkframe(iaxis,thpar,zpar,xpar, syscol)

*     .. PM
      x1pmt = -Rpmt
      x2pmt =  Rpmt
*     .. pmt frame
      iaxis = 3
      thpmt = -90.
      zpmt = Dbox
      xpmt = Rbox + (h2par-h1par)
      call mkframe(iaxis,thpmt,zpmt,xpmt, syspmt)
      
      print*, ' '
      print*, 'Distance to PM is', xpmt, ' cm from counter axis'

*     .. Cherenkov angle
      thetac = ACOS(1.D0/(1.D0 + rindex*1.D-6))
*     .. the number of photoelectrons
      pe = Dbox*PMconst*(SIN(thetac))**2
*     .. the number of photons
      phot = pe/PMconv
*     .. mean free path length between emitances
      alambda = Dbox/phot
      print*, ' '
      print*, 'Cherenkov angle', thetac, ' rad'
      print*, 'The average number of photoelectrons', pe
      print*, 'The average number of photons', phot
      print*, 'Mean free path length', alambda

      dummy = BAGINT(init)
      print*, ' '

      call FILELUN(ntlun,afile)
      print*, 'With unit', ntlun, ' connected file ', afile
      print*, ' '
      call WAIT('Begin the calculation')
      print*, 'Running ..'

      tlim = 10000.
      call TIMEST(tlim)
      
      do i=1,6
      	 ray(i) = 0.
      enddo
      
*
*--   MAIN PART
*
      npmt = 0
            
      idummy = 0
      nemit = 0
      nref  = 0
*------------------------------------- main loop begin
      DO np=1,nran
*     .. particle radiant in circle with radius Rbox
  100 continue
      xp = Rbox*2.*(.5 - RANBAG(idummy))
      idummy = idummy+1
      yp = Rbox*2.*(.5 - RANBAG(idummy))
      if (xp**2 + yp**2 .GT. Rbox**2) goto 100
      zp = Dbox
      idummy = idummy+1
*     .. particle direction
      idummy = idummy+1
      thp = pi - angrad*RANBAG(idummy)
      idummy = idummy+1
      php = 2.*pi*RANBAG(idummy)

      part(4) = zp
      part(5) = xp
      part(6) = yp
      part(1) = cos(thp)
      part(2) = sin(thp)*cos(php)
      part(3) = sin(thp)*sin(php)

*     .. consts for coordinates transformations
      sinthp = sin(thp)
      costhp = cos(thp)
      sinphp = sin(php)
      cosphp = cos(php)
      
*     .. initialize particle current position
      z = zp
      x = xp
      y = yp

 1000 continue      
*     .. current range
      idummy = idummy+1
      range = -alambda*LOG(RANBAG(idummy))
*     .. current position
      z = z + range*part(1)
      x = x + range*part(2)
      y = y + range*part(3)
      if (z.LE.0.) goto 1999
      
*     .. emit the photon in the particle frame
      nemit = nemit+1
*     .. clear reflector flag
      iref = 0      
*     .. clear registration flag
      ireg = 0

*     .. photon radiant
      ray(4) = z
      ray(5) = x
      ray(6) = y      
*     .. photon direction
      thg = thetac
      idummy = idummy+1
      phg = 2.*pi*RANBAG(idummy)

*     .. photon dir.cosines in particle system
      cz = cos(thg)
      cx = sin(thg)*cos(phg)
      cy = sin(thg)*sin(phg)	  
*      print*, 'photon dir.cos. in particle sys.:', cz,cx,cy

*     .. convert photon dir.cosines to lab
*     	 step 1: rotate round y by thp
*     	 step 2: rotate round z by php
*     	 Result:
      ray(1) = cz*costhp        - cx*sinthp
      ray(2) = cz*sinthp*cosphp + cx*costhp*cosphp - cy*sinphp
      ray(3) = cz*sinthp*sinphp + cx*costhp*sinphp + cy*cosphp
*     .. store initial photon vector in ray0
      do i=1,6
      	 ray0(i) = ray(i)
      enddo
*
*     Ray runs
*
*     .. reflection from mirror
      call toloc(sysmir,ray,rayloc)
      iermir = MISPH(R,rayloc,rayref)
*     .. return to lab
      call tolab(sysmir,rayref,raylab)
      if (iermir.GT.0) then
*        .. absorbed at the mirror.
         goto 2000
      endif

*     .. reflection from reflector
      call toloc(sysref,raylab,rayloc)
      ierref = MIRCIR(Rref,rayloc,rayref)
      if (ierref.NE.0) then
*        .. lost ray
      	 goto 2000
      endif
      call tolab(sysref,rayref,raylab)
      
      nref = nref+1
*     .. up reflector flag      
      iref = 1
      do i=1,6
      	 ray(i) = raylab(i)
      enddo
      
*     .. clear registration flag      
      ireg = 0
      nray = 1
      do i=1,6
      	 rayhis(i,nray) = ray(i)
      enddo
      
*--   Reflection from paraboloid collector
*     .. ray in collector local frame
      call toloc(syscol,rayhis(1,nray),rayhis(1,nray))
*      print*, 'In loc: rayhis, nray =', nray
*      print*, (rayhis(j,nray), j=1,6)
      irest = MAXRAY-nray
      call PARAB(fcol,h1par,h2par,irest,rayhis(1,nray),npro,iercol)
*     .. return to lab
      do i=nray,nray+npro
         call tolab(syscol,rayhis(1,i),rayhis(1,i))
      enddo
      nray = nray+npro
      if (iercol.GT.0) then
*        .. absorbed at the collector. Next photon
         goto 2000
      endif

*--   Intersection with PM
*     .. ray in PM local frame
      call toloc(syspmt,rayhis(1,nray),rayloc)
      call PMT(Rpmt,rayloc, rayref,ierpmt)
      if (ierpmt.GE.0) then
*        .. return to lab output ray
         call tolab(syspmt,rayref,raylab)
   
         if (ierpmt.EQ.0) then
*           .. absorbed at the PMT.
      	    ireg = 1
      	    npmt = npmt+1
*     	    .. coordinates in PM system      	    
      	    xypmt(1) = rayref(5)
      	    xypmt(2) = rayref(6)
      	 else
*     	    .. lost
      	    goto 2000
         endif
      endif
*     .. goto store      
      goto 3000
      
 2000 continue
*     .. clear hit flag and ray
      iref = 0
      ireg = 0
      do i=1,6
      	 ray(i) = 0.
      enddo
      xypmt(1) = 10.*Rpmt
      xypmt(2) = 10.*Rpmt

 3000 continue
*     .. store raw
      call HFNT(id)
      
*     .. for the next photon
      goto 1000
 1999 continue
      if (MOD(np,1000).EQ.0) then
      	 print*, 'processed', np, ' particles from', nran
      endif
      ENDDO
*------------------------------------- main loop end
*
*--   MAIN PART END
*

      print*, ' '
      print*, 'There are', nemit, ' photons were emitted'
      print*, 'There are', nref, ' photons were reflected'
      phpart = REAL(nemit)/REAL(nran)
      print*, 'The average number of emitted photons/particle', phpart

      eff = 100.*REAL(npmt)/REAL(nemit)
      print*, ' '
      print*, 'Hit PM', npmt, ' photons from', nref, ' reflected'
      print*, 'Hit PM', npmt, ' photons from', nemit, ' emitted'
      print*, 'Registration efficiency is', eff, '%'
      print*, '--------------------------------------'

*      dummy = RNOREP(nhit,miss)
      dummy = BAGREP(iused)
      
      print*, ' '
      call TIMEX(time)
      print*,'Job time', time, ' seconds'
      print*, ' '

      CALL HPRNT(id)
*
*-- write batch version of analisys routine to file
*
      lunb = LUNFREE(1)
      OPEN(lunb, FILE=bfile, STATUS='UNKNOWN')
      CALL HUWFUN(lunb, id, thisn(1:lename), 0, 'B')
*
*-- write Ntuple buffer to disk and close RZ file
*
      CALL HROUT(0,ICYCLE,' ')
*      print*, 'After HROUT(0,ICYCLE,'' '') ICYCLE =', ICYCLE

      CHPATH = ' '
      call HCDIR(CHPATH,'R')
      print*, ' '
      call HLDIR('//PAWC','T')
      print*, ' '
      print '('' Current directory is '', A)', CHPATH(1:LENOCC(CHPATH))
      call HLDIR(' ','T')

      CALL HREND(CHTOP)
      STOP
10000 print*, 'Ini file ', ifile(1:LENOCC(ifile)), ' did not found'
      END

*                                         @METAGS LOOK
*                                         03-26-98 01:30am
*--------------- LOOK ---------------
*
      SUBROUTINE LOOK(id)
      REAL F,E1,E2,THin,THout,Xin,Xout
      INTEGER pdf,PChit,PCch,PCcl, Thit,Tch,Tcl, err
*      CHARACTER distr*4
      COMMON/NTDATA/ F,E1,E2,THin,THout,Xin,Xout,
     &               pdf,PChit,PCch,PCcl,Thit,Tch,Tcl, err
*      COMMON /NTDATAC/ distr

 1000 I = 1
      print*, 'Enter a number of event to look (<CR>=Exit)'
      read 1, I
    1 format(I8)
      if (I.EQ.0) goto 2000
      CALL HGNT(id, I, IER)
      IF (IER .NE. 0) THEN
         PRINT *, 'Error reading row ', I
      ENDIF
      PRINT *,'Event',I,' E1,E2,PCch,pdf:',E1,E2,PCch,pdf
      goto 1000
 2000 continue
      END

*                                         @METAGS mkray
*                                         03/07/98 19:15
*--------------- mkray ---------------
*
      SUBROUTINE mkray(z,x,y,dthmr,phrad,ray)
      real z,x,y,dthmr,phrad,ray(6)
      double precision pi,thr,phr
      pi = ACOS(-1.d0)
      thr = pi/2.d0 - .001*dthmr
*      phr = pi      - .001*dphmr
      phr = phrad
*     .. spherical --> Cartesian
      ray(1) = sin(thr)*cos(phr)
      ray(2) = cos(thr)
      ray(3) = sin(thr)*sin(phr)
      ray(4) = z
      ray(5) = x
      ray(6) = y
      END

*                                         @METAGS mkframe
*                                         03/07/98 10:11
*--------------- mkframe ---------------
*
      SUBROUTINE mkframe(iaxis,thdeg,d1,d2,frame)
      real thdeg,d1,d2,frame(5)
      double precision pi,torad,thrad
      pi = ACOS(-1.d0)
      torad = pi/180.d0
      thrad = thdeg*torad
      frame(1) = iaxis
      frame(2) = COS(thrad)
      frame(3) = SIN(thrad)
      frame(4) = d1
      frame(5) = d2
      END

*                                         @METAGS toloc
*                                         16/06/98 09:40
*--------------- toloc ---------------
*
      SUBROUTINE toloc(frame,raylab, rayloc)
*
*     Rotate the right coord. system round the INT(axis3) axis.
*     frame = (axis3,cosA,sinA,d1,d2)
*     d1 - shift along the axis1
*     d2 - shift along the axis2
*
      real frame(5),raylab(6),rayloc(6),c(6),x(6)
      i3 = frame(1)
      if (i3.EQ.3) then
         i1 = 1
         i2 = 2
      elseif (i3.EQ.1) then
         i1 = 2
         i2 = 3
      elseif (i3.EQ.2) then
         i1 = 3
         i2 = 1
      else
         print*, 'ERROR toloc: Wrong axis number in frame(1):', frame(1)
         STOP
      endif

      do i=1,3
         c(i) = raylab(i)
         x(i) = raylab(3+i)
      enddo
*     .. directrix cosines
      rayloc(i1) =  c(i1)*frame(2) + c(i2)*frame(3)
      rayloc(i2) = -c(i1)*frame(3) + c(i2)*frame(2)
      rayloc(i3) =  c(i3)
*     .. point
      rayloc(3+i1)=  (x(i1)-frame(4))*frame(2)+(x(i2)-frame(5))*frame(3)
      rayloc(3+i2)= -(x(i1)-frame(4))*frame(3)+(x(i2)-frame(5))*frame(2)
      rayloc(3+i3)=  x(i3)
      END

*                                         @METAGS tolab
*                                         16/06/98 09:40
*--------------- tolab ---------------
*
      SUBROUTINE tolab(frame,rayloc, raylab)
*
*     Rotate the right coord. system round the INT(axis3) axis.
*     frame = (axis3,cosA,sinA,d1,d2)
*     d1 - shift along the axis1
*     d2 - shift along the axis2
*
      real frame(5),rayloc(6),raylab(6),c(6),x(6)
      i3 = frame(1)
      if (i3.EQ.3) then
         i1 = 1
         i2 = 2
      elseif (i3.EQ.1) then
         i1 = 2
         i2 = 3
      elseif (i3.EQ.2) then
         i1 = 3
         i2 = 1
      else
         print*, 'ERROR tolab: Wrong axis number in frame(1):', frame(1)
         STOP
      endif

      do i=1,3
         c(i) = rayloc(i)
         x(i) = rayloc(3+i)
      enddo
*     .. directrix cosines
      raylab(i1) = c(i1)*frame(2) - c(i2)*frame(3)
      raylab(i2) = c(i1)*frame(3) + c(i2)*frame(2)
      raylab(i3) = c(i3)
*     .. point
      raylab(3+i1) = x(i1)*frame(2) - x(i2)*frame(3) + frame(4)
      raylab(3+i2) = x(i1)*frame(3) + x(i2)*frame(2) + frame(5)
      raylab(3+i3) = x(i3)
      END

*                                         @METAGS raycop
*                                         08/06/98 20:16
*--------------- raycop ---------------
*
      SUBROUTINE raycop(ray1,ray2)
      real ray1(6),ray2(6)
      do i=1,6
         ray2(i) = ray1(i)
      enddo
      END

*                                         @METAGS DROOT
*                                         06-02-98 03:47am
*--------------- DROOT ---------------
*
      subroutine DROOT(a,b,c,x1,x2,ierr)
*
*     For square equation a*x**2 + b*x + c = 0
*     calculates REAL roots.F
*     For COMPLEX returns in x1 real part and in x2 imaginary part
*     of root -b/(2*a) + i*sqrt(-D)/(2*a)
*     To preserve computer accuracy use formula x1 = 2*c/(-b-sqrt(D))
*     Error flags:
*      0: normal
*      1: any roots: a=b=c=0
*      2: complex roots
*     -1: no roots: a=b=0, c<>0
*
      double precision a,b,c,x1,x2
      integer ierr
      double precision D

      x1 = 0.d0
      x2 = 0.d0
      ierr = 0

*      print*, a,b,c
*      call WAIT('DROOT: a,b,c')
      if (a .EQ. 0d0) then
         if (b.EQ.0d0) then
            if (c.EQ.0.d0) then
*               print*, 'DROOT: Any roots: a=b=c=0'
               ierr = 1
            else
               print*, 'DROOT ERROR: No roots: a=b=0, c<>0'
               ierr = -1
            endif
            RETURN
         endif
         x1 = -c/b
         x2 = x1
         RETURN
      endif

      D = b**2 - 4d0*a*c
*      print*, 'D =', D
*      call WAIT('DROOT: D')
      if (D .LT. 0.d0) then
*         print*, 'DROOT: Complex roots'
         x1 = -b/(2.d0*a)
         x2 = sqrt(-D)/(2.d0*a)
         ierr = 2
         RETURN
      elseif (D .EQ. 0.d0) then
         x1 = -b/(2.d0*a)
         x2 = x1
         RETURN
      endif

      if (b .GT. 0.d0) then
*        .. compute "minus" root
         x1 = 2.d0*c/(-b - sqrt(D))
      else
*        .. compute "plus" root
         x1 = 2.d0*c/(-b + sqrt(D))
      endif

*     .. use Viete formulae for the other root
      if (x1 .NE. 0.d0) then
*        .. to preserve computer accuracy
         x2 =  (c/a) / x1
      else
         x2 = -(b/a) - x1
      endif
      END

*                                         @METAGS MISPH
*                                         06-02-98 01:17pm
*--------------- MISPH ---------------
*
      integer function MISPH(R,ray,ray1)
*
*     Minimized routine: only one reflection
*
*     Reflection from spherical mirror  x**2 + y**2 + (z-R)**2 = R**2
*     Focus f=R/2
*     To simplify calculation we use transfer to center of sphera system
*
*     Note: ray,ray1 in frame of spheroid
*           (i.e. normal in sphera top is (1,0,0))
*
*  Intersection point ray1(4:6) of spheroid x**2 + y**2 + (z-R)**2 = R**2
*  and ray with radiant ray(4:6) and directrix cosines ray(1:3)
*  Notation: ray(cos1,cos2,cos3,x1,x2,x3) == ray(cosz,cosx,cosy,z,x,y)
*
*     Returns error code:
*      0: normal reflection
*      1: ray is absorbed at the back side of sphera
*     -1: ray have missed the sphera
*
*              x|          "Mirror" system
*               |      .   ---------------
*               |     /|
*           ray1(4:6)/ |   ray(1:3) ray(4:6)
*               |   *-------<--------*-----
*               |  /   |
*               | |    |
*              0| |zmin|zmax  
*            ---+-+----+-----*-----------------------
*              /| |    |      R                     z
*             / | |    |
*            /  |  \   |           
*           /   |   \  |
*          /    |    \ |
*        y/     |     \|
*               |      `
*
      IMPLICIT NONE
      real R,ray(6),ray1(6)
      double precision R2
      double precision c1,c2,c3,x1,x2,x3
      double precision a,b,c,t,t1,t2
      double precision o1,o2,o3,oo,oray
      double precision wavel2,vlen2
*     .. square of wave length (in cm) used as an epsilon
      parameter (wavel2 = 4.d-10)
      integer i,idroot

      do i=1,6
         ray1(i) = ray(i)
      enddo

      R2 = DBLE(R)**2

*     .. cosines of incident ray
      c1 = ray(1)
      c2 = ray(2)
      c3 = ray(3)
*     .. radiant of incident ray in center of sphera system (shift for x1)
      x1 = DBLE(ray(4)) - DBLE(R)
      x2 = ray(5)
      x3 = ray(6)
*      print*, 'rayhis(4,nray) in "mirror" system =', rayhis(4,nray)
*      print*, 'x1 in "sphera" system =', x1

*
*     Rewrite line equation (z-z0)/cosz=(x-x0)/cosx=(y-y0)/cosz
*     in parametric form:
*     z = z0 + cosz*t
*     x = x0 + cosx*t
*     y = y0 + cosy*t
*     here t - parameter. t=0 corresponds to radiant,
*     t > 0 corresponds going downstream the directrix vector
*     t < 0 - going in backward direction.
*     To find t solve the equation x**2 + y**2 + z**2 = R**2
*--
*     Perform the calculation in center of sphera system.
*     Here x**2 + y**2 + z**2 = R**2
*--
*     .. coeffs. of square equation
      a = c1**2 + c2**2 + c3**2
      b = 2.*(c1*x1 + c2*x2 + c3*x3)
      c = x1**2 + x2**2 + x3**2 - R2
      call DROOT(a,b,c,t1,t2,idroot)
*      print*, 'roots:', t1,t2
      if (idroot.NE.0) then
*         print*, 'MISPH: idroot =', idroot
*         call WAIT('MISPH: No roots of square equation')
         goto 1001

*         if (idroot.LT.0) call MESS('MISPH: No roots')
*         if (idroot.EQ.1) call MESS('MISPH: Any roots')
*         if (idroot.EQ.2) call MESS('MISPH: Complex roots')
*         goto 1001
      endif

*     .. start test from MIN root becouse both roots may be positive
      t = MIN(t1,t2)
*      print*, 'MIN root = ', t
*     .. test vector length before the sign: if t ~ 0 may be any sign     
      vlen2 = (c1**2 + c2**2 + c3**2)*t**2
      if (vlen2.LT.wavel2) then
*         print*, '--- vlen too small. Change to MAX root ---'
*         print*, 'vlen2 =', vlen2
         goto 100
      endif
      if (t .LT. 0.d0) then
*         print*, '--- Negative MIN root. Change to MAX root ---'
         goto 100
      endif
      goto 200

  100 continue
      t = MAX(t1,t2)
*      print*, 'MAX root =', t
*     .. test vector length before the sign: if t ~ 0 may be any sign     
      vlen2 = (c1**2 + c2**2 + c3**2)*t**2
      if (vlen2.LT.wavel2) then
*         print*, '--- MAX vlen too small. Exit ---'
*         print*, 'vlen2 =', vlen2
         goto 1001
      endif
      if (t .LT. 0.d0) then
*         print*, '--- Negative MAX root. Exit ---'
         goto 1001
      endif
  200 continue

*      print*, 't1,t2,t:'
*      print*, t1,t2,t
*      print*, 'eq.:', a*t**2 + b*t + c

*     .. adjust point on spheroid
      x2 = x2 + c2*t
      x3 = x3 + c3*t
      x1 = -1.d0*sqrt(R2 - x2**2 - x3**2)
*      print*, 'z in "sphera" system =', z
*     .. radiant of the output ray in mirror system (shift for z)
      ray1(4) = x1 + R
      ray1(5) = x2
      ray1(6) = x3
*
*  Gets directrix cosines o of nornal to spheroid in the point ray1(4:6)
*
*  Rewrite spheroid equation x**2 + y**2 + z**2 = R**2 in form
*  F(z,x,y)=0. Becouse spheroid reflect by its inner surface,
*  F(z,x,y) = R**2 - z**2 + x**2 + y**2
*  Surface normal in the point r=(z0,x0,y0) is grad(F(z,x,y)) in this point.
*  Direct normal vector along the z-axis: F = z**2 - (x**2 + y**2)/tg**2
*  Normalized normal vector o=((1/|o|)*dF/dz, (1/|o|)*dF/dx, (1/|o|)*dF/dy),
*  here |o| = sqrt((dF/dz)**2 + (dF/dx)**2 + (dF/dy)**2)
*
*     .. omit 2 in o(i) and 2*2 in oo
      o1 = -x1
      o2 = -x2
      o3 = -x3
      oo = o1**2 + o2**2 + o3**2
*
*     Directrix cosines of the reflected ray.
*     Vectorize reflection law: ray1 = ray - 2*(o*ray)*o
*                            or ray1 = ray - 2*(o*ray)*o/(o*o)
*     here ray,ray1 - input and output rays (directrix cosines),
*          o - surface normal, (o*ray) - scalar product. Note: (o*ray) < 0.
*
*     .. include factor 2 to oray
      oray = 2.d0*(o1*c1 + o2*c2 + o3*c3)/oo
*      print*, 'o:', o1,o2,o3
      if (oray .GE. 0.d0) then
*         print*, 'o:', o1,o2,o3
*         print*, '--- Positive oray =', oray
         goto 1002
      endif
*     .. directrix cosines of the output ray
      ray1(1) = c1 - oray*o1
      ray1(2) = c2 - oray*o2
      ray1(3) = c3 - oray*o3

*     .. normal reflection
      MISPH = 0
      goto 1000
*      RETURN

 1001 continue
*     .. ray have missed the spheroid
      MISPH = -1
      goto 1000
*      RETURN
 1002 continue
*     .. ray is absorbed at the back side of spheroid
      MISPH = 1
      goto 1000
*      RETURN
 1000 continue
*      print*, 'MISPH: ray,ray1:'
*      print*, ray
*      print*, ray1
      RETURN
      END

*                                         @METAGS MIRCIR
*                                         07-08-98 09:40pm
*--------------- MIRCIR ---------------
*
      integer function MIRCIR(R,ray,ray1)
*
*        x|
*         |
*         |R
*         O-------
*        /       z
*      y/
*
      IMPLICIT NONE
      real R,ray(6),ray1(6)
      real t
      integer i

      do i=1,6
         ray1(i) = ray(i)
      enddo
*
*     Intersection point of line and plane z=0
*     Line: direcrix cosines ray(3) and point aray(3)
*     Solve system of linear equations:
*     Ray: dir.cosines=(cosz,cosx,cosy), point aray=(z0,x0,y0)
*     eq.: (x-x0)/cosx = (y-y0)/cosy = (z-z0)/cosz
*     Plane: z=0
*
      if (ray(1) .NE. 0.) then
         t = -ray(4)/ray(1)
      else
*        .. ray is parallel to plane
         goto 1001
      endif
*      print*, 'MIRCIR: t=', t
      if (t .LT. 0.) then
*        .. ray have missed the plane
         goto 1001
      endif

*     .. intersection point ray1(3+i) = ray(3+i) + ray(i)*t
      ray1(5) = ray(5) + ray(2)*t
      ray1(6) = ray(6) + ray(3)*t
*     .. z precisely on plane
      ray1(4) = 0.

*     .. test to hit
      if (ray1(5)**2 + ray1(6)**2 .GT. R**2) then
*        .. miss
         goto 1001
      endif
*
*     Directrix cosines of the reflected ray.
*     Vectorize reflection law: ray1 = ray - 2*(o*ray)*o
*     here ray,ray1 - input and output rays (directrix cosines),
*          o - surface normal, (o*ray) - scalar product. Note: (o*ray) < 0.
*     In our frame o=(1,0,0) and (o*ray) = ray(1)
*
      if (ray(1) .GE. 0.) then
*        .. (o*ray) > 0: ray hitted the back side of plane
         goto 1002
      endif
*     .. to reflect enough invert z-component (already ray1(1:3)=ray(1:3))
      ray1(1) = -ray(1) 

*     .. normal reflection
      MIRCIR = 0
      goto 1000
*      RETURN
 1001 continue
*     .. ray have missed the plane
      MIRCIR = -1
      goto 1000
*      RETURN
 1002 continue
*     .. ray is absorbed at the back side of plane
      MIRCIR = 1
      goto 1000
*      RETURN
 1000 continue
*      print*, 'MIRCIR: ray,ray1:'
*      print*, ray
*      print*, ray1
      RETURN
      END


*                                         @METAGS PARABpar
*                                         07-10-98 04:25pm
*--------------- PARABpar ---------------
*
      SUBROUTINE PARABpar(Rfoc,Rmax, f,h)
*     .. Rfoc - paraboloid radius in the focus
      f = Rfoc/2.
      h = Rmax**2/(4.*f)
      END

*                                         @METAGS PARAB
*                                         06-02-98 01:17pm
*--------------- PARAB ---------------
*
      SUBROUTINE PARAB(f,zmin,zmax,maxray,rayhis,nray,ierr)
*
*     Commonly, p = q = 2f
*     Reflection from paraboloid mirror 2z = x**2/p + y**2/q
*
*     Note: ray,ray1 in frame of paraboloid
*           (i.e. normal in paraboloid top is (1,0,0))
*
*  Intersection point ray1(4:6) of paraboloid 2z = x**2/p + y**2/q
*  and ray with radiant ray(4:6) and directrix cosines ray(1:3)
*  Notation: ray(cos1,cos2,cos3,x1,x2,x3) == ray(cosz,cosx,cosy,z,x,y)
*
*     Returns error code:
*      0: normal reflection
*      1: ray is absorbed at the back side of paraboloid
*     -1: ray have missed the paraboloid
*
*              x|       
*               |      .
*               |     /|
*           ray1(4:6)/ |   ray(1:3) ray(4:6)
*               |   *-------<--------*-----
*               |  /   |
*               | |    |
*              0| |zmin|zmax  
*            ---+-+----+-----------------------------
*              /| |    |                            z
*             / | |    |
*            /  |  \   |           
*           /   |   \  |
*          /    |    \ |
*        y/     |     \|
*               |      `
*
      IMPLICIT NONE
      integer maxray,nray,ierr
      real f,zmin,zmax,rayhis(6,maxray+1)
      double precision p,q,z1,z2
      double precision c1,c2,c3,x1,x2,x3, z
      double precision a,b,c,t,t1,t2
      double precision o1,o2,o3,oo,oray
      double precision wavel2,vlen2
*     .. square of wave length (in cm) used as an epsilon
      parameter (wavel2 = 4.d-10)
      integer i,idroot

      ierr = 0
      nray = 1

      do i=1,6
         rayhis(i,nray+1) = rayhis(i,nray)
      enddo

      p = 2.d0*f
      q = p
*     .. paraboloid boundaries along the z
      z1 = zmin
      z2 = zmax
      if (zmax.EQ.0.) z2=f

*     .. cosines of incident ray
      c1 = rayhis(1,nray)
      c2 = rayhis(2,nray)
      c3 = rayhis(3,nray)
*     .. radiant of incident ray
      x1 = rayhis(4,nray)
      x2 = rayhis(5,nray)
      x3 = rayhis(6,nray)

*
*     Rewrite line equation (z-z0)/cosz=(x-x0)/cosx=(y-y0)/cosz
*     in parametric form:
*     z = z0 + cosz*t
*     x = x0 + cosx*t
*     y = y0 + cosy*t
*     here t - parameter. t=0 corresponds to radiant,
*     t > 0 corresponds going downstream the directrix vector
*     t < 0 - going in backward direction.
*     To find t solve the equation x**2/p + y**2/q - 2*z = 0
*
*     .. coeffs. of square equation
      a = c2**2/p + c3**2/q
      b = 2.*(c2*x2/p + c3*x3/q - c1)
      c = x2**2/p + x3**2/q - 2.*x1
      call DROOT(a,b,c,t1,t2,idroot)
*      print*, 'roots:', t1,t2
      if (idroot.NE.0) then
*         print*, 'PARAB: idroot =', idroot
*         call WAIT('PARAB: No roots of square equation')
*         goto 1001

         if (idroot.LT.0) call MESS('PARAB: No roots')
         if (idroot.EQ.1) call MESS('PARAB: Any roots')
         if (idroot.EQ.2) call MESS('PARAB: Complex roots')
         goto 1001
      endif

*     .. start test from MIN root becouse both roots may be positive
      t = MIN(t1,t2)
      z = x1 + c1*t
*      print*, 'MIN root = ', t, ', z =', z
      if ((z.LT.z1).OR.(z.GT.z2)) then
*        .. miss
*         print*, '--- Point out of paraboloid. Change to MAX root ---'
         goto 100
      endif
*     .. test vector length before the sign: if t ~ 0 may be any sign     
      vlen2 = (c1**2 + c2**2 + c3**2)*t**2
      if (vlen2.LT.wavel2) then
*         print*, '--- vlen too small. Change to MAX root ---'
*         print*, 'vlen2 =', vlen2
         goto 100
      endif
      if (t .LT. 0.d0) then
*         print*, '--- Negative MIN root. Change to MAX root ---'
         goto 100
      endif
      goto 200

  100 continue
      t = MAX(t1,t2)
      z = x1 + c1*t
*      print*, 'MAX root =', t, ', z =', z
      if ((z.LT.z1).OR.(z.GT.z2)) then
*        .. miss
*         print*, '--- MAX out of paraboloid. Exit ---'
         goto 1001
      endif
*     .. test vector length before the sign: if t ~ 0 may be any sign     
      vlen2 = (c1**2 + c2**2 + c3**2)*t**2
      if (vlen2.LT.wavel2) then
*         print*, '--- MAX vlen too small. Exit ---'
*         print*, 'vlen2 =', vlen2
         goto 1001
      endif
      if (t .LT. 0.d0) then
*         print*, '--- Negative MAX root. Exit ---'
         goto 1001
      endif
  200 continue

*      print*, 't1,t2,t:'
*      print*, t1,t2,t
*      print*, 'eq.:', a*t**2 + b*t + c

      nray = nray+1

*     .. adjust point on paraboloid
      x2 = x2 + c2*t
      x3 = x3 + c3*t
      x1 = (x2**2/p + x3**2/q)/2.d0
*      z = REAL(x1 + c1*t)
*     .. radiant of the output ray
      rayhis(4,nray) = x1
      rayhis(5,nray) = x2
      rayhis(6,nray) = x3
*
*  Gets directrix cosines o of nornal to paraboloid in the point ray1(4:6)
*
*  Rewrite paraboloid equation x**2/p + y**2/q - 2*z = 0 in form
*  F(z,x,y)=0. Becouse paraboloid reflect by its inner surface,
*  F(z,x,y) = 2*z - x**2/p + y**2/q
*  Surface normal in the point r=(z0,x0,y0) is grad(F(z,x,y)) in this point.
*  Direct normal vector along the z-axis: F = z**2 - (x**2 + y**2)/tg**2
*  Normalized normal vector o=((1/|o|)*dF/dz, (1/|o|)*dF/dx, (1/|o|)*dF/dy),
*  here |o| = sqrt((dF/dz)**2 + (dF/dx)**2 + (dF/dy)**2)
*
*     .. omit 2 in o(i) and 2*2 in oo
      o1 =  1.d0
      o2 = -x2/p
      o3 = -x3/q
      oo = o1**2 + o2**2 + o3**2
*
*     Directrix cosines of the reflected ray.
*     Vectorize reflection law: ray1 = ray - 2*(o*ray)*o
*                            or ray1 = ray - 2*(o*ray)*o/(o*o)
*     here ray,ray1 - input and output rays (directrix cosines),
*          o - surface normal, (o*ray) - scalar product. Note: (o*ray) < 0.
*
*     .. include factor 2 to oray
      oray = 2.d0*(o1*c1 + o2*c2 + o3*c3)/oo
      if (oray .GE. 0.d0) then
*         print*, '--- Positive oray =', oray
         goto 1002
      endif
*     .. directrix cosines of the output ray
      c1 = c1 - oray*o1
      c2 = c2 - oray*o2
      c3 = c3 - oray*o3
*     .. assign to output ray
      rayhis(1,nray) = c1
      rayhis(2,nray) = c2
      rayhis(3,nray) = c3

       if (nray.GT.maxray) then
*         .. history overflow
          goto 1001
       endif

*
*     Reflect reflected ray if possible
*
*     .. one root=0. Use Viete formula for the other root.
*     .. here c1,c2,c3 is cosines and x1,x2,x3 is radiant of reflected ray
*     .. coeffs. of square equation
      a = c2**2/p + c3**2/q
      b = 2.*(c2*x2/p + c3*x3/q - c1)
 
*      c = x2**2/p + x3**2/q - 2.*x1
*      call DROOT(a,b,c,t1,t2,idroot)
*      print*, 'Compare: roots:', t1,t2

      if (a .NE. 0.d0) then
         t1 = -b/a
*         print*, 'by Viete t1=', t1
*         print*, 'Compare: eq.:', a*t1**2 + b*t1 + c
      else
*        .. the other root=0 too
*           there is no intersection point downstream the ray
         nray = nray+1
         goto 1001
      endif
*     .. to use the same code
      t2 = t1

*     .. repeat
      goto 100

 1001 continue
*     .. ray have missed the paraboloid
      if (nray.EQ.1) ierr = -1
      nray = nray-1
      goto 1000
*      RETURN
 1002 continue
*     .. ray is absorbed at the back side of paraboloid
      ierr = 1
      goto 1000
*      RETURN
 1000 continue
*      print*, 'PARAB: nray =', nray, '   rays for nray-1, nray:'
*      print*, (rayhis(j,nray-1), j=1,6)
*      print*, (rayhis(j,nray), j=1,6)
      RETURN
      END


*                                         @METAGS CONEpar
*                                         07-10-98 04:28pm
*--------------- CONEpar ---------------
*
      SUBROUTINE CONEpar(r1,r2,dh, tgcone,h1,h2)
      tgcone = (r2-r1)/dh
      h1 = r1/tgcone
      h2 = r2/tgcone
      END

*                                         @METAGS CONE
*                                         06-02-98 01:17pm
*--------------- CONE ---------------
*
      SUBROUTINE CONE(tgcone,zmin,zmax,maxray,rayhis,nray,ierr)
*
*     Reflection from conical mirror (x**2 + y**2)/tgcone**2 - z**2 = 0
*
*     Note: ray,ray1 in frame of cone
*           (i.e. normal in the plane of cutted cone top is (1,0,0))
*
*  Intersection point ray1(4:6) of cone (x**2 + y**2)/tgcone**2 - z**2 = 0
*  and ray with radiant ray(4:6) and directrix cosines ray(1:3)
*  Notation: ray(cos1,cos2,cos3,x1,x2,x3) == ray(cosz,cosx,cosy,z,x,y)
*
*     Returns error code:
*      0: normal reflection
*      1: ray is absorbed at the back side of cone
*     -1: ray have missed the cone
*
*              x|       
*               |      .
*               |     /|
*           ray1(4:6)/ |   ray(1:3) ray(4:6)
*               |   *-------<--------*-----
*               |  /   |
*               | |    |
*              0| |zmin|zmax  
*            ---+-+----+-----------------------------
*              /| |    |                            z
*             / | |__  |
*            /  |  \)  | <-- tgcone
*           /   |   \  |
*          /    |    \ |
*        y/     |     \|
*               |      `
*
      IMPLICIT NONE
      integer maxray,nray,ierr
      real tgcone,zmin,zmax,rayhis(6,maxray+1)
      double precision tg2,z1,z2
      double precision c1,c2,c3,x1,x2,x3
      double precision a,b,c,t,t1,t2
      double precision o1,o2,o3,oo,oray
      double precision z,x,y
      double precision wavel2,vlen2
*     .. square of wave length (in cm) used as an epsilon
      parameter (wavel2 = 4.d-10)
      integer i,idroot

      ierr = 0
      nray = 1

      do i=1,6
         rayhis(i,nray+1) = rayhis(i,nray)
      enddo

      tg2 = tgcone**2
*     .. cone boundaries along the z
      z1 = zmin
      z2 = zmax

*     .. cosines of incident ray
      c1 = rayhis(1,nray)
      c2 = rayhis(2,nray)
      c3 = rayhis(3,nray)
*     .. radiant of incident ray
      x1 = rayhis(4,nray)
      x2 = rayhis(5,nray)
      x3 = rayhis(6,nray)

*
*     Rewrite line equation (z-z0)/cosz=(x-x0)/cosx=(y-y0)/cosz
*     in parametric form:
*     z = z0 + cosz*t
*     x = x0 + cosx*t
*     y = y0 + cosy*t
*     here t - parameter. t=0 corresponds to radiant,
*     t > 0 corresponds going downstream the directrix vector
*     t < 0 - going in backward direction.
*     To find t solve the equation (x**2 + y**2)/tgcone**2 - z**2 = 0
*
*     .. coeffs. of square equation
      a = (c2**2 + c3**2)/tg2 - c1**2
      b = 2.d0*((c2*x2 + c3*x3)/tg2 - c1*x1)
      c = (x2**2 + x3**2)/tg2 - x1**2
      call DROOT(a,b,c,t1,t2,idroot)
*      print*, a,b,c
*      print*, t1,t2
      if (idroot.NE.0) then
         if (idroot.EQ.2) then
            call MESS('CONE: Complex roots: no intersection')
            goto 1001
         endif
         if (idroot.LT.0) call MESS('ERROR CONE: No roots')
         if (idroot.EQ.1) call MESS('CONE: Any roots')
         if (x1.GT.z2) then
            z = 0.
         elseif (x1.LT.z1) then
            z = z1-z2
         else
            z = x1
         endif
         t = (z-x1)/c1
         rayhis(5,nray) = x2 + c2*t
         rayhis(6,nray) = x3 + c3*t
*        .. z precisely on cone
         rayhis(4,nray) = REAL(sqrt(x2**2 + x3**2)/tgcone)
         goto 1002
      endif

*     .. start test from MIN root becouse both roots may be positive
      t = MIN(t1,t2)
      z = x1 + c1*t
      if ((z.LT.z1).OR.(z.GT.z2)) then
*        .. miss
*         print*, '--- Point out of cone. Change to MAX root ---'
*         print*, z
         goto 100
      endif
*     .. test vector length before the sign: if t ~ 0 may be any sign     
      vlen2 = (c1**2 + c2**2 + c3**2)*t**2
      if (vlen2.LT.wavel2) then
*         print*, '--- vlen too small. Change to MAX root ---'
*         print*, 'vlen2 =', vlen2
         goto 100
      endif
      if (t .LT. 0.d0) then
*         print*, '--- Negative MIN root. Change to MAX root ---'
         goto 100
      endif
      goto 200

  100 continue
      t = MAX(t1,t2)
      z = x1 + c1*t
      if ((z.LT.z1).OR.(z.GT.z2)) then
*        .. miss
*         print*, '--- MAX out of cone. Exit ---'
         goto 1001
      endif
*     .. test vector length before the sign: if t ~ 0 may be any sign     
      vlen2 = (c1**2 + c2**2 + c3**2)*t**2
      if (vlen2.LT.wavel2) then
*         print*, '--- MAX vlen too small. Exit ---'
*         print*, 'vlen2 =', vlen2
         goto 1001
      endif
      if (t .LT. 0.d0) then
*         print*, '--- Negative MAX root. Exit ---'
         goto 1001
      endif
  200 continue

*      print*, 't1,t2,t:'
*      print*, t1,t2,t
*      print*, 'eq.:', a*t**2 + b*t + c

      nray = nray+1

*     .. adjust point on cone
      x = REAL(x2 + c2*t)
      y = REAL(x3 + c3*t)
      z = REAL(sqrt(x**2 + y**2)/DBLE(tgcone))
*     .. radiant of the output ray
      rayhis(4,nray) = z
      rayhis(5,nray) = x
      rayhis(6,nray) = y
*
*  Gets directrix cosines o of nornal to cone in the point ray1(4:6)
*
*  Rewrite cone equation (x**2 + y**2)/tg**2 - z**2 = 0 in form
*  F(z,x,y)=0. Becouse cone reflect by its inner surface,
*  F(z,x,y) = z**2 - (x**2 + y**2)/tg**2
*  Surface normal in the point r=(z0,x0,y0) is grad(F(z,x,y)) in this point.
*  Direct normal vector along the z-axis: F = z**2 - (x**2 + y**2)/tg**2
*  Normalized normal vector o=((1/|o|)*dF/dz, (1/|o|)*dF/dx, (1/|o|)*dF/dy),
*  here |o| = sqrt((dF/dz)**2 + (dF/dx)**2 + (dF/dy)**2)
*
*  For cone (o*o) = 4*(z**2 + (x**2+y**2)/tg**2) = 4*2*z**2
*  is only function of z.
*
*     .. omit 2 in o(i) and 2*2 in oo
      o1 =  z
      o2 = -x/tg2
      o3 = -y/tg2
      oo = o1**2 + o2**2 + o3**2
*
*     Directrix cosines of the reflected ray.
*     Vectorize reflection law: ray1 = ray - 2*(o*ray)*o
*                            or ray1 = ray - 2*(o*ray)*o/(o*o)
*     here ray,ray1 - input and output rays (directrix cosines),
*          o - surface normal, (o*ray) - scalar product. Note: (o*ray) < 0.
*
      oray = o1*c1 + o2*c2 + o3*c3
      if (oray .GE. 0.d0) then
*         print*, '--- Positive oray =', oray
         goto 1002
      endif
*     .. directrix cosines of the output ray
      rayhis(1,nray) = c1 - 2.d0*oray*o1/oo
      rayhis(2,nray) = c2 - 2.d0*oray*o2/oo
      rayhis(3,nray) = c3 - 2.d0*oray*o3/oo

       if (nray.GT.maxray) then
*         .. history overflow
          goto 1002
       endif

*
*     Reflect reflected ray if possible
*
*     .. one root=0. Use Viete formula for the other root.
*     .. cosines of incident ray
      c1 = rayhis(1,nray)
      c2 = rayhis(2,nray)
      c3 = rayhis(3,nray)
*     .. radiant of incident ray
      x1 = rayhis(4,nray)
      x2 = rayhis(5,nray)
      x3 = rayhis(6,nray)
*     .. coeffs. of square equation
      a = (c2**2 + c3**2)/tg2 - c1**2
      b = 2.d0*((c2*x2 + c3*x3)/tg2 - c1*x1)
*     .. becouse a.NE.0.
      t1 = -b/a
*     .. to use the same code
      t2 = t1
*     .. repeat
      goto 100

 1001 continue
*     .. ray have missed the cone
      if (nray.EQ.1) ierr = -1
      nray = nray-1
      goto 1000
*      RETURN
 1002 continue
*     .. ray is absorbed at the back side of cone
      ierr = 1
      goto 1000
*      RETURN
 1000 continue
*      print*, 'CONE: nray =', nray, '   rays for nray-1, nray:'
*      print*, (rayhis(j,nray-1), j=1,6)
*      print*, (rayhis(j,nray), j=1,6)
      RETURN
      END

*                                         @METAGS PMT
*                                         07-08-98 10:10pm
*--------------- PMT ---------------
*
      SUBROUTINE PMT(Rpmt,ray,ray1,ierr)
      IMPLICIT NONE
      real Rpmt,ray(6),ray1(6)
      integer ierr
      real t
      integer i

      do i=1,6
         ray1(i) = ray(i)
      enddo
*
*     Intersection point of line and plane z=0
*     Line: direcrix cosines ray(3) and point aray(3)
*     Solve system of linear equations:
*     Ray: dir.cosines=(cosz,cosx,cosy), point aray=(z0,x0,y0)
*     eq.: (x-x0)/cosx = (y-y0)/cosy = (z-z0)/cosz
*     Plane: z=0
*
      if (ray(1) .NE. 0.) then
         t = -ray(4)/ray(1)
      else
*        .. ray is parallel to PM plane
         goto 1001
      endif
*      print*, 'PMT: t=', t
      if (t .LT. 0.) then
*        .. ray have missed the plane
         goto 1001
      endif

*     .. intersection point ray1(3+i) = ray(3+i) + ray(i)*t
      ray1(5) = ray(5) + ray(2)*t
      ray1(6) = ray(6) + ray(3)*t
*     .. z precisely on plane
      ray1(4) = 0.

*     .. test to hit
      if (ray1(5)**2 + ray1(6)**2 .GT. Rpmt**2) then
*        .. miss
         goto 1001
      endif

*     .. hit the PM
      ierr = 0
      goto 1000
*      RETURN
 1001 continue
*     .. ray have missed the PM
      ierr = -1
      goto 1000
*      RETURN
 1002 continue
*     .. ray is absorbed at the back side of PM
      ierr = 1
      goto 1000
*      RETURN
 1000 continue
*      print*, 'PMT: ray,ray1:'
*      print*, ray
*      print*, ray1
      RETURN
      END

c                                      @METAGS RANBAG
c---------- RANBAG ----------
c
      real function ranbag(idummy)
c
c     Uniform random generator. Uses RANLUX.
c     Adjusting parameters:
c     LEN - bag size
c     LUX - RANLUX luxury level
c
      IMPLICIT NONE
      integer LEN,LUX
      parameter (LEN=10000)
      parameter (LUX=3)
      real bag(LEN)
      integer ibuf(25)
      integer label, iniran,irused,ncur,idummy
      logical first
      save bag,ibuf, first, iniran,irused, ncur
*     .. entry BAGINT
      real BAGINT
      integer init
*     .. entry BAGREP
      real BAGREP
      integer iused
*     .. entry RNOR
      real RNOR
      real amean,sigma
      real  S, T, A, B, R1, R2
      SAVE  S, T, A, B, R1, R2
      real U(2),V,X,Y,Q,DEVIAT
*     .. entry RNORINT
      real RNORINT
*     .. entry RNOREP
      real RNOREP
      integer nhit,miss, nhit1,miss1
      save    nhit,miss
      real ratio
c
      data iniran/1/
      data irused/0/
      data ncur/LEN/
      data first/.TRUE./
      data nhit,miss /0,0/
c
      assign 2000 to label
 1000 ncur = ncur+1
      if (ncur.GT.LEN) then
         if (first) then
c           .. initialization
            first=.FALSE.
            irused = 0
            call RLUXGO(LUX,iniran,0,0)
         else
c           .. restart from saved int vector
            call RLUXIN(ibuf)         
         endif
c        .. create new bag
         call RANLUX(bag,LEN)
c        .. save int vector to restart
         call RLUXUT(ibuf)
         ncur=1
      endif
      ranbag = bag(ncur)
      irused = irused+1
      goto label
 2000 RETURN
c
c---     @METAGS BAGINT_ENTRY
c
      ENTRY BAGINT(init)
*     .. Entry to initialize RANBAG
 3000 iniran = init
      ncur  = LEN
      first = .TRUE.
*     .. bring the random to initialize RANLUX
      assign 4000 to label
      goto 1000
 4000 BAGINT = 0.
      RETURN
c
c---     @METAGS BAGREP_ENTRY
c
      ENTRY BAGREP(iused)
*     .. Entry to report RANBAG
      print*, '*** RANBAG: There are', irused, ' randoms were used ***'
      iused = irused
      BAGREP = 0.
      RETURN
c                                   @METAGS RNOR_ENTRY
c------------- RNOR ---------------
c
      ENTRY RNOR(idummy,amean,sigma)
*     .. Entry to get normal random
*
*     Changed form of RNORMX to use RANBAG
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
      DATA  S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA  R1, R2/ 0.27597, 0.27846/
C         generate pair of uniform deviates
   49 continue
*     .. bring the random
      assign 50 to label
      goto 1000
   50 U(1) = bag(ncur)
*     .. bring the random
      assign 51 to label
      goto 1000
   51 U(2) = bag(ncur)
      V = 1.7156 * (U(2) - 0.5)
      X = U(1) - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)
C           accept P if inside inner ellipse
      IF (Q .LT. R1) THEN
         miss = miss+1
         GO TO 100
      ENDIF
C           reject P if outside outer ellipse
      IF (Q .GT. R2) THEN
         miss = miss+1
         GO TO 49
      ENDIF
C           reject P if outside acceptance region
      IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2)  GO TO 49
C           ratio of P's coordinates is normal deviate
  100 DEVIAT = V/U(1)
      nhit = nhit+1
      RNOR = amean + sigma*DEVIAT
      RETURN
*
*---     @METAGS RNORINT_ENTRY
*
      ENTRY RNORINT(init)
*     .. Entry to initialize RNOR
      nhit = 0
      miss = 0
      goto 3000
*
*---     @METAGS RNOREP
*
      ENTRY RNOREP(nhit1,miss1)
*     .. Entry to report RNOR
      ratio = REAL(miss)/REAL(nhit)
      print*, '*** RNOR report to genrated randoms:'
      print*, 'Generated ', nhit, ' points, missing shots =', miss
      print*, 'There are ', ratio, ' missings per hit'
      print*, '***'
      nhit1 = nhit
      miss1 = miss
      RNOREP = ratio
      END

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
         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      p =',LP
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
***      WRITE(6,'(5X,5I12)') ISDEXT
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
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
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
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
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


*                                         @METAGS MESS.COMIS
*                                         01-22-98 02:11pm
*--------------- MESS ---------------
*
      SUBROUTINE MESS(line)
      character line*(*)
      length = LENOCC(line)
      if (length.GT.0) print 1, line(1:length)
      RETURN
    1 FORMAT(A)
      END

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

