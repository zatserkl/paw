* 1000   2000   3000   nran    *** WARNING: FIRSTS LINES USED AS INI DATA! 
* 2      2       init INT
* 20.    Rbox, cm
* 500.   Dbox, cm
* 140    160.   f, cm
* 100.   Rref, cm
* 2.5    Rpmt, cm
* 34.9   rindex: index of refraction (n-1)*1e6
* 90.    PMconst
* 0.20   PMconv 
* 20.    angmax, mrad
* 25.    Ppi, GeV
* 20.    Rbeam
*
      program slope
*      IMPLICIT NONE
      PARAMETER (NWPAWC = 300000)
      PARAMETER (LRECL  = 1024)
      PARAMETER (id = 1)
      COMMON /PAWC/ IPAW(NWPAWC)

*     .. small mirror shift to prevent impact at the back surface
*        = 0 here
      parameter (smallshift=0.)

      parameter (MAXRAY=100)
      real rayhis(6,MAXRAY+1)
      
      real sysmir(5),syscol(5),syspmt(5), systst(5)
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
*                                   @METAGS NTcommon
      COMMON /OPTIC/ ireg,ip,part(6),ray0(6),raytst(6),intref,raypmt(6)
      real raylab(6),rayloc(6),rayref(6)

      character*32 thisf,thisn,ifile,bfile, afile
      double precision betan
      logical pion

      data CHTOP /'NTDIR'/
      data CHFILE/' '/
      data afile /'unknown'/

      call GETARG(0,thisf)
      call FPARSE(thisf,thisn,thisf,'x')
      lename = LENOCC(thisn)
      CHFILE = thisn(1:lename)//'.hbook'
      ifile  = thisn(1:lename)//'.f'
      bfile  = thisn(1:lename)//'_t.for'
**
**     Ntuple
**      
*      write(chtitl,10)
*   10 format ('A half of spherical mirror')
*      call MESS('Ntuple title:')
*      call MESS(chtitl)

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
      read (lun,*) c,f
      read (lun,*) c,Rref
      read (lun,*) c,Rpmt
      read (lun,*) c,rindex
      read (lun,*) c,PMconst
      read (lun,*) c,PMconv
      read (lun,*) c,angmax
      read (lun,*) c,Ppi
      read (lun,*) c,Rbeam

      close(lun)

      print*, 'Ini file data:'
      print*, 'N randoms/bin =', nran
      print*, 'Rbox  =', Rbox, ' cm'
      print*, 'Dbox =', Dbox, ' cm'
      print*, 'f =', f, ' cm'
      print*, 'Rref =', Rref, ' cm'
      print*, 'Rpmt =', Rpmt, ' cm'
      print*, 'rindex =', rindex, ' (n-1)*1e6'
      print*, 'PMconst =', PMconst, ' 1/cm'
      print*, 'PMconv =', PMconv
      print*, 'angmax  =', angmax, ' mrad'
      print*, 'Ppi  =', Ppi, ' GeV'
      print*, 'Rbeam  =', Rbeam, ' cm'

c     .. convert
      pi = ACOS(-1.)
      torad = pi/180.
      todeg = 180./pi
      
      angrad = angmax/1000.
      
*                                         @METAGS Geometry
*       Geometry
*
*                     _/\
*                   _/   \ cone+PM
*                 _/      \
*          x|    /       / 
*           |    \ Rref |    
*           |___  \    /____________ 
*  sloped   |   _  \  |             |
*  mirror   |   /|  \/              |
*  -------->|  /                    |
*           | /                     |
*          0|_ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _
*          /|\                      |        z
*         / | \    sloped           |
*        /  |  \   mirror           |
*      y/   |   \<--------          |
*           |____\__________________|
*

*     .. Rref is z-position of cone input center
*     .. collector and PM angle, rad
      alpha = atan2(Rref,Dbox)
*                                         @METAGS mirror_frame
*     .. mirror with Rmir = 2f
*      Rmir = 2.*Dbox + 1.
      Rmir = 2.*f
*     .. mirror frame
      iaxis = 3
*     .. mirror slope = alpha/2
      thmir = todeg*alpha/2.
      zmir = 0.-smallshift
      xmir = 0.
      call mkframe(iaxis,thmir,zmir,xmir, sysmir)

*                                         @METAGS tst_frame
*     .. test plane
      Rtst = 2.*Rbox
      if (Rtst.EQ.0.) Rtst=2.*10.
      iaxis = 3
      thtst = 180. + todeg*alpha
*      ztst = Dbox-smallshift
      ztst = Rref
      xtst = Rbox
      call mkframe(iaxis,thtst,ztst,xtst, systst)

*                                         @METAGS collector_frame
*     .. collector
      Rcol1 = Rpmt
      Rcol2 = 12.
*     .. cone height
      dhcol = 60.
      print*, 'Current Rcol2,dhcol:', Rcol2,dhcol
      print*, 'Enter you settings (''/''=untouched)'
      read *, Rcol2,dhcol
      print*, 'Calculate for Rcol2,dhcol:', Rcol2,dhcol

      call CONEpar(Rcol1,Rcol2,dhcol, tgcone,h1col,h2col)
*     .. collector frame
      iaxis = 3
      thcol = 180. + todeg*alpha
*      zcol = Dbox-smallshift + h2col*cos(alpha)
      zcol = Rref + h2col*cos(alpha)
      xcol = Rbox + h2col*sin(alpha)
      call mkframe(iaxis,thcol,zcol,xcol, syscol)

*                                         @METAGS pmt_frame
*     .. PM
*     .. pmt frame
      iaxis = 3
      thpmt = 180. + todeg*alpha
*      zpmt = Dbox-smallshift + (h2col-h1col)*cos(alpha)
      zpmt = Rref + (h2col-h1col)*cos(alpha)
      xpmt = Rbox + (h2col-h1col)*sin(alpha)
      call mkframe(iaxis,thpmt,zpmt,xpmt, syspmt)
      
*      print*, ' '
      print*, 'Distance to PM is', zpmt, ' cm from counter box back'
      print*, 'Distance to PM is', xpmt, ' cm from counter axis'

*--   .. electron
*     .. Cherenkov angle
      thetae = ACOS(1.D0/(1.D0 + rindex*1.D-6))
*     .. the number of photoelectrons
      pee = Dbox*PMconst*(SIN(thetae))**2
*     .. the number of photons
      phote = pee/PMconv
*     .. mean free path length between emitances
      alame = Dbox/phote
*--   .. pion
      betapi = 1.D0/(1.D0 + (0.140D0/Ppi)**2)
      betan = betapi*(1.D0 + rindex*1.D-6)
      if (betan .GT. 1.D0) then
         thetapi = ACOS(1.D0/(betan))
*        .. the number of photoelectrons
         pepi = Dbox*PMconst*(SIN(thetapi))**2
*        .. the number of photons
         photpi = pepi/PMconv
*        .. mean free path length between emitances
         alampi = Dbox/photpi
      else
         pepi = 0.
         photpi = 0.
         alampi = 2.*Dbox
      endif
      print*, ' '
      print*, 'electron Cherenkov angle', thetae, ' rad'
      print*, 'electron the average number of photoelectrons', pee
      print*, 'electron the average number of photons', phote
      print*, 'electron mean free path length', alame
      print*, 'pion     Cherenkov angle', thetapi, ' rad'
      print*, 'pion     the average number of photoelectrons', pepi
      print*, 'pion     the average number of photons', photpi
      print*, 'pion     mean free path length', alampi

      dummy = BAGINT(init)

*
*     Ntuple
*      
      write(chtitl,20) Rcol2,dhcol, f
   20 format ('Sloped mirror. Rcol2,dhcol=', 2F4.0, ' f=', F4.0, ' cm')
      call MESS('Ntuple title:')
      call MESS(chtitl)

      print*, ' '
      call WAIT('Begin the calculation')

      CALL HLIMIT(NWPAWC)
*     .. open a new RZ file
      ntlun = LUNFREE(1)
      CALL HROPEN(ntlun,CHTOP,CHFILE,'N',LRECL,ISTAT)
*     .. book Ntuple
      CALL HBNT(id,chtitl,' ')
*     .. define Ntuple
*                                         @METAGS HBNAME
      CALL HBNAME(id,'OPTIC', ireg,
     &            'ireg[-1,1]:I, ip[0,1]:U,
     '             part(6):R,ray0(6):R,raytst(6):R,
     '             intref:U,raypmt(6):R')

      call FILELUN(ntlun,afile)
      print*, 'With unit', ntlun, ' connected file ', afile
      print*, 'Running ..'

      tlim = 10000.
      call TIMEST(tlim)
      
*
*--   MAIN PART
*
      idummy = 0
      nemite = 0
      nrege  = 0
      nemitpi = 0
      nregpi  = 0
      print*, 'Total particles', nran
*      print*, ' '

*                                      @METAGS MAINloop
*------------------------------------- main loop begin
      DO np=1,nran
*     .. electron flag
      pion = .FALSE.
      ip = 0
      theta = thetae
      alam = alame
*     .. particle radiant in circle with radius Rbox
  100 continue
      xp = Rbeam*2.*(.5 - RANBAG(idummy))
      idummy = idummy+1
      yp = Rbeam*2.*(.5 - RANBAG(idummy))
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
      
 1001 continue
*     .. initialize particle current position
      z = zp
      x = xp
      y = yp

 1000 continue      
*     .. current range
      idummy = idummy+1
      range = -alam*LOG(RANBAG(idummy))
*     .. current position
      z = z + range*part(1)
      x = x + range*part(2)
      y = y + range*part(3)
      if (z.LE.0.) goto 1999
      
*     .. emit the photon in the particle frame
      if (pion) then
         nemitpi = nemitpi+1
      else
         nemite = nemite+1
      endif
*     .. clear reflector flag
      iref = 0      

*     .. photon radiant
      ray0(4) = z
      ray0(5) = x
      ray0(6) = y      
*     .. photon direction
      thg = theta
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
      ray0(1) = cz*costhp        - cx*sinthp
      ray0(2) = cz*sinthp*cosphp + cx*costhp*cosphp - cy*sinphp
      ray0(3) = cz*sinthp*sinphp + cx*costhp*sinphp + cy*cosphp

*                                         @METAGS RunRay
*     Ray runs
*
*     .. clear registration flag
      ireg = 0
*     .. clear the number of reflection in the collector
      intref = 0
*     .. clear arrays
      do i=1,6
         raytst(i) = 100.
         raypmt(i) = 100.
      enddo
*                                         @METAGS mirror_ref
*     .. reflection from mirror
      call toloc(sysmir,ray0,rayloc)
      iermir = MISPH(Rmir,rayloc,rayref)
*     .. return to lab
      call tolab(sysmir,rayref,raylab)
      if (iermir.NE.0) then
*        .. lost photon at mirror
**        .. do not save event, get next photon
*         goto 1000
*        .. save event with ireg=-1
         ireg = -1
         goto 2000
      endif

*                                         @METAGS tst_ref
*     .. look at test plane
      call toloc(systst,raylab,rayloc)
*      iertst = MIRCIR(Rtst,rayloc,raytst)
*     .. use routine PMT for hit test plane
      call PMT(Rtst,rayloc,raytst,iertst)
      if (iertst.NE.0) then
*        .. lost ray
         print*, 'Error in hit of test plane'
      endif
      
*                                         @METAGS collector_ref
*--   Reflection from collector
*     .. ray in collector local frame
      nray = 1
      call toloc(syscol,raylab,rayhis(1,nray))
      irest = MAXRAY-nray
      call CONE(tgcone,h1col,h2col,irest,rayhis(1,nray),npro,iercol)
*     .. return to lab
      if (iercol.EQ.0) ncol=ncol+1
      intref = npro
      do i=nray,nray+npro
         call tolab(syscol,rayhis(1,i),rayhis(1,i))
      enddo
      nray = nray+npro
      if (iercol.GT.0) then
*        .. absorbed at the collector. Next photon
         goto 2000
      endif

*                                         @METAGS pmt_ref
*--   Intersection with PM
*     .. ray in PM local frame
      call toloc(syspmt,rayhis(1,nray),rayloc)
      call PMT(Rpmt,rayloc, raypmt,ierpmt)
      if (ierpmt.GE.0) then
         if (ierpmt.EQ.0) then
*           .. absorbed at the PMT.
            ireg = 1
            npmt = npmt+1

            if (pion) then
               nregpi = nregpi+1
            else
               nrege = nrege+1
            endif
          else
*     	    .. lost
             print*, 'Ray is absorbed at PM back'
         endif
      endif

*     .. goto store      
      goto 3000
      
 2000 continue
*      ireg = 0
*      do i=1,6
*          raypmt(i) = 100.
*      enddo

 3000 continue
*     .. store raw
      call HFNT(id)
      
*     .. for the next photon
      goto 1000

 1999 continue
      if ((.NOT.pion) .AND. (photpi .GT. 0.)) then
*        .. pion flag
         pion = .TRUE.
         ip = 1
         theta = thetapi
         alam = alampi
*        .. run pion by the same way
         goto 1001
      endif
      if (MOD(np,100).EQ.0) then
          call IMESS('processed particles ', np)
      endif
      ENDDO
*------------------------------------- main loop end
*
*--   MAIN PART END
*

      print*, ' '
      print*, '--- electrons:'
      print*, 'There are', nemite, ' photons were emitted'
      print*, 'There are', nrege, ' photons were registrated'
      phpart = REAL(nemite)/REAL(nran)
      print*, 'The average number of emitted photons/particle', phpart
      percent = 100.*REAL(nrege)/REAL(nemite)
      print*, 'Registration efficiency is', percent, '%'
      print*, '--- pions:'
      print*, 'There are', nemitpi, ' photons were emitted'
      print*, 'There are', nregpi, ' photons were registrated'
      phpart = REAL(nemitpi)/REAL(nran)
      print*, 'The average number of emitted photons/particle', phpart
      percent = 100.*REAL(nregpi)/REAL(nemitpi)
      print*, 'Registration efficiency is', percent, '%'

      dummy = BAGREP(iused)
      
      print*, ' '
      call TIMEX(time)
      print*,'Job time', time, ' seconds'
      call WAIT('End of calculation. Close files')

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

      goto 99999
      CHPATH = ' '
      call HCDIR(CHPATH,'R')
      print*, ' '
      call HLDIR('//PAWC','T')
      print*, ' '
      print '('' Current directory is '', A)', CHPATH(1:LENOCC(CHPATH))
      call HLDIR(' ','T')

99999 CALL HREND(CHTOP)
      STOP
10000 print*, 'Ini file ', ifile(1:LENOCC(ifile)), ' did not found'
      END
