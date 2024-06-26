*     .. Autogenerated function
      REAL FUNCTION allow(x)
*TYPE Fit function EP
*EMIN    142.308
*EMAX    300.000
      real PAR( 4)
      data PAR/
     #  .0000000E+00    ,
     #  8.580540        ,
     # -.1905297E-01    ,
     #  22.41860        /
      ipower = PAR(1)
      e = 0.
      arg = PAR(2) + x*PAR(3)
      if (ABS(arg) .LT. 50.) e=exp(arg)
 
      p = 0.
*     .. the last non-polinom parameter
      last = 1 + 2
      do n=1, ipower+1
         p = p + PAR(last+n)*x**(n-1)
      enddo
      result = e+p
 
      allow = result
      END
*
*--
*
*     .. Autogenerated function
      REAL FUNCTION alhigh(x)
*TYPE Fit function EX
*EMIN    280.000
*EMAX    410.000
      real PAR( 3)
      data PAR/
     #  .0000000E+00    ,
     #  10.14016        ,
     # -.2226839E-01    /
      idummy = PAR(1)
      e = 0.
      arg = PAR(2) + x*PAR(3)
      if (ABS(arg) .LT. 50.) e=exp(arg)
      result = e
 
      alhigh = result
      END
*
*--
*
      REAL FUNCTION al(x)
      parameter (Eprim=282.6)
      if (x.LE.Eprim) then
         result = allow(x)
      else
         result = alhigh(x)
      endif
      al = result
      END
*
*--
*
*     .. Autogenerated function
      REAL FUNCTION d2p(x)
*     .. Only low D2+
*TYPE Fit function EX
*EMIN    144.824
*EMAX    290.000
      real PAR( 3)
      data PAR/
     #  .0000000E+00    ,
     #  8.993238        ,
     # -.1820010E-01    /
      idummy = PAR(1)
      e = 0.
      arg = PAR(2) + x*PAR(3)
      if (ABS(arg) .LT. 50.) e=exp(arg)
      result = e
 
      d2p = result
      END
*
*--
*
*     .. Autogenerated function
      REAL FUNCTION alp(x)
*TYPE Fit function EP
*EMIN    142.000
*EMAX    410.000
      real PAR( 4)
      data PAR/
     #  .0000000E+00    ,
     #  8.075356        ,
     # -.1845684E-01    ,
     #  .9864309        /
      ipower = PAR(1)
      e = 0.
      arg = PAR(2) + x*PAR(3)
      if (ABS(arg) .LT. 50.) e=exp(arg)
 
      p = 0.
*     .. the last non-polinom parameter
      last = 1 + 2
      do n=1, ipower+1
         p = p + PAR(last+n)*x**(n-1)
      enddo
      result = e+p
 
      alp = result
      END
*
*--
*
      REAL FUNCTION diffp(x)
      result = 0.
      result = result+d2p(x)
      result = result-alp(x)
      diffp = result
      END
*
*--
*
      REAL FUNCTION fonlow(x)
      result = 0.
      result = result+al(x)
      result = result+diffp(x)
      fonlow = result
      END
 
*     .. Autogenerated function
      REAL FUNCTION fonhigh(x)
*TYPE Fit function EX
*EMIN    290.000
*EMAX    404.929
      real PAR( 3)
      data PAR/
     #  .0000000E+00    ,
     #  11.51250        ,
     # -.2532100E-01    /
      idummy = PAR(1)
      e = 0.
      arg = PAR(2) + x*PAR(3)
      if (ABS(arg) .LT. 50.) e=exp(arg)
      result = e
 
      fonhigh = result
      END
 
*                                         @METAGS fon
*                                         02-05-99 05:25pm
*--------------- fon ---------------
*
      real FUNCTION TheFUN(x)
      parameter (Eprim=290.3)
      if (x.LE.Eprim) then
         result = fonlow(x)
      else
         result = fonhigh(x)
      endif
      TheFUN = result
      END
