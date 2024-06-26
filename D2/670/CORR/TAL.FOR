*     .. Autogenerated function
      REAL FUNCTION tallow(x)
*     .. Titov ALlow 9.03-18.01 up to 320 MeV
*TYPE Fit function EX
*EMIN    148.667
*EMAX    320.000
      real PAR( 3)
      data PAR/
     #  .0000000E+00    ,
     #  7.929532        ,
     # -.1369659E-01    /
      idummy = PAR(1)
      e = 0.
      arg = PAR(2) + x*PAR(3)
      if (ABS(arg) .LT. 50.) e=exp(arg)
      result = e
 
      tallow = result
      END
*
*--
*
*     .. Autogenerated function
      REAL FUNCTION talhigh(x)
*     .. Titov ALhigh
*TYPE Fit function EX
*EMIN    270.000
*EMAX    420.000
      real PAR( 3)
      data PAR/
     #  .0000000E+00    ,
     #  9.349513        ,
     # -.1843024E-01    /
      idummy = PAR(1)
      e = 0.
      arg = PAR(2) + x*PAR(3)
      if (ABS(arg) .LT. 50.) e=exp(arg)
      result = e
 
      talhigh = result
      END
*
*--
*
      REAL FUNCTION tal(x)
      parameter (Eprim=299.95)
      if (x.LT.Eprim) then
         result = tallow(x)
      else
         result = talhigh(x)
      endif
      tal = result
      END
