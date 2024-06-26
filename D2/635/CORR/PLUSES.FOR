*     .. Autogenerated function
      REAL FUNCTION d2p(x)
*TYPE Fit function EP
*EMIN    142.000
*EMAX    410.000
      real PAR( 4)
      data PAR/
     #  .0000000E+00    ,
     #  8.766567        ,
     # -.1778479E-01    ,
     #  2.263293        /
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
      REAL FUNCTION pluses(x)
      result = 0.
      result = result+d2p(x)
      result = result-alp(x)
      pluses = result
      END
