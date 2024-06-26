*     .. Autogenerated function
      REAL FUNCTION alhigh(x)
*TYPE Fit function EP
*EMIN    270.000
*EMAX    380.000
      real PAR( 4)
      data PAR/
     #  .0000000E+00    ,
     #  8.360350        ,
     # -.1529690E-01    ,
     # -9.872491        /
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
 
      alhigh = result
      END
