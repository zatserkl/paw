*     .. Autogenerated function
      REAL FUNCTION fiteff(x)
*TYPE Fit function EP
*EMIN    275.000
*EMAX    360.000
      real PAR( 4)
      data PAR/
     #  .0000000E+00    ,
     #  24.24819        ,
     # -.6136904E-01    ,
     #  6.795246        /
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
 
      fiteff = result
      END
