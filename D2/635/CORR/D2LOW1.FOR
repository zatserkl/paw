*     .. Autogenerated function
      REAL FUNCTION d2low1(x)
*     .. For effs. ID=2
*TYPE Fit function POLI
*EMIN    170.000
*EMAX    260.000
      real PAR( 6)
      data PAR/
     #  4.000000        ,
     #  126647.5        ,
     # -2251.453        ,
     #  15.30001        ,
     # -.4684978E-01    ,
     #  .5446061E-04    /
      ipower = PAR(1)
      p = 0.
*     .. the last non-polinom parameter
      last = 1
      do n=1, ipower+1
         p = p + PAR(last+n)*x**(n-1)
      enddo
      result = p
 
      d2low1 = result
      END
