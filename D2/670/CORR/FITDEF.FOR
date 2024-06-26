*     .. Autogenerated function
      REAL FUNCTION fitdef(x)
*TYPE Fit function POLI
*EMIN    266.677
*EMAX    313.577
      real PAR( 2)
      data PAR/
     #  .0000000E+00    ,
     #  69.93081        /
      ipower = PAR(1)                      
      p = 0.                               
*     .. the last non-polinom parameter  
      last = 1                             
      do n=1, ipower+1                     
         p = p + PAR(last+n)*x**(n-1)      
      enddo                                
      result = p                           
      
      fitdef = result
      END
