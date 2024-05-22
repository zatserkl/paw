      call limits
      end
*                                         @METAGS LIMITS
*                                         05/06/98 13:15
*--------------- LIMITS ---------------
*
      SUBROUTINE LIMITS()
      real REALEPS,REALMI,REALMA
      double precision DOUBEPS,DOUBMI,DOUBMA

      real a,ainf,ainfinf
      double precision d,dinf,dinfinf
      a = 1.
  100 a = a/2.
      if (1.+a .GT. 1.) goto 100
      REALEPS = a
      print*, 'REAL epsilon =', REALEPS

      a = 1.
  200 a = a/2.
      if (a/2. .GT. 0.) goto 200
      REALMI = a
      print*, 'REALMI =', REALMI

      a = 1.
  300 a = a*2.
      ainf = a*2.
      ainfinf = ainf*2.
      if (ainfinf .GT. ainf) goto 300
      REALMA = a
      print*, 'REALMA =', REALMA

      d = 1.d0
  400 d = d/2.d0
      if (1.d0+d .GT. 1.d0) goto 400
      DOUBEPS = d
      print*, 'DOUBLE epsilon =', DOUBEPS

      d = 1.d0
  500 d = d/2.d0
      if (d/2.d0 .GT. 0.d0) goto 500
      DOUBMI = d
      print*, 'DOUBMI =', DOUBMI

      d = 1.d0
  600 d = d*2.d0
      dinf = d*2.d0
      dinfinf = dinf*2.d0
      if (dinfinf .GT. dinf) goto 600
      DOUBMA = d
      print*, 'DOUBMA =', DOUBMA
      END
